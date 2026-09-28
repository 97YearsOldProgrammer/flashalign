// DNA mapping quality: minimap2's mm_set_mapq2 evaluated over whole-query
// colinear chains, with the vote evidence folded in. This is the only producer
// of a DNA MAPQ. The inputs are the committed chain's score and anchor count
// (minimap2's score0 / cnt), each catalogue rival's chain score and vote
// evidence, the rival count n_sub and, with CIGAR output, the realized DP
// scores dp1 / dp2. Branches, the C (int) truncations and the n_sub term
// follow minimap2 (q_coef 40, min_chain_score 40). The differences:
//   - x is the largest per-rival product of vote ratio, chain ratio and DP
//     ratio. The vote ratio discounts a rival only while its chain scores
//     below the winner's; from there the chain ratio is a floor.
//   - The winner's vote is an absolute term: it floors pen_cm and, below
//     kDnaChainMapqVoteCapMinVote votes, caps the Phred.
//   - A rival on the winner's own locus (a shadow, such as an SV flank) is not
//     a rival. A chained rival at another locus that ties the winner exactly
//     gives x = 1.
//   - minimap2's seed-repetitiveness factor (uniq_ratio) is not applied.
//   - Arithmetic is double rather than float, so a value on an integer
//     boundary may truncate differently, and n_sub counts a bounded rival set
//     rather than every secondary.
//
// The formula runs once per block-owning record of the committed family: each
// record is scored on its own candidate's chain, with its own dp1 and
// identity. The divergence contrast (kDnaChainMapqDivergentRecord*) is a
// family-level cap applied beside it.
//
// Two shadow refinements are carried on the evidence and on by default: a
// chained shadow's vote counts for the winner (shadow_vote_credit), and the
// chained shadow slack grows with read length (chain_shadow_read_slack).
//
// HiFi margin rule (hifi_margin_rule, HiFi presets only). At HiFi accuracy the
// ratio form 1 - x*x cannot tell a rival one substitution worse from a tie, so
// the raw ksw2 margin over a realized rival is read against the scoring row:
//   R1 The seat. The realized rivals are the competing alternative, the
//      winner's realized sibling chain (single-block families) and the block's
//      realized catalogue rival (multi-block families). Each has a raw margin
//      dp1_raw - dp2_raw and, where the two records' query spans intersect,
//      a margin m_in over the intersection. The effective margin m_eff is
//      min(raw, m_in) for a whole-read hypothesis and m_in alone for a block
//      rival, which is realized clipped to the block's span. The rival with
//      the smallest m_eff takes the dp2 seat; the others lost by more and
//      enter neither x nor x_other.
//   R2 The guard. The margin decides when m_eff >= A + B and the seat is
//      chained, not a tie, has fewer anchors than the winner, and every other
//      competing rival is weaker than the winner (x_other < 1).
//   R3 The ladder. q = (int)(6.02 * m_eff / a + 0.499), capped at
//      kDnaChainMapqHifiSingleEventMapq when one gap event explains the
//      margin, scaled by min(pen_s1, pen_cm) without the vote floor, then the
//      n_sub term and the clamp.
//   R4 The never-Q0 promotion needs a seat beaten by at least A + B and
//      x_other < 1.
//   R5 Otherwise minimap2's ratio form, with the seat's own score pair as
//      dp1 / dp2 and mapq_alt bounded by m_eff.
//   R6 On a block record of a multi-block family the MAPQ is capped at
//      kDnaChainMapqHifiSplitCapMapq when (a) there is no seat and a competing
//      rival reaches kDnaChainMapqHifiSplitRivalStrength, or (b) the record's
//      own divergence is substitution-dominated: HiFi errors are mostly
//      indels, a wrong repeat copy's differences are not.
// The sibling and block rivals are realized like the alternative into a
// local outcome that only the MAPQ reads; nothing is emitted from them.
//
// The kDnaChainMapqStudy* bits, on under the HiFi presets only, adapt these
// rules to real reads, which carry sample variants and SVs:
//   2   R6(b) fires only when the block's mismatch rate also exceeds the
//       family's cleanest record's, by the divergence contrast's z test.
//   8   The committed record vouches for a chained rival when both ends of
//       the rival's chain lie on the record's alignment (record_shadow).
//   16  A bridged unit belongs to the block whose candidate has the best
//       chain.
//   32  A realized block rival the shadow rule calls the winner's own locus
//       does not take the seat.
//   64  Any record of the family, or two records jointly (an SV's flanks),
//       can vouch as in bit 8.
//   128 The block-rival path skips a candidate the own-locus verdict accepts,
//       before realizing it.
//   256 Each block of a multi-block family also chains its own top rivals, so
//       a smaller flank is compared against the rivals of its own span.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include DNA chain MAPQ"
#endif

#include "context.h" // kDnaMinChainScore

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

// minimap2 constants.
inline constexpr double kDnaChainMapqCoef = 40.0; // q_coef
// minimap2's mask_level: a rival overlapping the committed hypothesis over
// more than half the shorter query span competes with it; less overlap means
// another segment of the same read.
inline constexpr double kDnaChainMapqMaskLevel = 0.5;
// The vote floor under which pen_cm is discounted linearly, and the vote cap
// q <= (int)(4.343 * w / 3).
inline constexpr double kDnaChainMapqVoteFloor = 40.0;
inline constexpr double kDnaChainMapqVoteCapCoef = 4.343;
inline constexpr double kDnaChainMapqVoteCapDiv = 3.0;
// The vote cap applies only while the winner has fewer votes than this; above
// it the chain's own score and anchors carry the absolute evidence.
inline constexpr int kDnaChainMapqVoteCapMinVote = 10;
// Divergence contrast (record_divergence.h): a record whose event divergence
// (PAF de:f) exceeds the family minimum by z > kDnaChainMapqDivergentRecordZ
// is likely on the wrong repeat copy and has its MAPQ capped at
// kDnaChainMapqDivergentRecordMapq. Only records that clear the emission
// floor take part, and the primary is judged like any other record. The event
// denominator is block_len + ambiguities - gap_bases + gap_opens.
inline constexpr double kDnaChainMapqDivergentRecordZ = 8.0;
inline constexpr int kDnaChainMapqDivergentRecordMapq = 3;
// HiFi margin rule (R1-R6 above).
inline constexpr int kDnaChainMapqHifiRuleVersion = 2;
// R2: every competing rival other than the seat must be strictly below this
// strength.
inline constexpr double kDnaChainMapqHifiRivalStrengthMax = 1.0;
// R3: the cap for a margin above one substitution (A + B) and below a
// two-base gap (O1 + 2 E1), which one gap event (the HiFi homopolymer error)
// explains.
inline constexpr int kDnaChainMapqHifiSingleEventMapq = 9;
// R6: the split-read cap; (a) the rival strength that triggers it; (b) the
// minimum event count (mismatches + indel events) and mismatch share. (b)
// applies with or without a seat, since a realized block rival need not be
// the true copy.
inline constexpr int kDnaChainMapqHifiSplitCapMapq = 3;
inline constexpr double kDnaChainMapqHifiSplitRivalStrength = 0.4;
inline constexpr int kDnaChainMapqHifiSubstitutionEvents = 20;
inline constexpr double kDnaChainMapqHifiSubstitutionShare = 0.45;
// Shadow rule, chained rivals: a chained rival is the winner's own locus when
// its chain overlaps the winner's reference span over more than half the
// shorter span, lies within that span extended by the slack on each side, and
// starts and ends within the slack of the winner's diagonals. A repeat copy
// one period away sits on a diagonal of its own. The default slack is the
// interval pad of the whole-query scan window.
inline constexpr int kDnaChainMapqDefaultShadowSlackBp = 512;
// With chain_shadow_read_slack the slack is max(shadow_slack_bp, floor,
// read_len / fraction): a rival chain harvests its own anchors and can reach a
// few kb past the winner's span at an SV flank. The fraction is wider than the
// peak window below because the chain's containment must also absorb that
// overhang.
inline constexpr int kDnaChainMapqChainShadowFloorBp = 2000;
inline constexpr int kDnaChainMapqChainShadowFraction = 5;
// Shadow rule, un-chained rivals: a vote peak on the winner's contig and strand
// within max(floor, read_len / fraction) of the winner's diagonal, with no
// more vote than the winner, is the winner's own locus. The same window as
// the alternative hypothesis's shadow test.
inline constexpr int kDnaChainMapqPeakShadowFloorBp = 2000;
inline constexpr int kDnaChainMapqPeakShadowFraction = 10;
// The winner's sibling chain is realized only when
// max(sib_f2, kDnaMinChainScore) / f1 reaches this; below it the
// ratio form is flat and the DP would change nothing.
inline constexpr double kDnaChainMapqSiblingRealizeMin = 0.5;
// At most this many block-rival realizations per read, blocks taken by query
// length descending, then owner id; the rest record Capped.
inline constexpr int kDnaChainMapqBlockRivalMax = 3;
// The rule bits (header comment).
inline constexpr int kDnaChainMapqStudyRelativeShareCap = 2;
inline constexpr int kDnaChainMapqStudyRecordShadow = 8;
inline constexpr int kDnaChainMapqStudyBridgeOwner = 16;
inline constexpr int kDnaChainMapqStudyNoShadowSeat = 32;
inline constexpr int kDnaChainMapqStudyFamilyVouch = 64;
inline constexpr int kDnaChainMapqStudyLaneVerdict = 128;
inline constexpr int kDnaChainMapqStudyBlockRivals = 256;

// The bits in force under the HiFi presets.
inline constexpr int kDnaChainMapqStudy =
    kDnaChainMapqStudyRelativeShareCap | kDnaChainMapqStudyRecordShadow |
    kDnaChainMapqStudyBridgeOwner | kDnaChainMapqStudyNoShadowSeat |
    kDnaChainMapqStudyFamilyVouch | kDnaChainMapqStudyLaneVerdict |
    kDnaChainMapqStudyBlockRivals;
static_assert(kDnaChainMapqStudy == 506, "the HiFi MAPQ rule set");

// The rule bits for a run: kDnaChainMapqStudy under the HiFi presets, else 0.
inline constexpr int dna_chain_mapq_study_bits(bool hifi_margin_rule) noexcept {
  return hifi_margin_rule ? kDnaChainMapqStudy : 0;
}

// One catalogue rival of the committed hypothesis. A chained rival carries its
// whole-query chain; an un-chained one only its vote evidence and screening
// query span, enough for the vote factor but not for f2 / n_sub.
struct DnaChainMapqRival {
  // candidate and status are not read by the formula.
  int candidate = -1;
  int status = 0;
  int chain_score = 0;
  int chain_anchors = 0;
  int q_begin = -1;
  int q_end = -1;
  int vote_evidence = 0;
  bool chained = false;
  // The evidence's dp2 is this rival's DP score, as minimap2 attributes
  // dp_max2 to one rival; n_sub then counts it once under either rule.
  bool owns_dp2 = false;
  // Locus, for the shadow rule: the vote peak's contig and strand, the
  // chain's contig-local half-open reference span (chained rivals only), and
  // the vote peak's reference-start diagonal.
  int contig = -1;
  bool reverse = false;
  int ref_begin = -1;
  int ref_end = -1;
  std::int64_t peak_diagonal = 0;
  // Bit 8 / 64: the committed alignment maps both ends of this rival's chain
  // to within the chained shadow slack of its reference span, so the rival is
  // part of this placement. The formula treats it as a shadow.
  bool record_shadow = false;
  // Which record vouched: 0 none, 1 the block's own record (bit 8), 2 another
  // record or the family jointly (bit 64). Not read by the formula.
  int vouch = 0;
};

// Why a rival was not realized; None means it was. (B) codes are the sibling
// path's, (A) codes the block-rival path's.
enum class DnaChainMapqRealizeRefusal : std::uint8_t {
  None = 0,
  NotAttempted = 1,
  NoSibling = 2,         // (B) the winner chain binds no rival sibling
  WeakSibling = 3,       // (B) the sibling floor is under the realize floor
  NoRival = 4,           // (A) no catalogue candidate outside the family
  Unchained = 5,         // (A) rivals exist, none has a whole-query chain
  Inadmissible = 6,      // (A) chained rivals are all disjoint / shadow /
                         //     incredible for this block
  Capped = 7,            // (A) an admissible rival lost to the per-read cap
  RealizationFailed = 8, // the fixpoint or the realizer refused it
  OwnLocus = 9,          // (A) judged the read's own locus (bit 128)
};

// One realized rival: its DP result, geometry, divergence and its
// interval-matched score against the record it rivals. R1 reads realized,
// dp2_raw, chain_score, chain_anchors and the im_* fields. Fields are -1 / 0
// when not realized.
struct DnaChainMapqRealizedRival {
  bool realized = false;
  DnaChainMapqRealizeRefusal refusal = DnaChainMapqRealizeRefusal::NotAttempted;
  int candidate = -1; // catalogue id of the realized rival
  // The realization's raw ksw2 decision score, and the realized primary
  // record's own score, which differs by end accounting.
  double dp2_raw = 0.0;
  int rec_score = -1;
  // The realized primary record's forward-query span, contig, strand and
  // contig-local reference span.
  int q_begin = -1;
  int q_end = -1;
  int contig = -1;
  bool reverse = false;
  int ref_begin = -1;
  int ref_end = -1;
  // The realized rival chain's score and anchors.
  int chain_score = 0;
  int chain_anchors = 0;
  // Start diagonal of the rival chain minus that of the rivaled record, in
  // the rival's strand frame; 0 across contigs or strands.
  std::int64_t diag_shift = 0;
  // The realized record's mismatches (-1 without accounting) and its CIGAR
  // event and aligned-column counts.
  int mismatches = -1;
  int ins_events = -1;
  int del_events = -1;
  int aligned_bases = -1;
  // Dual-affine scores (affine_score_over_query_interval) over [im_lo, im_hi),
  // the forward-query intersection of the two records' spans (im_hi < im_lo
  // means disjoint): im_a_* for the rivaled record inside / outside, im_b_*
  // for the rival.
  int im_lo = -1;
  int im_hi = -1;
  long long im_a_inside = -1;
  long long im_a_outside = -1;
  long long im_b_inside = -1;
  long long im_b_outside = -1;
  // Block rivals: 1 = realized clipped to the block record's query span; 2 =
  // the clipped realization failed and the whole-query one stands in. 0
  // otherwise.
  int clip_arm = 0;
  // What the realization cost, accepted or not; both attempts when clip_arm is 2.
  int ksw2_attempts = 0;
  long long estimated_cells = 0;
};

struct DnaChainMapqEvidence {
  int f1 = 0;  // committed whole-query chain score      (mm2 score0 / score)
  int cnt = 0; // committed whole-query chain anchors    (mm2 cnt)
  // The committed candidate's sibling f2 from its own chain partition, a free
  // lower bound on subsc.
  int sib_f2 = 0;
  int winner_q_begin = -1;
  int winner_q_end = -1;
  int winner_vote = 0;
  int read_len = 0;
  // The winner's locus, for the shadow rule; see DnaChainMapqRival.
  int winner_contig = -1;
  bool winner_reverse = false;
  int winner_ref_begin = -1;
  int winner_ref_end = -1;
  std::int64_t winner_peak_diagonal = 0;
  // The chained shadow slack: the whole-query scan window's interval pad.
  int shadow_slack_bp = kDnaChainMapqDefaultShadowSlackBp;
  std::vector<DnaChainMapqRival> rivals;
  // CIGAR output only; 0 otherwise.
  double dp1 = 0.0;
  double dp2 = 0.0;
  // mlen/blen of the realized primary; 1.0 when there is no base-level
  // accounting (map-only), which makes the identity factor inert.
  double identity = 1.0;
  int match_sc = 1; // minimap2 opt->a
  // minimap2's opt->a * 2 + opt->b: a realized rival within this of dp1
  // counts toward n_sub (mm_set_parent).
  int sub_diff = 0;
  // mm_set_parent's `ri->cnt >= rp->cnt` clause: a chained rival with at
  // least the winner's anchors counts toward n_sub.
  bool n_sub_cnt_clause = true;
  // Shadow refinements, set from ResolvedDnaOptions (on by default) and
  // off on a bare struct.
  //   shadow_vote_credit: a chained shadow's vote counts for the winner, in
  //   the absolute term and the strength denominators, since one locus's
  //   vote can be split between twin candidates.
  //   chain_shadow_read_slack: see dna_chain_mapq_chain_shadow_slack.
  bool shadow_vote_credit = false;
  bool chain_shadow_read_slack = false;
  // HiFi margin rule (header comment), HiFi presets only. dp1_raw / dp2_raw
  // are the raw ksw2 decision scores (dp1 / dp2 above may be post-DP
  // prices), read against the scoring row below.
  bool hifi_margin_rule = false;
  double dp1_raw = 0.0;
  double dp2_raw = 0.0;
  int substitution_cost = 0; // A + B
  int gap_open1 = 0;         // O1
  int gap_extend1 = 0;       // E1
  // R1: the alternative's interval margin, winner minus rival after any
  // promotion swap. Valid only when the primaries' query spans intersect.
  bool alternative_im_valid = false;
  double alternative_im_margin = 0.0;
  // R6: the family's block-owning records, and this record's mismatches (-1
  // without accounting) and indel events (-1 without a CIGAR).
  int family_blocks = 1;
  int record_mismatches = -1;
  int record_ins_events = -1;
  int record_del_events = -1;
  // The rule bits (dna_chain_mapq_study_bits); the formula reads only bit 2.
  // For it: this record's aligned bases, and the (mismatches, aligned bases)
  // of the family's cleanest other accounted record, -1 when there is none.
  int study_bits = 0;
  int record_aligned_bases = -1;
  int family_min_mismatches = -1;
  int family_min_aligned = -1;
  // The realized sibling (single-block families) and block rival
  // (multi-block families), read by R1 under hifi_margin_rule only.
  DnaChainMapqRealizedRival sibling;
  DnaChainMapqRealizedRival block_rival;
};

// The formula's verdict on one rival, for the caller.
struct DnaChainMapqRivalVerdict {
  // minimap2's mask_level test on the query spans.
  bool admissible = false;
  // The rival is the winner's own locus.
  bool shadow = false;
  // Not read by the formula: the half-overlap test without containment.
  bool shadow_half_only = false;
  // Chained exact tie at another locus.
  bool tie = false;
  // The rival's contribution to x (0 when it contributes nothing).
  double strength = 0.0;
  // The rival counted toward n_sub. The standalone verdict applies only the
  // anchor rule; dna_chain_mapq adds minimap2's DP rule for the dp2 owner.
  bool counted = false;
};

// Every intermediate the formula forms. Filled only when the caller passes a
// non-null pointer.
struct DnaChainMapqBreakdown {
  int rivals_total = 0;
  int rivals_admissible = 0; // minimap2's mask_level test, shadows included
  int rivals_shadow = 0;     // admissible rivals the shadow rule excluded
  int rivals_tie = 0;
  int f2 = 0;
  int subsc = 0;
  int n_sub = 0;
  double vote_ratio = 0.0; // largest vote_i / w over the competing rivals
  double x_chain = 0.0;    // minimap2's subsc / f1 over the competing rivals
  double x_used = 0.0;     // the x the emitted branch consumed
  double pen_cm = 0.0;     // after every discount, as consumed
  int cap = 60;            // the absolute vote cap as applied
  bool dp_branch = false;
  // HiFi margin rule: the largest strength among competing rivals other than
  // the seat, whether the rule decided the MAPQ, and its Phred before the
  // tail (-1 when it did not).
  double x_other = 0.0;
  bool hifi_margin_applied = false;
  int hifi_margin_mapq = -1;
  // The seat (0 none, 1 alternative, 2 sibling, 3 block rival) and its R1
  // margins, the pen the margin rule uses, and the split-read cap applied
  // (0 none, 1 rival strength, 2 substitution share, 3 both).
  int seat_kind = 0;
  double seat_raw_margin = 0.0;
  bool seat_in_valid = false;
  double seat_in_margin = 0.0;
  double seat_m_eff = 0.0;
  double hifi_pen = 0.0;
  int split_cap = 0;
};

// minimap2's mask_level test with uncov_len = 0: a rival competes when it
// overlaps the committed chain over more than half the shorter query span. A
// missing span counts as unknown, so the rival stays admissible.
bool dna_chain_mapq_admissible(const DnaChainMapqEvidence& evidence,
                               const DnaChainMapqRival& rival) noexcept;

// The chained shadow slack in bp: max(shadow_slack_bp,
// kDnaChainMapqChainShadowFloorBp, read_len / kDnaChainMapqChainShadowFraction)
// with chain_shadow_read_slack, else shadow_slack_bp.
std::int64_t dna_chain_mapq_chain_shadow_slack(
    const DnaChainMapqEvidence& evidence) noexcept;

// The chain part of the shadow rule as geometry: the winner's chain against a
// rival's, each a contig-local reference span, a forward-query span and a
// strand, with `slack` bp of tolerance. The caller checks contig and strand.
// Without query spans on_diagonal stays false; without reference spans
// nothing is set. `shadow` is half && contained && on_diagonal.
struct DnaChainMapqChainShadow {
  bool half = false;
  bool contained = false;
  bool on_diagonal = false;
  bool shadow = false;
};

DnaChainMapqChainShadow dna_chain_mapq_chain_shadow(
    int winner_ref_begin, int winner_ref_end, int winner_q_begin,
    int winner_q_end, bool winner_reverse, int rival_ref_begin,
    int rival_ref_end, int rival_q_begin, int rival_q_end, bool rival_reverse,
    int read_len, std::int64_t slack) noexcept;

// The mask_level test on two forward-query spans, in integer form: they
// compete when they overlap over more than half the shorter. Unlike
// dna_chain_mapq_admissible, a missing span never competes.
bool dna_chain_mapq_spans_compete(int a_begin, int a_end, int b_begin,
                                  int b_end) noexcept;

// The per-rival verdict, without the parts that depend on the whole
// evidence: the DP factor of `strength`, the DP rule of `counted` and the
// shadow vote credit, which dna_chain_mapq adds.
DnaChainMapqRivalVerdict
dna_chain_mapq_rival_verdict(const DnaChainMapqEvidence& evidence,
                             const DnaChainMapqRival& rival) noexcept;

// minimap2's MAPQ for the committed hypothesis, in [0,60]. Returns 0 when
// there is no committed whole-query chain (f1 <= 0). `verdicts`, when
// non-null, receives one entry per evidence.rivals element, in order.
int dna_chain_mapq(const DnaChainMapqEvidence& evidence,
                   DnaChainMapqBreakdown* breakdown = nullptr,
                   std::vector<DnaChainMapqRivalVerdict>* verdicts = nullptr);

} // namespace fa::cpu::lr
