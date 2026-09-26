// The RNA mapping quality: minimap2's mm_set_mapq2 branch structure evaluated over the
// committed fine chain, the coarse catalogue's rivals and, with CIGAR output, the
// realized primary's DP score and identity. A pure function of its evidence.
//
//   map-only (dp1 == 0):  q = pen * 40 * (1 - x) * ln(f1)                          [A]
//   CIGAR, dp2 == 0:      q = identity * pen * 40 * (1 - x) * ln(dp1/a)            [B]
//   CIGAR, dp2 > 0:       q = identity * pen * 40 * (1 - x_c^2) * ln(dp1/a)        [C]
//   pen = min(min(1, 0.01*f1), cnt > 10 ? 1 : 0.1*cnt, clip((qcov/read_len)/tau, 0, 1))
//   x   = max(max(sib_f2, 20)/f1, max_i s_i, v)
//   s_i = (rs_i/rs1) * (chained_i ? max(f_i, 20)/f1 : 1), 1 on an exact tie, over the
//         rivals that pass the mask_level test
//   v   = min(1, damp * rs2/rs1), rs2 the largest catalogue vote at any other locus
//   C only: rho = dp2/dp1, x_c = max(rho * x, x_floor),
//           x_floor = max(max(sib_f2, 20)/f1, max_i {s_i : rival i not realized})
//   C only: q = min(q, (int)(6.02 * 0.6 * identity^2 * (dp1 - dp2)/a + 0.499))
//   q -= (int)(4.343 * ln(n_sub + 1) + 0.499), clamped to [0, 60]
//
// n_sub counts the admissible chained rivals with cnt_i >= cnt, plus the dp2 owner when
// dp1 - dp2 <= sub_diff (minimap2's DP rule).
//
// Differences from minimap2:
//   - min_chain_score is the spliced chain's own min_sc = 20, not 40.
//   - The query-coverage brake replaces the rep_len discount in pen.
//   - v has no minimap2 analogue. It reads the coarse catalogue, so a second locus the
//     read's k-mers vote for counts even if it was never chained; it is damped because
//     nothing tested it. The catalogue is ordered by vote, so rs2 bounds every rival's.
//   - Only a few rivals are realized, so an unrealized rival keeps its undamped
//     strength through x_floor, while rho bounds the realized ones. An
//     unrealized exact tie gives x_c = 1 and q = 0.
//   - The near-tie cap's 6.02 is scaled by kRnaChainMapqBwaCapScale, and the final
//     "mapq 0 -> 1 when dp1 > dp2" floor is not applied.
//   - Rival votes come from the deduplicated catalogue, whose strand-blind window dedup
//     keeps the winner's own opposite-strand mirror out of the rivals.
//   - double and std::log instead of float and logf, so a value on an integer boundary
//     can truncate either way; n_sub is bounded by the catalogue size.
#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include RNA chain MAPQ"
#endif

#include "anchoring/query_span.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr::rna {

// minimap2's q_coef.
inline constexpr double kRnaChainMapqCoef = 40.0;
// The spliced chain's min_sc, standing in for minimap2's min_chain_score in the subsc
// floor.
inline constexpr int kRnaChainMapqMinChainScore = 20;
// Default tau of the query-coverage brake: pen is discounted linearly below an
// anchor-union coverage of tau * read_len. Presets set RnaChainMapqEvidence::qcov_tau.
inline constexpr double kRnaChainMapqQcovTau = 0.7;
// Scale on minimap2's 6.02 in the near-tie cap; the cap's shape, operands and order
// stay minimap2's.
inline constexpr double kRnaChainMapqBwaCapScale = 0.6;

// One catalogue rival of the committed hypothesis. A chained rival carries its fine-chain
// measurement; an unchained one only its vote and screening span (censored evidence is
// not zero evidence).
struct RnaChainMapqRival {
  // Catalogue index; not used by the formula.
  int candidate = -1;
  // Fine-chain measurement (chained rivals only).
  int chain_score = 0;
  int chain_anchors = 0;
  // The rival's chain query span in its own mapping frame (see `reverse`); -1 when
  // censored.
  int q_begin = -1;
  int q_end = -1;
  // The coarse catalogue's rank_score (the RNA vote evidence).
  std::int64_t rank_score = 0;
  bool chained = false;
  // The strand the rival's chain was taken on, the frame of q_begin/q_end: the
  // post-repair locus's `reverse` (a repaired rival's locus is already flipped, so do not
  // XOR again). Not read by rna_chain_mapq.
  bool reverse = false;
  // This rival's realized DP score is evidence.dp2; lets n_sub keep minimap2's per-rival
  // OR.
  bool owns_dp2 = false;
  // Realized with a positive DP maximum, so rho already bounds this rival and it stays
  // out of the undamped x_floor. Several rivals may be realized while one owns dp2; the
  // formula reads `realized_dp || owns_dp2`.
  bool realized_dp = false;
  // This rival's fine chain is the committed chain's own placement seen from a
  // neighbouring catalogue window (same contig and strand, overlapping or continuing
  // reference spans; rna_same_place). The opposite-lane postings of the winner's envelope
  // can form a small window of their own whose harvest re-chains the committed anchors;
  // priced as a rival that is an exact realized tie and a false Q0. A shadow stays in the
  // lifecycle but never supplies dp2 and is inadmissible to the formula (x, x_floor, f2,
  // n_sub, tie). Set in backend.cpp against the committed and the elected chain.
  bool shadow = false;
};

struct RnaChainMapqEvidence {
  int f1 = 0;  // committed fine-chain score            (mm2 score0 / score)
  int cnt = 0; // committed fine-chain anchors          (mm2 cnt)
  // The committed window's second-best chain from the same recurrence, a free lower
  // bound on subsc; <= 0 when absent.
  int sib_f2 = 0;
  int winner_q_begin = -1;
  int winner_q_end = -1;
  // The committed locus's coarse rank_score (rs1), the vote denominator.
  std::int64_t rank1_score = 0;
  // rs2: the largest rank_score over catalogue loci other than the committed one (on an
  // adopted read, the demoted incumbent's). Taken from the catalogue rather than
  // `rivals` so the floor v counts a locus that was never chained. 0 when the catalogue
  // held no other locus.
  std::int64_t rank2_score = 0;
  int read_len = 0;
  // The committed hypothesis's strand, with RnaChainMapqRival::reverse's semantics; the
  // frame of winner_q_begin/winner_q_end. Not read by the formula.
  bool winner_reverse = false;
  // Anchor-union query bases covered by the committed fine chain, the brake's
  // numerator (the map-only PAF matches field).
  std::uint64_t selected_query_covered_bases = 0;
  // The brake's denominator: 0 means read_len (the primary); a second family passes its
  // realized hull's width.
  int qcov_bases = 0;
  std::vector<RnaChainMapqRival> rivals;
  // Realized DP maxima with CIGAR output; 0 for map-only.
  double dp1 = 0.0;
  double dp2 = 0.0;
  // matches/block_len of the realized primary; 1.0 for map-only output, which makes the
  // identity factor inert.
  double identity = 1.0;
  int match_sc = 1; // minimap2 opt->a (cigar_dp_match; 1 in both presets)
  // minimap2 opt->a * 2 + opt->b: 4 (splice) / 6 (splice:hq). The threshold
  // under which the realized rival still counts toward n_sub.
  int sub_diff = 0;
  // With CIGAR output, which chains every catalogue rival, an unchained (vote-only)
  // rival contributes x strength 0 instead of its vote ratio. False for map-only
  // output, which chains no rival. x is deliberately not clamped at 1: x > 1 from a
  // demoted incumbent brakes an adopted read.
  bool unchained_zero = false;
  // The preset calibration values travel with the evidence so the formula stays a pure
  // function of this struct; backend.cpp fills them from ResolvedOptions.
  //
  // Coverage brake threshold, in (0, 1]: pen's third limb is
  // clip((qcov/read_len)/qcov_tau, 0, 1).
  double qcov_tau = kRnaChainMapqQcovTau;
  // Coefficient of the damped disjoint-vote floor v, in [0, 1]; 0 disables it with no
  // float operation.
  double disjoint_vote_damp = 0.0;
};

// The formula's verdict on one rival, for the caller.
struct RnaChainMapqRivalVerdict {
  // minimap2's mask_level test on the query spans.
  bool admissible = false;
  // Chained exact tie (f_i == f1 && cnt_i == cnt): undecidable, x = 1.
  bool tie = false;
  // The rival's contribution to x (0 when it contributes nothing).
  double strength = 0.0;
  // The rival counted toward n_sub. The standalone verdict applies only the
  // anchor rule; rna_chain_mapq adds minimap2's DP rule for the dp2 owner.
  bool counted = false;
};

// Every intermediate the formula forms. Filled only when the caller passes a
// non-null pointer.
struct RnaChainMapqBreakdown {
  int rivals_total = 0;
  int rivals_admissible = 0;
  int rivals_tie = 0;
  int f2 = 0; // best rival or sibling chain score
  int n_sub = 0;
  double vote_ratio = 0.0; // largest rs_i / rs1 over the admissible rivals
  double x_used = 0.0;     // base x, before branch C's damping
  double pen = 0.0;        // after every discount, as consumed
  double qcov_brake = 1.0; // the brake term alone, before the min
  double identity = 1.0;   // as consumed (1.0 for map-only)
  bool dp_branch = false;  // the emitted q read ln(dp1) instead of ln(f1)
  // `branch_c`: the emitted q used (1 - x_c^2). rho = dp2/dp1, x_floor the undamped
  // floor over the unrealized rivals (sibling included), x_c = max(rho * x_used,
  // x_floor). Off branch C: rho = 0, x_floor = the sibling share, x_c = x_used.
  bool branch_c = false;
  double rho = 0.0;
  double x_floor = 0.0;
  double x_c = 0.0;
  // minimap2's BWA-MEM-like near-tie cap, or -1 when not evaluated (no realized dp1/dp2
  // pair); `bwa_cap_bound` says whether it lowered q.
  int bwa_cap = -1;
  bool bwa_cap_bound = false;
};

// minimap2's mask_level test (uncov_len = 0) on the raw oriented spans: a rival competes
// when it overlaps the committed hypothesis over more than half the shorter span. A
// missing span on either side is censored geometry and stays admissible; a shadow never
// is. The forward-frame helpers in query_span.h serve the lifecycle, the chimeric floor
// and rival pricing.
bool rna_chain_mapq_admissible(const RnaChainMapqEvidence& evidence,
                               const RnaChainMapqRival& rival) noexcept;

// The per-rival verdict, except `counted`'s DP rule, which rna_chain_mapq fills.
RnaChainMapqRivalVerdict
rna_chain_mapq_rival_verdict(const RnaChainMapqEvidence& evidence,
                             const RnaChainMapqRival& rival) noexcept;

// The RNA MAPQ for the committed hypothesis, in [0,60]. Returns 0 when there
// is no committed fine chain (f1 <= 0). `verdicts`, when non-null, receives
// one entry per evidence.rivals element, in order.
int rna_chain_mapq(const RnaChainMapqEvidence& evidence,
                   RnaChainMapqBreakdown* breakdown = nullptr,
                   std::vector<RnaChainMapqRivalVerdict>* verdicts = nullptr);

} // namespace fa::cpu::lr::rna
