// DNA mapping quality: minimap2's mm_set_mapq2 over whole-query chains, each record on its own block, written onto the emitted family.
#include "dna.h"

#include "../core/core.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

/**** chain_mapq ****/
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
//   - The margin terms (mapq_alt, and the HiFi rule's guard, ladder and R4
//     below) read dp1 and dp2 before the controller's rank rescale
//     (minimap2's dp_max0, which it prints as ms:i). The ratio, the n_sub DP
//     rule and the dp1 > dp2 never-Q0 promotion read the rescaled pair.
//     minimap2 reads every term, mapq_alt included, on the rescaled pair.
//   - Arithmetic is double rather than float, so a value on an integer
//     boundary may truncate differently, and n_sub counts a bounded rival set
//     rather than every secondary.
//
// The formula runs once per block-owning record of the committed family: each
// record is scored on its own candidate's chain, with its own dp1 and
// identity. The divergence contrast (kDnaChainMapqDivergentRecord*) is a
// family-level cap applied beside it.
//
// The rivals are the catalogue candidates that own no block and the candidate
// dp2 came from, each on its whole-query chain (whole_query_chain). With CIGAR
// output every record the controller realized is a rival too, the secondaries
// it drops included, each on its own region's chain, and dp2 is owned by the
// realized record it came from.
//
// Two shadow refinements are carried on the evidence and on by default: a
// chained shadow's vote counts for the winner (shadow_vote_credit), and the
// chained shadow slack grows with read length (chain_shadow_read_slack).
//
// HiFi margin rule (hifi_margin_rule, HiFi presets only). At HiFi accuracy the
// ratio form 1 - x*x cannot tell a rival one substitution worse from a tie, so
// the DP margin over the rival that owns dp2 is read against the scoring row:
//   R1 The seat: the rival that owns dp2, when it competes. Its margin m_eff
//      is dp1_raw - dp2_raw, the pair before the controller's rank rescale,
//      on which one substitution costs A + B.
//   R2 The guard. The margin decides when m_eff >= A + B and the seat is
//      chained, not a tie, has fewer anchors than the winner, and every other
//      competing rival is weaker than the winner (x_other < 1).
//   R3 The ladder. q = (int)(6.02 * m_eff / a + 0.499), capped at
//      kDnaChainMapqHifiSingleEventMapq when one gap event explains the
//      margin, scaled by min(pen_s1, pen_cm) without the vote floor, then the
//      n_sub term and the clamp.
//   R4 The never-Q0 promotion needs a seat beaten by at least A + B and
//      x_other < 1.
//   R5 Otherwise minimap2's ratio form on the rescaled dp1 / dp2, with its
//      mapq_alt on m_eff.
//   R6 On a block record of a multi-block family the MAPQ is capped at
//      kDnaChainMapqHifiSplitCapMapq when (a) there is no seat and a competing
//      rival reaches kDnaChainMapqHifiSplitRivalStrength, or (b) the record's
//      own divergence is substitution-dominated: HiFi errors are mostly
//      indels, a wrong repeat copy's differences are not.
//
// The kDnaChainMapqStudy* bits, on under the HiFi presets only, adapt these
// rules to real reads, which carry sample variants and SVs:
//   2   R6(b) fires only when the block's mismatch rate also exceeds the
//       family's cleanest record's, by the divergence contrast's z test.
//   8   The committed record vouches for a chained rival when both ends of
//       the rival's chain lie on the record's alignment (record_shadow).
//   64  Any record of the family, or two records jointly (an SV's flanks),
//       can vouch as in bit 8.
//   256 Each block of a multi-block family also chains its own top rivals, so
//       a smaller flank is compared against the rivals of its own span.

namespace fa::cpu::lr {

namespace {

// A rival scoring under min_chain_score counts as min_chain_score, as
// minimap2's subsc.
double chain_ratio(const DnaChainMapqEvidence& evidence,
                   const DnaChainMapqRival& rival) noexcept {
  return static_cast<double>(std::max(rival.chain_score, evidence.min_chain_score)) /
         static_cast<double>(evidence.f1);
}

bool same_locus(const DnaChainMapqEvidence& evidence,
                const DnaChainMapqRival& rival) noexcept {
  return evidence.winner_contig >= 0 &&
         rival.contig == evidence.winner_contig &&
         rival.reverse == evidence.winner_reverse;
}

// The start and end diagonals (reference minus oriented query) of a chain
// given in forward-query coordinates.
void chain_diagonals(int ref_begin, int ref_end, int q_begin, int q_end,
                     bool reverse, int read_len, std::int64_t& start,
                     std::int64_t& end) noexcept {
  const std::int64_t oriented_begin = reverse ? read_len - q_end : q_begin;
  const std::int64_t oriented_end = reverse ? read_len - q_begin : q_end;
  start = static_cast<std::int64_t>(ref_begin) - oriented_begin;
  end = static_cast<std::int64_t>(ref_end) - oriented_end;
}

// The shadow rule. A chained rival is judged by its whole-query chain's
// reference span and diagonals; an un-chained one by its vote peak's diagonal.
void shadow_verdict(const DnaChainMapqEvidence& evidence,
                    const DnaChainMapqRival& rival,
                    DnaChainMapqRivalVerdict& verdict) noexcept {
  verdict.shadow = false;
  verdict.shadow_half_only = false;
  if (!same_locus(evidence, rival))
    return;
  // The committed alignment explains the rival's chain (bits 8 / 64), even
  // when it lies off the winner chain's diagonals, as an SV flank does.
  if (rival.record_shadow) {
    verdict.shadow = true;
    verdict.shadow_half_only = true;
    return;
  }
  const bool rival_span =
      rival.ref_begin >= 0 && rival.ref_end > rival.ref_begin;
  const bool winner_span = evidence.winner_ref_begin >= 0 &&
                           evidence.winner_ref_end > evidence.winner_ref_begin;
  if (rival.chained && rival_span && winner_span) {
    const DnaChainMapqChainShadow geometry = dna_chain_mapq_chain_shadow(
        evidence.winner_ref_begin, evidence.winner_ref_end,
        evidence.winner_q_begin, evidence.winner_q_end,
        evidence.winner_reverse, rival.ref_begin, rival.ref_end, rival.q_begin,
        rival.q_end, rival.reverse, evidence.read_len,
        dna_chain_mapq_chain_shadow_slack(evidence));
    verdict.shadow_half_only = geometry.half;
    verdict.shadow = geometry.shadow;
    return;
  }
  if (rival.chained)
    return; // a chain without a span is censored, not shadow
  // Un-chained peak: a shadow when within the window of the winner's diagonal
  // and with no more vote than the winner. A nearby peak that out-votes the
  // winner stays a rival and drives x >= 1.
  const std::int64_t window = std::max<std::int64_t>(
      kDnaChainMapqPeakShadowFloorBp,
      std::max(0, evidence.read_len) / kDnaChainMapqPeakShadowFraction);
  const std::int64_t offset =
      std::llabs(rival.peak_diagonal - evidence.winner_peak_diagonal);
  if (offset <= window && rival.vote_evidence <= evidence.winner_vote) {
    verdict.shadow = true;
    verdict.shadow_half_only = true;
  }
}

// A non-tie rival's strength against a winner with `winner_vote` votes.
double rival_strength(const DnaChainMapqEvidence& evidence,
                      const DnaChainMapqRival& rival,
                      int winner_vote) noexcept {
  const double denominator = static_cast<double>(std::max(1, winner_vote));
  double strength = static_cast<double>(rival.vote_evidence) / denominator;
  if (rival.chained) {
    const double ratio = chain_ratio(evidence, rival);
    strength *= ratio;
    // The vote discounts a rival only while its whole-query chain scores below
    // the winner's; from the winner's score up, the chain ratio is a floor.
    if (ratio >= 1.0)
      strength = std::max(strength, ratio);
  }
  return strength;
}

// The rival holding the dp2 seat (R1), or none.
struct DnaChainMapqSeat {
  int kind = 0;         // 0 none, 1 the dp2 owner
  std::size_t slot = 0; // index into evidence.rivals when kind is 1
  bool chained = false;
  int chain_anchors = 0;
  bool tie = false;
};

} // namespace

std::int64_t dna_chain_mapq_chain_shadow_slack(
    const DnaChainMapqEvidence& evidence) noexcept {
  std::int64_t slack = std::max(0, evidence.shadow_slack_bp);
  if (evidence.chain_shadow_read_slack)
    slack = std::max(
        {slack, static_cast<std::int64_t>(kDnaChainMapqChainShadowFloorBp),
         static_cast<std::int64_t>(std::max(0, evidence.read_len) /
                                   kDnaChainMapqChainShadowFraction)});
  return slack;
}

DnaChainMapqChainShadow dna_chain_mapq_chain_shadow(
    int winner_ref_begin, int winner_ref_end, int winner_q_begin,
    int winner_q_end, bool winner_reverse, int rival_ref_begin,
    int rival_ref_end, int rival_q_begin, int rival_q_end, bool rival_reverse,
    int read_len, std::int64_t slack) noexcept {
  DnaChainMapqChainShadow out;
  if (rival_ref_begin < 0 || rival_ref_end <= rival_ref_begin ||
      winner_ref_begin < 0 || winner_ref_end <= winner_ref_begin)
    return out;
  const std::int64_t overlap =
      static_cast<std::int64_t>(std::min(rival_ref_end, winner_ref_end)) -
      static_cast<std::int64_t>(std::max(rival_ref_begin, winner_ref_begin));
  const std::int64_t shorter =
      std::min(static_cast<std::int64_t>(rival_ref_end - rival_ref_begin),
               static_cast<std::int64_t>(winner_ref_end - winner_ref_begin));
  out.half = overlap > 0 && static_cast<double>(overlap) >
                                kDnaChainMapqMaskLevel *
                                    static_cast<double>(shorter);
  out.contained = rival_ref_begin >= winner_ref_begin - slack &&
                  rival_ref_end <= winner_ref_end + slack;
  // A shadow starts and ends on the winner's diagonals; a local repeat copy
  // has a diagonal of its own. Without query spans there are no diagonals.
  const bool rival_q = rival_q_begin >= 0 && rival_q_end > rival_q_begin;
  const bool winner_q = winner_q_begin >= 0 && winner_q_end > winner_q_begin;
  if (rival_q && winner_q) {
    std::int64_t rival_start = 0, rival_end = 0, winner_start = 0,
                 winner_end = 0;
    chain_diagonals(rival_ref_begin, rival_ref_end, rival_q_begin, rival_q_end,
                    rival_reverse, read_len, rival_start, rival_end);
    chain_diagonals(winner_ref_begin, winner_ref_end, winner_q_begin,
                    winner_q_end, winner_reverse, read_len, winner_start,
                    winner_end);
    const std::int64_t low = std::min(winner_start, winner_end) - slack;
    const std::int64_t high = std::max(winner_start, winner_end) + slack;
    out.on_diagonal = rival_start >= low && rival_start <= high &&
                      rival_end >= low && rival_end <= high;
  }
  out.shadow = out.half && out.contained && out.on_diagonal;
  return out;
}

namespace {

// minimap2's mask_level test with uncov_len = 0: a rival competes when it
// overlaps the committed chain over more than half the shorter query span. A
// missing span counts as unknown, so the rival stays admissible.
bool dna_chain_mapq_admissible(const DnaChainMapqEvidence& evidence,
                               const DnaChainMapqRival& rival) noexcept {
  const bool rival_span = rival.q_begin >= 0 && rival.q_end > rival.q_begin;
  const bool winner_span = evidence.winner_q_begin >= 0 &&
                           evidence.winner_q_end > evidence.winner_q_begin;
  if (!rival_span || !winner_span)
    return true;
  const int overlap = std::min(rival.q_end, evidence.winner_q_end) -
                      std::max(rival.q_begin, evidence.winner_q_begin);
  if (overlap <= 0)
    return false;
  const int shorter = std::min(rival.q_end - rival.q_begin,
                               evidence.winner_q_end - evidence.winner_q_begin);
  return static_cast<double>(overlap) >
         kDnaChainMapqMaskLevel * static_cast<double>(shorter);
}

} // namespace

DnaChainMapqRivalVerdict
dna_chain_mapq_rival_verdict(const DnaChainMapqEvidence& evidence,
                             const DnaChainMapqRival& rival) noexcept {
  DnaChainMapqRivalVerdict verdict;
  verdict.admissible = dna_chain_mapq_admissible(evidence, rival);
  shadow_verdict(evidence, rival, verdict);
  const bool competing = verdict.admissible && !verdict.shadow;
  if (!competing || evidence.f1 <= 0)
    return verdict;
  verdict.tie = rival.chained && rival.chain_score == evidence.f1 &&
                rival.chain_anchors == evidence.cnt;
  if (verdict.tie) {
    verdict.strength = 1.0;
  } else {
    verdict.strength = rival_strength(evidence, rival, evidence.winner_vote);
  }
  verdict.counted = rival.chained && evidence.n_sub_cnt_clause &&
                    rival.chain_anchors >= evidence.cnt;
  return verdict;
}

int dna_chain_mapq(const DnaChainMapqEvidence& evidence,
                   DnaChainMapqBreakdown* breakdown,
                   std::vector<DnaChainMapqRivalVerdict>* verdicts) {
  if (breakdown != nullptr)
    *breakdown = DnaChainMapqBreakdown();
  if (verdicts != nullptr) {
    verdicts->clear();
    verdicts->reserve(evidence.rivals.size());
  }
  if (evidence.f1 <= 0) {
    if (verdicts != nullptr)
      for (const DnaChainMapqRival& rival : evidence.rivals)
        verdicts->push_back(dna_chain_mapq_rival_verdict(evidence, rival));
    return 0;
  }

  const double dp1 = evidence.dp1;
  const bool hifi = evidence.hifi_margin_rule;

  // First pass: per-rival verdicts. dp2 counts only when its owner competes
  // (admissible and not a shadow), as minimap2 takes dp_max2 from a hit that
  // survived mm_set_parent.
  std::vector<DnaChainMapqRivalVerdict> local;
  std::vector<DnaChainMapqRivalVerdict>& judged =
      verdicts != nullptr ? *verdicts : local;
  bool alternative_competes = false;
  std::size_t alternative_slot = 0;
  for (std::size_t index = 0; index < evidence.rivals.size(); ++index) {
    const DnaChainMapqRival& rival = evidence.rivals[index];
    judged.push_back(dna_chain_mapq_rival_verdict(evidence, rival));
    const DnaChainMapqRivalVerdict& verdict = judged.back();
    if (verdict.admissible && !verdict.shadow && rival.owns_dp2) {
      alternative_competes = true;
      alternative_slot = index;
    }
  }
  alternative_competes =
      alternative_competes && dp1 > 0.0 && evidence.dp2 > 0.0;

  // shadow_vote_credit: the winner's vote is the largest among itself and its
  // chained or record-explained shadows, since one locus's vote can split
  // between twin candidates. Un-chained shadows are not credited.
  int effective_vote = evidence.winner_vote;
  if (evidence.shadow_vote_credit) {
    for (std::size_t index = 0; index < evidence.rivals.size(); ++index)
      if (judged[index].shadow && (evidence.rivals[index].chained ||
                                   evidence.rivals[index].record_shadow))
        effective_vote =
            std::max(effective_vote, evidence.rivals[index].vote_evidence);
    if (effective_vote != evidence.winner_vote)
      for (std::size_t index = 0; index < evidence.rivals.size(); ++index)
        if (!judged[index].tie)
          judged[index].strength =
              rival_strength(evidence, evidence.rivals[index], effective_vote);
  }

  // R1: the seat is the dp2 owner.
  const double f1 = static_cast<double>(evidence.f1);
  DnaChainMapqSeat seat;
  if (alternative_competes) {
    const DnaChainMapqRival& rival = evidence.rivals[alternative_slot];
    seat.kind = 1;
    seat.slot = alternative_slot;
    seat.chained = rival.chained;
    seat.chain_anchors = rival.chain_anchors;
    seat.tie = judged[alternative_slot].tie;
  }
  const bool dp_branch = seat.kind != 0;
  // dp2 is the seat's, or nothing.
  const double dp2 = alternative_competes ? evidence.dp2 : 0.0;
  // The seat's margin on the rescaled pair, which the n_sub DP rule reads.
  const double m_eff = dp_branch ? dp1 - dp2 : 0.0;
  // R1: the seat's margin on the pair before the rank rescale, which the
  // guard, the ladder and mapq_alt read. x, the n_sub DP rule and the
  // dp1 > dp2 promotion keep the rescaled pair.
  const double dp1_margin = evidence.dp1_raw;
  const double dp2_margin = alternative_competes ? evidence.dp2_raw : 0.0;
  const double m_margin = dp_branch ? dp1_margin - dp2_margin : 0.0;

  // Second pass: f2, n_sub and x over the competing rivals. As minimap2's
  // mm_set_parent, a rival counts toward n_sub once if it meets the anchor
  // rule or the DP rule; the DP rule applies to the seat, on its own pair.
  int f2 = evidence.sib_f2;
  // A record priced on its block part starts from the part's own n_sub.
  int n_sub = evidence.part_priced ? evidence.part_n_sub : 0;
  int admissible_rivals = 0;
  int shadow_rivals = 0;
  int tie_rivals = 0;
  double vote_ratio = 0.0;
  double max_competing_strength = 0.0;
  const double winner_vote = static_cast<double>(std::max(1, effective_vote));
  // The winner's sibling path seeds x with its chain ratio.
  const double sibling_floor =
      static_cast<double>(std::max(evidence.sib_f2, evidence.min_chain_score)) / f1;
  double x = sibling_floor;
  // x over the competing rivals other than the seat.
  double x_other = x;
  for (std::size_t index = 0; index < evidence.rivals.size(); ++index) {
    const DnaChainMapqRival& rival = evidence.rivals[index];
    DnaChainMapqRivalVerdict& verdict = judged[index];
    if (!verdict.admissible)
      continue;
    ++admissible_rivals;
    if (verdict.shadow) {
      ++shadow_rivals;
      continue;
    }
    if (verdict.tie)
      ++tie_rivals;
    vote_ratio = std::max(vote_ratio, static_cast<double>(rival.vote_evidence) /
                                          winner_vote);
    max_competing_strength = std::max(max_competing_strength, verdict.strength);
    const bool is_seat = seat.kind != 0 && seat.slot == index;
    // R1: a dp2 owner that is not the seat enters neither x nor x_other.
    const bool lost_seat = !is_seat && rival.owns_dp2 && alternative_competes;
    if (rival.chained) {
      f2 = std::max(f2, rival.chain_score);
      bool counts =
          evidence.n_sub_cnt_clause && rival.chain_anchors >= evidence.cnt;
      if (is_seat && m_eff <= static_cast<double>(evidence.sub_diff))
        counts = true;
      verdict.counted = counts;
      if (counts)
        ++n_sub;
    }
    if (is_seat) {
      // The DP ratio discounts the seat once; a tie is not discounted.
      if (!verdict.tie)
        verdict.strength *= dp2 / dp1;
      x = std::max(x, verdict.strength);
    } else if (!lost_seat) {
      x = std::max(x, verdict.strength);
      x_other = std::max(x_other, verdict.strength);
    }
  }
  // minimap2's own x (chain ratio alone), for the breakdown.
  const int subsc = std::max(f2, evidence.min_chain_score);
  const double x_chain = static_cast<double>(subsc) / f1;

  // minimap2's pen_s1, without its uniq_ratio factor.
  const double pen_s1 =
      evidence.f1 > 100 ? 1.0 : 0.01 * static_cast<double>(evidence.f1);
  double pen_cm =
      evidence.cnt > 10 ? 1.0 : 0.1 * static_cast<double>(evidence.cnt);
  pen_cm = std::min(pen_s1, pen_cm);
  // A record priced on its block part: minimap2's pen on the part's score and
  // anchors, as mm_set_mapq2 prices a split chain.
  if (evidence.part_priced)
    pen_cm = std::min(
        evidence.part_score > 100
            ? 1.0
            : 0.01 * static_cast<double>(evidence.part_score),
        evidence.part_anchors > 10
            ? 1.0
            : 0.1 * static_cast<double>(evidence.part_anchors));
  // R3: the margin rule's pen has no vote floor.
  const double pen_margin = pen_cm;
  const double w_abs = static_cast<double>(std::max(0, effective_vote));
  // The vote floor: a winner the vote barely supports is discounted as
  // minimap2 discounts a short chain.
  pen_cm = std::min(pen_cm, std::min(1.0, w_abs / kDnaChainMapqVoteFloor));

  const double match_sc = static_cast<double>(std::max(1, evidence.match_sc));
  // R2. `margin_beats_rival` also gates the never-Q0 promotion below.
  const bool margin_beats_rival =
      hifi && dp_branch && dp1_margin > 0.0 && dp2_margin > 0.0 &&
      evidence.substitution_cost > 0 &&
      m_margin >= static_cast<double>(evidence.substitution_cost);
  const bool margin_rule = margin_beats_rival && seat.chained && !seat.tie &&
                           seat.chain_anchors < evidence.cnt &&
                           x_other < kDnaChainMapqHifiRivalStrengthMax;
  int hifi_margin_mapq = -1;
  int mapq;
  if (margin_rule) {
    // R3: BWA-MEM's margin Phred (minimap2's mapq_alt).
    const double two_base_gap =
        static_cast<double>(evidence.gap_open1 + 2 * evidence.gap_extend1);
    int margin_q = truncate_to_int(6.02 * m_margin / match_sc + 0.499);
    if (m_margin > static_cast<double>(evidence.substitution_cost) &&
        m_margin < two_base_gap)
      margin_q = std::min(margin_q, kDnaChainMapqHifiSingleEventMapq);
    hifi_margin_mapq = margin_q;
    mapq = truncate_to_int(pen_margin * static_cast<double>(margin_q));
  } else if (dp_branch) {
    // minimap2's DP branch; the DP ratio is already inside the seat's
    // strength.
    mapq = truncate_to_int(evidence.identity * pen_cm * kDnaChainMapqCoef *
                           (1.0 - x * x) * std::log(dp1 / match_sc));
    // mapq_alt on the seat's margin (R5).
    const int mapq_alt = truncate_to_int(
        6.02 * evidence.identity * evidence.identity * m_margin / match_sc +
        0.499);
    mapq = std::min(mapq, mapq_alt);
  } else if (dp1 > 0.0) {
    // Aligned, but no rival DP score.
    mapq = truncate_to_int(evidence.identity * pen_cm * kDnaChainMapqCoef *
                           (1.0 - x) * std::log(dp1 / match_sc));
  } else if (evidence.part_priced) {
    // Map-only, on the block part's own score.
    mapq = truncate_to_int(pen_cm * kDnaChainMapqCoef * (1.0 - x) *
                           std::log(static_cast<double>(evidence.part_score)));
  } else {
    // Map-only: scored on the chain alone.
    mapq =
        truncate_to_int(pen_cm * kDnaChainMapqCoef * (1.0 - x) * std::log(f1));
  }
  mapq -= truncate_to_int(4.343 * std::log(static_cast<double>(n_sub) + 1.0) +
                          0.499);
  mapq = std::max(0, mapq);
  mapq = std::min(60, mapq);
  // An aligned hit that beat its rival is never Q0; under the HiFi rule (R4)
  // the seat must have lost by at least one substitution and x_other < 1.
  const bool beat_rival =
      hifi ? (margin_beats_rival && x_other < 1.0) : (dp1 > 0.0 && dp1 > dp2);
  if (beat_rival && mapq == 0)
    mapq = 1;
  // The vote cap, only while the vote is below kDnaChainMapqVoteCapMinVote.
  const int cap =
      w_abs < static_cast<double>(kDnaChainMapqVoteCapMinVote)
          ? std::max(0, truncate_to_int(kDnaChainMapqVoteCapCoef * w_abs /
                                        kDnaChainMapqVoteCapDiv))
          : 60;
  mapq = std::min(mapq, cap);
  // R6: (a) no seat and a strong competing rival; (b) a
  // substitution-dominated record, seat or not.
  int split_cap = 0;
  if (hifi && evidence.family_blocks >= 2) {
    if (!dp_branch &&
        max_competing_strength >= kDnaChainMapqHifiSplitRivalStrength)
      split_cap |= 1;
    if (evidence.record_mismatches >= 0 && evidence.record_ins_events >= 0 &&
        evidence.record_del_events >= 0) {
      const int events = evidence.record_mismatches +
                         evidence.record_ins_events +
                         evidence.record_del_events;
      bool share_dominated = events >= kDnaChainMapqHifiSubstitutionEvents &&
                             static_cast<double>(evidence.record_mismatches) /
                                     static_cast<double>(events) >=
                                 kDnaChainMapqHifiSubstitutionShare;
      // Bit 2: also require the block's mismatch rate to exceed the family's
      // cleanest record's by the divergence contrast's two-sample test.
      // Without such a record there is no cap.
      if (share_dominated &&
          (evidence.study_bits & kDnaChainMapqStudyRelativeShareCap) != 0) {
        share_dominated = false;
        if (evidence.record_aligned_bases > 0 &&
            evidence.family_min_aligned > 0 &&
            evidence.family_min_mismatches >= 0) {
          const double length =
              static_cast<double>(evidence.record_aligned_bases);
          const double length_min =
              static_cast<double>(evidence.family_min_aligned);
          const double rate =
              static_cast<double>(evidence.record_mismatches) / length;
          const double rate_min =
              static_cast<double>(evidence.family_min_mismatches) / length_min;
          const double pooled =
              static_cast<double>(evidence.record_mismatches +
                                  evidence.family_min_mismatches) /
              (length + length_min);
          const double variance = std::max(pooled * (1.0 - pooled), 1e-12) *
                                  (1.0 / length + 1.0 / length_min);
          share_dominated = (rate - rate_min) / std::sqrt(variance) >
                            kDnaChainMapqDivergentRecordZ;
        }
      }
      if (share_dominated)
        split_cap |= 2;
    }
    if (split_cap != 0)
      mapq = std::min(mapq, kDnaChainMapqHifiSplitCapMapq);
  }

  if (breakdown != nullptr) {
    breakdown->rivals_total = static_cast<int>(evidence.rivals.size());
    breakdown->rivals_admissible = admissible_rivals;
    breakdown->rivals_shadow = shadow_rivals;
    breakdown->rivals_tie = tie_rivals;
    breakdown->f2 = f2;
    breakdown->subsc = subsc;
    breakdown->n_sub = n_sub;
    breakdown->vote_ratio = vote_ratio;
    breakdown->x_chain = x_chain;
    breakdown->x_used = x;
    breakdown->pen_cm = pen_cm;
    breakdown->cap = cap;
    breakdown->dp_branch = dp_branch;
    breakdown->x_other = x_other;
    breakdown->hifi_margin_applied = margin_rule;
    breakdown->hifi_margin_mapq = hifi_margin_mapq;
    breakdown->seat_kind = seat.kind;
    breakdown->seat_m_eff = m_eff;
    breakdown->hifi_pen = pen_margin;
    breakdown->split_cap = split_cap;
  }
  return mapq;
}

} // namespace fa::cpu::lr

/**** record_divergence ****/

namespace fa::cpu::lr {
namespace dna {

void dna_divergence_contrast(DnaFamilyDivergence& family,
                             const std::vector<double>& denominators) {
  family.accounted = 0;
  std::size_t lowest = family.records.size();
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    if (!family.records[index].accounted)
      continue;
    if (lowest == family.records.size() ||
        family.records[index].divergence < family.records[lowest].divergence)
      lowest = index;
    ++family.accounted;
  }
  if (family.accounted < 2)
    return;
  const double de_min = family.records[lowest].divergence;
  const double length_min = denominators[lowest];
  const double events_min = de_min * length_min;
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    DnaRecordDivergence& record = family.records[index];
    if (!record.accounted || index == lowest)
      continue;
    const double length = denominators[index];
    const double events = record.divergence * length;
    const double pooled = (events + events_min) / (length + length_min);
    const double variance = std::max(pooled * (1.0 - pooled), 1e-12) *
                            (1.0 / length + 1.0 / length_min);
    record.z = (record.divergence - de_min) / std::sqrt(variance);
    record.divergent = record.z > kDnaChainMapqDivergentRecordZ;
  }
}

} // namespace dna
} // namespace fa::cpu::lr

/**** mapq_routing ****/
// Writes the scored MAPQs onto the records of the emitted DNA family.

namespace fa::cpu::lr {
namespace {

bool inversion_record(const AlignResult& record) noexcept {
  return record.origin == AlignmentOrigin::DnaLocalInversion;
}

// As minimap2's mm_set_inv_mapq: records are ordered by (contig, position)
// and an inversion record between two others takes the smaller of their
// MAPQs. Secondaries take no part.
void apply_inversion_mapq(dna::Result& realized) {
  std::vector<AlignResult*> ordered;
  ordered.reserve(1 + realized.supplementary.size());
  bool has_inversion = false;
  ordered.push_back(static_cast<AlignResult*>(&realized));
  for (AlignResult& supplementary : realized.supplementary) {
    has_inversion = has_inversion || inversion_record(supplementary);
    ordered.push_back(&supplementary);
  }
  if (!has_inversion || ordered.size() < 3)
    return;
  std::stable_sort(ordered.begin(), ordered.end(),
                   [](const AlignResult* left, const AlignResult* right) {
                     return std::tie(left->chromosome, left->pos) <
                            std::tie(right->chromosome, right->pos);
                   });
  for (std::size_t index = 1; index + 1 < ordered.size(); ++index) {
    if (!inversion_record(*ordered[index]))
      continue;
    ordered[index]->mapq =
        std::min(ordered[index - 1]->mapq, ordered[index + 1]->mapq);
  }
}

// `supplementary_mapq` is parallel to realized.supplementary: an entry >= 0 is
// that record's own MAPQ, -1 inherits the primary's (records that own no
// block). Secondaries get 0, and an inversion middle gets min of its flanks'
// MAPQs, as minimap2's mm_set_inv_mapq.
void route_dna_mapq(const DnaPlacementFamily& catalogue,
                    ::fa::cpu::voting::CandidateId primary_candidate, int mapq,
                    const std::vector<int>& supplementary_mapq,
                    dna::Result& realized) {
  if (!realized.mapped())
    return;

  realized.mapq = mapq;

  const DnaPlacementCandidate* candidate = catalogue.find(primary_candidate);
  realized.median_occurrence =
      candidate == nullptr ? 0 : candidate->peak.anchor.median_occurrence;
  for (std::size_t index = 0; index < realized.supplementary.size(); ++index) {
    AlignResult& supplementary = realized.supplementary[index];
    const int own =
        index < supplementary_mapq.size() ? supplementary_mapq[index] : -1;
    supplementary.mapq =
        inversion_record(supplementary) ? 0 : own >= 0 ? own : mapq;
  }
  for (AlignResult& secondary : realized.secondary) {
    secondary.mapq = 0;
    for (AlignResult& supplementary : secondary.supplementary)
      supplementary.mapq = 0;
  }
  // Last, as in minimap2's mm_set_mapq2.
  apply_inversion_mapq(realized);
}

} // namespace

} // namespace fa::cpu::lr

/**** evidence ****/
// Each record's MAPQ evidence: its chain and the catalogue's rivals, the
// committed family's diagonal maps for the shadow rules, its census for the
// split-read caps and the divergence contrast's inputs.

namespace fa::cpu::lr {
namespace dna {

namespace {

// The chain MAPQ rule bits in force for this run.
int chain_mapq_study_bits(const DnaContext& dctx) noexcept {
  return dna_chain_mapq_study_bits(dctx.opts.chain_mapq_hifi_margin);
}

// One catalogue rival's MAPQ evidence (see whole_query_chain; `realized` for
// the dp_max2 owner). A rival with no accepted whole-query chain keeps only
// its vote and its record's query span.
DnaChainMapqRival
chain_mapq_rival_from_placement(const DnaPlacementFamily& catalogue,
                                const DnaPlacementChainingResult& placement,
                                ::fa::cpu::voting::CandidateId id,
                                bool realized = false) {
  DnaChainMapqRival rival;
  rival.candidate = id;
  if (const DnaPlacementCandidate* candidate = catalogue.find(id)) {
    rival.vote_evidence = candidate->vote_evidence;
    rival.contig = candidate->peak.chr;
    rival.reverse = candidate->peak.is_rc;
    rival.peak_diagonal = candidate->peak.raw_ref_start;
  }
  const DnaPlacementCandidateChain* whole =
      placement.whole_query_chain(id, realized);
  const DnaPlacementCandidateChain* record = placement.find(id);
  if (whole != nullptr)
    rival.status = static_cast<int>(whole->status);
  else if (record != nullptr)
    rival.status = static_cast<int>(record->status);
  if (whole != nullptr && whole->status == DnaPlacementChainStatus::Accepted &&
      whole->chain_anchors > 0) {
    rival.chained = true;
    rival.chain_score = whole->chain_score;
    rival.chain_anchors = whole->chain_anchors;
    rival.q_begin = whole->forward_query_begin;
    rival.q_end = whole->forward_query_end;
    rival.ref_begin = whole->reference_begin;
    rival.ref_end = whole->reference_end;
  } else if (record != nullptr) {
    rival.q_begin = record->forward_query_begin;
    rival.q_end = record->forward_query_end;
  }
  return rival;
}

// The MAPQ evidence for one block owner of the committed family: the primary,
// or a supplementary scored on its own candidate's chain. Every owner sees the
// same rivals: candidates that own no selected block and are not the winner.
// Only the primary has a dp2 owner.
DnaChainMapqEvidence build_chain_mapq_evidence(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement, bool placement_ran,
    const DnaPlacementCandidateChain* winner_chain,
    ::fa::cpu::voting::CandidateId winner_candidate,
    ::fa::cpu::voting::CandidateId dp2_owner, double dp1, double dp2,
    double identity, int read_len) {
  DnaChainMapqEvidence evidence;
  evidence.min_chain_score = dctx.opts.min_chain_score;
  evidence.read_len = read_len;
  evidence.match_sc = dctx.opts.cigar_dp_match;
  // As minimap2: sub_diff = opt->a * 2 + opt->b.
  evidence.sub_diff =
      2 * dctx.opts.cigar_dp_match + dctx.opts.cigar_dp_mismatch;
  evidence.dp1 = dp1;
  evidence.dp2 = dp2;
  evidence.identity = identity;
  // HiFi margin rule: the scoring row its margin is read against.
  evidence.hifi_margin_rule = dctx.opts.chain_mapq_hifi_margin;
  evidence.substitution_cost =
      dctx.opts.cigar_dp_match + dctx.opts.cigar_dp_mismatch;
  evidence.gap_open1 = dctx.opts.cigar_dp_gap_open1;
  evidence.gap_extend1 = dctx.opts.cigar_dp_gap_extend1;
  evidence.shadow_slack_bp = dctx.opts.cigar_local_interval_anchor_interval_pad;
  // The shadow refinements (dna.h), on by default.
  static_assert(ResolvedDnaOptions{}.dna_chain_mapq_shadow_vote_credit &&
                    ResolvedDnaOptions{}.dna_chain_mapq_chain_shadow_read_slack,
                "the dense-chain MAPQ locus rules are on by default");
  evidence.shadow_vote_credit = dctx.opts.dna_chain_mapq_shadow_vote_credit;
  evidence.chain_shadow_read_slack =
      dctx.opts.dna_chain_mapq_chain_shadow_read_slack;
  evidence.study_bits = chain_mapq_study_bits(dctx);
  if (winner_chain != nullptr) {
    evidence.f1 = winner_chain->chain_score;
    evidence.cnt = winner_chain->chain_anchors;
    evidence.sib_f2 = winner_chain->rival_chain_score;
    evidence.winner_q_begin = winner_chain->forward_query_begin;
    evidence.winner_q_end = winner_chain->forward_query_end;
    evidence.winner_ref_begin = winner_chain->reference_begin;
    evidence.winner_ref_end = winner_chain->reference_end;
  }
  if (const DnaPlacementCandidate* winner = catalogue.find(winner_candidate)) {
    evidence.winner_vote = winner->vote_evidence;
    evidence.winner_contig = winner->peak.chr;
    evidence.winner_reverse = winner->peak.is_rc;
    evidence.winner_peak_diagonal = winner->peak.raw_ref_start;
  }
  if (!placement_ran)
    return evidence;

  std::vector<::fa::cpu::voting::CandidateId> ids;
  ids.reserve(catalogue.candidates.size());
  for (const DnaPlacementCandidate& candidate : catalogue.candidates) {
    if (candidate.id == winner_candidate)
      continue;
    bool owns_selected_block = false;
    for (const auto& block : catalogue.partition.selected.blocks) {
      if (block.candidate == candidate.id) {
        owns_selected_block = true;
        break;
      }
    }
    if (owns_selected_block)
      continue;
    ids.push_back(candidate.id);
  }
  // The dp2 owner always competes, even when it owns a block: the candidate
  // of the record dp_max2 came from.
  if (dp2_owner != ::fa::cpu::voting::kNullCandidate &&
      dp2_owner != winner_candidate &&
      std::find(ids.begin(), ids.end(), dp2_owner) == ids.end())
    ids.push_back(dp2_owner);

  evidence.rivals.reserve(ids.size());
  for (const ::fa::cpu::voting::CandidateId id : ids) {
    DnaChainMapqRival rival = chain_mapq_rival_from_placement(
        catalogue, placement, id, id == dp2_owner);
    rival.owns_dp2 = id == dp2_owner;
    evidence.rivals.push_back(rival);
  }
  // Strongest vote first, then catalogue id; the formula does not depend on
  // the order.
  std::stable_sort(evidence.rivals.begin(), evidence.rivals.end(),
                   [](const DnaChainMapqRival& a, const DnaChainMapqRival& b) {
                     return std::tie(b.vote_evidence, a.candidate) <
                            std::tie(a.vote_evidence, b.candidate);
                   });
  return evidence;
}

// The realized secondaries of other candidates join the rivals, each on its
// own region's chain, and the dp2 seat moves to the one dp_max2 came from. A
// path -1 secondary is its candidate's whole-query chain: it chains that
// candidate's entry when the entry has no chain and is not added twice; the
// segments of one region are one rival. The record's own candidate's
// secondaries stay in sib_f2 and n_sub.
void admit_realized_rivals(const DnaPlacementFamily& catalogue,
                           const std::vector<DnaRealizedRival>& realized,
                           const DnaRegionPrice& price,
                           DnaChainMapqEvidence& evidence) {
  for (DnaChainMapqRival& rival : evidence.rivals) rival.owns_dp2 = false;
  // The entry each region took, and each candidate's first entry.
  std::unordered_map<int, std::size_t> taken;
  std::unordered_map<int, std::size_t> first;
  for (std::size_t at = 0; at < evidence.rivals.size(); ++at)
    first.emplace(evidence.rivals[at].candidate, at);
  for (std::size_t index = 0; index < realized.size(); ++index) {
    const DnaRealizedRival& source = realized[index];
    if (source.candidate == price.candidate) continue;
    std::size_t slot = evidence.rivals.size();
    if (const auto region = taken.find(source.region); region != taken.end())
      slot = region->second;
    else if (source.path < 0)
      if (const auto entry = first.find(source.candidate); entry != first.end())
        slot = entry->second;
    if (slot == evidence.rivals.size()) {
      DnaChainMapqRival rival;
      rival.candidate = source.candidate;
      if (const DnaPlacementCandidate* candidate = catalogue.find(
              static_cast<::fa::cpu::voting::CandidateId>(source.candidate))) {
        rival.vote_evidence = candidate->vote_evidence;
        rival.contig = candidate->peak.chr;
        rival.reverse = candidate->peak.is_rc;
        rival.peak_diagonal = candidate->peak.raw_ref_start;
      }
      first.emplace(rival.candidate, slot);
      evidence.rivals.push_back(rival);
    }
    taken.emplace(source.region, slot);
    DnaChainMapqRival& entry = evidence.rivals[slot];
    if (!entry.chained) {
      entry.status = static_cast<int>(DnaPlacementChainStatus::Accepted);
      entry.chained = true;
      entry.chain_score = source.score0;
      entry.chain_anchors = source.anchors0;
      entry.q_begin = source.chain_q_begin;
      entry.q_end = source.chain_q_end;
      entry.ref_begin = source.chain_ref_begin;
      entry.ref_end = source.chain_ref_end;
    }
    if (static_cast<int>(index) == price.dp_max2_rival) entry.owns_dp2 = true;
  }
  std::stable_sort(evidence.rivals.begin(), evidence.rivals.end(),
                   [](const DnaChainMapqRival& a, const DnaChainMapqRival& b) {
                     return std::tie(b.vote_evidence, a.candidate) <
                            std::tie(a.vote_evidence, b.candidate);
                   });
}

// A record printed from family.block_parts[part] is priced on that part
// (dna.h part_priced); its rivals stay.
void price_block_part(const DnaPlacementFamily& family, int part,
                      DnaChainMapqEvidence& evidence) {
  if (part < 0 || static_cast<std::size_t>(part) >= family.block_parts.size() ||
      evidence.f1 <= 0)
    return;
  const DnaBlockPart& block_part =
      family.block_parts[static_cast<std::size_t>(part)];
  evidence.f1 = block_part.item_score;
  evidence.cnt = block_part.item_anchors;
  // As minimap2 rounds a split piece's score.
  evidence.sib_f2 = static_cast<int>(block_part.pool_subsc + .499);
  evidence.part_priced = true;
  evidence.part_score = block_part.score;
  evidence.part_anchors = block_part.anchors;
  evidence.part_n_sub = block_part.pool_n_sub;
}

// A record's gap-compressed divergence and its event denominator, as the PAF
// de:f tag computes them. False without base-level accounting.
bool record_event_divergence(const AlignResult& record, double& divergence,
                             double& denominator) {
  if (!record.alignment_accounting_valid || record.edit_distance < 0)
    return false;
  int gap_bases = 0;
  int gap_opens = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    if (op.second == 'I' || op.second == 'D') {
      gap_bases += op.first;
      ++gap_opens;
    }
  }
  const int events = record.block_len + std::max(0, record.ambiguities) -
                     gap_bases + gap_opens;
  if (events <= 0)
    return false;
  denominator = static_cast<double>(events);
  divergence = 1.0 - static_cast<double>(record.matches) / denominator;
  return true;
}

// minimap2's emission floor (mm_filter_regs): no record with matched bases
// under min_chain_score or DP score under min_dp_max (-s). min_dp_max <= 0
// disables the DP half.
bool clears_dna_emission_floor(const AlignResult& record, int min_chain_score,
                               int min_dp_max) {
  if (record.matches < min_chain_score)
    return false;
  return min_dp_max <= 0 || record.score >= min_dp_max;
}

// Collects each record's divergence for the divergence contrast (index 0 the
// primary, i + 1 the i-th supplementary) and runs it. Inversion middles and
// records below the emission floor take no part.
DnaFamilyDivergence divergence_contrast(const dna::Result& realized,
                                        int min_chain_score, int min_dp_max) {
  DnaFamilyDivergence family;
  family.records.assign(1 + realized.supplementary.size(),
                        DnaRecordDivergence());
  std::vector<double> denominators(family.records.size(), 0.0);
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    const AlignResult& record = index == 0
                                    ? static_cast<const AlignResult&>(realized)
                                    : realized.supplementary[index - 1];
    if (record.origin == AlignmentOrigin::DnaLocalInversion)
      continue;
    if (!clears_dna_emission_floor(record, min_chain_score, min_dp_max))
      continue;
    if (!record_event_divergence(record, family.records[index].divergence,
                                 denominators[index]))
      continue;
    family.records[index].accounted = true;
  }
  dna_divergence_contrast(family, denominators);
  return family;
}

// A record's divergence structure: mismatches and ambiguous bases from its
// replay accounting, indel events and bases and aligned columns from its
// CIGAR. All -1 without a CIGAR; the first two -1 without accounting.
struct DnaChainMapqRecordCensus {
  int mismatches = -1;
  int ambiguous = -1;
  int ins_events = -1;
  int ins_bases = -1;
  int del_events = -1;
  int del_bases = -1;
  int aligned_bases = -1;
};

DnaChainMapqRecordCensus chain_mapq_record_census(const AlignResult& record) {
  DnaChainMapqRecordCensus census;
  if (record.cigar.empty())
    return census;
  census.ins_events = 0;
  census.ins_bases = 0;
  census.del_events = 0;
  census.del_bases = 0;
  census.aligned_bases = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    if (op.first <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X':
      census.aligned_bases += op.first;
      break;
    case 'I':
      ++census.ins_events;
      census.ins_bases += op.first;
      break;
    case 'D':
      ++census.del_events;
      census.del_bases += op.first;
      break;
    default:
      break;
    }
  }
  if (record.alignment_accounting_valid) {
    census.mismatches = record.mismatches;
    census.ambiguous = record.ambiguities;
  }
  return census;
}

// The contig index of a name, or -1 when unknown.
int chain_mapq_contig_index(const DnaContext& dctx, const std::string& name) {
  if (dctx.ref.names == nullptr)
    return -1;
  const std::vector<std::string>& names = *dctx.ref.names;
  for (std::size_t index = 0; index < names.size(); ++index)
    if (names[index] == name)
      return static_cast<int>(index);
  return -1;
}

// A committed record's alignment as a map from oriented query position to
// diagonal, one entry per M/=/X run of its CIGAR.
struct DnaChainMapqRecordDiagonals {
  struct Run {
    int q_begin = 0; // oriented query, half-open
    int q_end = 0;
    std::int64_t diagonal = 0;
  };
  std::vector<Run> runs;
};

DnaChainMapqRecordDiagonals
chain_mapq_record_diagonals(const AlignResult& record) {
  DnaChainMapqRecordDiagonals out;
  if (record.cigar.empty() || record.pos < 0)
    return out;
  int qoff = 0;
  std::int64_t toff = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    const int len = op.first;
    if (len <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X':
      out.runs.push_back({qoff, qoff + len,
                          static_cast<std::int64_t>(record.pos) + toff - qoff});
      qoff += len;
      toff += len;
      break;
    case 'I':
    case 'S':
      qoff += len;
      break;
    case 'D':
    case 'N':
      toff += len;
      break;
    default: // H, P
      break;
    }
  }
  return out;
}

// The record's diagonal at oriented query position `oq`: the run holding it,
// else the nearest run boundary within `slack` (an insertion run between two
// aligned runs, or a clipped end). False when the record does not reach it.
bool chain_mapq_record_diagonal_at(const DnaChainMapqRecordDiagonals& map,
                                   std::int64_t oq, std::int64_t slack,
                                   std::int64_t& diagonal) {
  std::int64_t best_distance = -1;
  for (const DnaChainMapqRecordDiagonals::Run& run : map.runs) {
    std::int64_t distance = 0;
    if (oq < run.q_begin)
      distance = run.q_begin - oq;
    else if (oq >= run.q_end)
      distance = oq - (run.q_end - 1);
    if (best_distance < 0 || distance < best_distance) {
      best_distance = distance;
      diagonal = run.diagonal;
    }
    if (distance == 0)
      break;
  }
  return best_distance >= 0 && best_distance <= slack;
}

// A whole-query chain's contig, strand, forward-query span and reference
// span. Without a chain or either span it is invalid: it explains nothing and
// nothing explains it.
struct DnaChainMapqChainGeometry {
  int contig = -1;
  bool reverse = false;
  int q_begin = -1;
  int q_end = -1;
  int ref_begin = -1;
  int ref_end = -1;

  bool valid() const {
    return contig >= 0 && q_begin >= 0 && q_end > q_begin && ref_begin >= 0 &&
           ref_end > ref_begin;
  }
};

DnaChainMapqChainGeometry
chain_mapq_rival_geometry(const DnaChainMapqRival& rival) {
  DnaChainMapqChainGeometry chain;
  if (!rival.chained)
    return chain;
  chain.contig = rival.contig;
  chain.reverse = rival.reverse;
  chain.q_begin = rival.q_begin;
  chain.q_end = rival.q_end;
  chain.ref_begin = rival.ref_begin;
  chain.ref_end = rival.ref_end;
  return chain;
}

// Whether one committed record explains a chain: both chain ends lie on the
// record's diagonal there, within `slack`, on the same contig and strand.
bool chain_mapq_record_explains(const DnaChainMapqRecordDiagonals& map,
                                int record_contig, bool record_reverse,
                                const DnaChainMapqChainGeometry& chain,
                                int read_len, std::int64_t slack) {
  if (!chain.valid() || chain.contig != record_contig ||
      chain.reverse != record_reverse)
    return false;
  const std::int64_t oriented_begin =
      chain.reverse ? read_len - chain.q_end : chain.q_begin;
  const std::int64_t oriented_end =
      chain.reverse ? read_len - chain.q_begin : chain.q_end;
  const std::int64_t start =
      static_cast<std::int64_t>(chain.ref_begin) - oriented_begin;
  const std::int64_t end =
      static_cast<std::int64_t>(chain.ref_end) - oriented_end;
  std::int64_t at_start = 0, at_end = 0;
  if (!chain_mapq_record_diagonal_at(map, oriented_begin, slack, at_start) ||
      !chain_mapq_record_diagonal_at(map, oriented_end - 1, slack, at_end))
    return false;
  return std::llabs(at_start - start) <= slack &&
         std::llabs(at_end - end) <= slack;
}

// The committed family's diagonal maps, built once per read (index 0 the
// primary, i + 1 the i-th supplementary). A record without CIGAR, position,
// known contig or aligned run explains nothing.
struct DnaChainMapqFamilyMap {
  struct Record {
    DnaChainMapqRecordDiagonals diagonals;
    int contig = -1;
    bool reverse = false;
    int q_begin = -1; // forward query, half-open
    int q_end = -1;
    bool inversion = false;
    // The candidate this record came from, or null for an inversion middle.
    // A record does not vouch for its own owner.
    ::fa::cpu::voting::CandidateId owner = ::fa::cpu::voting::kNullCandidate;
  };
  std::vector<Record> records;

  // Bit 8: the block's own record explains the chain.
  bool explains_by_record(std::size_t index,
                          const DnaChainMapqChainGeometry& chain, int read_len,
                          std::int64_t slack) const {
    if (index >= records.size())
      return false;
    const Record& record = records[index];
    return chain_mapq_record_explains(record.diagonals, record.contig,
                                      record.reverse, chain, read_len, slack);
  }

  // Bit 64: the chain's start lies on some record's diagonal and its end on
  // some record's, the same or another, as for a chain spanning a deletion
  // between two flank records. Inversion middles and the record `candidate`
  // owns (which trivially explains its own chain) take no part.
  bool explains_by_family(const DnaChainMapqChainGeometry& chain, int read_len,
                          std::int64_t slack,
                          ::fa::cpu::voting::CandidateId candidate) const {
    if (!chain.valid())
      return false;
    const std::int64_t oriented_begin =
        chain.reverse ? read_len - chain.q_end : chain.q_begin;
    const std::int64_t oriented_end =
        chain.reverse ? read_len - chain.q_begin : chain.q_end;
    const std::int64_t start =
        static_cast<std::int64_t>(chain.ref_begin) - oriented_begin;
    const std::int64_t end =
        static_cast<std::int64_t>(chain.ref_end) - oriented_end;
    bool start_held = false;
    bool end_held = false;
    for (const Record& record : records) {
      if (record.inversion || record.contig != chain.contig ||
          record.reverse != chain.reverse)
        continue;
      if (candidate != ::fa::cpu::voting::kNullCandidate &&
          record.owner == candidate)
        continue;
      std::int64_t at = 0;
      if (!start_held &&
          chain_mapq_record_diagonal_at(record.diagonals, oriented_begin, slack,
                                        at) &&
          std::llabs(at - start) <= slack)
        start_held = true;
      if (!end_held &&
          chain_mapq_record_diagonal_at(record.diagonals, oriented_end - 1,
                                        slack, at) &&
          std::llabs(at - end) <= slack)
        end_held = true;
      if (start_held && end_held)
        return true;
    }
    return false;
  }
};

DnaChainMapqFamilyMap chain_mapq_family_map(
    const DnaContext& dctx, const dna::Result& realized,
    ::fa::cpu::voting::CandidateId primary_candidate) {
  DnaChainMapqFamilyMap map;
  map.records.resize(1 + realized.supplementary.size());
  for (std::size_t index = 0; index < map.records.size(); ++index) {
    const AlignResult& member =
        index == 0 ? static_cast<const AlignResult&>(realized)
                   : realized.supplementary[index - 1];
    map.records[index].owner =
        index == 0 ? primary_candidate
        : index - 1 < realized.supplementary_candidates.size()
            ? realized.supplementary_candidates[index - 1]
            : ::fa::cpu::voting::kNullCandidate;
    if (member.cigar.empty() || member.pos < 0)
      continue;
    const int contig = chain_mapq_contig_index(dctx, member.chromosome);
    if (contig < 0)
      continue;
    DnaChainMapqFamilyMap::Record& record = map.records[index];
    record.diagonals = chain_mapq_record_diagonals(member);
    if (record.diagonals.runs.empty())
      continue;
    record.contig = contig;
    record.reverse = member.is_reverse;
    record.q_begin = member.query_start;
    record.q_end = member.query_end;
    record.inversion = member.origin == AlignmentOrigin::DnaLocalInversion;
  }
  return map;
}

// Bits 8 / 64: marks every chained rival the committed record (or family)
// explains as a record shadow, and adds the block owners it explains as
// rivals, so their vote is credited to the winner. `record_index` is the
// block's own record in the family map.
void mark_chain_mapq_record_shadows(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement,
    const DnaChainMapqFamilyMap& family_map, std::size_t record_index,
    ::fa::cpu::voting::CandidateId winner_candidate, int read_len,
    DnaChainMapqEvidence& evidence) {
  const int bits = chain_mapq_study_bits(dctx);
  if ((bits & (kDnaChainMapqStudyRecordShadow |
               kDnaChainMapqStudyFamilyVouch)) == 0)
    return;
  // The block's own record is asked first, then the family under bit 64.
  const bool family_arm = (bits & kDnaChainMapqStudyFamilyVouch) != 0;
  const std::int64_t slack = dna_chain_mapq_chain_shadow_slack(evidence);
  const auto vouch_for = [&](const DnaChainMapqRival& rival) {
    const DnaChainMapqChainGeometry chain = chain_mapq_rival_geometry(rival);
    if (family_map.explains_by_record(record_index, chain, read_len, slack))
      return 1;
    if (family_arm && family_map.explains_by_family(chain, read_len, slack,
                                                    rival.candidate))
      return 2;
    return 0;
  };
  for (DnaChainMapqRival& rival : evidence.rivals) {
    const int vouch = vouch_for(rival);
    rival.vouch = vouch;
    rival.record_shadow = vouch != 0;
  }
  bool appended = false;
  for (const DnaPlacementCandidate& candidate : catalogue.candidates) {
    if (candidate.id == winner_candidate)
      continue;
    bool owns_selected_block = false;
    for (const auto& block : catalogue.partition.selected.blocks)
      if (block.candidate == candidate.id) {
        owns_selected_block = true;
        break;
      }
    if (!owns_selected_block)
      continue;
    bool present = false;
    for (const DnaChainMapqRival& rival : evidence.rivals)
      if (rival.candidate == candidate.id) {
        present = true;
        break;
      }
    if (present)
      continue;
    DnaChainMapqRival rival =
        chain_mapq_rival_from_placement(catalogue, placement, candidate.id);
    const int vouch = vouch_for(rival);
    if (vouch == 0)
      continue;
    rival.record_shadow = true;
    rival.vouch = vouch;
    evidence.rivals.push_back(rival);
    appended = true;
  }
  if (appended)
    std::stable_sort(
        evidence.rivals.begin(), evidence.rivals.end(),
        [](const DnaChainMapqRival& a, const DnaChainMapqRival& b) {
          return std::tie(b.vote_evidence, a.candidate) <
                 std::tie(a.vote_evidence, b.candidate);
        });
}

} // namespace

} // namespace dna
} // namespace fa::cpu::lr

/**** pricing ****/

namespace fa::cpu::lr {
namespace dna {

// The MAPQ of the records realization.cpp emits, through the router:
// each primary record on its own chain with dp1 its dp_max and dp2 its
// dp_max2 (owned by the record that gave it, when that is another catalogue
// candidate's; every realized secondary of another candidate is a rival,
// admit_realized_rivals). The region stands in for the block part
// (price_block_part): f1 and cnt its whole chain's, score and cnt its own,
// sib_f2 and n_sub mm_set_parent's over the secondaries from its own
// candidate's chains. An inversion middle takes its flanks' MAPQ.
void price_dna_regions(const DnaContext& dctx,
                       const DnaPlacementFamily& family,
                       const DnaPlacementChainingResult& placement,
                       const std::vector<DnaRegionPrice>& prices,
                       const std::vector<DnaRealizedRival>& rivals,
                       int read_len, dna::Result& realized) {
  const std::size_t count = 1 + realized.supplementary.size();
  if (prices.size() != count)
    return;
  const auto member = [&realized](std::size_t index) -> AlignResult& {
    return index == 0 ? static_cast<AlignResult&>(realized)
                      : realized.supplementary[index - 1];
  };
  const bool family_map_wanted = dctx.opts.chain_mapq_hifi_margin;
  const ::fa::cpu::voting::CandidateId primary_candidate =
      prices[0].candidate >= 0
          ? static_cast<::fa::cpu::voting::CandidateId>(prices[0].candidate)
          : ::fa::cpu::voting::kNullCandidate;
  DnaChainMapqFamilyMap family_map;
  if (family_map_wanted)
    family_map = chain_mapq_family_map(dctx, realized, primary_candidate);
  int family_blocks = 0;
  for (const DnaRegionPrice& price : prices)
    if (!price.inversion)
      ++family_blocks;
  const bool census_wanted =
      dctx.opts.chain_mapq_hifi_margin && family_blocks >= 2;
  std::vector<std::pair<int, int>> family_mismatch_pairs;
  if ((chain_mapq_study_bits(dctx) & kDnaChainMapqStudyRelativeShareCap) !=
          0 &&
      census_wanted) {
    family_mismatch_pairs.assign(count, std::pair<int, int>(-1, -1));
    for (std::size_t index = 0; index < count; ++index) {
      const AlignResult& record = member(index);
      if (record.origin == AlignmentOrigin::DnaLocalInversion ||
          !record.alignment_accounting_valid ||
          !clears_dna_emission_floor(record, dctx.opts.min_chain_score,
                                     dctx.opts.cigar_dp_min_dp_max))
        continue;
      const DnaChainMapqRecordCensus census = chain_mapq_record_census(record);
      if (census.mismatches < 0 || census.aligned_bases <= 0)
        continue;
      family_mismatch_pairs[index] = {census.mismatches, census.aligned_bases};
    }
  }
  std::vector<int> mapqs(count, -1);
  for (std::size_t index = 0; index < count; ++index) {
    const DnaRegionPrice& price = prices[index];
    if (price.inversion)
      continue;
    AlignResult& record = member(index);
    const ::fa::cpu::voting::CandidateId owner =
        price.candidate >= 0
            ? static_cast<::fa::cpu::voting::CandidateId>(price.candidate)
            : ::fa::cpu::voting::kNullCandidate;
    const ::fa::cpu::voting::CandidateId dp2_owner =
        price.dp_max2_candidate >= 0 &&
                price.dp_max2_candidate != price.candidate
            ? static_cast<::fa::cpu::voting::CandidateId>(
                  price.dp_max2_candidate)
            : ::fa::cpu::voting::kNullCandidate;
    const double identity =
        record.alignment_accounting_valid && record.block_len > 0
            ? record.identity()
            : 1.0;
    const double dp1 = static_cast<double>(price.dp_max);
    const double dp2 = static_cast<double>(price.dp_max2);
    DnaChainMapqEvidence evidence = build_chain_mapq_evidence(
        dctx, family, placement, true,
        owner != ::fa::cpu::voting::kNullCandidate ? placement.find(owner)
                                                   : nullptr,
        owner, dp2_owner, dp1, dp2, identity, read_len);
    evidence.winner_q_begin = price.chain_q_begin;
    evidence.winner_q_end = price.chain_q_end;
    evidence.winner_ref_begin = price.chain_ref_begin;
    evidence.winner_ref_end = price.chain_ref_end;
    evidence.dp1_raw = static_cast<double>(price.dp_max0);
    evidence.dp2_raw = static_cast<double>(price.dp_max2_0);
    admit_realized_rivals(family, rivals, price, evidence);
    if (family_map_wanted)
      mark_chain_mapq_record_shadows(dctx, family, placement, family_map, index,
                                     owner, read_len, evidence);
    evidence.f1 = price.score0;
    evidence.cnt = price.anchors0;
    evidence.part_priced = true;
    evidence.part_score = price.score;
    evidence.part_anchors = price.cnt;
    evidence.sib_f2 = static_cast<int>(price.pool_subsc + .499);
    evidence.part_n_sub = price.pool_n_sub;
    const DnaChainMapqRecordCensus census =
        census_wanted ? chain_mapq_record_census(record)
                      : DnaChainMapqRecordCensus();
    evidence.family_blocks = family_blocks;
    evidence.record_mismatches = census.mismatches;
    evidence.record_ins_events = census.ins_events;
    evidence.record_del_events = census.del_events;
    evidence.record_aligned_bases = census.aligned_bases;
    for (std::size_t other = 0; other < family_mismatch_pairs.size(); ++other) {
      if (other == index || family_mismatch_pairs[other].second <= 0)
        continue;
      const bool cleaner =
          evidence.family_min_aligned <= 0 ||
          static_cast<double>(family_mismatch_pairs[other].first) *
                  static_cast<double>(evidence.family_min_aligned) <
              static_cast<double>(evidence.family_min_mismatches) *
                  static_cast<double>(family_mismatch_pairs[other].second);
      if (cleaner) {
        evidence.family_min_mismatches = family_mismatch_pairs[other].first;
        evidence.family_min_aligned = family_mismatch_pairs[other].second;
      }
    }
    DnaChainMapqBreakdown breakdown;
    mapqs[index] = dna_chain_mapq(evidence, &breakdown, nullptr);
    if (evidence.f1 > 0) {
      record.chain_anchors = price.cnt;
      record.chain_score = price.score;
      record.secondary_chain_score = breakdown.f2;
    }
  }
  const DnaFamilyDivergence divergence = divergence_contrast(
      realized, dctx.opts.min_chain_score, dctx.opts.cigar_dp_min_dp_max);
  for (std::size_t index = 0; index < count; ++index) {
    if (!divergence.records[index].divergent)
      continue;
    const int carried = mapqs[index] >= 0 ? mapqs[index] : mapqs[0];
    mapqs[index] = std::min(carried, kDnaChainMapqDivergentRecordMapq);
  }
  const int mapq = std::max(0, mapqs[0]);
  mapqs.erase(mapqs.begin());
  route_dna_mapq(family, primary_candidate, mapq, mapqs, realized);
}

// The MAPQ of the records map-only emits, through the router: the primary and
// every supplementary that owns a selected block, each on its own candidate's
// whole-query chain and priced on the block part it prints from
// (price_block_part), then the divergence contrast cap. A record that owns no
// block inherits the primary's MAPQ.
void price_dna_blocks(const DnaContext& dctx,
                      const DnaPlacementFamily& placement_family,
                      const DnaPlacementChainingResult& placement_chaining,
                      bool placement_chaining_ran,
                      ::fa::cpu::voting::CandidateId primary_candidate,
                      int read_len, dna::Result& realized) {
  // Confidence belongs to the committed hypothesis, so MAPQ cannot affect
  // placement or record-family decisions.
  const DnaPlacementCandidateChain* winner_chain =
      placement_chaining_ran ? placement_chaining.find(primary_candidate)
                             : nullptr;
  // Confidence is minimap2's formula over the whole-query chains, scored
  // once per block-owning record (dna.h), after everything else is
  // committed.
  const double identity =
      realized.alignment_accounting_valid && realized.block_len > 0
          ? realized.identity()
          : 1.0;
  // The committed family's diagonal map, built once for the shadow rules.
  // HiFi preset only.
  const bool family_map_wanted = dctx.opts.chain_mapq_hifi_margin;
  DnaChainMapqFamilyMap family_map;
  if (family_map_wanted)
    family_map = chain_mapq_family_map(dctx, realized, primary_candidate);
  // The committed family's block-owning records. The split-read caps
  // (dna.h R6) read this and each block record's divergence, so
  // the per-record counts are taken only where the caps can fire.
  int family_blocks = 1;
  for (const ::fa::cpu::voting::CandidateId owner :
       realized.supplementary_candidates)
    if (owner != ::fa::cpu::voting::kNullCandidate)
      ++family_blocks;
  const bool census_wanted =
      dctx.opts.chain_mapq_hifi_margin && family_blocks >= 2;
  const auto block_census = [census_wanted](const AlignResult& record) {
    return census_wanted ? chain_mapq_record_census(record)
                         : DnaChainMapqRecordCensus();
  };
  // For kDnaChainMapqStudyRelativeShareCap: each record's (mismatches,
  // aligned bases), over the records the divergence contrast collects, so
  // each block's evidence can carry the family's cleanest other record.
  // Index 0 is the primary.
  std::vector<std::pair<int, int>> family_mismatch_pairs;
  if ((chain_mapq_study_bits(dctx) & kDnaChainMapqStudyRelativeShareCap) !=
          0 &&
      census_wanted) {
    family_mismatch_pairs.assign(1 + realized.supplementary.size(),
                                 std::pair<int, int>(-1, -1));
    for (std::size_t index = 0; index < family_mismatch_pairs.size();
         ++index) {
      const AlignResult& member =
          index == 0 ? static_cast<const AlignResult&>(realized)
                     : realized.supplementary[index - 1];
      if (member.origin == AlignmentOrigin::DnaLocalInversion ||
          !member.alignment_accounting_valid ||
          !clears_dna_emission_floor(member, dctx.opts.min_chain_score,
                                     dctx.opts.cigar_dp_min_dp_max))
        continue;
      const DnaChainMapqRecordCensus member_census =
          chain_mapq_record_census(member);
      if (member_census.mismatches < 0 || member_census.aligned_bases <= 0)
        continue;
      family_mismatch_pairs[index] = {member_census.mismatches,
                                      member_census.aligned_bases};
    }
  }
  const auto rule_fields = [family_blocks, &family_mismatch_pairs](
                               const DnaChainMapqRecordCensus& census,
                               DnaChainMapqEvidence& into,
                               std::size_t record_index) {
    into.family_blocks = family_blocks;
    into.record_mismatches = census.mismatches;
    into.record_ins_events = census.ins_events;
    into.record_del_events = census.del_events;
    into.record_aligned_bases = census.aligned_bases;
    // The cleanest other accounted record by mismatch rate.
    for (std::size_t other = 0; other < family_mismatch_pairs.size();
         ++other) {
      if (other == record_index || family_mismatch_pairs[other].second <= 0)
        continue;
      const bool cleaner =
          into.family_min_aligned <= 0 ||
          static_cast<double>(family_mismatch_pairs[other].first) *
                  static_cast<double>(into.family_min_aligned) <
              static_cast<double>(into.family_min_mismatches) *
                  static_cast<double>(family_mismatch_pairs[other].second);
      if (cleaner) {
        into.family_min_mismatches = family_mismatch_pairs[other].first;
        into.family_min_aligned = family_mismatch_pairs[other].second;
      }
    }
  };
  DnaChainMapqEvidence evidence = build_chain_mapq_evidence(
      dctx, placement_family, placement_chaining, placement_chaining_ran,
      winner_chain, primary_candidate, ::fa::cpu::voting::kNullCandidate,
      0.0, 0.0, identity, read_len);
  mark_chain_mapq_record_shadows(dctx, placement_family,
                                 placement_chaining, family_map, 0,
                                 primary_candidate, read_len, evidence);
  price_block_part(placement_family, realized.primary_part, evidence);
  const DnaChainMapqRecordCensus primary_census =
      block_census(static_cast<const AlignResult&>(realized));
  rule_fields(primary_census, evidence, 0);
  DnaChainMapqBreakdown breakdown;
  std::vector<DnaChainMapqRivalVerdict> verdicts;
  int mapq = dna_chain_mapq(evidence, &breakdown, &verdicts);
  // cm:i and s1:i: the block part the record prints from; s2:i: the
  // competing chain this record's MAPQ weighed. With f1 <= 0 there is no
  // chain, the breakdown is default-constructed and "no evidence" is not
  // "no rival", so the record keeps -1 and writes none of the three.
  if (evidence.f1 > 0) {
    const DnaBlockPart& part = placement_family.block_parts
        [static_cast<std::size_t>(realized.primary_part)];
    realized.chain_anchors = part.anchors;
    realized.chain_score = part.score;
    realized.secondary_chain_score = breakdown.f2;
  }
  // Every supplementary that owns a selected block is scored on its own
  // candidate's whole-query chain, against the same rivals. One that owns
  // no block keeps -1 and inherits the primary's MAPQ during routing.
  std::vector<int> supplementary_mapq(realized.supplementary.size(), -1);
  for (std::size_t index = 0; index < realized.supplementary.size();
       ++index) {
    const ::fa::cpu::voting::CandidateId owner =
        index < realized.supplementary_candidates.size()
            ? realized.supplementary_candidates[index]
            : ::fa::cpu::voting::kNullCandidate;
    if (owner == ::fa::cpu::voting::kNullCandidate ||
        !placement_chaining_ran)
      continue;
    const DnaPlacementCandidateChain* block_chain =
        placement_chaining.find(owner);
    if (block_chain == nullptr)
      continue;
    const AlignResult& record = realized.supplementary[index];
    const double block_identity =
        record.alignment_accounting_valid && record.block_len > 0
            ? record.identity()
            : 1.0;
    DnaChainMapqEvidence block_evidence = build_chain_mapq_evidence(
        dctx, placement_family, placement_chaining, placement_chaining_ran,
        block_chain, owner, ::fa::cpu::voting::kNullCandidate, 0.0, 0.0,
        block_identity, read_len);
    mark_chain_mapq_record_shadows(dctx, placement_family,
                                   placement_chaining, family_map,
                                   index + 1, owner, read_len,
                                   block_evidence);
    price_block_part(placement_family,
                     index < realized.supplementary_parts.size()
                         ? realized.supplementary_parts[index]
                         : -1,
                     block_evidence);
    const DnaChainMapqRecordCensus census = block_census(record);
    rule_fields(census, block_evidence, index + 1);
    DnaChainMapqBreakdown block_breakdown;
    std::vector<DnaChainMapqRivalVerdict> block_verdicts;
    supplementary_mapq[index] =
        dna_chain_mapq(block_evidence, &block_breakdown, &block_verdicts);
    // cm:i, s1:i and s2:i for this block's record, as for the primary.
    // Written through the vector because `record` is a const reference.
    if (block_evidence.f1 > 0) {
      AlignResult& stamped = realized.supplementary[index];
      const DnaBlockPart& part = placement_family.block_parts
          [static_cast<std::size_t>(realized.supplementary_parts[index])];
      stamped.chain_anchors = part.anchors;
      stamped.chain_score = part.score;
      stamped.secondary_chain_score = block_breakdown.f2;
    }
  }
  // Divergence contrast cap (dna.h): a record whose event
  // divergence sits far above the family's minimum is the wrong copy,
  // whatever its chain says. The cap applies to the MAPQ each record would
  // otherwise carry, inherited ones included.
  const DnaFamilyDivergence divergence =
      divergence_contrast(realized, dctx.opts.min_chain_score,
                          dctx.opts.cigar_dp_min_dp_max);
  if (divergence.records[0].divergent)
    mapq = std::min(mapq, kDnaChainMapqDivergentRecordMapq);
  for (std::size_t index = 0; index < supplementary_mapq.size(); ++index) {
    if (!divergence.records[index + 1].divergent)
      continue;
    const int carried =
        supplementary_mapq[index] >= 0 ? supplementary_mapq[index] : mapq;
    supplementary_mapq[index] =
        std::min(carried, kDnaChainMapqDivergentRecordMapq);
  }
  route_dna_mapq(placement_family, primary_candidate, mapq,
                 supplementary_mapq, realized);
}

} // namespace dna
} // namespace fa::cpu::lr
