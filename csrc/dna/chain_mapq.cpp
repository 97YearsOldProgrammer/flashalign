#include "chain_mapq.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace fa::cpu::lr {

namespace {

// C's (int) conversion, truncating toward zero, saturated so an out-of-range
// value is not undefined behaviour.
int truncate_to_int(double value) noexcept {
  if (!(value == value))
    return 0; // NaN
  if (value >= 2147483647.0)
    return 2147483647;
  if (value <= -2147483648.0)
    return -2147483648;
  return static_cast<int>(value);
}

// A rival scoring under min_chain_score counts as min_chain_score, as
// minimap2's subsc.
double chain_ratio(const DnaChainMapqEvidence& evidence,
                   const DnaChainMapqRival& rival) noexcept {
  return static_cast<double>(std::max(rival.chain_score, kDnaMinChainScore)) /
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

// The realized rival holding the dp2 seat (R1), or none. The guard and the
// ladder read m_eff; the ratio form reads the seat's own score pair
// (pair_dp1 / pair_dp2): dp1 / dp2 for the alternative, the raw pair for the
// sibling, and for the block rival the two inside scores of the contested
// interval, or the raw pair when there is none or the block's own is not
// positive.
struct DnaChainMapqSeat {
  int kind = 0; // 0 none, 1 alternative, 2 sibling, 3 block rival
  bool has_slot = false;
  std::size_t slot = 0; // index into evidence.rivals when has_slot
  double dp2_raw = 0.0;
  double raw_margin = 0.0;
  bool in_valid = false;
  double in_margin = 0.0;
  double m_eff = 0.0;
  bool chained = false;
  int chain_score = 0;
  int chain_anchors = 0;
  bool tie = false;
  double pair_dp1 = 0.0;
  double pair_dp2 = 0.0;
  double strength = 0.0; // before the DP factor; 1.0 on a tie
};

// R1's margins for one realized rival: raw, interval-matched where the spans
// intersect, and effective: min(raw, in) for a whole-read hypothesis, the
// interval margin alone for a block rival (`interval_only`), raw when there is
// no interval.
void seat_margins(double dp1_raw, double dp2_raw, bool in_valid,
                  double in_margin, bool interval_only,
                  DnaChainMapqSeat& seat) noexcept {
  seat.dp2_raw = dp2_raw;
  seat.raw_margin = dp1_raw - dp2_raw;
  seat.in_valid = in_valid;
  seat.in_margin = in_valid ? in_margin : 0.0;
  seat.m_eff = !in_valid       ? seat.raw_margin
               : interval_only ? in_margin
                               : std::min(seat.raw_margin, in_margin);
}

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

bool dna_chain_mapq_spans_compete(int a_begin, int a_end, int b_begin,
                                  int b_end) noexcept {
  if (a_begin < 0 || a_end <= a_begin || b_begin < 0 || b_end <= b_begin)
    return false;
  const int overlap = std::min(a_end, b_end) - std::max(a_begin, b_begin);
  if (overlap <= 0)
    return false;
  const int shorter = std::min(a_end - a_begin, b_end - b_begin);
  return static_cast<std::int64_t>(overlap) * 2 >
         static_cast<std::int64_t>(shorter);
}

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

  // R1: the seat goes to the realized rival with the smallest effective
  // margin; ties go to the sibling, then the alternative, then the block
  // rival.
  const double f1 = static_cast<double>(evidence.f1);
  const bool sibling_realized = hifi && dp1 > 0.0 &&
                                evidence.sibling.realized &&
                                evidence.sibling.dp2_raw > 0.0;
  const bool block_rival_realized = hifi && dp1 > 0.0 &&
                                    evidence.block_rival.realized &&
                                    evidence.block_rival.dp2_raw > 0.0;
  DnaChainMapqSeat seat;
  const auto consider = [&seat](const DnaChainMapqSeat& candidate) {
    if (seat.kind == 0 || candidate.m_eff < seat.m_eff)
      seat = candidate;
  };
  const auto realized_seat = [&](const DnaChainMapqRealizedRival& rival,
                                 int kind) {
    DnaChainMapqSeat candidate;
    candidate.kind = kind;
    const bool block_rival = kind == 3;
    const bool in_valid = rival.im_lo >= 0 && rival.im_hi > rival.im_lo;
    // A block rival is realized clipped to the block's query span, so it is
    // compared on the contested interval alone.
    seat_margins(evidence.dp1_raw, rival.dp2_raw, in_valid,
                 static_cast<double>(rival.im_a_inside - rival.im_b_inside),
                 /*interval_only=*/block_rival, candidate);
    candidate.chained = true;
    candidate.chain_score = rival.chain_score;
    candidate.chain_anchors = rival.chain_anchors;
    candidate.tie =
        rival.chain_score == evidence.f1 && rival.chain_anchors == evidence.cnt;
    candidate.pair_dp1 = evidence.dp1_raw;
    candidate.pair_dp2 = rival.dp2_raw;
    // R5: the block rival's pair is the two inside scores, while the block's
    // own is positive (the ratio needs a positive denominator).
    if (block_rival && in_valid && rival.im_a_inside > 0) {
      candidate.pair_dp1 = static_cast<double>(rival.im_a_inside);
      candidate.pair_dp2 = static_cast<double>(rival.im_b_inside);
    }
    return candidate;
  };
  if (sibling_realized) {
    // Its strength for the ratio form is the sibling floor,
    // max(sib_f2, min_chain_score) / f1.
    DnaChainMapqSeat candidate = realized_seat(evidence.sibling, 2);
    candidate.strength =
        candidate.tie ? 1.0
                      : static_cast<double>(
                            std::max(evidence.sib_f2, kDnaMinChainScore)) /
                            f1;
    consider(candidate);
  }
  if (alternative_competes) {
    const DnaChainMapqRival& rival = evidence.rivals[alternative_slot];
    const DnaChainMapqRivalVerdict& verdict = judged[alternative_slot];
    DnaChainMapqSeat candidate;
    candidate.kind = 1;
    candidate.has_slot = true;
    candidate.slot = alternative_slot;
    seat_margins(evidence.dp1_raw, evidence.dp2_raw,
                 hifi && evidence.alternative_im_valid,
                 evidence.alternative_im_margin, /*interval_only=*/false,
                 candidate);
    candidate.chained = rival.chained;
    candidate.chain_score = rival.chain_score;
    candidate.chain_anchors = rival.chain_anchors;
    candidate.tie = verdict.tie;
    candidate.pair_dp1 = dp1;
    candidate.pair_dp2 = evidence.dp2;
    candidate.strength = verdict.strength;
    consider(candidate);
  }
  if (block_rival_realized) {
    DnaChainMapqSeat candidate = realized_seat(evidence.block_rival, 3);
    // The strength is the rival's catalogue slot verdict, or, when the slot
    // is inadmissible or a shadow, the realized chain's ratio.
    candidate.strength =
        static_cast<double>(
            std::max(evidence.block_rival.chain_score, kDnaMinChainScore)) /
        f1;
    bool slot_shadow = false;
    for (std::size_t index = 0; index < evidence.rivals.size(); ++index) {
      if (evidence.rivals[index].candidate != evidence.block_rival.candidate)
        continue;
      slot_shadow = judged[index].shadow;
      if (judged[index].admissible && !judged[index].shadow) {
        candidate.has_slot = true;
        candidate.slot = index;
        candidate.strength = judged[index].strength;
      }
      break;
    }
    if (candidate.tie)
      candidate.strength = 1.0;
    // Bit 32: a block rival the shadow rule calls the winner's own locus
    // cannot take the seat.
    const bool shadow_seat =
        slot_shadow &&
        (evidence.study_bits & kDnaChainMapqStudyNoShadowSeat) != 0;
    if (!shadow_seat)
      consider(candidate);
  }
  const bool dp_branch = seat.kind != 0;
  // Outside the HiFi rule dp2 is the alternative's, or nothing.
  const double dp2 = alternative_competes ? evidence.dp2 : 0.0;

  // Second pass: f2, n_sub and x over the competing rivals. As minimap2's
  // mm_set_parent, a rival counts toward n_sub once if it meets the anchor
  // rule or the DP rule; the DP rule applies to the seat, on its own pair.
  int f2 = evidence.sib_f2;
  int n_sub = 0;
  int admissible_rivals = 0;
  int shadow_rivals = 0;
  int tie_rivals = 0;
  double vote_ratio = 0.0;
  double max_competing_strength = 0.0;
  const double winner_vote = static_cast<double>(std::max(1, effective_vote));
  // The winner's sibling path seeds x with its chain ratio, unless it was
  // realized and so is judged as a realized rival instead.
  const double sibling_floor =
      static_cast<double>(std::max(evidence.sib_f2, kDnaMinChainScore)) / f1;
  double x = sibling_realized ? 0.0 : sibling_floor;
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
    const bool is_seat = seat.has_slot && seat.slot == index;
    // R1: a realized rival that lost the seat enters neither x nor x_other.
    const bool lost_seat =
        !is_seat && ((rival.owns_dp2 && alternative_competes) ||
                     (block_rival_realized &&
                      rival.candidate == evidence.block_rival.candidate));
    if (rival.chained) {
      f2 = std::max(f2, rival.chain_score);
      bool counts =
          evidence.n_sub_cnt_clause && rival.chain_anchors >= evidence.cnt;
      if (is_seat && seat.pair_dp1 - seat.pair_dp2 <=
                         static_cast<double>(evidence.sub_diff))
        counts = true;
      verdict.counted = counts;
      if (counts)
        ++n_sub;
    }
    if (is_seat) {
      // The DP ratio discounts the seat once; a tie is not discounted.
      if (!verdict.tie)
        verdict.strength *= seat.pair_dp2 / seat.pair_dp1;
      x = std::max(x, verdict.strength);
    } else if (!lost_seat) {
      x = std::max(x, verdict.strength);
      x_other = std::max(x_other, verdict.strength);
    }
  }
  // A seat without a catalogue slot enters x the same way.
  if (dp_branch && !seat.has_slot)
    x = std::max(x, seat.tie ? 1.0
                             : seat.strength * seat.pair_dp2 / seat.pair_dp1);

  // minimap2's own x (chain ratio alone), for the breakdown.
  const int subsc = std::max(f2, kDnaMinChainScore);
  const double x_chain = static_cast<double>(subsc) / f1;

  // minimap2's pen_s1, without its uniq_ratio factor.
  const double pen_s1 =
      evidence.f1 > 100 ? 1.0 : 0.01 * static_cast<double>(evidence.f1);
  double pen_cm =
      evidence.cnt > 10 ? 1.0 : 0.1 * static_cast<double>(evidence.cnt);
  pen_cm = std::min(pen_s1, pen_cm);
  // R3: the margin rule's pen has no vote floor.
  const double pen_margin = pen_cm;
  const double w_abs = static_cast<double>(std::max(0, effective_vote));
  // The vote floor: a winner the vote barely supports is discounted as
  // minimap2 discounts a short chain.
  pen_cm = std::min(pen_cm, std::min(1.0, w_abs / kDnaChainMapqVoteFloor));

  const double match_sc = static_cast<double>(std::max(1, evidence.match_sc));
  // R2. `margin_beats_rival` also gates the never-Q0 promotion below.
  const bool margin_beats_rival =
      hifi && dp_branch && evidence.dp1_raw > 0.0 && seat.dp2_raw > 0.0 &&
      evidence.substitution_cost > 0 &&
      seat.m_eff >= static_cast<double>(evidence.substitution_cost);
  const bool margin_rule = margin_beats_rival && seat.chained && !seat.tie &&
                           seat.chain_anchors < evidence.cnt &&
                           x_other < kDnaChainMapqHifiRivalStrengthMax;
  int hifi_margin_mapq = -1;
  int mapq;
  if (margin_rule) {
    // R3: BWA-MEM's margin Phred on the raw scores (minimap2's mapq_alt).
    const double two_base_gap =
        static_cast<double>(evidence.gap_open1 + 2 * evidence.gap_extend1);
    int margin_q = truncate_to_int(6.02 * seat.m_eff / match_sc + 0.499);
    if (seat.m_eff > static_cast<double>(evidence.substitution_cost) &&
        seat.m_eff < two_base_gap)
      margin_q = std::min(margin_q, kDnaChainMapqHifiSingleEventMapq);
    hifi_margin_mapq = margin_q;
    mapq = truncate_to_int(pen_margin * static_cast<double>(margin_q));
  } else if (dp_branch) {
    // minimap2's DP branch; the DP ratio is already inside the seat's
    // strength.
    mapq = truncate_to_int(evidence.identity * pen_cm * kDnaChainMapqCoef *
                           (1.0 - x * x) * std::log(dp1 / match_sc));
    // mapq_alt on the seat's pair margin, bounded by m_eff under the HiFi
    // rule (R5).
    const double pair_margin = seat.pair_dp1 - seat.pair_dp2;
    const int mapq_alt = truncate_to_int(
        6.02 * evidence.identity * evidence.identity *
            (hifi ? std::min(pair_margin, seat.m_eff) : pair_margin) /
            match_sc +
        0.499);
    mapq = std::min(mapq, mapq_alt);
  } else if (dp1 > 0.0) {
    // Aligned, but no rival DP score.
    mapq = truncate_to_int(evidence.identity * pen_cm * kDnaChainMapqCoef *
                           (1.0 - x) * std::log(dp1 / match_sc));
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
    breakdown->seat_raw_margin = seat.raw_margin;
    breakdown->seat_in_valid = seat.in_valid;
    breakdown->seat_in_margin = seat.in_margin;
    breakdown->seat_m_eff = seat.m_eff;
    breakdown->hifi_pen = pen_margin;
    breakdown->split_cap = split_cap;
  }
  return mapq;
}

} // namespace fa::cpu::lr
