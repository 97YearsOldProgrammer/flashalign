#include "chain_mapq.h"

#include <algorithm>
#include <cmath>

namespace fa::cpu::lr::rna {

namespace {

// C's (int) truncates toward zero. Saturate first so a negative or huge product stays
// defined.
int truncate_to_int(double value) noexcept {
  if (!(value == value))
    return 0; // NaN
  if (value >= 2147483647.0)
    return 2147483647;
  if (value <= -2147483648.0)
    return -2147483648;
  return static_cast<int>(value);
}

// minimap2's per-rival chain-score floor in subsc, at min_sc.
double chain_ratio(const RnaChainMapqEvidence& evidence,
                   const RnaChainMapqRival& rival) noexcept {
  return static_cast<double>(
             std::max(rival.chain_score, kRnaChainMapqMinChainScore)) /
         static_cast<double>(evidence.f1);
}

} // namespace

bool rna_chain_mapq_admissible(const RnaChainMapqEvidence& evidence,
                               const RnaChainMapqRival& rival) noexcept {
  // A shadow is the committed placement itself.
  if (rival.shadow)
    return false;
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
         kRnaChainMapqMaskLevel * static_cast<double>(shorter);
}

RnaChainMapqRivalVerdict
rna_chain_mapq_rival_verdict(const RnaChainMapqEvidence& evidence,
                             const RnaChainMapqRival& rival) noexcept {
  RnaChainMapqRivalVerdict verdict;
  verdict.admissible = rna_chain_mapq_admissible(evidence, rival);
  if (!verdict.admissible || evidence.f1 <= 0)
    return verdict;
  verdict.tie = rival.chained && rival.chain_score == evidence.f1 &&
                rival.chain_anchors == evidence.cnt;
  if (verdict.tie) {
    verdict.strength = 1.0;
  } else {
    const double winner_vote =
        static_cast<double>(std::max<std::int64_t>(1, evidence.rank1_score));
    verdict.strength =
        static_cast<double>(std::max<std::int64_t>(0, rival.rank_score)) /
        winner_vote;
    if (rival.chained)
      verdict.strength *= chain_ratio(evidence, rival);
    else if (evidence.unchained_zero)
      verdict.strength = 0.0;
  }
  verdict.counted = rival.chained && rival.chain_anchors >= evidence.cnt;
  return verdict;
}

int rna_chain_mapq(const RnaChainMapqEvidence& evidence,
                   RnaChainMapqBreakdown* breakdown,
                   std::vector<RnaChainMapqRivalVerdict>* verdicts) {
  if (breakdown != nullptr)
    *breakdown = RnaChainMapqBreakdown();
  if (verdicts != nullptr) {
    verdicts->clear();
    verdicts->reserve(evidence.rivals.size());
  }
  if (evidence.f1 <= 0) {
    if (verdicts != nullptr)
      for (const RnaChainMapqRival& rival : evidence.rivals)
        verdicts->push_back(rna_chain_mapq_rival_verdict(evidence, rival));
    return 0;
  }

  const double f1 = static_cast<double>(evidence.f1);
  const double dp1 = evidence.dp1;
  const double dp2 = evidence.dp2;

  // The committed window's sibling chain enters x with the same floor. It is never
  // realized, so it also seeds x_floor.
  const double sib_limb = static_cast<double>(std::max(
                              evidence.sib_f2, kRnaChainMapqMinChainScore)) /
                          f1;
  double x = sib_limb;
  double x_floor = sib_limb;
  int f2 = std::max(evidence.sib_f2, 0);
  int n_sub = 0;
  int admissible_rivals = 0;
  int tie_rivals = 0;
  double vote_ratio = 0.0;
  const double winner_vote =
      static_cast<double>(std::max<std::int64_t>(1, evidence.rank1_score));
  // The damped disjoint-vote floor v, from the catalogue's demoted locus rather than
  // the rivals below, which often do not include it. It reaches x only: not f2, n_sub,
  // x_floor or vote_ratio. Guarded so damp = 0 performs no float operation.
  if (evidence.disjoint_vote_damp > 0.0) {
    const double demoted_vote =
        evidence.disjoint_vote_damp *
        static_cast<double>(std::max<std::int64_t>(0, evidence.rank2_score)) /
        winner_vote;
    x = std::max(x, std::min(1.0, demoted_vote));
  }

  std::vector<RnaChainMapqRivalVerdict> local;
  std::vector<RnaChainMapqRivalVerdict>& judged =
      verdicts != nullptr ? *verdicts : local;
  for (const RnaChainMapqRival& rival : evidence.rivals) {
    judged.push_back(rna_chain_mapq_rival_verdict(evidence, rival));
    RnaChainMapqRivalVerdict& verdict = judged.back();
    if (!verdict.admissible)
      continue;
    ++admissible_rivals;
    if (verdict.tie)
      ++tie_rivals;
    vote_ratio = std::max(
        vote_ratio,
        static_cast<double>(std::max<std::int64_t>(0, rival.rank_score)) /
            winner_vote);
    if (rival.chained) {
      f2 = std::max(f2, rival.chain_score);
      bool counts = rival.chain_anchors >= evidence.cnt;
      // minimap2's DP rule for the realized rival: within sub_diff of the primary it
      // counts whatever its anchor count; a per-rival OR, as in mm_set_parent.
      if (rival.owns_dp2 && dp1 > 0.0 && dp2 > 0.0 &&
          dp1 - dp2 <= static_cast<double>(evidence.sub_diff))
        counts = true;
      verdict.counted = counts;
      if (counts)
        ++n_sub;
    }
    x = std::max(x, verdict.strength);
    // A rival that was never realized is invisible to rho and keeps its undamped
    // strength in x_floor; a realized one is already bounded by rho (dp2 is the
    // largest realized non-elected maximum).
    if (!(rival.realized_dp || rival.owns_dp2))
      x_floor = std::max(x_floor, verdict.strength);
  }
  // minimap2's pen, with the query-coverage brake as a third limb in place of the
  // rep_len discount.
  const double pen_s1 = evidence.f1 > 100 ? 1.0 : 0.01 * f1;
  double pen =
      evidence.cnt > 10 ? 1.0 : 0.1 * static_cast<double>(evidence.cnt);
  pen = std::min(pen_s1, pen);
  // The denominator is read_len unless the evidence names its own (a second family's
  // realized hull).
  const double qcov =
      static_cast<double>(evidence.selected_query_covered_bases) /
      static_cast<double>(std::max(1, evidence.qcov_bases > 0
                                          ? evidence.qcov_bases
                                          : evidence.read_len));
  const double qcov_brake =
      std::min(1.0, std::max(0.0, qcov / evidence.qcov_tau));
  pen = std::min(pen, qcov_brake);

  const double match_sc = static_cast<double>(std::max(1, evidence.match_sc));
  const bool dp_branch = dp1 > 0.0;
  // Branch C fires when a rival was realized: rho = dp2/dp1 damps x by how close that
  // rival came at base level, and x_floor keeps the rivals rho cannot describe at their
  // undamped strength.
  const bool branch_c = dp_branch && dp2 > 0.0;
  double rho = 0.0;
  double x_c = x;
  if (branch_c) {
    rho = dp2 / dp1;
    x_c = std::max(rho * x, x_floor);
  }
  int mapq;
  if (branch_c) {
    // Branch C: the committed hypothesis is scored on the square of the damped
    // ambiguity.
    mapq = truncate_to_int(evidence.identity * pen * kRnaChainMapqCoef *
                           (1.0 - x_c * x_c) * std::log(dp1 / match_sc));
  } else if (dp_branch) {
    // Branch B: aligned, but no realized rival to damp against.
    mapq = truncate_to_int(evidence.identity * pen * kRnaChainMapqCoef *
                           (1.0 - x) * std::log(dp1 / match_sc));
  } else {
    // Branch A: map-only, scored on the chain alone.
    mapq = truncate_to_int(pen * kRnaChainMapqCoef * (1.0 - x) * std::log(f1));
  }
  // minimap2's BWA-MEM-like near-tie cap, in minimap2's order: inside the realized
  // branch, before the n_sub subtraction, and only when both realized maxima exist. A
  // min, so it can only demote.
  int bwa_cap = -1;
  bool bwa_cap_bound = false;
  if (branch_c) {
    bwa_cap =
        truncate_to_int(6.02 * kRnaChainMapqBwaCapScale * evidence.identity *
                            evidence.identity * (dp1 - dp2) / match_sc +
                        0.499);
    if (bwa_cap < mapq) {
      mapq = bwa_cap;
      bwa_cap_bound = true;
    }
  }
  // minimap2's n_sub penalty.
  mapq -= truncate_to_int(4.343 * std::log(static_cast<double>(n_sub) + 1.0) +
                          0.499);
  mapq = std::max(0, mapq);
  mapq = std::min(60, mapq);

  if (breakdown != nullptr) {
    breakdown->rivals_total = static_cast<int>(evidence.rivals.size());
    breakdown->rivals_admissible = admissible_rivals;
    breakdown->rivals_tie = tie_rivals;
    breakdown->f2 = f2;
    breakdown->n_sub = n_sub;
    breakdown->vote_ratio = vote_ratio;
    breakdown->x_used = x;
    breakdown->pen = pen;
    breakdown->qcov_brake = qcov_brake;
    breakdown->identity = dp_branch ? evidence.identity : 1.0;
    breakdown->dp_branch = dp_branch;
    breakdown->branch_c = branch_c;
    breakdown->rho = rho;
    breakdown->x_floor = x_floor;
    breakdown->x_c = x_c;
    breakdown->bwa_cap = bwa_cap;
    breakdown->bwa_cap_bound = bwa_cap_bound;
  }
  return mapq;
}

} // namespace fa::cpu::lr::rna
