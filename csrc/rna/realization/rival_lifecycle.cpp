#include "rival_lifecycle.h"

#include <algorithm>
#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

RnaRivalClass rna_classify_rival(bool spans_compete) noexcept {
  return spans_compete ? RnaRivalClass::Competitor : RnaRivalClass::CoPrimary;
}

bool rna_competitor_admitted(int chain_score, int primary_chain_score,
                             const RnaRivalLifecycleConfig& config) noexcept {
  if (chain_score <= 0)
    return false;
  // The ratio band in double; the exact integer band below carries the near ties.
  if (static_cast<double>(chain_score) >=
      config.pri_ratio * static_cast<double>(primary_chain_score))
    return true;
  return chain_score + config.min_diff >= primary_chain_score;
}

bool rna_coprimary_admitted(int chain_score, int chain_anchors,
                            int uncovered_query_bases) noexcept {
  return chain_score >= kRnaCoprimaryMinChainScore &&
         chain_anchors >= kRnaCoprimaryMinAnchors &&
         uncovered_query_bases >= kRnaCoprimaryMinQueryBases;
}

bool rna_chimera_query_disjoint(RnaQuerySpan a, RnaQuerySpan b) noexcept {
  const auto measured = [](RnaQuerySpan span) noexcept {
    return span.begin >= 0 && span.end > span.begin;
  };
  // Censored geometry is not disjointness; an unmeasurable hull is not emitted.
  if (!measured(a) || !measured(b))
    return false;
  const int overlap = std::min(a.end, b.end) - std::max(a.begin, b.begin);
  if (overlap <= 0)
    return true;
  const int shorter = std::min(a.end - a.begin, b.end - b.begin);
  // Integer floor: overlap at most shorter * numerator / denominator.
  const std::int64_t threshold =
      (static_cast<std::int64_t>(shorter) *
       static_cast<std::int64_t>(kRnaChimeraOverlapNumerator)) /
      static_cast<std::int64_t>(kRnaChimeraOverlapDenominator);
  return static_cast<std::int64_t>(overlap) <= threshold;
}

bool rna_chimera_family_admitted(std::optional<int> dp_segment0,
                                 RnaQuerySpan hull,
                                 int minimum_dp_maximum) noexcept {
  if (!dp_segment0.has_value() || *dp_segment0 < minimum_dp_maximum)
    return false;
  if (hull.begin < 0 || hull.end <= hull.begin)
    return false;
  return hull.end - hull.begin >= kRnaChimeraMinQueryBases;
}

bool rna_rival_queue_before(const RnaRealizedHypothesis& a,
                            const RnaRealizedHypothesis& b) noexcept {
  if (a.chain_score != b.chain_score)
    return a.chain_score > b.chain_score;
  if (a.chain_anchors != b.chain_anchors)
    return a.chain_anchors > b.chain_anchors;
  return a.catalogue_index < b.catalogue_index;
}

namespace {

// minimap2's mg_log2, copied verbatim so the value is identical on every platform
// (std::log2's last bit is unspecified, and this feeds a placement decision). The
// argument is 1 + a CIGAR run length, and the realizer refuses zero-length runs.
float mg_log2(float x) noexcept { // NB: this doesn't work when x<2
  union {
    float f;
    std::uint32_t i;
  } z = {x};
  float log_2 = static_cast<float>(((z.i >> 23) & 255) - 128);
  z.i &= ~(255u << 23);
  z.i += 127u << 23;
  log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
  return log_2;
}

// Insertion or deletion in the realizer's encoding, minimap2's `op == 1 || op == 2`. An
// intron (op 3) is not a gap here, so it costs nothing.
bool is_gap_operation(std::uint32_t encoded) noexcept {
  const std::uint32_t operation = encoded & 0xfu;
  return operation == 1u || operation == 2u;
}

// A borrowed CIGAR fit to walk: present, non-empty and without zero-length runs.
bool recal_cigar_walkable(const RnaRankRecalAccounting& accounting) noexcept {
  if (accounting.cigar == nullptr || accounting.cigar_len == 0)
    return false;
  for (std::size_t i = 0; i < accounting.cigar_len; ++i)
    if ((accounting.cigar[i] >> 4) == 0u)
      return false;
  return true;
}

} // namespace

double rna_event_identity(const RnaRankRecalAccounting& accounting) noexcept {
  if (!recal_cigar_walkable(accounting))
    return -1.0;
  if (accounting.matches < 0 || accounting.mismatches < 0 ||
      accounting.ambiguities < 0)
    return -1.0;
  // n_gap counts gap bases, n_gapo gap operations.
  std::int64_t n_gap = 0;
  std::int64_t n_gapo = 0;
  for (std::size_t i = 0; i < accounting.cigar_len; ++i) {
    if (!is_gap_operation(accounting.cigar[i]))
      continue;
    n_gap += static_cast<std::int64_t>(accounting.cigar[i] >> 4);
    ++n_gapo;
  }
  // minimap2's blen, rebuilt from the counters; it equals the controller's block length.
  const std::int64_t blen = static_cast<std::int64_t>(accounting.matches) +
                            static_cast<std::int64_t>(accounting.mismatches) +
                            n_gap;
  const std::int64_t denominator =
      blen + static_cast<std::int64_t>(accounting.ambiguities) - n_gap + n_gapo;
  if (denominator <= 0)
    return -1.0;
  return static_cast<double>(accounting.matches) /
         static_cast<double>(denominator);
}

double rna_rank_recal_b2(double top_event_identity, int match_sc,
                         int mismatch_sc) noexcept {
  double divergence = 1.0 - top_event_identity;
  if (divergence < kRnaRankRecalMinDivergence)
    divergence = kRnaRankRecalMinDivergence;
  double b2 = 0.5 / divergence;
  // minimap2's line, as written (see the header). The guard is ours: a non-positive pair
  // cannot form the ratio.
  if (match_sc > 0 && mismatch_sc > 0 &&
      b2 * static_cast<double>(match_sc) < static_cast<double>(mismatch_sc))
    b2 = static_cast<double>(match_sc) / static_cast<double>(mismatch_sc);
  return b2;
}

int rna_recal_max_dp(const RnaRankRecalAccounting& accounting, double b2,
                     int match_sc) noexcept {
  if (!recal_cigar_walkable(accounting))
    return -1;
  // A gap costs b2 + log2(1 + len) per operation; introns never enter the sum.
  double gap_cost = 0.0;
  for (std::size_t i = 0; i < accounting.cigar_len; ++i) {
    if (!is_gap_operation(accounting.cigar[i]))
      continue;
    const std::uint32_t length = accounting.cigar[i] >> 4;
    gap_cost +=
        b2 + static_cast<double>(mg_log2(1.0f + static_cast<float>(length)));
  }
  // minimap2's n_mis is blen + n_ambi - mlen - n_gap, which under these counters is
  // exactly mismatches + ambiguities: the ambiguous gap bases subtracted from
  // `mismatches` are added back. The realizer's clamp of `mismatches` at 0 can only
  // charge more, never promote.
  const double n_mis = static_cast<double>(accounting.mismatches) +
                       static_cast<double>(accounting.ambiguities);
  // minimap2 truncates match_sc * (...) + .499 toward zero, as static_cast<int> does.
  return static_cast<int>(
      static_cast<double>(match_sc) *
          (static_cast<double>(accounting.matches) - b2 * n_mis - gap_cost) +
      0.499);
}

bool rna_rank_recal_gate(int read_len, int top_aq_begin, int top_aq_end,
                         int top_dp, int second_dp) noexcept {
  if (read_len < kRnaRankRecalMinReadLen)
    return false;
  // Both prices must be real evidence, in the same sense the election uses.
  if (top_dp <= 0 || second_dp <= 0)
    return false;
  // An unmeasured hull (-1/-1) is not a span that covers the read.
  if (top_aq_begin < 0 || top_aq_end < top_aq_begin)
    return false;
  // minimap2's two floating-point tests with the sense inverted. Each is one multiply of
  // exact operands and one comparison, so equality at the boundary passes on every
  // platform, as in minimap2.
  const double span = static_cast<double>(top_aq_end - top_aq_begin);
  if (span < kRnaRankRecalFrac * static_cast<double>(read_len))
    return false;
  return static_cast<double>(second_dp) >=
         kRnaRankRecalFrac * static_cast<double>(top_dp);
}

namespace {

// A hypothesis's realized DP maximum as election evidence: present, realized,
// and strictly positive. Anything else is "no evidence", which the election
// and the dp2 selection both treat as absence rather than as a low score.
bool election_evidence(const RnaRealizedHypothesis& hypothesis,
                       int& score) noexcept {
  if (!hypothesis.realized || !hypothesis.dp_maximum.has_value())
    return false;
  if (*hypothesis.dp_maximum <= 0)
    return false;
  score = *hypothesis.dp_maximum;
  return true;
}

bool competes_for_primary(const RnaRealizedHypothesis& hypothesis) noexcept {
  return hypothesis.klass != RnaRivalClass::CoPrimary;
}

} // namespace

std::size_t rna_elect_primary_hypothesis(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept {
  if (hypotheses.empty())
    return 0;
  int incumbent_dp = 0;
  // Fail closed: without the incumbent's own realized maximum there is no
  // scale to beat, so nothing can be promoted against it.
  if (!election_evidence(hypotheses.front(), incumbent_dp))
    return 0;
  std::size_t elected = 0;
  int best_dp = incumbent_dp;
  int best_catalogue = hypotheses.front().catalogue_index;
  for (std::size_t i = 1; i < hypotheses.size(); ++i) {
    const RnaRealizedHypothesis& hypothesis = hypotheses[i];
    if (!competes_for_primary(hypothesis))
      continue;
    int dp = 0;
    if (!election_evidence(hypothesis, dp))
      continue;
    // Strictly greater promotes; an exact tie keeps whoever holds the slot,
    // and since the incumbent starts holding it a tie retains the incumbent.
    // Between two competitors the lower catalogue index wins the tie.
    if (dp > best_dp || (dp == best_dp && elected != 0 &&
                         hypothesis.catalogue_index < best_catalogue)) {
      elected = i;
      best_dp = dp;
      best_catalogue = hypothesis.catalogue_index;
    }
  }
  return elected;
}

std::size_t rna_elect_primary_rank_recal(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept {
  // The raw verdict is both the fallback for every failure of evidence and
  // one side of the veto below, so it is computed first and exactly once.
  const std::size_t raw_winner = rna_elect_primary_hypothesis(hypotheses);
  if (hypotheses.empty())
    return raw_winner;
  // Raw plus recalibrated, both required. Repricing is atomic, so a missing dp_recal is
  // missing evidence.
  const auto sum_price = [](const RnaRealizedHypothesis& hypothesis,
                            std::int64_t& price) noexcept {
    int raw = 0;
    if (!election_evidence(hypothesis, raw) || !hypothesis.dp_recal.has_value())
      return false;
    price = static_cast<std::int64_t>(raw) +
            static_cast<std::int64_t>(*hypothesis.dp_recal);
    return true;
  };
  std::int64_t incumbent_sum = 0;
  if (!sum_price(hypotheses.front(), incumbent_sum))
    return raw_winner;
  // Same discipline as the raw election: strictly greater promotes, a tie keeps the
  // holder, and between competitors the lower catalogue index wins.
  std::size_t sum_winner = 0;
  std::int64_t best_sum = incumbent_sum;
  int best_catalogue = hypotheses.front().catalogue_index;
  for (std::size_t i = 1; i < hypotheses.size(); ++i) {
    const RnaRealizedHypothesis& hypothesis = hypotheses[i];
    if (!competes_for_primary(hypothesis))
      continue;
    std::int64_t price = 0;
    if (!sum_price(hypothesis, price))
      continue;
    if (price > best_sum || (price == best_sum && sum_winner != 0 &&
                             hypothesis.catalogue_index < best_catalogue)) {
      sum_winner = i;
      best_sum = price;
      best_catalogue = hypothesis.catalogue_index;
    }
  }
  if (sum_winner == raw_winner)
    return sum_winner;
  // The veto: the recalibrated half may overturn at most kRnaRankRecalRawVetoMargin of
  // raw margin (inclusive). Any missing price keeps the raw verdict.
  int raw_of_raw_winner = 0;
  int raw_of_sum_winner = 0;
  if (!election_evidence(hypotheses[raw_winner], raw_of_raw_winner) ||
      !election_evidence(hypotheses[sum_winner], raw_of_sum_winner))
    return raw_winner;
  return raw_of_raw_winner - raw_of_sum_winner <= kRnaRankRecalRawVetoMargin
             ? sum_winner
             : raw_winner;
}

std::optional<std::size_t>
rna_select_dp2_owner(const std::vector<RnaRealizedHypothesis>& hypotheses,
                     std::size_t elected) noexcept {
  std::optional<std::size_t> owner;
  int best_dp = 0;
  int best_catalogue = 0;
  for (std::size_t i = 0; i < hypotheses.size(); ++i) {
    if (i == elected)
      continue;
    const RnaRealizedHypothesis& hypothesis = hypotheses[i];
    if (!competes_for_primary(hypothesis))
      continue;
    // A shadow of the incumbent (RnaRealizedHypothesis::shadow) carries the
    // committed placement's own maximum; it never supplies dp2.
    if (hypothesis.shadow)
      continue;
    int dp = 0;
    if (!election_evidence(hypothesis, dp))
      continue;
    if (!owner.has_value() || dp > best_dp ||
        (dp == best_dp && hypothesis.catalogue_index < best_catalogue)) {
      owner = i;
      best_dp = dp;
      best_catalogue = hypothesis.catalogue_index;
    }
  }
  return owner;
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
