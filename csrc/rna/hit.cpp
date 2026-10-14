// Rivals, minimap2's parent rule and the RNA MAPQ: the rival lifecycle, mm_set_parent over the chain-stage hits, mm_update_dp_max, the pricing of a realized hypothesis, mm_set_mapq2 over the elected chain and its admissible rivals.
#include "rna.h"

#include "../core/core.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>

/**** rival_lifecycle ****/

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

/**** rna_parent ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

void rna_parent_set(std::vector<RnaParentHit>& r, float mask_level,
                    int sub_diff) {
  const int n = static_cast<int>(r.size());
  if (n <= 0)
    return;
  std::vector<std::pair<int, int>> cov;
  cov.reserve(static_cast<std::size_t>(n));
  std::vector<int> w(static_cast<std::size_t>(n));
  w[0] = 0;
  r[0].parent = 0;
  int k = 1;
  for (int i = 1; i < n; ++i) {
    RnaParentHit& ri = r[i];
    const int si = ri.qs;
    const int ei = ri.qe;
    cov.clear();
    for (int j = 0; j < k; ++j) {
      const RnaParentHit& rp = r[w[j]];
      int sj = rp.qs;
      int ej = rp.qe;
      if (ej <= si || sj >= ei)
        continue;
      if (sj < si)
        sj = si;
      if (ej > ei)
        ej = ei;
      cov.emplace_back(sj, ej);
    }
    int j = k;
    if (!cov.empty()) {
      // The query length no existing primary covers.
      int uncov_len = 0;
      std::sort(cov.begin(), cov.end());
      int x = si;
      for (const auto& [cs, ce] : cov) {
        if (cs > x)
          uncov_len += cs - x;
        x = ce > x ? ce : x;
      }
      if (ei > x)
        uncov_len += ei - x;
      for (j = 0; j < k; ++j) {
        RnaParentHit& rp = r[w[j]];
        const int sj = rp.qs;
        const int ej = rp.qe;
        if (ej <= si || sj >= ei)
          continue;
        const int min = ej - sj < ei - si ? ej - sj : ei - si;
        const int max = ej - sj > ei - si ? ej - sj : ei - si;
        const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                               : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
        // minimap2 divides by min unguarded; a zero min gives NaN there and fails the
        // test, as the guard does here.
        if (min > 0 &&
            static_cast<float>(ol) / static_cast<float>(min) -
                    static_cast<float>(uncov_len) / static_cast<float>(max) >
                mask_level) {
          bool cnt_sub = false;
          ri.parent = rp.parent;
          rp.subsc = rp.subsc > ri.score ? rp.subsc : ri.score;
          if (ri.cnt >= rp.cnt)
            cnt_sub = true;
          // The last condition excludes identical hits after DP.
          if (rp.has_dp && ri.has_dp &&
              (rp.rid != ri.rid || rp.rs != ri.rs || rp.re != ri.re ||
               ol != min)) {
            rp.dp_max2 = rp.dp_max2 > ri.dp_max ? rp.dp_max2 : ri.dp_max;
            if (rp.dp_max - ri.dp_max <= sub_diff)
              cnt_sub = true;
          }
          if (cnt_sub)
            ++rp.n_sub;
          break;
        }
      }
    }
    if (j == k) {
      w[static_cast<std::size_t>(k++)] = i;
      ri.parent = i;
      ri.n_sub = 0;
    }
  }
}

bool rna_parent_update_dp_max(int qlen, int min_len,
                              std::vector<RnaParentHit>& r, float frac, int a,
                              int b) {
  if (qlen < min_len || r.size() < 2)
    return false;
  int max = -1;
  int max2 = -1;
  int max_i = -1;
  for (std::size_t i = 0; i < r.size(); ++i) {
    if (!r[i].has_dp)
      continue;
    if (r[i].dp_max > max) {
      max2 = max;
      max = r[i].dp_max;
      max_i = static_cast<int>(i);
    } else if (r[i].dp_max > max2) {
      max2 = r[i].dp_max;
    }
  }
  if (max_i < 0 || max < 0 || max2 < 0)
    return false;
  const RnaParentHit& top = r[static_cast<std::size_t>(max_i)];
  if (top.qe - top.qs < static_cast<double>(qlen) * frac)
    return false;
  if (max2 < static_cast<double>(max) * frac)
    return false;
  for (const RnaParentHit& h : r)
    if (h.has_dp && !h.accounting_ok)
      return false;
  const double identity = rna_event_identity(top.accounting);
  if (identity < 0.0)
    return false;
  const double b2 = rna_rank_recal_b2(identity, a, b);
  for (RnaParentHit& h : r) {
    if (!h.has_dp)
      continue;
    const int dp = rna_recal_max_dp(h.accounting, b2, a);
    h.dp_max = dp < 0 ? 0 : dp;
  }
  return true;
}

void rna_parent_hit_sort(std::vector<RnaParentHit>& r) {
  std::stable_sort(r.begin(), r.end(),
                   [](const RnaParentHit& x, const RnaParentHit& y) {
                     const int kx = x.has_dp ? x.dp_max : x.score;
                     const int ky = y.has_dp ? y.dp_max : y.score;
                     if (kx != ky)
                       return kx > ky;
                     if (x.catalogue != y.catalogue)
                       return x.catalogue < y.catalogue;
                     return x.segment < y.segment;
                   });
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** rank2_arbitration ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** rival_pricing ****/

namespace fa::cpu::lr::rna {

namespace {

RnaQuerySpan segment_span(const SpliceRegionEvidence& evidence) noexcept {
  RnaQuerySpan span;
  span.begin = evidence.query_begin;
  span.end = evidence.query_end;
  return span;
}

bool span_valid(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

} // namespace

const RnaSpliceHypothesisResult*
rna_winning_hypothesis(const RnaSpliceRealizationResult& result) noexcept {
  if (result.refused || result.segments.empty())
    return nullptr;
  for (const RnaSpliceHypothesisResult& hypothesis : result.hypotheses) {
    if (hypothesis.orientation == result.winning_hypothesis &&
        !hypothesis.refused && !hypothesis.segments.empty())
      return &hypothesis;
  }
  return nullptr;
}

std::optional<int>
rna_segment0_dp_maximum(const RnaSpliceRealizationResult& result) noexcept {
  const RnaSpliceHypothesisResult* winner = rna_winning_hypothesis(result);
  if (winner == nullptr)
    return std::nullopt;
  return winner->primary_dp_maximum;
}

std::optional<RnaQuerySpan>
rna_realized_query_span(const RnaSpliceRealizationResult& result) noexcept {
  const RnaSpliceHypothesisResult* winner = rna_winning_hypothesis(result);
  if (winner == nullptr)
    return std::nullopt;
  std::optional<RnaQuerySpan> hull;
  for (const SpliceRegionEvidence& evidence : winner->region_evidence) {
    const RnaQuerySpan span = segment_span(evidence);
    if (!span_valid(span))
      continue;
    if (!hull) {
      hull = span;
    } else {
      hull->begin = std::min(hull->begin, span.begin);
      hull->end = std::max(hull->end, span.end);
    }
  }
  return hull;
}

} // namespace fa::cpu::lr::rna

/**** chain_mapq ****/
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

namespace fa::cpu::lr::rna {

namespace {

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
  if (evidence.parent_rule && rival.chained)
    return rival.parented;
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

namespace {

// The per-rival verdict, except `counted`'s DP rule, which rna_chain_mapq fills.
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

} // namespace

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
  if (evidence.parent_rule)
    n_sub = evidence.parent_n_sub;
  // minimap2's pen, with the query-coverage brake as a third limb in place of the
  // rep_len discount.
  const double pen_s1 = evidence.f1 > 100 ? 1.0 : 0.01 * f1;
  double pen =
      evidence.cnt > 10 ? 1.0 : 0.1 * static_cast<double>(evidence.cnt);
  pen = std::min(pen_s1, pen);
  // The denominator is read_len unless the evidence names its own (a piece's region).
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

/**** plain_election ****/

namespace fa::cpu::lr::rna {

namespace {

// The chain's junctions: consecutive chained anchors whose reference gap exceeds their
// query gap by at least min_intron. `path` must be a bundle summarize_exact_anchor_path
// accepted, so its selected indices are in range and in chain order.
int rna_plain_elect_junctions(const ExactAnchorPath& path,
                              int min_intron) noexcept {
  int junctions = 0;
  for (std::size_t i = 1; i < path.selected_raw_indices.size(); ++i) {
    const SelectedPoolAnchor& left =
        path.selected_pool_anchors[path.selected_raw_indices[i - 1]];
    const SelectedPoolAnchor& right =
        path.selected_pool_anchors[path.selected_raw_indices[i]];
    const std::int64_t reference_gap =
        static_cast<std::int64_t>(right.reference_begin) - left.reference_begin;
    const std::int64_t query_gap =
        static_cast<std::int64_t>(right.query_begin) - left.query_begin;
    if (reference_gap - query_gap >= min_intron)
      ++junctions;
  }
  return junctions;
}

// Anchor count plus the junction credit, in thousandths of an anchor.
std::int64_t rna_plain_elect_score(const ExactAnchorPath& path,
                                   int min_intron) noexcept {
  return static_cast<std::int64_t>(path.anchor_path_anchor_count) *
             kRnaPlainElectScale +
         rna_plain_elect_junctions(path, min_intron) *
             kRnaPlainElectJunctionCredit;
}

} // namespace

// The election (rna.h): the incumbent and every chained rival are scored,
// and the best eligible rival takes the primary when it scores strictly higher; a
// tie between rivals goes to the lower catalogue index. An adoption swaps the roles
// as the CIGAR lane's lifecycle does: the elected chain owns the MAPQ evidence and
// the record, the demoted incumbent joins the rivals in catalogue order, rs2 is
// rank 1's vote, and every shadow verdict is re-read against the elected chain.
void rna_plain_election(const ReadFacts& facts, const ReadState& rs,
                        const AttemptState& as, const ExactAnchorPath& anchor_path,
                        const std::vector<ProductionRival>& production_rivals,
                        RnaChainMapqEvidence& chain_evidence,
                        const placement::CoarseLocus*& committed_locus,
                        const ExactAnchorPathSummary*& committed_summary,
                        std::uint64_t& committed_qcov) {
  const auto& min_intron = facts.min_intron;
  const auto& coarse_pass = rs.coarse_pass;
  const auto& catalogue_rank1_score = rs.catalogue_rank1_score;
  const auto& read_len = facts.read_len;
  const auto& selected_path_summary = as.selected_path_summary;
  const auto& selected_locus = as.selected_locus;
  const auto& selected_reverse = as.selected_reverse;
  const auto& committed_geometry = as.committed_geometry;
  const auto& max_intron = facts.max_intron;
  const std::int64_t incumbent_score =
      rna::rna_plain_elect_score(anchor_path, min_intron);
  const std::size_t rivals = production_rivals.size();
  std::vector<std::int64_t> rival_scores(rivals, 0);
  std::vector<std::uint64_t> rival_qcov(rivals, 0);
  std::optional<std::size_t> elected;
  for (std::size_t i = 0; i < rivals; ++i) {
    const ProductionRival& entry = production_rivals[i];
    if (!entry.rival.chained)
      continue;
    const int qspan = entry.rival.q_end - entry.rival.q_begin;
    rival_scores[i] =
        rna::rna_plain_elect_score(entry.chained_path, min_intron);
    // Eligible on the incumbent's terms: an anchor union the record can report,
    // and on a second-pass catalogue the rescue floor.
    const bool eligible =
        rna::anchor_union_query_coverage(entry.chained_path,
                                         rival_qcov[i]) &&
        (coarse_pass != 2 || qspan >= rna::kRnaChimeraMinQueryBases);
    if (eligible &&
        (!elected || rival_scores[i] > rival_scores[*elected]))
      elected = i;
  }
  if (elected && rival_scores[*elected] > incumbent_score) {
    const ProductionRival& winner = production_rivals[*elected];
    committed_locus = &winner.locus;
    committed_summary = &winner.summary;
    committed_qcov = rival_qcov[*elected];
    chain_evidence.f1 = winner.rival.chain_score;
    chain_evidence.cnt = winner.rival.chain_anchors;
    chain_evidence.sib_f2 = winner.siblings.sib_score;
    chain_evidence.winner_q_begin = winner.rival.q_begin;
    chain_evidence.winner_q_end = winner.rival.q_end;
    chain_evidence.rank1_score = winner.locus.rank_score;
    chain_evidence.rank2_score = catalogue_rank1_score;
    chain_evidence.winner_reverse = winner.rival.reverse;
    chain_evidence.selected_query_covered_bases = committed_qcov;
    // The shadow rule's geometry, as chain_one_rival builds it.
    const auto chain_geometry = [&](const ProductionRival& entry) {
      rna::RnaSegmentGeometry geometry;
      geometry.reference_id = entry.reference_id;
      geometry.reverse = entry.rival.reverse;
      geometry.reference_begin = static_cast<std::uint64_t>(
          std::max(0, entry.summary.first_reference_begin));
      geometry.reference_end = static_cast<std::uint64_t>(
          std::max(0, entry.summary.last_reference_end));
      geometry.forward_span = rna::rna_forward_query_span(
          entry.rival.q_begin, entry.rival.q_end, entry.rival.reverse,
          read_len);
      return geometry;
    };
    const rna::RnaSegmentGeometry winner_geometry =
        chain_geometry(winner);
    chain_evidence.rivals.clear();
    rna::RnaChainMapqRival demoted;
    demoted.candidate = 0;
    demoted.chained = true;
    demoted.chain_score = anchor_path.anchor_path_score;
    demoted.chain_anchors =
        static_cast<int>(anchor_path.anchor_path_anchor_count);
    demoted.q_begin = selected_path_summary.first_query_begin;
    demoted.q_end = selected_path_summary.last_query_end;
    demoted.rank_score = selected_locus->rank_score;
    demoted.reverse = selected_reverse;
    // committed_geometry is still the incumbent's chain.
    demoted.shadow = rna::rna_same_place(winner_geometry,
                                         committed_geometry, max_intron);
    chain_evidence.rivals.push_back(demoted);
    for (std::size_t i = 0; i < rivals; ++i) {
      if (i == *elected)
        continue;
      rna::RnaChainMapqRival rival = production_rivals[i].rival;
      rival.shadow =
          rival.chained &&
          rna::rna_same_place(winner_geometry,
                              chain_geometry(production_rivals[i]),
                              max_intron);
      chain_evidence.rivals.push_back(rival);
    }
  }
}

} // namespace fa::cpu::lr::rna

/**** parent_rule ****/
// The lifecycle's helpers: binding and pricing a realization, mm_set_parent over a
// stage's hits, the shadow rule's geometry, a hypothesis's realized regions, and the
// supplementary pieces.

namespace fa::cpu::lr::rna {

namespace {

// minimap2's mask_level for the parent rule.
constexpr float kParentMaskLevel = 0.5f;

} // namespace

// Binds a realization request to contig reference_id: its sequence, its name and
// the known junctions on it. False when the contig will not bind.
bool bind_realization_reference(const ReadFacts& facts,
                                RnaSpliceRealizationRequest& target,
                                int reference_id) {
  const auto& rctx = facts.rctx;
  return rna::bind_splice_realization_reference(
      target, *rctx.ref.encoded, *rctx.ref.names,
      rctx.opts.known_junctions.get(), reference_id);
}

// Prices one realization for the lifecycle: the election price is segment 0's
// dp_maximum, and the realized hull feeds the recalibration gate and the emission
// floor.
void price_hypothesis(RnaRealizedHypothesis& h,
                      const RnaSpliceRealizationResult& result) {
  h.dp_segment0 = rna::rna_segment0_dp_maximum(result);
  if (const auto hull = rna::rna_realized_query_span(result)) {
    h.aq_begin = hull->begin;
    h.aq_end = hull->end;
  }
  h.dp_maximum = h.dp_segment0;
}

// mm_set_parent over hits already in their stage's order, less the hits whose owner
// `left_out` marks.
std::vector<RnaParentHit> parent_regions(const Lifecycle& lc,
                                         const std::vector<RnaParentHit>& ordered,
                                         const std::vector<char>* left_out) {
  const auto& parent_sub_diff = lc.parent_sub_diff;
  std::vector<rna::RnaParentHit> regions;
  regions.reserve(ordered.size());
  for (const rna::RnaParentHit& hit : ordered) {
    const std::size_t owner = static_cast<std::size_t>(hit.owner);
    if (left_out != nullptr && hit.owner >= 0 && owner < left_out->size() &&
        (*left_out)[owner] != 0)
      continue;
    regions.push_back(hit);
  }
  rna::rna_parent_set(regions, kParentMaskLevel, parent_sub_diff);
  return regions;
}

// A chained hypothesis's place for the shadow rule; false without a chain.
bool shadow_geometry(const ReadFacts& facts, const RnaRealizedHypothesis& h,
                     RnaSegmentGeometry& g) {
  const auto& read_len = facts.read_len;
  if (!h.chained || h.summary == nullptr ||
      h.summary->last_reference_end <= h.summary->first_reference_begin)
    return false;
  g.reference_id = h.reference_id;
  g.reverse = h.reverse;
  g.reference_begin = static_cast<std::uint64_t>(
      std::max(0, h.summary->first_reference_begin));
  g.reference_end = static_cast<std::uint64_t>(
      std::max(0, h.summary->last_reference_end));
  g.forward_span = rna::rna_forward_query_span(h.q_begin, h.q_end,
                                               h.reverse, read_len);
  return true;
}

// Whether hypotheses a and b are chained at the same place (rna_same_place); false
// when either has no chain.
bool same_place(const ReadFacts& facts, const Lifecycle& lc, std::size_t a,
                std::size_t b) {
  const auto& hypotheses = lc.hypotheses;
  const auto& max_intron = facts.max_intron;
  rna::RnaSegmentGeometry ga;
  rna::RnaSegmentGeometry gb;
  return shadow_geometry(facts, hypotheses[a], ga) &&
         shadow_geometry(facts, hypotheses[b], gb) &&
         rna::rna_same_place(ga, gb, max_intron);
}

// Hypothesis i's realized regions, one per segment of `result`, as the parent rule
// reads them.
void append_dp_hits(const Lifecycle& lc, std::size_t i,
                    const RnaSpliceRealizationResult& result,
                    std::vector<RnaParentHit>& hits) {
  const auto& hypotheses = lc.hypotheses;
  const rna::RnaRealizedHypothesis& h = hypotheses[i];
  const rna::RnaSpliceHypothesisResult* winner =
      rna::rna_winning_hypothesis(result);
  if (winner == nullptr)
    return;
  for (std::size_t s = 0; s < winner->segments.size(); ++s) {
    const rna::SpliceRegionEvidence& evidence = winner->region_evidence[s];
    const AlignResult& segment = winner->segments[s];
    const std::vector<std::uint32_t>& cigar = winner->segment_cigars[s];
    rna::RnaParentHit hit;
    hit.qs = evidence.query_begin;
    hit.qe = evidence.query_end;
    hit.score = evidence.chain_score;
    hit.cnt = evidence.region_anchor_count;
    hit.rid = h.reference_id;
    hit.rs = segment.pos;
    hit.re = segment.target_end;
    hit.has_dp = true;
    hit.dp_max = evidence.dp_maximum;
    hit.accounting.matches = segment.matches;
    hit.accounting.mismatches = segment.mismatches;
    hit.accounting.ambiguities = segment.ambiguities;
    hit.accounting.cigar = cigar.data();
    hit.accounting.cigar_len = cigar.size();
    hit.accounting_ok = segment.alignment_accounting_valid && !cigar.empty();
    hit.owner = static_cast<int>(i);
    hit.catalogue = h.catalogue_index;
    hit.segment = static_cast<int>(s);
    hits.push_back(hit);
  }
}

// The pieces from the after-election regions `dp` and the chain stage's `chain`,
// rebuilt from scratch: every after-DP primary but the elected family's own regions, in
// hit order, unless its hypothesis is at the place of a printed record. Its children
// less its shadows are counted as mm_set_parent counts them.
void decide_pieces(const ReadFacts& facts, Lifecycle& lc,
                   const std::vector<RnaParentHit>& dp,
                   const std::vector<RnaParentHit>& chain) {
  auto& pieces = lc.pieces;
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& parent_sub_diff = lc.parent_sub_diff;
  auto& hypotheses = lc.hypotheses;
  pieces.clear();
  for (std::size_t j = 1; j < dp.size(); ++j) {
    const rna::RnaParentHit& region = dp[j];
    const std::size_t owner = static_cast<std::size_t>(region.owner);
    if (region.parent != static_cast<int>(j) || owner == lifecycle_elected ||
        same_place(facts, lc, owner, lifecycle_elected))
      continue;
    if (std::any_of(pieces.begin(), pieces.end(),
                    [&](const PieceRegion& printed) {
                      const std::size_t other =
                          static_cast<std::size_t>(printed.region.owner);
                      return other != owner && same_place(facts, lc, other, owner);
                    }))
      continue;
    PieceRegion piece;
    piece.region = region;
    piece.region.n_sub = 0;
    piece.region.dp_max2 = 0;
    const auto add_child = [&](std::size_t child_owner) {
      if (child_owner != owner &&
          std::find(piece.children.begin(), piece.children.end(),
                    child_owner) == piece.children.end())
        piece.children.push_back(child_owner);
    };
    for (std::size_t c = j + 1; c < dp.size(); ++c) {
      const rna::RnaParentHit& child = dp[c];
      const std::size_t child_owner = static_cast<std::size_t>(child.owner);
      if (child.parent != static_cast<int>(j) ||
          (child_owner != owner && same_place(facts, lc, owner, child_owner)))
        continue;
      const int si = child.qs;
      const int ei = child.qe;
      const int sj = region.qs;
      const int ej = region.qe;
      const int min = std::min(ej - sj, ei - si);
      const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                             : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
      bool cnt_sub = child.cnt >= region.cnt;
      if (region.has_dp && child.has_dp &&
          (region.rid != child.rid || region.rs != child.rs ||
           region.re != child.re || ol != min)) {
        piece.region.dp_max2 = std::max(piece.region.dp_max2, child.dp_max);
        if (region.dp_max - child.dp_max <= parent_sub_diff)
          cnt_sub = true;
      }
      if (cnt_sub)
        ++piece.region.n_sub;
      add_child(child_owner);
    }
    for (std::size_t c = j + 1; c < dp.size() && !piece.dp2_owner; ++c) {
      const rna::RnaParentHit& child = dp[c];
      const std::size_t child_owner = static_cast<std::size_t>(child.owner);
      if (child.parent == static_cast<int>(j) && piece.region.dp_max2 > 0 &&
          child.dp_max == piece.region.dp_max2 && child_owner != owner &&
          !same_place(facts, lc, owner, child_owner))
        piece.dp2_owner = child_owner;
    }
    // Its chain's chain-stage children, when that chain was a primary there.
    for (std::size_t e = 0; e < chain.size(); ++e) {
      if (static_cast<std::size_t>(chain[e].owner) != owner)
        continue;
      if (chain[e].parent == static_cast<int>(e))
        for (std::size_t c = e + 1; c < chain.size(); ++c) {
          const std::size_t child_owner =
              static_cast<std::size_t>(chain[c].owner);
          if (chain[c].parent == static_cast<int>(e) &&
              !same_place(facts, lc, owner, child_owner))
            add_child(child_owner);
        }
      break;
    }
    piece.place.push_back(owner);
    for (std::size_t i = 0; i < hypotheses.size(); ++i)
      if (i != owner && same_place(facts, lc, owner, i))
        piece.place.push_back(i);
    pieces.push_back(std::move(piece));
  }
}

} // namespace fa::cpu::lr::rna

/**** lifecycle ****/
// Stage 5's rival lifecycle, in minimap2's order: the hypotheses, the realization
// queue, the rank recalibration, the election and the pieces. map.cpp's run_lifecycle
// calls them around the rival chains, the realizations and the adoption.

namespace fa::cpu::lr::rna {

// The lifecycle's stages (2) and (3): the hypothesis vector, the incumbent first and
// then one hypothesis per production rival, classified against the incumbent and
// admitted.
void build_hypotheses(const ReadFacts& facts, const ReadState& rs,
                      const AttemptState& as, Lifecycle& lc,
                      const ExactAnchorPath& anchor_path,
                      const RnaRivalLifecycleConfig& lifecycle_cfg,
                      const RnaSpliceRealizationResult& realized) {
  auto& hypotheses = lc.hypotheses;
  auto& production_rivals = lc.production_rivals;
  const auto& chr = as.chr;
  const auto& selected_locus = as.selected_locus;
  const auto& selected_path_summary = as.selected_path_summary;
  const auto& rank1_siblings = as.rank1_siblings;
  const auto& selected_reverse = as.selected_reverse;
  const auto& repaired_orientation = as.repaired_orientation;
  const auto& selected_query_covered_bases = as.selected_query_covered_bases;
  const auto& path_coverage_ok = as.path_coverage_ok;
  const auto& read_len = facts.read_len;
  const auto& loci = rs.loci;
  // (2) The hypothesis vector, incumbent first. Every pointer borrows storage that
  // outlives it (the rank-1 chain, or a rival chain in `production_rivals`, which no
  // longer grows), so reordering `hypotheses` moves no bundle.
  hypotheses.reserve(production_rivals.size() + 1);
  {
    rna::RnaRealizedHypothesis incumbent;
    incumbent.catalogue_index = 0;
    incumbent.reference_id = chr;
    incumbent.locus = selected_locus;
    incumbent.bundle = &anchor_path;
    incumbent.summary = &selected_path_summary;
    incumbent.siblings = rank1_siblings;
    incumbent.rank_score = selected_locus->rank_score;
    incumbent.chained = true;
    incumbent.chain_score = anchor_path.anchor_path_score;
    incumbent.chain_anchors =
        static_cast<int>(anchor_path.anchor_path_anchor_count);
    incumbent.q_begin = selected_path_summary.first_query_begin;
    incumbent.q_end = selected_path_summary.last_query_end;
    incumbent.reverse = selected_reverse;
    incumbent.repaired = repaired_orientation;
    incumbent.query_covered_bases = selected_query_covered_bases;
    incumbent.coverage_ok = path_coverage_ok;
    incumbent.klass = rna::RnaRivalClass::Incumbent;
    incumbent.admitted = true;
    // Rank 1 was realized above into `realized`; the hypothesis carries only its
    // price.
    incumbent.realized = true;
    price_hypothesis(incumbent, realized);
    hypotheses.push_back(std::move(incumbent));
  }

  // (3) Classification and admission, in the read's forward frame: chains on opposite
  // strands can carry identical oriented spans over opposite ends of the read.
  // `reverse` is the post-repair strand on both sides.
  const rna::RnaQuerySpan incumbent_span = rna::rna_forward_query_span(
      selected_path_summary.first_query_begin,
      selected_path_summary.last_query_end, selected_reverse, read_len);
  for (const ProductionRival& entry : production_rivals) {
    const std::size_t cat =
        static_cast<std::size_t>(std::max(0, entry.rival.candidate));
    rna::RnaRealizedHypothesis h;
    h.catalogue_index = entry.rival.candidate;
    h.reference_id = entry.reference_id;
    // A chained rival owns its window (the catalogue entry or its adopted mirror); an
    // unchained one has only the catalogue's nomination.
    h.locus = entry.rival.chained || entry.repaired
                  ? &entry.locus
                  : (cat < loci.size() ? &loci[cat] : nullptr);
    h.rank_score = entry.rival.rank_score;
    h.bundle = entry.owns_path ? &entry.chained_path : nullptr;
    h.summary = &entry.summary;
    h.siblings = entry.siblings;
    h.chained = entry.rival.chained;
    h.chain_score = entry.rival.chain_score;
    h.chain_anchors = entry.rival.chain_anchors;
    h.q_begin = entry.rival.q_begin;
    h.q_end = entry.rival.q_end;
    h.reverse = entry.rival.reverse;
    h.repaired = entry.repaired;
    h.shadow = entry.rival.shadow;
    const rna::RnaQuerySpan rival_span =
        rna::rna_forward_query_span(entry.rival.q_begin, entry.rival.q_end,
                                    entry.rival.reverse, read_len);
    h.klass = rna::rna_classify_rival(
        rna::rna_query_spans_compete(incumbent_span, rival_span));
    if (h.chained && h.bundle != nullptr) {
      h.coverage_ok =
          rna::anchor_union_query_coverage(*h.bundle, h.query_covered_bases);
      if (h.klass == rna::RnaRivalClass::CoPrimary) {
        // A co-primary explains different query bases, so it is neither an election
        // candidate nor dp2/n_sub evidence; its admission is read by the traces only.
        h.admitted = rna::rna_coprimary_admitted(
            h.chain_score, h.chain_anchors,
            rna::rna_uncovered_query_bases(incumbent_span, rival_span));
      } else {
        h.admitted =
            h.coverage_ok &&
            rna::rna_competitor_admitted(
                h.chain_score, anchor_path.anchor_path_score, lifecycle_cfg);
      }
    }
    hypotheses.push_back(std::move(h));
  }
}

// The lifecycle's stage (4): the realization queue, the admitted competitors within
// the budget and then the chain-stage co-primaries, whose chain-stage hits it leaves
// in lc.parent_chain_hits.
std::vector<std::size_t> queue_realizations(const ReadFacts& facts, ReadState& rs,
                                            Lifecycle& lc,
                                            const RnaRivalLifecycleConfig& lifecycle_cfg) {
  auto& hypotheses = lc.hypotheses;
  auto& parent_chain_hits = lc.parent_chain_hits;
  const auto& read_len = facts.read_len;
  // (4) The realization budget, spent on admitted competitors in a total order (chain
  // score, anchor count, catalogue index) so the choice does not depend on the sort.
  std::vector<std::size_t> queue;
  for (std::size_t i = 1; i < hypotheses.size(); ++i) {
    const rna::RnaRealizedHypothesis& candidate = hypotheses[i];
    if (candidate.admitted &&
        candidate.klass == rna::RnaRivalClass::Competitor &&
        candidate.bundle != nullptr)
      queue.push_back(i);
  }
  std::sort(queue.begin(), queue.end(),
            [&hypotheses](std::size_t a, std::size_t b) {
              return rna::rna_rival_queue_before(hypotheses[a],
                                                 hypotheses[b]);
            });
  const std::size_t realize_budget =
      static_cast<std::size_t>(std::max(0, lifecycle_cfg.realize_max));
  if (queue.size() > realize_budget)
    queue.resize(realize_budget);
  // The chain-stage co-primaries: mm_set_parent over the incumbent and every chained
  // rival with a reportable anchor union, shadows left out, in the queue's order. Each
  // primary the queue does not hold is realized after it, outside its budget, and
  // classed a co-primary: it takes no part in the election, the recalibration, the
  // dp2 owner, the MAPQ's prices or the runner-ups.
  {
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < hypotheses.size(); ++i) {
      const rna::RnaRealizedHypothesis& h = hypotheses[i];
      if (i == 0 || (h.chained && h.bundle != nullptr && h.coverage_ok &&
                     !h.shadow))
        order.push_back(i);
    }
    std::stable_sort(order.begin(), order.end(),
                     [&hypotheses](std::size_t a, std::size_t b) {
                       return rna::rna_rival_queue_before(hypotheses[a],
                                                          hypotheses[b]);
                     });
    parent_chain_hits.reserve(order.size());
    for (const std::size_t i : order) {
      const rna::RnaRealizedHypothesis& h = hypotheses[i];
      rna::RnaParentHit hit;
      const rna::RnaQuerySpan span =
          rna::rna_forward_query_span(h.q_begin, h.q_end, h.reverse, read_len);
      hit.qs = span.begin;
      hit.qe = span.end;
      hit.score = h.chain_score;
      hit.cnt = h.chain_anchors;
      hit.rid = h.reference_id;
      hit.rs = h.summary->first_reference_begin;
      hit.re = h.summary->last_reference_end;
      hit.owner = static_cast<int>(i);
      hit.catalogue = h.catalogue_index;
      parent_chain_hits.push_back(hit);
    }
    const std::vector<rna::RnaParentHit> regions =
        parent_regions(lc, parent_chain_hits, nullptr);
    for (std::size_t j = 0; j < regions.size(); ++j) {
      const std::size_t i = static_cast<std::size_t>(regions[j].owner);
      if (regions[j].parent != static_cast<int>(j) || i == 0 ||
          std::find(queue.begin(), queue.end(), i) != queue.end())
        continue;
      hypotheses[i].klass = rna::RnaRivalClass::CoPrimary;
      queue.push_back(i);
    }
  }
  if (!queue.empty()) {
  }
  return queue;
}

// Rank recalibration, as minimap2's mm_update_dp_max. An intronless retrocopy of a
// spliced gene can win the raw-DP election, because the spliced hit pays a long-gap
// open at every junction. minimap2 reprices every aligned hit on an intron-free,
// divergence-amplified scale and elects on that alone; here the election combines
// the raw and recalibrated prices (rna_elect_primary_rank_recal), since the
// recalibrated scale alone is blind to the aligned-mass differences the raw one
// carries. dp_maximum stays the raw segment-0 price, the currency the MAPQ is
// calibrated on, so when a flip leaves the demoted incumbent's raw dp2 above the
// winner's dp1, branch C sees rho > 1 and lowers the MAPQ. Everything here fails
// open: a gate that does not fire, a missing divergence estimate or missing
// accounting leaves the read on raw prices.
void rank_recalibrate(const ReadFacts& facts, Lifecycle& lc,
                      const RnaSpliceRealizationResult& realized) {
  auto& hypotheses = lc.hypotheses;
  const auto& read_len = facts.read_len;
  const auto& rctx = facts.rctx;
  auto& rank_recal_fired = lc.rank_recal_fired;
  // One hypothesis's primary-segment accounting. The incumbent's realization is in
  // `realized` (element 0's `result` is empty by the borrowed-storage contract). The
  // numeric CIGAR is the winning transcript hypothesis's segment_cigars entry, found
  // by orientation, rather than a re-parse of the clip-padded public string.
  const auto recal_accounting =
      [&](std::size_t index,
          rna::RnaRankRecalAccounting& accounting) -> bool {
    const rna::RnaSpliceRealizationResult& result =
        index == 0 ? realized : hypotheses[index].result;
    if (result.refused || result.segments.empty())
      return false;
    const rna::RnaSpliceHypothesisResult* winner = nullptr;
    for (const rna::RnaSpliceHypothesisResult& candidate :
         result.hypotheses)
      if (candidate.orientation == result.winning_hypothesis) {
        winner = &candidate;
        break;
      }
    if (winner == nullptr ||
        winner->segment_cigars.size() != winner->segments.size() ||
        winner->segment_cigars.empty())
      return false;
    const AlignResult& primary = result.segments.front();
    if (!primary.alignment_accounting_valid)
      return false;
    accounting.matches = primary.matches;
    accounting.mismatches = primary.mismatches;
    accounting.ambiguities = primary.ambiguities;
    // Borrowed, never copied: the vector lives in the realization result,
    // which outlives this block.
    accounting.cigar = winner->segment_cigars.front().data();
    accounting.cigar_len = winner->segment_cigars.front().size();
    return accounting.cigar_len != 0;
  };

  // The dp leader and runner-up, as the election orders them: argmax on the realized
  // maximum, ties to the lower catalogue index. Only the realized competitors count,
  // not the co-primaries, so this is minimap2's aligned-hit set.
  const auto off_scale = [](const rna::RnaRealizedHypothesis& h) {
    return !h.realized || !h.dp_maximum.has_value() || *h.dp_maximum <= 0 ||
           h.klass == rna::RnaRivalClass::CoPrimary;
  };
  std::size_t leader = hypotheses.size();
  for (std::size_t i = 0; i < hypotheses.size(); ++i) {
    const rna::RnaRealizedHypothesis& h = hypotheses[i];
    if (off_scale(h))
      continue;
    if (leader == hypotheses.size() ||
        *h.dp_maximum > *hypotheses[leader].dp_maximum ||
        (*h.dp_maximum == *hypotheses[leader].dp_maximum &&
         h.catalogue_index < hypotheses[leader].catalogue_index))
      leader = i;
  }
  int second_dp = 0;
  if (leader != hypotheses.size())
    for (std::size_t i = 0; i < hypotheses.size(); ++i) {
      if (i == leader)
        continue;
      const rna::RnaRealizedHypothesis& h = hypotheses[i];
      if (off_scale(h))
        continue;
      second_dp = std::max(second_dp, *h.dp_maximum);
    }
  // Fewer than two aligned hits is minimap2's max2 == 0: nothing to recalibrate
  // against.
  if (leader != hypotheses.size() && second_dp > 0) {
    const rna::RnaRealizedHypothesis& top = hypotheses[leader];
    rna::RnaRankRecalAccounting top_accounting;
    double identity = -1.0;
    if (rna::rna_rank_recal_gate(read_len, top.aq_begin, top.aq_end,
                                 *top.dp_maximum, second_dp) &&
        recal_accounting(leader, top_accounting))
      identity = rna::rna_event_identity(top_accounting);
    // A divergence the top hit cannot supply is no scale to reprice on.
    if (identity >= 0.0) {
      const int match_sc = rctx.opts.splice_controller.match;
      const double b2 = rna::rna_rank_recal_b2(
          identity, match_sc, rctx.opts.splice_controller.mismatch);
      // Repricing is atomic per read: a hypothesis left on its raw price beside
      // repriced rivals would be compared across currencies.
      std::vector<std::pair<std::size_t, int>> repriced;
      repriced.reserve(hypotheses.size());
      bool recal_complete = true;
      for (std::size_t i = 0; i < hypotheses.size(); ++i) {
        const rna::RnaRealizedHypothesis& h = hypotheses[i];
        if (off_scale(h))
          continue;
        rna::RnaRankRecalAccounting accounting;
        if (!recal_accounting(i, accounting)) {
          recal_complete = false; // fail open for the whole read
          break;
        }
        // A negative price floors to 0, as in minimap2.
        const int recal = rna::rna_recal_max_dp(accounting, b2, match_sc);
        repriced.emplace_back(i, recal < 0 ? 0 : recal);
      }
      if (recal_complete) {
        // dp_maximum keeps the raw price; the recalibrated one lands beside it.
        for (const auto& [index, price] : repriced)
          hypotheses[index].dp_recal = price;
        rank_recal_fired = true;
      }
    }
  }
}

// The lifecycle's stage (5) up to the adoption: the election, the shadow rule against
// the elected hypothesis, the dp2 owner, and the pieces after the election.
void elect_and_decide_pieces(const ReadFacts& facts, Lifecycle& lc,
                             const RnaSpliceRealizationResult& realized) {
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& rank_recal_fired = lc.rank_recal_fired;
  auto& hypotheses = lc.hypotheses;
  const auto& max_intron = facts.max_intron;
  auto& lifecycle_dp2_owner = lc.lifecycle_dp2_owner;
  const auto& read_len = facts.read_len;
  const auto& rctx = facts.rctx;
  auto& pieces_chain = lc.pieces_chain;
  auto& parent_chain_hits = lc.parent_chain_hits;
  auto& pieces_dp_input = lc.pieces_dp_input;
  auto& pieces_left_out = lc.pieces_left_out;
  // (5) Election and adoption. The raw election fails closed (a tie or a missing
  // maximum keeps the incumbent); a repriced read runs the combined rule with a raw
  // veto band (rna.h).
  lifecycle_elected = rank_recal_fired
                          ? rna::rna_elect_primary_rank_recal(hypotheses)
                          : rna::rna_elect_primary_hypothesis(hypotheses);
  // The shadow rule again, against the elected hypothesis: the winner's shadows may
  // include the demoted incumbent or a tie rival. The election is not touched.
  rna::RnaSegmentGeometry winner_geometry;
  if (shadow_geometry(facts, hypotheses[lifecycle_elected], winner_geometry)) {
    for (std::size_t i = 0; i < hypotheses.size(); ++i) {
      if (i == lifecycle_elected)
        continue;
      rna::RnaSegmentGeometry other_geometry;
      hypotheses[i].shadow =
          shadow_geometry(facts, hypotheses[i], other_geometry) &&
          rna::rna_same_place(winner_geometry, other_geometry, max_intron);
    }
  }
  lifecycle_dp2_owner =
      rna::rna_select_dp2_owner(hypotheses, lifecycle_elected);
  // The pieces, after the election: every realized region, repriced by
  // mm_update_dp_max and in mm_hit_sort's order with the elected family's regions
  // first (they own the query they cover), through mm_set_parent less the elected
  // family's shadows. A region outside the elected family stays only if it passes
  // minimap2's region filter, before and after the repricing (mm_filter_regs in
  // mm_align_skeleton). The election and the MAPQ do not read it.
  {
    const int elected_owner = static_cast<int>(lifecycle_elected);
    // Every segment of every realized hypothesis; a refused or empty realization has
    // no region.
    std::vector<rna::RnaParentHit> dp_hits;
    for (std::size_t i = 0; i < hypotheses.size(); ++i)
      if (i == 0 || hypotheses[i].realized)
        append_dp_hits(lc, i, i == 0 ? realized : hypotheses[i].result, dp_hits);
    const auto filter_regions = [&]() {
      const auto passes = [&](const rna::RnaParentHit& hit) {
        if (hit.owner == elected_owner)
          return true;
        const std::size_t i = static_cast<std::size_t>(hit.owner);
        const rna::RnaSpliceHypothesisResult* winner =
            rna::rna_winning_hypothesis(i == 0 ? realized : hypotheses[i].result);
        rna::SpliceRegionEvidence evidence =
            winner->region_evidence[static_cast<std::size_t>(hit.segment)];
        evidence.dp_maximum = hit.dp_max;
        return rna::classify_splice_region(evidence, read_len,
                                           rctx.opts.splice_controller) ==
               rna::SpliceRegionReject::None;
      };
      dp_hits.erase(
          std::stable_partition(dp_hits.begin(), dp_hits.end(), passes),
          dp_hits.end());
    };
    filter_regions();
    rna::rna_parent_update_dp_max(
        read_len, rna::kRnaRankRecalMinReadLen, dp_hits,
        static_cast<float>(rna::kRnaRankRecalFrac),
        rctx.opts.splice_controller.match,
        rctx.opts.splice_controller.mismatch);
    filter_regions();
    rna::rna_parent_hit_sort(dp_hits);
    std::stable_partition(dp_hits.begin(), dp_hits.end(),
                          [elected_owner](const rna::RnaParentHit& hit) {
                            return hit.owner == elected_owner;
                          });
    std::vector<char> left_out(hypotheses.size(), 0);
    for (std::size_t i = 0; i < hypotheses.size(); ++i)
      left_out[i] = i != lifecycle_elected && hypotheses[i].shadow ? 1 : 0;
    pieces_chain = parent_regions(lc, parent_chain_hits, &left_out);
    decide_pieces(facts, lc, parent_regions(lc, dp_hits, &left_out), pieces_chain);
    pieces_dp_input = std::move(dp_hits);
    pieces_left_out = std::move(left_out);
  }
}

} // namespace fa::cpu::lr::rna

/**** committed_evidence ****/
// The CIGAR lane's chain MAPQ evidence for the committed record, after the lifecycle.

namespace fa::cpu::lr::rna {

// The committed record's chain MAPQ evidence in the committed frame: its chain, vote,
// coverage, identity and realized DP maxima, and every other hypothesis as a rival.
void committed_chain_evidence(const ReadFacts& facts, const ReadState& rs,
                              const Lifecycle& lc, const Committed& cf,
                              const std::optional<int>& committed_dp_max,
                              const std::optional<int>& rival_dp_max,
                              const RnaSpliceRealizationResult& realized,
                              RnaChainMapqEvidence& chain_evidence) {
  const auto& committed_path = cf.committed_path;
  const auto& committed_siblings = cf.committed_siblings;
  const auto& committed_path_summary = cf.committed_path_summary;
  const auto& committed_locus = cf.committed_locus;
  const auto& committed_rank = cf.committed_rank;
  const auto& catalogue_rank1_score = rs.catalogue_rank1_score;
  const auto& catalogue_rank2_score = rs.catalogue_rank2_score;
  const auto& read_len = facts.read_len;
  const auto& committed_reverse = cf.committed_reverse;
  const auto& commit_qcov = cf.commit_qcov;
  const auto& rctx = facts.rctx;
  const auto& hypotheses = lc.hypotheses;
  const auto& lifecycle_elected = lc.lifecycle_elected;
  const auto& lifecycle_dp2_owner = lc.lifecycle_dp2_owner;
  chain_evidence.f1 = committed_path->anchor_path_score;
  chain_evidence.cnt =
      static_cast<int>(committed_path->anchor_path_anchor_count);
  chain_evidence.sib_f2 = committed_siblings->sib_score;
  chain_evidence.winner_q_begin = committed_path_summary->first_query_begin;
  chain_evidence.winner_q_end = committed_path_summary->last_query_end;
  chain_evidence.rank1_score = committed_locus->rank_score;
  // The largest vote at any locus but the committed one: rank 2's without adoption,
  // rank 1's on any adoption (the catalogue is ordered).
  chain_evidence.rank2_score =
      committed_rank >= 2 ? catalogue_rank1_score : catalogue_rank2_score;
  chain_evidence.read_len = read_len;
  // As for map-only, in the committed frame: an adopted rival carries its own
  // post-repair strand.
  chain_evidence.winner_reverse = committed_reverse;
  chain_evidence.selected_query_covered_bases = commit_qcov;
  // Every catalogue rival is chained here, so a vote-only rival contributes x strength
  // 0 (see unchained_zero).
  chain_evidence.unchained_zero = true;
  if (committed_dp_max && *committed_dp_max > 0)
    chain_evidence.dp1 = static_cast<double>(*committed_dp_max);
  if (chain_evidence.dp1 > 0.0 && rival_dp_max && *rival_dp_max > 0)
    chain_evidence.dp2 = static_cast<double>(*rival_dp_max);
  if (!realized.segments.empty()) {
    const AlignResult& primary = realized.segments.front();
    chain_evidence.identity =
        primary.alignment_accounting_valid && primary.block_len > 0
            ? static_cast<double>(primary.matches) /
                  static_cast<double>(primary.block_len)
            : 1.0;
  }
  chain_evidence.match_sc = rctx.opts.splice_controller.match;
  // minimap2 map.c: sub_diff = 2a + b, per preset (4 splice / 6 splice:hq),
  // following any -A/-B override.
  chain_evidence.sub_diff = 2 * rctx.opts.splice_controller.match +
                            rctx.opts.splice_controller.mismatch;
  // The preset MAPQ calibration values ride the evidence.
  chain_evidence.qcov_tau = rctx.opts.mapq_qcov_tau;
  chain_evidence.disjoint_vote_damp = rctx.opts.mapq_disjoint_vote_damp;
  // The rivals are assembled from the hypothesis vector rather than chained again, in
  // catalogue order, with a demoted incumbent among them. owns_dp2 lands on the
  // lifecycle's dp2 owner, never on a co-primary.
  for (std::size_t i = 0; i < hypotheses.size(); ++i) {
    if (i == lifecycle_elected)
      continue;
    const rna::RnaRealizedHypothesis& other = hypotheses[i];
    rna::RnaChainMapqRival rival;
    rival.candidate = other.catalogue_index;
    rival.chained = other.chained;
    rival.chain_score = other.chain_score;
    rival.chain_anchors = other.chain_anchors;
    rival.q_begin = other.q_begin;
    rival.q_end = other.q_end;
    rival.rank_score = other.rank_score;
    rival.reverse = other.reverse;
    // The shadow verdict, re-read against the elected winner.
    rival.shadow = other.shadow;
    rival.owns_dp2 = chain_evidence.dp2 > 0.0 &&
                     lifecycle_dp2_owner.has_value() &&
                     *lifecycle_dp2_owner == i;
    // A realized rival is already bounded by rho = dp2/dp1 and stays out of x_floor;
    // a chained-only or censored rival enters it, and so does a co-primary, realized
    // for the pieces rather than as a price.
    rival.realized_dp = other.realized && other.dp_maximum.has_value() &&
                        *other.dp_maximum > 0 &&
                        other.klass != rna::RnaRivalClass::CoPrimary;
    chain_evidence.rivals.push_back(rival);
  }
}

} // namespace fa::cpu::lr::rna

/**** pieces ****/
// The supplementary pieces of minimap2's parent rule as records, after the leftover.

namespace fa::cpu::lr::rna {

// The pieces' records: each piece's region segment priced by the chain MAPQ on its own
// region against its children and stored as one flat supplementary family; the
// places the pieces print go on the wire.
void price_pieces(const ReadFacts& facts, Lifecycle& lc,
                  const RnaSpliceRealizationResult& realized,
                  const RnaChainMapqEvidence& chain_evidence,
                  std::vector<RnaSpliceRealizationResult>& chimeric_results,
                  std::vector<RnaChimericFamily>& chimeric_families) {
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& demoted_incumbent = lc.demoted_incumbent;
  auto& hypotheses = lc.hypotheses;
  auto& pieces = lc.pieces;
  const auto& read_len = facts.read_len;
  auto& locus_on_wire = lc.locus_on_wire;
  const auto realization_of =
      [&](std::size_t i) -> const rna::RnaSpliceRealizationResult& {
    if (i == lifecycle_elected)
      return realized;
    if (i == 0)
      return demoted_incumbent;
    return hypotheses[i].result;
  };
  // Reserved once: `chimeric_families` borrows into `chimeric_results`.
  chimeric_results.reserve(pieces.size());
  chimeric_families.reserve(pieces.size());
  for (const PieceRegion& piece : pieces) {
    const rna::RnaParentHit& region = piece.region;
    const std::size_t owner = static_cast<std::size_t>(region.owner);
    const rna::RnaSpliceRealizationResult& result = realization_of(owner);
    const std::size_t s = static_cast<std::size_t>(region.segment);
    if (result.refused || s >= result.segments.size())
      continue;
    const AlignResult& segment = result.segments[s];
    const rna::RnaRealizedHypothesis& candidate = hypotheses[owner];
    rna::RnaChainMapqEvidence piece_evidence;
    piece_evidence.f1 = candidate.chain_score;
    piece_evidence.cnt = candidate.chain_anchors;
    piece_evidence.sib_f2 = candidate.siblings.sib_score;
    piece_evidence.winner_q_begin = candidate.q_begin;
    piece_evidence.winner_q_end = candidate.q_end;
    piece_evidence.rank1_score = candidate.rank_score;
    // The demoted vote: the largest vote among its children.
    for (const std::size_t c : piece.children)
      piece_evidence.rank2_score =
          std::max(piece_evidence.rank2_score, hypotheses[c].rank_score);
    piece_evidence.read_len = read_len;
    piece_evidence.winner_reverse = candidate.reverse;
    piece_evidence.selected_query_covered_bases =
        candidate.query_covered_bases;
    piece_evidence.qcov_bases = region.qe - region.qs;
    piece_evidence.unchained_zero = true;
    if (region.dp_max > 0)
      piece_evidence.dp1 = static_cast<double>(region.dp_max);
    if (piece_evidence.dp1 > 0.0 && region.dp_max2 > 0)
      piece_evidence.dp2 = static_cast<double>(region.dp_max2);
    piece_evidence.identity =
        segment.alignment_accounting_valid && segment.block_len > 0
            ? static_cast<double>(segment.matches) /
                  static_cast<double>(segment.block_len)
            : 1.0;
    piece_evidence.match_sc = chain_evidence.match_sc;
    piece_evidence.sub_diff = chain_evidence.sub_diff;
    piece_evidence.qcov_tau = chain_evidence.qcov_tau;
    piece_evidence.disjoint_vote_damp = chain_evidence.disjoint_vote_damp;
    piece_evidence.rivals.reserve(piece.children.size());
    for (const std::size_t c : piece.children) {
      const rna::RnaRealizedHypothesis& other = hypotheses[c];
      rna::RnaChainMapqRival rival;
      rival.candidate = other.catalogue_index;
      rival.chained = other.chained;
      rival.chain_score = other.chain_score;
      rival.chain_anchors = other.chain_anchors;
      rival.q_begin = other.q_begin;
      rival.q_end = other.q_end;
      rival.rank_score = other.rank_score;
      rival.reverse = other.reverse;
      rival.owns_dp2 = piece_evidence.dp2 > 0.0 && piece.dp2_owner == c;
      rival.realized_dp = other.realized && other.dp_maximum.has_value() &&
                          *other.dp_maximum > 0;
      rival.parented = true;
      piece_evidence.rivals.push_back(rival);
    }
    piece_evidence.parent_rule = true;
    piece_evidence.parent_n_sub = region.n_sub;
    rna::RnaChainMapqBreakdown piece_breakdown;
    const int piece_mapq =
        rna::rna_chain_mapq(piece_evidence, &piece_breakdown);
    rna::RnaSpliceRealizationResult piece_result;
    piece_result.segments.push_back(segment);
    chimeric_results.push_back(std::move(piece_result));
    rna::RnaChimericFamily family;
    family.result = &chimeric_results.back();
    family.mapq = piece_mapq;
    if (piece_evidence.f1 > 0) {
      family.chain_tags.chain_anchors = piece_evidence.cnt;
      family.chain_tags.chain_score = piece_evidence.f1;
      family.chain_tags.secondary_chain_score = piece_breakdown.f2;
    }
    chimeric_families.push_back(family);
    for (const std::size_t i : piece.place)
      locus_on_wire[i] = 1;
  }
}

} // namespace fa::cpu::lr::rna

/**** runner_ups ****/
// The runner-ups beside the committed record, after the commit.

namespace fa::cpu::lr::rna {

// The runner-ups: one realized family per surviving competitor of the committed record,
// the best priced first and one per place, moved out of the hypotheses.
void attach_runner_ups(const ReadFacts& facts, ReadState& rs, Lifecycle& lc) {
  auto& locus_on_wire = lc.locus_on_wire;
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& hypotheses = lc.hypotheses;
  const auto& read_len = facts.read_len;
  auto& demoted_incumbent = lc.demoted_incumbent;
  auto& out = rs.out;
  const auto& max_intron = facts.max_intron;
  // One survivor: what to move, where it sits, and the keys the printed
  // order uses.
  struct RunnerUp {
    rna::RnaSpliceRealizationResult* result;
    rna::RnaSegmentGeometry place;
    int dp_maximum;
    int catalogue_index;
    // cm:i and s1:i from the hypothesis's own chain, as a chimeric family's.
    rna::RnaChainTags chain_tags;
  };
  std::vector<RunnerUp> runner_ups;
  const auto on_wire = [&](std::size_t j) {
    return j < locus_on_wire.size() && locus_on_wire[j] != 0;
  };
  // The realized place, in the geometry every duplicate test here
  // reads: the segments' reference hull on their own strand, with the
  // realized forward-frame query hull the lifecycle already computed.
  const auto collect = [&](rna::RnaSpliceRealizationResult& result,
                           const rna::RnaRealizedHypothesis& h) {
    if (result.refused || result.segments.empty() ||
        !h.dp_maximum.has_value() || *h.dp_maximum <= 0)
      return;
    rna::RnaSegmentGeometry place;
    place.reference_id = h.reference_id;
    place.reverse = result.segments.front().is_reverse;
    place.forward_span = rna::RnaQuerySpan{h.aq_begin, h.aq_end};
    for (std::size_t s = 0; s < result.segments.size(); ++s) {
      const AlignResult& segment = result.segments[s];
      const std::uint64_t begin = static_cast<std::uint64_t>(
          std::max<std::int64_t>(0, static_cast<std::int64_t>(segment.pos)));
      const std::uint64_t end =
          static_cast<std::uint64_t>(std::max<std::int64_t>(
              0, static_cast<std::int64_t>(segment.target_end)));
      place.reference_begin =
          s == 0 ? begin : std::min(place.reference_begin, begin);
      place.reference_end = s == 0 ? end : std::max(place.reference_end, end);
    }
    rna::RnaChainTags chain_tags;
    if (h.chain_score > 0) {
      chain_tags.chain_anchors = h.chain_anchors;
      chain_tags.chain_score = h.chain_score;
    }
    runner_ups.push_back(
        {&result, place, *h.dp_maximum, h.catalogue_index, chain_tags});
  };
  // Whether hypothesis i is an alternative placement of the bases the committed
  // record explains: not a shadow, and competing with the committed chain's span at
  // mask_level in the forward frame. `klass` was read against the incumbent, so after
  // an adoption it cannot answer this; rna_chain_mapq_admissible compares raw oriented
  // spans and would drop an opposite-strand alternative. With no lifecycle the span
  // stays censored and the vector is empty.
  rna::RnaQuerySpan committed_span;
  if (lifecycle_elected < hypotheses.size()) {
    const rna::RnaRealizedHypothesis& committed =
        hypotheses[lifecycle_elected];
    committed_span = rna::rna_forward_query_span(
        committed.q_begin, committed.q_end, committed.reverse, read_len);
  }
  const auto weighed = [&](std::size_t i) {
    const rna::RnaRealizedHypothesis& h = hypotheses[i];
    return !h.shadow &&
           rna::rna_query_spans_compete(
               committed_span, rna::rna_forward_query_span(
                                   h.q_begin, h.q_end, h.reverse, read_len));
  };
  for (std::size_t i = 1; i < hypotheses.size(); ++i) {
    rna::RnaRealizedHypothesis& h = hypotheses[i];
    // Rivals of the committed record only: not the elected hypothesis (it is the
    // record), not a co-primary, not a shadow.
    if (i == lifecycle_elected || h.klass != rna::RnaRivalClass::Competitor ||
        !h.realized || on_wire(i) || !weighed(i))
      continue;
    collect(h.result, h);
  }
  // The demoted incumbent, from its own local; it was among the MAPQ's rivals, so the
  // same test applies.
  if (lifecycle_elected != 0 && !on_wire(0) && weighed(0))
    collect(demoted_incumbent, hypotheses[0]);
  if (!runner_ups.empty()) {
    // Richest realized maximum first, ties to the lower catalogue index
    // (rna_select_dp2_owner's rule), so the first alternative is the dp2 the MAPQ was
    // priced against.
    std::stable_sort(runner_ups.begin(), runner_ups.end(),
                     [](const RunnerUp& a, const RunnerUp& b) {
                       if (a.dp_maximum != b.dp_maximum)
                         return a.dp_maximum > b.dp_maximum;
                       return a.catalogue_index < b.catalogue_index;
                     });
    // Two catalogue windows can re-find one locus; the sort makes the kept copy the
    // better priced one.
    out.secondary.reserve(runner_ups.size());
    std::vector<std::size_t> kept;
    kept.reserve(runner_ups.size());
    for (std::size_t i = 0; i < runner_ups.size(); ++i) {
      bool duplicate = false;
      for (const std::size_t j : kept) {
        if (rna::rna_same_place(runner_ups[j].place, runner_ups[i].place,
                                max_intron)) {
          duplicate = true;
          break;
        }
      }
      if (duplicate)
        continue;
      kept.push_back(i);
      rna::attach_runner_up_family(*runner_ups[i].result,
                                   runner_ups[i].chain_tags, out);
    }
  }
}

} // namespace fa::cpu::lr::rna
