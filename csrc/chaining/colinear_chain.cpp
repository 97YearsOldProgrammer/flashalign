// Colinear chaining DP, traceback and best-chain choice, as minimap2's mg_lchain_dp. DNA and
// spliced chaining share one compile-time driver.
#include "colinear_chain.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace chaining {
namespace {

constexpr int32_t kNegInf = INT32_MIN;
// DNA anchors here are sparse whole-read representatives: past this step the linear gap
// penalty scales down with the step, charging diagonal drift as a rate over the spacing.
constexpr int32_t kDnaLinearPenaltyPitch = 256;

float mg_log2(float x) {
  union {
    float f;
    uint32_t i;
  } z = {x};
  float log_2 = static_cast<float>(((z.i >> 23) & 255)) - 128.0f;
  z.i &= ~(255u << 23);
  z.i += 127u << 23;
  log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
  return log_2;
}

// The scalars the recurrence reads, passed by value: a store into `marks` could alias an
// int32_t field of the caller's params and force reloads in the inner loop.
struct ScoreParams {
  int32_t max_dist_x;
  int32_t max_dist_y;
  int32_t bw;
  float chn_pen_gap;
  float chn_pen_skip;
};

template <class Params>
ScoreParams score_params_of(const Params& params) {
  return ScoreParams{params.max_dist_x, params.max_dist_y, params.bw,
                     params.chn_pen_gap, params.chn_pen_skip};
}

template <bool Spliced, bool ScaleSparseLongSteps>
int32_t compute_score(const Anchor& ai, const Anchor& aj, ScoreParams params) {
  const int32_t dq = ai.q - aj.q;
  if (dq <= 0 || dq > params.max_dist_x) return kNegInf;
  const int32_t dr = ai.r - aj.r;
  if (dr == 0 || dq > params.max_dist_y) return kNegInf;
  const int32_t dd = dr > dq ? dr - dq : dq - dr;
  if (dd > params.bw) return kNegInf;
  const int32_t dg = dr < dq ? dr : dq;
  const int32_t q_span = Spliced ? aj.span
                                 : (aj.span > 255 ? 255 : aj.span);
  int32_t score = q_span < dg ? q_span : dg;
  if (dd) {
    float gap_penalty = params.chn_pen_gap;
    if constexpr (!Spliced && ScaleSparseLongSteps) {
      if (dg > kDnaLinearPenaltyPitch) {
        gap_penalty *=
            static_cast<float>(kDnaLinearPenaltyPitch) / static_cast<float>(dg);
      }
    }
    const float linear_penalty = gap_penalty * static_cast<float>(dd) +
                                 params.chn_pen_skip * static_cast<float>(dg);
    const float log_penalty = mg_log2(static_cast<float>(dd + 1));
    if constexpr (Spliced) {
      if (dr > dq) {
        score -= static_cast<int32_t>(
            linear_penalty < log_penalty ? linear_penalty : log_penalty);
      } else {
        score -= static_cast<int32_t>(linear_penalty + 0.5f * log_penalty);
      }
    } else {
      score -= static_cast<int32_t>(linear_penalty + 0.5f * log_penalty);
    }
  } else if (dg > q_span) {
    // With dd == 0 only the skip charge remains; this equals the full expression.
    score -= static_cast<int32_t>(params.chn_pen_skip *
                                  static_cast<float>(dg));
  }
  return score;
}

int64_t backtrack_end(int32_t max_drop, int32_t score_end,
                      const std::vector<int32_t>& scores,
                      const std::vector<int32_t>& predecessors,
                      std::vector<int32_t>& marks, int64_t start) {
  int64_t i = start;
  int64_t end_i = -1;
  int64_t max_i = i;
  int32_t max_score = 0;
  if (i < 0 || marks[i] != 0) return i;
  do {
    marks[i] = 2;
    end_i = i = predecessors[i];
    const int32_t score = i < 0 ? score_end : score_end - scores[i];
    if (score > max_score) {
      max_score = score;
      max_i = i;
    } else if (max_score - score > max_drop) {
      break;
    }
  } while (i >= 0 && marks[i] == 0);
  for (i = start; i >= 0 && i != end_i; i = predecessors[i]) marks[i] = 0;
  return max_i;
}

template <bool Spliced, bool ScaleSparseLongSteps, class Params>
ChainResult chain_driver(std::vector<Anchor> anchors, Params params) {
  ChainResult result;
  const int64_t count = static_cast<int64_t>(anchors.size());
  if (count == 0) return result;

  if (params.max_dist_x < params.bw) params.max_dist_x = params.bw;
  if constexpr (!Spliced) {
    if (params.max_dist_y < params.bw) params.max_dist_y = params.bw;
  }

  // Callers usually pass anchors strictly sorted by (r, q), so sort only when an O(n) check
  // fails. Equal pairs also sort, since an unstable sort may permute them.
  const auto rq_strictly_less = [](const Anchor& left, const Anchor& right) {
    if (left.r != right.r) return left.r < right.r;
    return left.q < right.q;
  };
  const bool rq_strictly_sorted =
      std::adjacent_find(anchors.begin(), anchors.end(),
                         [&rq_strictly_less](const Anchor& left,
                                             const Anchor& right) {
                           return !rq_strictly_less(left, right);
                         }) == anchors.end();
  if (!rq_strictly_sorted) {
    std::sort(anchors.begin(), anchors.end(), rq_strictly_less);
  }
  result.anchors = std::move(anchors);
  const std::vector<Anchor>& a = result.anchors;

  // Hoisted out of the j-loop for the aliasing reason on ScoreParams.
  const ScoreParams score_params = score_params_of(params);
  const int32_t max_dist_x = params.max_dist_x;
  const int32_t max_iter = params.max_iter;
  const int32_t max_skip = params.max_skip;

  std::vector<int32_t> scores(static_cast<size_t>(count));
  std::vector<int32_t> predecessors(static_cast<size_t>(count));
  std::vector<int32_t> marks(static_cast<size_t>(count), 0);
  const int32_t max_drop = Spliced ? INT32_MAX : params.bw;
  // upper[j] = max over j' <= j of scores[j'] + j''s capped span. With non-negative
  // penalties compute_score(ai, aj) never exceeds aj's capped span, so once
  // upper[j] <= best_score no remaining predecessor can strictly improve the score, and
  // stopping there gives the same result as the full scan.
  std::vector<int32_t> upper(static_cast<size_t>(count));
  const bool prune_by_upper_bound =
      params.chn_pen_gap >= 0.0f && params.chn_pen_skip >= 0.0f;

  int64_t start = 0;
  int64_t max_index = -1;
  for (int64_t i = 0; i < count; ++i) {
    const Anchor ai = a[i];
    int32_t best_predecessor = -1;
    int64_t end_j;
    const int32_t span_i =
        Spliced ? ai.span : (ai.span > 255 ? 255 : ai.span);
    int32_t best_score = span_i;
    int32_t skipped = 0;
    while (start < i && ai.r > a[start].r + max_dist_x) {
      ++start;
    }
    if (i - start > max_iter) start = i - max_iter;
    int64_t j;
    for (j = i - 1; j >= start; --j) {
      if (prune_by_upper_bound && upper[j] <= best_score) break;
      int32_t score = compute_score<Spliced, ScaleSparseLongSteps>(
          ai, a[j], score_params);
      if (score == kNegInf) continue;
      score += scores[j];
      if (score > best_score) {
        best_score = score;
        best_predecessor = static_cast<int32_t>(j);
        if (skipped > 0) --skipped;
      } else if (marks[j] == static_cast<int32_t>(i)) {
        if (++skipped > max_skip) break;
      }
      if (predecessors[j] >= 0)
        marks[predecessors[j]] = static_cast<int32_t>(i);
    }
    end_j = j;
    if (max_index < 0 ||
        ai.r - a[max_index].r > static_cast<int64_t>(max_dist_x)) {
      int32_t maximum = kNegInf;
      max_index = -1;
      for (j = i - 1; j >= start; --j) {
        if (maximum < scores[j]) {
          maximum = scores[j];
          max_index = j;
        }
      }
    }
    if (max_index >= 0 && max_index < end_j) {
      const int32_t extra =
          compute_score<Spliced, ScaleSparseLongSteps>(
              ai, a[max_index], score_params);
      if (extra != kNegInf && best_score < extra + scores[max_index]) {
        best_score = extra + scores[max_index];
        best_predecessor = static_cast<int32_t>(max_index);
      }
    }
    scores[i] = best_score;
    predecessors[i] = best_predecessor;
    const int32_t upper_i = best_score + span_i;
    upper[i] = (i > 0 && upper[i - 1] > upper_i) ? upper[i - 1] : upper_i;
    if (max_index < 0 ||
        (ai.r - a[max_index].r <= static_cast<int64_t>(max_dist_x) &&
         scores[max_index] < scores[i])) {
      max_index = i;
    }
  }

  std::vector<std::pair<int32_t, int64_t>> ranked_ends;
  ranked_ends.reserve(static_cast<size_t>(count));
  for (int64_t i = 0; i < count; ++i) {
    if (scores[i] >= params.min_sc) ranked_ends.emplace_back(scores[i], i);
  }
  if (ranked_ends.empty()) return result;
  // Pairs are appended in ascending index order, so sorting by (score, index) equals a
  // stable sort by score.
  std::sort(ranked_ends.begin(), ranked_ends.end(),
            [](const std::pair<int32_t, int64_t>& left,
               const std::pair<int32_t, int64_t>& right) {
              if (left.first != right.first) return left.first < right.first;
              return left.second < right.second;
            });

  std::fill(marks.begin(), marks.end(), 0);
  for (int64_t k = static_cast<int64_t>(ranked_ends.size()) - 1; k >= 0; --k) {
    const int64_t chain_start = ranked_ends[static_cast<size_t>(k)].second;
    if (marks[chain_start] != 0) continue;
    const int32_t score_end = ranked_ends[static_cast<size_t>(k)].first;
    const int64_t end_i = backtrack_end(
        max_drop, score_end, scores, predecessors, marks, chain_start);
    // Count the chain first so the indices fill back to front in one allocation.
    int64_t length = 0;
    for (int64_t walk = chain_start; walk != end_i; walk = predecessors[walk]) {
      ++length;
    }
    std::vector<int32_t> indices(static_cast<size_t>(length));
    size_t fill = static_cast<size_t>(length);
    int64_t i = chain_start;
    for (; i != end_i; i = predecessors[i]) {
      indices[--fill] = static_cast<int32_t>(i);
      marks[i] = 1;
    }
    const int32_t score = i < 0 ? score_end : score_end - scores[i];
    if (score >= params.min_sc &&
        static_cast<int32_t>(indices.size()) >= params.min_cnt) {
      result.chains.push_back(Chain{std::move(indices), score});
    }
  }
  return result;
}

}  // namespace

ChainResult chain_colinear(std::vector<Anchor> anchors,
                           ColinearChainParams params) {
  return chain_driver<false, true>(std::move(anchors), params);
}

ChainResult chain_spliced_colinear(std::vector<Anchor> anchors,
                                   SplicedChainParams params) {
  return chain_driver<true, false>(std::move(anchors), params);
}

int best_chain_index(const ChainResult& result) {
  int best = -1;
  int32_t best_score = INT32_MIN;
  for (size_t i = 0; i < result.chains.size(); ++i) {
    if (result.chains[i].score > best_score) {
      best_score = result.chains[i].score;
      best = static_cast<int>(i);
    }
  }
  return best;
}

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
