// Peak-to-peak transition scoring for the coarse colinear recurrence
// (best_colinear_score); query_partition.cpp's staircases chain peaks with the same rule.
#pragma once

#include "coarse_chain.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

namespace detail {

inline uint32_t ceil_log2_u64(uint64_t value) {
  uint32_t r = 0;
  uint64_t v = 1;
  while (v < value) {
    v <<= 1;
    ++r;
  }
  return r;
}

enum class TransitionKind : uint8_t { kSameExon, kSplice, kInvalid };

struct Transition {
  TransitionKind kind = TransitionKind::kInvalid;
  int64_t penalty = 0;
};

// Asymmetric transition score left -> right. The caller guarantees same contig and strand
// and reference-monotone order, so this scores geometry only.
inline Transition score_transition(const CoarseDiagonalPeak& left,
                                   const CoarseDiagonalPeak& right,
                                   const CoarseLocusOptions& o) {
  Transition t;
  if (right.reference_begin <= left.reference_begin)
    return t; // non-monotone
  const int64_t dq = static_cast<int64_t>(right.oriented_query_begin) -
                     static_cast<int64_t>(left.oriented_query_end);
  const int64_t dr = static_cast<int64_t>(right.reference_begin) -
                     static_cast<int64_t>(left.reference_end);
  if (dr <= 0)
    return t;
  if (dq < -static_cast<int64_t>(o.max_query_overlap) ||
      dq > static_cast<int64_t>(o.max_query_gap)) {
    return t;
  }
  const int64_t overlap = std::max<int64_t>(0, -dq);
  const int64_t query_gap = std::max<int64_t>(0, dq);
  const int64_t diagonal_delta = std::llabs(dr - dq);
  if (diagonal_delta <= static_cast<int64_t>(o.diagonal_band)) {
    t.kind = TransitionKind::kSameExon;
    t.penalty = static_cast<int64_t>(o.diagonal_drift_scale) * diagonal_delta +
                static_cast<int64_t>(o.overlap_scale) * overlap;
    return t;
  }
  const int64_t intron = dr - dq;
  if (intron < static_cast<int64_t>(o.min_intron) ||
      intron > static_cast<int64_t>(o.max_intron)) {
    return t; // kInvalid: an over-long / impossible jump, never free-chained
  }
  // Sub-linear in intron length: a zero-cost ref jump over-chains noise, a
  // linear penalty erases real long introns.
  t.kind = TransitionKind::kSplice;
  t.penalty = static_cast<int64_t>(o.splice_open) +
              static_cast<int64_t>(o.splice_log_scale) *
                  static_cast<int64_t>(
                      ceil_log2_u64(static_cast<uint64_t>(intron) + 1)) +
              static_cast<int64_t>(o.query_gap_scale) * query_gap +
              static_cast<int64_t>(o.overlap_scale) * overlap;
  return t;
}

inline int64_t incremental_query_reward(const CoarseDiagonalPeak& left,
                                        const CoarseDiagonalPeak& right) {
  const uint32_t start =
      std::max(left.oriented_query_end, right.oriented_query_begin);
  return right.oriented_query_end > start
             ? static_cast<int64_t>(right.oriented_query_end - start)
             : 0;
}

} // namespace detail

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
