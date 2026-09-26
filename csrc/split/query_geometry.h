// Query-overlap test that decides whether two placements are rivals for the same read span or
// separate pieces of the read: minimap2's mask_level test (mm_set_parent).
#pragma once

#include "coordinates.h"

#include <algorithm>
#include <cstdint>

namespace fa {
namespace cpu {
namespace split {

// Relationship between two query intervals for split-mapping selection.
enum class QueryRelation : uint8_t {
  DisjointSegment, // separately alignable piece of the read (supplementary)
  Rival            // alternative placement of the same span (secondary/f2)
};

// Overlap fraction of the shorter interval at or above which the pair are rivals. The
// threshold is rounded down.
struct QueryMaskPolicy {
  uint32_t numerator = 1;
  uint32_t denominator = 2;
};

inline constexpr QueryMaskPolicy kHalfFloor{1, 2};

// Rival when the overlap reaches the policy fraction of the shorter interval; no overlap or
// an empty interval is always DisjointSegment.
inline QueryRelation classify_query_relation(QueryInterval a, QueryInterval b,
                                             QueryMaskPolicy p = {}) noexcept {
  const int overlap =
      std::min(a.end, b.end) - std::max(a.begin, b.begin);
  if (overlap <= 0)
    return QueryRelation::DisjointSegment;
  const int shorter = std::min(a.length(), b.length());
  if (shorter <= 0)
    return QueryRelation::DisjointSegment;
  const int64_t scaled =
      static_cast<int64_t>(shorter) * static_cast<int64_t>(p.numerator);
  const int64_t den = static_cast<int64_t>(p.denominator);
  const int threshold = static_cast<int>(scaled / den);
  return overlap >= threshold ? QueryRelation::Rival
                              : QueryRelation::DisjointSegment;
}

} // namespace split
} // namespace cpu
} // namespace fa
