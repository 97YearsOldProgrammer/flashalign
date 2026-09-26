// Rank-1 search regions. Coordinates are zero-based, half-open, chromosome-local and in
// the selected strand's query frame.
#pragma once

#include "refusal.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class SearchRegionKind : uint8_t { kCore = 0 };

struct SkeletonNode {
  uint32_t node_id = 0;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int64_t reference_begin = 0;
  int64_t reference_end = 0;
  int64_t median_diagonal = 0;
  uint32_t support = 0;
};

struct SearchRegion {
  uint32_t region_id = 0;
  SearchRegionKind kind = SearchRegionKind::kCore;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int64_t reference_begin = 0;
  int64_t reference_end = 0;
  int64_t diagonal_low = 0;
  int64_t diagonal_high = 0;  // inclusive
  uint32_t per_region_budget = 0;
};

struct SkeletonRegionParams {
  // Sanity limit on one gene's genomic footprint; posting enumeration is restricted to
  // this interval.
  int64_t maximum_locus_reference_span = 4000000;
  uint32_t locus_budget = 262144;
};

struct SearchRegionBuildResult {
  std::vector<SearchRegion> regions;
  uint64_t total_volume = 0;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

// One region covering the whole query against the selected locus envelope, with no
// diagonal restriction.
SearchRegionBuildResult derive_locus_search_region(
    int32_t query_length, int64_t reference_begin, int64_t reference_end,
    int64_t reference_length, const SkeletonRegionParams& params);

bool search_region_overlaps_anchor(const SearchRegion& region,
                                   int32_t query_begin,
                                   int32_t reference_begin, int32_t span);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
