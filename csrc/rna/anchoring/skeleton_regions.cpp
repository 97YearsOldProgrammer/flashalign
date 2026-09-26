#include "skeleton_regions.h"

#include "../../core/checked_range.h"

#include <algorithm>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

int64_t clamp_reference(int64_t value, int64_t reference_length) {
  return std::max<int64_t>(
      0, std::min<int64_t>(value, std::max<int64_t>(0, reference_length)));
}

uint64_t saturating_volume(int64_t query_span, int64_t reference_span) {
  if (query_span <= 0 || reference_span <= 0) return 0;
  const uint64_t q = static_cast<uint64_t>(query_span);
  const uint64_t r = static_cast<uint64_t>(reference_span);
  uint64_t product = 0;
  if (!mapping::checked_u64_product(q, r, product))
    return std::numeric_limits<uint64_t>::max();
  return product;
}

}  // namespace

SearchRegionBuildResult derive_locus_search_region(
    int32_t query_length, int64_t reference_begin, int64_t reference_end,
    int64_t reference_length, const SkeletonRegionParams& params) {
  SearchRegionBuildResult out;
  if (reference_length > mapping::kSignedCoordinateMaximum) {
    out.refused = true;
    out.refusal = AnchoringRefusal::RegionCoordinateDomain;
    return out;
  }
  if (query_length <= 0 || reference_length <= 0 ||
      reference_end <= reference_begin) {
    out.refused = true;
    out.refusal = AnchoringRefusal::EmptyOrInvalidInput;
    return out;
  }
  const int64_t lo = clamp_reference(reference_begin, reference_length);
  const int64_t hi = clamp_reference(reference_end, reference_length);
  if (hi <= lo) {
    out.refused = true;
    out.refusal = AnchoringRefusal::EmptyOrInvalidInput;
    return out;
  }
  if (hi - lo > params.maximum_locus_reference_span) {
    out.refused = true;
    out.refusal = AnchoringRefusal::RegionReferenceSpanOverBudget;
    return out;
  }

  SearchRegion region;
  region.region_id = 0;
  region.kind = SearchRegionKind::kCore;
  region.query_begin = 0;
  region.query_end = query_length;
  region.reference_begin = lo;
  region.reference_end = hi;
  // Unrestricted: every diagonal an anchor inside the window can occupy is
  // admissible. The spliced chain, not the region, decides which are colinear.
  region.diagonal_low = -mapping::kSignedCoordinateMaximum;
  region.diagonal_high = mapping::kSignedCoordinateMaximum;
  region.per_region_budget = params.locus_budget;
  out.regions.push_back(region);
  out.total_volume = saturating_volume(query_length, hi - lo);
  return out;
}

bool search_region_overlaps_anchor(const SearchRegion& region,
                                   int32_t query_begin,
                                   int32_t reference_begin, int32_t span) {
  if (span <= 0) return false;
  const int64_t query_end = static_cast<int64_t>(query_begin) + span;
  const int64_t reference_end = static_cast<int64_t>(reference_begin) + span;
  const int64_t diagonal =
      static_cast<int64_t>(reference_begin) - query_begin;
  return query_end > region.query_begin && query_begin < region.query_end &&
         reference_end > region.reference_begin &&
         reference_begin < region.reference_end &&
         diagonal >= region.diagonal_low && diagonal <= region.diagonal_high;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
