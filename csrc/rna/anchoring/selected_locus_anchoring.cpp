#include "selected_locus_anchoring.h"

#include "../../core/checked_range.h"
#include "skeleton_regions.h"

#include <algorithm>
#include <cstdint>
#include <utility>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

Rank1HarvestResult harvest_rank1(
    const SeedIndex& index, const std::vector<QuerySeed>& strand_seeds,
    const placement::CoarseLocus& rank1_locus,
    const std::vector<placement::CoarseDiagonalPeak>& retained_peaks,
    int read_length, int reference_length, const Rank1HarvestOptions& options,
    SeedPostingMemo* posting_memo) {
  Rank1HarvestResult result;
  // Skeleton nodes are anchor priority, not a search restriction: anchors on a voted
  // peak's diagonal win the deterministic caps first.
  std::vector<SkeletonNode> skeleton;
  int skeleton_count = 0;
  if (!mapping::checked_size_to_int(rank1_locus.peak_indices.size(),
                                    skeleton_count)) {
    result.refused = true;
    result.refusal = AnchoringRefusal::Rank1SkeletonCountDomain;
    return result;
  }
  const size_t skeleton_nodes = rank1_locus.peak_indices.size();
  skeleton.reserve(skeleton_nodes);
  for (size_t ordinal = 0; ordinal < skeleton_nodes; ++ordinal) {
    const uint32_t peak_index = rank1_locus.peak_indices[ordinal];
    if (peak_index >= retained_peaks.size()) {
      result.refused = true;
      result.refusal = AnchoringRefusal::InvalidRank1Skeleton;
      return result;
    }
    const placement::CoarseDiagonalPeak& peak = retained_peaks[peak_index];
    int query_begin = 0;
    int query_end = 0;
    int reference_begin = 0;
    int reference_end = 0;
    if (!mapping::checked_nonnegative_to_int(
            static_cast<int64_t>(peak.oriented_query_begin), query_begin) ||
        !mapping::checked_nonnegative_to_int(
            static_cast<int64_t>(peak.oriented_query_end), query_end) ||
        !mapping::checked_nonnegative_to_int(
            static_cast<int64_t>(peak.reference_begin), reference_begin) ||
        !mapping::checked_nonnegative_to_int(
            static_cast<int64_t>(peak.reference_end), reference_end) ||
        query_begin >= query_end || query_end > read_length ||
        reference_begin >= reference_end || reference_end > reference_length) {
      result.refused = true;
      result.refusal = AnchoringRefusal::Rank1SkeletonCoordinateDomain;
      return result;
    }
    SkeletonNode node;
    node.node_id = static_cast<uint32_t>(ordinal);
    node.query_begin = query_begin;
    node.query_end = query_end;
    node.reference_begin = reference_begin;
    node.reference_end = reference_end;
    node.median_diagonal = peak.median_diagonal;
    node.support = peak.distinct_seed_support;
    skeleton.push_back(node);
  }

  SkeletonRegionParams region_params;

  // Pad the envelope by a full maximum intron on each side: a read flank that produced no
  // peak must still reach its exon, and the maximum intron is the largest gap the chain
  // may cross. The pad is symmetric and ignores the query, so mirroring the locus leaves
  // the window unchanged and the mirrored query frame can be projected out of this pool
  // (mirror_candidate_pool_query_frame).
  const int64_t pad = std::max<int64_t>(1, options.max_intron);
  int64_t window_begin =
      static_cast<int64_t>(rank1_locus.reference_begin) - pad;
  int64_t window_end = static_cast<int64_t>(rank1_locus.reference_end) + pad;
  if (window_begin < 0)
    window_begin = 0;
  // Shrink symmetrically rather than refuse when the padded window exceeds the
  // single-locus span bound.
  if (window_end - window_begin > region_params.maximum_locus_reference_span) {
    const int64_t excess = (window_end - window_begin) -
                           region_params.maximum_locus_reference_span;
    const int64_t trim_left = std::min<int64_t>(
        excess / 2,
        std::max<int64_t>(0, static_cast<int64_t>(rank1_locus.reference_begin) -
                                 window_begin));
    window_begin += trim_left;
    window_end -= (excess - trim_left);
  }

  SearchRegionBuildResult region_result = derive_locus_search_region(
      read_length, window_begin, window_end, reference_length, region_params);
  if (region_result.refused) {
    result.refused = true;
    result.refusal = region_result.refusal;
    return result;
  }

  SkeletonHarvestParams harvest_params;
  harvest_params.occurrence_cap =
      static_cast<uint32_t>(std::max(1, options.global_occurrence_cap));
  // Bound the pool by the read's own evidence: in the right locus a seed has about one
  // true placement, so a small multiple of the seed count holds every real anchor, and
  // the rest is repeat multiplicity the O(n * max_iter) chain cannot use.
  const uint64_t seed_count = static_cast<uint64_t>(strand_seeds.size());
  harvest_params.global_budget = static_cast<uint32_t>(std::min<uint64_t>(
      harvest_params.global_budget, std::max<uint64_t>(1024, 8 * seed_count)));
  // The stream's frame lane is the locus's own orientation.
  harvest_params.frame_reverse_lane = rank1_locus.reverse;
  SkeletonCandidatePools pools = harvest_skeleton_candidate_pools(
      index, strand_seeds, rank1_locus.reference_id, read_length, skeleton,
      region_result.regions, harvest_params, posting_memo);
  result.bounded_pool = std::move(pools.bounded);
  result.nominated_pool = std::move(pools.nominated);
  result.opposite_pool = std::move(pools.opposite);
  result.routed = pools.routed;
  result.routing = pools.routing;
  result.work = pools.work;
  result.skipped_blocks = std::move(pools.skipped_blocks);
  result.opposite_skipped_blocks = std::move(pools.opposite_skipped_blocks);
  if (!result.routed) {
    // An unrouted result is always a refused harvest: carry its reason onto the routed
    // pools rather than fall back to the unfiltered pool.
    result.nominated_pool.refused = true;
    result.opposite_pool.refused = true;
    result.nominated_pool.refusal = result.opposite_pool.refusal =
        result.bounded_pool.refusal == AnchoringRefusal::None
            ? AnchoringRefusal::PoolRefused
            : result.bounded_pool.refusal;
  }
  result.regions = std::move(region_result.regions);
  result.harvest_params = harvest_params;

  result.anchor_path_params.bw = std::max(1, options.max_intron);
  result.anchor_path_params.max_dist_x = std::max(
      result.anchor_path_params.bw,
      options.max_intron > 0 ? options.max_intron : options.chain_query_gap);
  result.anchor_path_params.max_dist_y = std::max(1, options.chain_query_gap);
  result.anchor_path_params.min_cnt = std::max(1, options.chain_min_count);
  result.anchor_path_params.min_sc = std::max(1, options.chain_min_score);
  result.anchor_path_params.chn_pen_gap =
      options.chain_gap_penalty > 0.0f ? options.chain_gap_penalty
                                       : 0.008f * static_cast<float>(index.k());
  result.anchor_path_params.chn_pen_skip = 0.0f;
  return result;
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
