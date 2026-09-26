// Candidate pool harvest around the selected (rank-1) coarse locus.
#pragma once

#include "../../chaining/colinear_params.h"
#include "../../index/index.h"
#include "../../index/seed.h"
#include "../placement/types.h"
#include "skeleton_harvest.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct Rank1HarvestOptions {
  int global_occurrence_cap = 200;
  int chain_query_gap = 5000;
  int chain_min_count = 2;
  int chain_min_score = 20;
  float chain_gap_penalty = 0.0f;
  // Also the reference pad the harvest adds on each side of the coarse envelope.
  int max_intron = 200000;
};

struct Rank1HarvestResult {
  CandidatePool bounded_pool;
  // `nominated_pool` is the locus's own strand lane and `opposite_pool` the other lane,
  // already in its query coordinates; each is capped on its own. A refused harvest never
  // sets `routed`.
  CandidatePool nominated_pool;
  CandidatePool opposite_pool;
  bool routed = false;
  RoutingCounters routing;
  chaining::SplicedChainParams anchor_path_params;
  // Harvest work counters.
  HarvestWorkCounters work;
  // What the pool was built against, kept so mirror_candidate_pool_query_frame can
  // project the opposite query frame out of this pool instead of harvesting again. Set
  // whenever the region build succeeded.
  std::vector<SearchRegion> regions;
  SkeletonHarvestParams harvest_params;
  std::vector<SkippedBlock> skipped_blocks;
  // The same blocks in the opposite lane's coordinates.
  std::vector<SkippedBlock> opposite_skipped_blocks;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

inline placement::CoarseLocus
mirror_locus_query_orientation(const placement::CoarseLocus& locus,
                               std::uint32_t query_length) {
  placement::CoarseLocus mirror = locus;
  mirror.reverse = !locus.reverse;
  mirror.oriented_query_begin = locus.oriented_query_end <= query_length
                                    ? query_length - locus.oriented_query_end
                                    : 0u;
  mirror.oriented_query_end = locus.oriented_query_begin <= query_length
                                  ? query_length - locus.oriented_query_begin
                                  : query_length;
  return mirror;
}

// `posting_memo` is the caller's per-read memo, shared by both strand streams. It only
// saves posting probes across the read's harvests; the result is the same without it.
Rank1HarvestResult harvest_rank1(
    const SeedIndex& index, const std::vector<QuerySeed>& strand_seeds,
    const placement::CoarseLocus& rank1_locus,
    const std::vector<placement::CoarseDiagonalPeak>& retained_peaks,
    int read_length, int reference_length, const Rank1HarvestOptions& options,
    SeedPostingMemo* posting_memo = nullptr);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
