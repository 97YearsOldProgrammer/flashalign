// The candidate catalogue: ranked vote peaks deduplicated and admitted per lane.
#pragma once

#include "query_tiles.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::voting {

using CandidateId = std::int16_t;
inline constexpr CandidateId kNullCandidate = -1;

struct CandidateInput {
  std::uint64_t equivalence_key = 0;
  int lane = 0;
  int catalogue_rank = 0;
  int vote_evidence = 0;
  int chain_evidence = 0;
  QueryTileMask support;
  // Ratio admission with a CandidateMaskSource: the source slot of this input's tile mask,
  // and a coarse forward-tile range known to contain it (invalid when
  // coarse_tile_lo > coarse_tile_hi, and the mask is then always materialized).
  int mask_slot = -1;
  int coarse_tile_lo = 0;
  int coarse_tile_hi = -1;
};

// Materializes tile masks on demand for ratio admission, so a mask is built only when a
// decision or an admitted candidate needs it. Implementations memoize per slot.
class CandidateMaskSource {
 public:
  virtual ~CandidateMaskSource() = default;
  virtual const QueryTileMask& mask(int slot) = 0;
};

struct QueryCandidate {
  CandidateId id = kNullCandidate;
  std::uint64_t equivalence_key = 0;
  int lane = 0;
  int catalogue_rank = 0;
  int vote_evidence = 0;
  int chain_evidence = 0;
  QueryTileMask support;
};

struct CandidateCatalogue {
  std::vector<QueryCandidate> candidates;
};

// Inputs are ranked vote peaks. Equivalent keys are deduplicated deterministically, then
// at most max_per_lane entries are kept per lane. admission_ratio > 0 also requires
// vote_evidence >= ceil(ratio * competing_max), competing_max being the best vote among
// inputs whose tile masks overlap this one's by at least half of the smaller mask (an empty
// mask competes with everything); max_per_lane is then only a cost ceiling. With
// mask_source, inputs carry a mask_slot and masks are built lazily; without it,
// input.support is used as given. Masks and coarse ranges are on a grid of
// tile_count tiles. max_per_lane is at most kMaxCatalogueLaneBound, or
// kAllChainsLaneBound for an `unpartitioned` catalogue, which no solver takes
// (query_partition.h).
CandidateCatalogue build_candidate_catalogue(
    std::vector<CandidateInput> inputs, int max_per_lane,
    double admission_ratio = 0.0, CandidateMaskSource* mask_source = nullptr,
    int tile_count = kFixedQueryTiles, bool unpartitioned = false);

}  // namespace fa::cpu::voting
