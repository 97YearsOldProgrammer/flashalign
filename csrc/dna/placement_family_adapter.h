#pragma once

#include "context.h"
#include "../seeding/scratch.h"
#include "../seeding/types.h"
#include "../voting/query_partition.h"

#include <cstdint>
#include <optional>
#include <vector>

namespace fa::cpu::lr {

struct DnaPlacementCandidate {
  ::fa::cpu::voting::CandidateId id =
      ::fa::cpu::voting::kNullCandidate;
  VotePeak peak;
  ::fa::cpu::voting::QueryTileMask support;
  std::uint64_t equivalence_key = 0;
  int lane = 0;
  int catalogue_rank = 0;
  int vote_evidence = 0;
  // Chain score over the screening-pass anchors.
  int screening_chain_score = 0;
  // A recovered candidate has no whole-read vote behind it, so it may compete
  // for supplementary blocks but must never be promoted to primary.
  bool residue_admitted = false;
};

struct DnaPlacementFamily {
  int read_length = 0;
  int seed_length = 0;
  // The partition's tiles (--tiles); masks, blocks and assignments use them.
  int tile_count = ::fa::cpu::voting::kQueryTileCount;
  // Set only for a restricted alternative-hypothesis catalogue, whose solver
  // ids are positional: the id of the candidate in the original catalogue.
  std::optional<::fa::cpu::voting::CandidateId> original_candidate_id;
  int forward_candidates = 0;
  int reverse_candidates = 0;
  std::uint64_t exact_posting_tests = 0;
  std::vector<DnaPlacementCandidate> candidates;
  ::fa::cpu::voting::QueryPartitionResult partition;
  // Set when build_dna_placement_family left `partition` unsolved for
  // build_dna_placement_chains to solve.
  bool partition_deferred = false;
  bool valid = false;

  const DnaPlacementCandidate* find(
      ::fa::cpu::voting::CandidateId id) const noexcept;
};

// Maps an oriented seed position to its forward-read query tile of
// tile_count. The generic catalogue and partition code only sees forward-read
// tiles.
int dna_forward_query_tile(
    int oriented_seed_position, bool reverse, int read_length, int seed_length,
    int tile_count) noexcept;

// Builds the candidate catalogue from the whole-read vote peaks and, for a
// non-empty catalogue, defers the query partition solve (see
// partition_deferred). A candidate's support is the set of forward query tiles
// holding a retained seed with a posting compatible with its geometry.
DnaPlacementFamily build_dna_placement_family(
    const DnaContext& context, const std::vector<std::uint8_t>& forward_query,
    const std::vector<VotePeak>& whole_read_peaks,
    const ChainWindowPeakScratch& forward_scratch,
    const ChainWindowPeakScratch& reverse_scratch);

}  // namespace fa::cpu::lr
