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
  // The lowest and highest diagonal (oriented r - q) of the screening chain's
  // anchors, when the screening pass held a chain.
  bool screening_diagonals = false;
  int screening_diagonal_low = 0;
  int screening_diagonal_high = 0;
};

// The chain a block prints from and its bounds in forward-read bp
// (dna/chain_ownership.h). path -1 is the candidate's primary, j its
// sibling_paths[j]; [anchor_begin, anchor_end) in chain order. score and
// anchors are the block's own; item_score and item_anchors its whole chain's;
// pool_subsc and pool_n_sub minimap2's subsc and n_sub over the chains of the
// same candidate attached to the block's owner.
struct DnaBlockPart {
  int path = -1;
  int anchor_begin = 0;
  int anchor_end = 0;
  int forward_begin = 0;
  int forward_end = 0;
  int score = 0;
  int anchors = 0;
  int item_score = 0;
  int item_anchors = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
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
  // Parallel to partition.selected.blocks once the owners are decided on the
  // dense chains' anchors; empty in a family placed on tiles alone (the
  // retained alternative, ranks 2..n, a rival placement).
  std::vector<DnaBlockPart> block_parts;
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
