// Placement chaining. The whole-read vote catalogue is the only source of
// candidates; this stage screens each candidate with a cheap chain, solves the
// partition on the screening chains, chains its owners over the whole query
// and decides the read's blocks on those chains' anchors. Both map-only
// projection and CIGAR realization consume the resulting family and anchor
// paths. No DP runs here.
#pragma once

#include "context.h"
#include "retained_seed_density.h"
#include "placement_family_adapter.h"
#include "alternative_hypothesis.h"
#include "../chaining/anchor.h"
#include "../chaining/colinear_params.h"
#include "../chaining/dense_chain.h"
#include "../chaining/partition.h"
#include "../index/seed.h"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>

namespace fa::cpu::lr {

enum class DnaPlacementChainStatus : std::uint8_t {
  NotSelected,
  Sparse,
  Accepted,
  InvalidReference,
  MissingSeeds,
  NoIntervalHits,
  NoChain,
};

struct DnaPlacementCandidateChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  ::fa::cpu::voting::QueryTileMask sparse_support;
  ::fa::cpu::voting::QueryTileMask dense_support;
  std::vector<chaining::Anchor> primary;
  // Sibling paths from the same exact DP. Realization and the map-only
  // projection may use one when the primary path has no anchors in a selected
  // block (dna_selected_sibling_path).
  std::vector<std::vector<chaining::Anchor>> sibling_paths;
  // Parallel to sibling_paths: each sibling's chain score. rival_sibling is
  // the index of the sibling rival_chain_score came from
  // (ChainPartition::f2_index), or -1 when there is none or it was not
  // materialized. Whole-query passes only; read by dna_sibling_rival_chain,
  // and by the map-only projection for a record projected from a sibling.
  std::vector<int> sibling_scores;
  int rival_sibling = -1;
  DnaPlacementChainStatus status = DnaPlacementChainStatus::NotSelected;
  std::int64_t interval_hits = 0;
  // Undercounts: the harvest gate skips slicing the seeds it drops.
  int filtered_hits = 0;
  int chain_score = 0;
  int rival_chain_score = 0;
  int overlapping_rivals = 0;
  int chain_anchors = 0;
  // The screening pass's result, kept because the exact restore later
  // overwrites the chain fields of selected candidates.
  int screening_chain_score = 0;
  int screening_chain_anchors = 0;
  int screening_forward_query_begin = -1;
  int screening_forward_query_end = -1;
  int oriented_query_begin = -1;
  int oriented_query_end = -1;
  int forward_query_begin = -1;
  int forward_query_end = -1;
  int reference_begin = -1;
  int reference_end = -1;
  int initial_selected_tiles = 0;
  int final_selected_tiles = 0;
  bool exact = false;
  // Harvest counters. sparse_anchors: anchors built inside the harvest loop.
  // deferred_anchors: postings of the slices deferred as metadata, within
  // this candidate's window. rescued_anchors: anchors built by the deferred
  // restore. The pool is sparse + rescued. The whole-query pass defers every
  // slice, so there sparse_anchors is 0 and rescued_anchors is the pool; the
  // screening pass defers nothing.
  std::uint64_t sparse_anchors = 0;
  std::uint64_t deferred_anchors = 0;
  std::uint64_t rescued_anchors = 0;
  // Runs produced by the whole-query pass's collapse; 0 on the screening pass.
  std::int64_t dense_runs = 0;
};

// Catalogue rivals chained over the whole query for the MAPQ, per ranking.
// mm_set_mapq2 reads one best rival and a count; two let the best rival by
// vote lose to the runner-up on chain score. Under
// kDnaChainMapqStudyBlockRivals each block also ranks its own.
inline constexpr int kDnaMapqRivalChains = 2;
// Only rivals with at least 1/kDnaMapqRivalVoteDenominator of the strongest
// owner's vote are chained; the rest still enter the MAPQ through their vote.
inline constexpr int kDnaMapqRivalVoteDenominator = 4;

// A catalogue rival chained over the whole query for the MAPQ; a candidate's
// own chain partition only sees rivals inside its window. Evidence only: kept
// out of `candidates`, which is parallel to family.candidates and drives
// selection.
struct DnaRivalExactChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  DnaPlacementCandidateChain chain;
  // Taken from a record that already had a whole-query chain instead of
  // chained again.
  bool reused = false;
};

// A rank 2..n alternative (-N n on a DNA preset) with its exact whole-query
// chain. Output only: no selection, MAPQ or rival lookup reads it.
struct DnaRankedAlternative {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  DnaPlacementCandidateChain exact;
};

struct DnaPlacementChainingResult {
  DnaPlacementFamily family;
  std::vector<DnaPlacementCandidateChain> candidates;
  DnaAlternativeSelection alternative;
  DnaPlacementCandidateChain alternative_exact;
  // Ranks 2..n of the alternative ranking whose exact restore succeeded, in
  // rank order. Empty unless n >= 2.
  std::vector<DnaRankedAlternative> ranked_alternatives;
  // The MAPQ's cross-locus rivals; empty when none is worth chaining.
  std::vector<DnaRivalExactChain> rival_exact;
  // Whole-query chains built for the MAPQ rivals, reused ones excluded.
  int mapq_rival_chains = 0;
  // The read's seed density, kept past placement when realization gates the
  // late inversion probe on it (inv_local_chain.h); null otherwise.
  std::shared_ptr<const RetainedSeedDensity> inversion_gate_seeds;
  std::int64_t initial_score = 0;
  int initial_blocks = 0;
  bool selection_changed = false;
  bool accepted = false;
  // Not accepted because the partition's owners have no whole-query chain
  // that qualifies to own a block (dna/chain_ownership.h).
  bool no_owner_chain = false;

  const DnaPlacementCandidateChain*
  find(::fa::cpu::voting::CandidateId candidate) const noexcept;

  // A candidate's whole-query exact chain: its MAPQ rival chain (any status),
  // else the retained alternative's exact restore, else its own record if it
  // was exact-restored. Null when there is none. Callers that need an
  // accepted chain check the status.
  const DnaPlacementCandidateChain*
  whole_query_chain(::fa::cpu::voting::CandidateId candidate) const noexcept;
};

// The winner's rival sibling as a chain record of its own, for
// build_dna_rival_placement: the sibling path as `primary`, exact and
// Accepted, with its own score, anchor count, dense_support (on tile_count
// tiles) and spans; other fields are left empty. Empty (NotSelected) when the
// winner has no rival sibling.
DnaPlacementCandidateChain
dna_sibling_rival_chain(const DnaPlacementCandidateChain& winner, bool reverse,
                        int read_length, int seed_length, int tile_count);

// The whole-query chain of the committed hypothesis. A promoted alternative
// uses its restricted rerun, not the stable run's record for that candidate.
const DnaPlacementCandidateChain* dna_committed_winner_chain(
    const DnaPlacementChainingResult* stable,
    const DnaPlacementChainingResult* promoted_alternative,
    ::fa::cpu::voting::CandidateId primary_candidate,
    bool primary_is_alternative) noexcept;

// Chain parameters for the candidate paths under diagonal band `band`.
chaining::ColinearChainParams
dna_candidate_chain_params(const DnaContext& context, int band,
                           int seed_length, int read_length);

// The same parameters for the dense run chain: the scoring is unchanged, and
// max_iter and max_skip have no counterpart. The diagonal-keyed search runs
// when diag_min_runs >= 0 && runs >= diag_min_runs, the linear scan
// otherwise; both are exact.
chaining::DenseChainParams
dna_dense_chain_params(const chaining::ColinearChainParams& params,
                       int seed_length, int diag_min_runs);

// Runs the restore/stabilize pipeline over {retained alternative, null}. The
// returned family uses solver id 0 and records the catalogue id in
// family.original_candidate_id.
DnaPlacementChainingResult
build_dna_alternative_placement(const DnaContext& context,
                                const DnaPlacementFamily& stable_family,
                                const DnaPlacementChainingResult& stable,
                                const std::vector<std::uint8_t>& forward_query,
                                const std::vector<std::uint8_t>& reverse_query);

// The same over {one catalogue candidate, null}, from a given exact
// whole-query chain.
DnaPlacementChainingResult
build_dna_rival_placement(const DnaContext& context,
                          const DnaPlacementFamily& stable_family,
                          ::fa::cpu::voting::CandidateId original,
                          const DnaPlacementCandidateChain& exact_chain,
                          const std::vector<std::uint8_t>& forward_query,
                          const std::vector<std::uint8_t>& reverse_query);

DnaPlacementChainingResult build_dna_placement_chains(
    const DnaContext& context, DnaPlacementFamily family,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    const std::vector<ChainWindowRetainedSeed>* forward_seeds,
    const std::vector<ChainWindowRetainedSeed>* reverse_seeds,
    const std::vector<QuerySeed>* fine_forward_seeds,
    const std::vector<QuerySeed>* fine_reverse_seeds,
    ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots);

} // namespace fa::cpu::lr
