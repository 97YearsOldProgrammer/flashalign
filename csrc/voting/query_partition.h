// Exact partition of a read's query tiles among catalogue candidates: each tile is owned
// by one candidate or by none, under an additive tile objective.
#pragma once

#include "candidate_catalogue.h"

#include <cstdint>
#include <limits>
#include <vector>

namespace fa::cpu::voting {

inline constexpr int kQueryTileCount = 128;
// Negative disables the same-candidate resume discount (the ordinary objective).
inline constexpr int kNoCandidateResumeDiscount = -1;
// Catalogue lane bound for the catalogue builder, the solver and residue admission. Ratio
// admission (vote_admission_ratio > 0) widens a lane up to this bound as a cost guard.
inline constexpr int kCatalogueLaneBound = 16;
// Lane bound of count admission (vote_admission_ratio == 0, --vote-ratio 0).
inline constexpr int kCountAdmissionLaneBound = 4;
static_assert(kCountAdmissionLaneBound <= kCatalogueLaneBound);
// Each of the four objective terms may contribute on every one of the 128 tiles, so the
// terms may use the full int range while the accumulated objective fits int64_t.
inline constexpr int kMaxQueryTileObjectiveTerm =
    std::numeric_limits<int>::max();
static_assert(static_cast<std::int64_t>(kMaxQueryTileObjectiveTerm) * 4 *
                  kQueryTileCount <
              std::numeric_limits<std::int64_t>::max());

// Preset-owned integer objective of the query partition. The DNA presets set
// {4, 12, 0, 1, 2}; RNA keeps these defaults.
struct QueryPartitionParameters {
  int supported_tile_reward = 1;
  int block_open_cost = 1;
  int unsupported_ownership_cost = 0;
  int null_tile_cost = 0;
  int minimum_supported_tiles_per_non_null_block = 2;
  // Negative: every non-null block pays block_open_cost. Otherwise a block that opens for
  // the candidate owning the run just before the ending one pays this instead: in A-B-A,
  // the interruption B would otherwise bill A the full open twice. Null tiles are a run
  // like any other, so A-null-A is a resumption and A-null-B-null-A is not. A residue-
  // admitted candidate also opens at this cost. The supporting-tile floor is unchanged.
  int same_candidate_resume_cost = kNoCandidateResumeDiscount;
};

// What the caller wants besides the winner; the runner-up costs extra and is opt-in.
enum class QueryPartitionRival : std::uint8_t {
  // Winner only. `runner_up` is absent and the first-difference window is -1.
  WinnerOnly,
  // Additionally return the exact best complete path whose assignment differs
  // from the winner's, under the same deterministic comparator. Costs a second
  // path slot per state.
  ExactSecondSolution,
};

struct QueryPartitionProblem {
  int tile_count = kQueryTileCount;
  QueryTileMask valid_tiles = valid_query_tiles(kQueryTileCount);
  CandidateCatalogue catalogue;
  QueryPartitionParameters parameters;
  QueryPartitionRival rival = QueryPartitionRival::WinnerOnly;
};

struct QueryBlock {
  CandidateId candidate = kNullCandidate;
  int query_tile_begin = 0;
  int query_tile_end = 0;
  int supporting_tiles = 0;

  friend bool operator==(const QueryBlock& left,
                         const QueryBlock& right) noexcept {
    return left.candidate == right.candidate &&
           left.query_tile_begin == right.query_tile_begin &&
           left.query_tile_end == right.query_tile_end &&
           left.supporting_tiles == right.supporting_tiles;
  }
};

struct QueryPartitionPath {
  std::int64_t score = 0;
  std::vector<CandidateId> assignment;
  std::vector<QueryBlock> blocks;
  int non_null_blocks = 0;
  int candidate_changes = 0;
  std::int64_t chain_evidence = 0;
  std::int64_t vote_evidence = 0;
  std::int64_t vote_rank_sum = 0;
};

struct QueryPartitionResult {
  QueryPartitionPath selected;
  // Present only when the problem requested ExactSecondSolution.
  QueryPartitionPath runner_up;
  bool has_runner_up = false;
  int first_difference_begin = -1;
  int first_difference_end = -1;
  // Work performed: (tile, state) cells advanced and edges relaxed.
  std::uint64_t dp_cells = 0;
  std::uint64_t dp_transitions = 0;
  // Exact ties of the additive key resolved by the lexicographic tail, and the tiles those
  // resolutions walked.
  std::uint64_t tail_comparisons = 0;
  std::uint64_t tail_steps = 0;
};

// Exact DP solver. The objective is additive over tiles plus a per-block open cost, and a
// non-null block may close only after the configured number of supporting tiles, so a
// left-to-right automaton over (candidate, supporting tiles in the open block, saturated at the
// minimum) is exact, and every tie-break key is additive over the same steps. The resume
// discount adds the previous run's candidate to the state, a factor of (candidates + 1), paid
// only when enabled. One comparable key per state, O(candidates) edges per tile, and
// assignments rebuilt at the end; exact key ties compare only the divergent backpointer
// suffixes.
QueryPartitionResult
solve_query_partition(const QueryPartitionProblem& problem);

} // namespace fa::cpu::voting
