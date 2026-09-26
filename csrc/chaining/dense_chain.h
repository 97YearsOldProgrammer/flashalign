// Dense-pool colinear chaining: collapse an anchor pool to same-diagonal runs, run the
// minimap2 recurrence over the runs and expand the run chains back to anchors. With one
// anchor per run the recurrence is minimap2's anchor chain.
#pragma once

#include "anchor.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace chaining {

// A run is one gap-free exact match: a same-diagonal anchor joins the open run when its
// k-mer overlaps or touches the run's last one (an r-step of at most k), so a run's weight
// is its length. Every chain of the collapsed DP then has an equal-score twin at the
// anchor level (its own expanded anchor path), so the collapsed optimum never exceeds the
// anchor-level one. The runs are the perfect chains of Rizzo, Caceres & Makinen (ISBRA 2025).

// Pools of at least this many runs take the exact arm below; smaller ones take the linear
// scan under the prefix-max bound. Both are exact searches of the recurrence, so the
// threshold chooses a cost, not an approximation.
inline constexpr int32_t kDenseDiagMinRuns = 512;

// The exact arm: branch and bound over a static k-d tree of the runs (median splits
// alternating on q_end and the end diagonal), whose nodes aggregate over their scored runs
// the max score, max index, coordinate minima, max q_end and diagonal interval. A node's
// bound is its max score plus the largest gain any member could be credited (at the
// smallest overlap its minima allow) minus the smallest drift any could be charged (at the
// nearest diagonal); it assumes non-negative penalties. A node is dropped when its bound
// cannot beat the running best, or can only tie it from an index at or below the elected
// predecessor. That is the linear scan's election rule, so both arms return the same scores
// and predecessors. The search is seeded with the last run on B's own begin diagonal, where
// the drift is 0, and descends into the better child first.
//
// The wide implementation runs the same search on a 16-ary tree whose nodes and leaf blocks
// hold sixteen lanes (AVX-512F on x86-64, NEON on AArch64), bounding a node's children or
// pricing a leaf's runs in one vector pass. Its bound also charges the drift's log term,
// after checking per (gap, bw) that the drift is non-decreasing in dd, and computes the
// trimmed gain in double precision, exact for runs shorter than 2^26 with coordinates below
// 2^30; a pool outside that domain takes the portable arm. dense_chain.cpp is compiled with
// -ffp-contract=off so the vector floats equal dense_step's.

// Runs per candidate above which the chain is refused (or the satellite skip runs). It
// bounds the collapsed pool, not the anchor count: a huge anchor pool that collapses to few
// runs is the case the collapse exists for.
inline constexpr int64_t kDenseRunCap = 65536;
// Memory ceiling on the anchor pool the collapse consumes; the run cap is the operative bound.
inline constexpr size_t kDenseRawAnchorCeiling = 1000000;
// Query-tile width, in oriented-query bp, of the pool's spatial cut; a run's tile is that of
// its query midpoint. The DNA pool trim shares its per-read budget over these tiles and the
// satellite skip classifies them, so a satellite stretch gets a share by the query it covers
// rather than by how many runs it piles up. Far below max_dist_x, far above a seed length.
inline constexpr int32_t kDenseAdmitTileBp = 1024;
// Runs in one query tile at or above which the tile is dense. A unique kilobase collapses
// to a handful of runs and a few-copy duplication to tens; hundreds means a tandem ladder.
inline constexpr int64_t kDenseSkipTileRuns = 256;
// Share of a dense tile's weight its heaviest diagonal must hold for the tile to be coherent
// (one densely supported alignment, kept at that diagonal) rather than a ladder of rungs of
// comparable weight (dropped).
inline constexpr double kDenseSkipCoherentShare = 0.5;

// A maximal same-diagonal anchor run, cut on a diagonal change or an r-step above `span`.
// `w` is its length, span plus the r-steps: the score the anchor-level DP gives it.
struct DenseRun {
  int64_t d = 0;
  int32_t q_begin = 0;
  int32_t q_end = 0;
  int32_t r_begin = 0;
  int32_t r_end = 0;
  int32_t w = 0;
  int32_t raw_count = 0;
  // Pool indices of the run's first and last anchor; with the builder's `links` they
  // expand the run to its anchors in ascending r (and q).
  int32_t first = -1;
  int32_t last = -1;

  int32_t length() const { return q_end - q_begin; }
};

struct DenseChainParams {
  // The seed length k: every anchor's span, and the run rule's reach.
  int32_t span = 21;
  int32_t max_dist_x = 20000;
  int32_t max_dist_y = 20000;
  int32_t bw = 20000;
  int32_t min_sc = 1;
  int32_t min_cnt = 1;
  float chn_pen_gap = 0.168f;
  float chn_pen_skip = 0.0f;
  int32_t diag_min_runs = kDenseDiagMinRuns;
  // Exact-arm implementation: 1 uses the wide arm where the CPU and the pool allow it and
  // the portable one otherwise; 2 forces the portable one. Both give the same result.
  int32_t exact = 1;
  // Runs above which the candidate is refused before the DP runs; 0 disables the cap.
  int64_t run_cap = kDenseRunCap;
  // Over the run cap, drop the satellite query tiles (dense_skip_satellite_tiles) and
  // chain the rest instead of refusing. Nothing reads it below the cap.
  bool skip_satellites = false;
};

// Builds the runs in one pass: a same-diagonal anchor joins its diagonal's open run when
// its r-step from the run's last anchor is at most `span`, else it opens a new run. `pool`
// must be ascending in r. links[i] is the next anchor in i's run, or -1. The runs come back
// sorted by (r_begin, q_begin, q_end).
std::vector<DenseRun>
dense_collapse_runs_stream(const std::vector<Anchor>& pool, int32_t span,
                           std::vector<int32_t>& links);

struct DenseRunChain {
  std::vector<int32_t> idx; // run indices, query order
  int32_t score = 0;
};

struct DenseRunResult {
  std::vector<DenseRunChain> chains;
  std::vector<int32_t> scores;
  // Exact-arm work, 0 when it did not run: successors searched, bounds computed (tree
  // nodes, or lanes in the wide arm) and dense_step calls. exact_steps can differ between
  // standard libraries, since the run order inside a wide leaf block is unspecified.
  int64_t exact_runs = 0;
  int64_t exact_nodes = 0;
  int64_t exact_steps = 0;
  // The wide implementation ran.
  bool exact_wide = false;
};

// The run DP and chain assembly. `runs` must be in ascending r_begin order.
DenseRunResult dense_run_chain(const std::vector<DenseRun>& runs,
                               const DenseChainParams& params);

// The satellite skip for a pool over the run cap: classifies each query tile and returns
// the kept runs' indices in ascending order (every index when no tile is dense):
//   satellite  kDenseSkipTileRuns or more runs and no diagonal holding
//              kDenseSkipCoherentShare of the weight: dropped whole, and the DP crosses
//              the hole as one colinear jump;
//   coherent   dense, with such a diagonal: only that diagonal's runs are kept;
//   sparse     everything else, kept whole.
// If the kept set still exceeds `run_cap` (> 0), the run_cap heaviest runs are kept.
// `skipped_tiles` and `coherent_tiles` receive the tile counts when non-null. The kept set
// is a pure function of the pool.
std::vector<int32_t>
dense_skip_satellite_tiles(const std::vector<DenseRun>& runs, int64_t run_cap,
                           int64_t* skipped_tiles, int64_t* coherent_tiles);

struct DenseChainStats {
  int64_t runs = 0;
  int64_t anchors = 0;
  // The DP input had diag_min_runs or more runs, so the exact arm ran.
  bool used_exact_arm = false;
  // Over the run cap with skip_satellites off: no DP ran and the result has no chains.
  bool run_cap_refused = false;
  // DenseRunResult's exact-arm counters.
  int64_t exact_runs = 0;
  int64_t exact_nodes = 0;
  int64_t exact_steps = 0;
  // The satellite skip, run instead of the refusal: tiles dropped, dense tiles kept at
  // their dominant diagonal, runs handed to the DP and runs dropped. The counts are -1
  // when the skip did not run; `runs` and `anchors` describe the full collapse.
  bool pool_skipped = false;
  int64_t skipped_tiles = -1;
  int64_t skipped_runs = -1;
  int64_t kept_runs = -1;
  int64_t coherent_tiles = -1;
};

// Pool -> runs -> run DP -> anchor chains. The result owns the pool; every chain's anchor
// indices are strictly ascending in both q and r.
ChainResult chain_dense_colinear(std::vector<Anchor> anchors,
                                 DenseChainParams params,
                                 DenseChainStats* stats = nullptr);

} // namespace chaining
} // namespace cpu
} // namespace fa
