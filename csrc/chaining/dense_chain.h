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


struct DenseChainStats {
  int64_t runs = 0;
  int64_t anchors = 0;
  // The DP input had diag_min_runs or more runs, so the exact arm ran.
  bool used_exact_arm = false;
  // DenseRunResult's exact-arm counters.
  int64_t exact_runs = 0;
  int64_t exact_nodes = 0;
  int64_t exact_steps = 0;
};

// Pool -> runs -> run DP -> anchor chains. The result owns the pool; every chain's anchor
// indices are strictly ascending in both q and r.
ChainResult chain_dense_colinear(std::vector<Anchor> anchors,
                                 DenseChainParams params,
                                 DenseChainStats* stats = nullptr);

} // namespace chaining
} // namespace cpu
} // namespace fa
