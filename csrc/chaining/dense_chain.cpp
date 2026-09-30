#include "dense_chain.h"
#include "../core/flat_int64_map.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// The wide exact arm: AVX-512F lanes on x86-64 (a per-function target attribute, entered
// only when the CPU reports the feature) or NEON on AArch64. Other targets use ExactArm.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define FA_DENSE_EXACT_WIDE 1
#define FA_DENSE_EXACT_WIDE_NEON 0
#elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
#include <arm_neon.h>
#define FA_DENSE_EXACT_WIDE 1
#define FA_DENSE_EXACT_WIDE_NEON 1
#else
#define FA_DENSE_EXACT_WIDE 0
#define FA_DENSE_EXACT_WIDE_NEON 0
#endif

namespace fa {
namespace cpu {
namespace chaining {
namespace {

constexpr int32_t kNegInf = INT32_MIN;

// minimap2's fast log2, as in colinear_chain.cpp.
float mg_log2(float x) {
  union {
    float f;
    uint32_t i;
  } z = {x};
  float log_2 = static_cast<float>(((z.i >> 23) & 255)) - 128.0f;
  z.i &= ~(255u << 23);
  z.i += 127u << 23;
  log_2 += (-0.34484843f * z.f + 2.02466578f) * z.f - 0.67487759f;
  return log_2;
}

bool run_order_less(const DenseRun& left, const DenseRun& right) {
  if (left.r_begin != right.r_begin)
    return left.r_begin < right.r_begin;
  if (left.q_begin != right.q_begin)
    return left.q_begin < right.q_begin;
  return left.q_end < right.q_end;
}

// The pairwise recurrence over runs, A before B. Admissibility is compute_score's test on
// the boundary pair (A's last anchor, B's first), plus a frontier rule: B must advance both
// q_end and r_end, since a run nested inside A has no anchor past A's last one. The gain
// trims B's weight by the overlap and caps the first anchor's loss at the span:
//     o        = max(q-overlap, r-overlap) of A's end into B
//     W_trim   = floor(W_B * (L_q(B) - o) / L_q(B))
//     deficit  = max(0, -min(g_q, g_r))
//     gain     = max(0, min(W_trim, W_B - min(deficit, span)))
// With one anchor per run this is compute_score's min(q_span, dg). RunA and RunB are
// DenseRun or HotRun.
template <class RunB, class RunA>
int32_t dense_step(const RunB& b, const RunA& a,
                   const DenseChainParams& params) {
  const int64_t g_q = static_cast<int64_t>(b.q_begin) - a.q_end;
  const int64_t g_r = static_cast<int64_t>(b.r_begin) - a.r_end;
  const int64_t dq = g_q + params.span;
  const int64_t dr = g_r + params.span;
  // Forward progress on the run starts (compute_score's dq <= 0 / dr == 0 test for single
  // anchors); an overlapping successor is priced by the trim instead of rejected.
  if (b.q_begin <= a.q_begin || b.r_begin <= a.r_begin)
    return kNegInf;
  // Forward progress on the ends too, so A's end is the chain's frontier and the trim exact.
  if (b.q_end <= a.q_end || b.r_end <= a.r_end)
    return kNegInf;
  // The second reach test reads dq, not dr; that is minimap2's own text.
  if (dq > params.max_dist_x || dq > params.max_dist_y)
    return kNegInf;
  const int64_t dd = dr > dq ? dr - dq : dq - dr;
  if (dd > params.bw)
    return kNegInf;
  const int64_t dg = dr < dq ? dr : dq;
  const int64_t length = b.q_end - b.q_begin;
  const int64_t overlap = std::max<int64_t>(
      0, std::max<int64_t>(static_cast<int64_t>(a.q_end) - b.q_begin,
                           static_cast<int64_t>(a.r_end) - b.r_begin));
  int64_t trimmed = 0;
  // No overlap is the common case, and there the trimmed weight is W_B without a division.
  if (overlap == 0 && length > 0)
    trimmed = b.w;
  else if (length > 0 && overlap < length)
    trimmed = (static_cast<int64_t>(b.w) * (length - overlap)) / length;
  const int64_t g_min = g_q < g_r ? g_q : g_r;
  // The head loss is the first anchor's lost contribution, bounded by the span as in
  // compute_score's min(q_span, dg).
  const int64_t deficit = g_min < 0 ? -g_min : 0;
  const int64_t head_loss = deficit < params.span ? deficit : params.span;
  int64_t gain =
      std::min<int64_t>(trimmed, static_cast<int64_t>(b.w) - head_loss);
  if (gain < 0)
    gain = 0;
  int32_t drift = 0;
  if (dd) {
    // An overlapping successor (dg <= 0) skips no reference and pays no skip term;
    // otherwise a positive chn_pen_skip would reward overlap and break both search bounds,
    // which need drift >= 0.
    const float skip = dg > 0
                           ? params.chn_pen_skip * static_cast<float>(dg)
                           : 0.0f;
    const float linear = params.chn_pen_gap * static_cast<float>(dd) + skip;
    const float log_penalty = mg_log2(static_cast<float>(dd + 1));
    drift = static_cast<int32_t>(linear + 0.5f * log_penalty);
  } else if (dg > params.span) {
    drift = static_cast<int32_t>(params.chn_pen_skip * static_cast<float>(dg));
  }
  return static_cast<int32_t>(gain - drift);
}

// colinear_chain.cpp's backtrack_end, over runs.
int64_t backtrack_end(int32_t max_drop, int32_t score_end,
                      const std::vector<int32_t>& scores,
                      const std::vector<int32_t>& predecessors,
                      std::vector<int32_t>& marks, int64_t start) {
  int64_t i = start;
  int64_t end_i = -1;
  int64_t max_i = i;
  int32_t max_score = 0;
  if (i < 0 || marks[static_cast<size_t>(i)] != 0)
    return i;
  do {
    marks[static_cast<size_t>(i)] = 2;
    end_i = i = predecessors[static_cast<size_t>(i)];
    const int32_t score =
        i < 0 ? score_end : score_end - scores[static_cast<size_t>(i)];
    if (score > max_score) {
      max_score = score;
      max_i = i;
    } else if (max_score - score > max_drop) {
      break;
    }
  } while (i >= 0 && marks[static_cast<size_t>(i)] == 0);
  for (i = start; i >= 0 && i != end_i;
       i = predecessors[static_cast<size_t>(i)])
    marks[static_cast<size_t>(i)] = 0;
  return max_i;
}

// The DP's compact view of a run: the coordinates dense_step reads, the run's score and
// `upper`, the prefix maximum of scores that the linear scan prunes on.
struct HotRun {
  int32_t q_begin;
  int32_t q_end;
  int32_t r_begin;
  int32_t r_end;
  int32_t w;
  int32_t score;
  int32_t upper;
};

// Exact-arm tree node (see dense_chain.h), one per run (a leaf) and one per split. The
// aggregates cover the node's scored runs. 64 bytes, one cache line.
struct ExactNode {
  int32_t left;
  int32_t right;
  int32_t parent;
  int32_t run; // the leaf's run; -1 at a split
  int32_t max_index;
  int32_t score;
  int32_t min_q_begin;
  int32_t min_r_begin;
  int32_t min_q_end;
  int32_t min_r_end;
  int32_t max_q_end;
  int64_t min_d;
  int64_t max_d;
};

// Search stack entry: a node and the bound its parent computed for it.
struct ExactVisit {
  int64_t bound;
  int32_t node;
};

// The build's explicit stack entry: the index range still to split, the node
// that owns it and the depth that picks the split key.
struct ExactBuild {
  int32_t lo;
  int32_t hi;
  int32_t node;
  int32_t depth;
};

// The wide exact arm's records: 16-ary tree nodes that hold their children's aggregates in
// SoA lanes, and leaf blocks of up to 16 runs, so one vector pass bounds a node's children
// or prices a block's runs.
#if FA_DENSE_EXACT_WIDE

constexpr int32_t kWideLanes = 16;
// |coordinate| < 2^30 keeps every lane difference inside int32; a run
// shorter than 2^26 keeps w * length under 2^53, where the trimmed gain's
// double-precision division is exact.
constexpr int32_t kWideCoordinateLimit = 1 << 30;
constexpr int32_t kWideLengthLimit = 1 << 26;

// One internal node: its CHILDREN's aggregates, lane c = child c, each array
// one cache line. The children themselves are contiguous from first_child,
// in tree_ when leaf_children is 0 and in leaves_ when it is 1.
struct alignas(64) WideNode {
  int32_t max_index[kWideLanes];
  int32_t score[kWideLanes];
  int32_t min_q_begin[kWideLanes];
  int32_t min_r_begin[kWideLanes];
  int32_t min_q_end[kWideLanes];
  int32_t min_r_end[kWideLanes];
  int32_t max_q_end[kWideLanes];
  int32_t min_d[kWideLanes];
  int32_t max_d[kWideLanes];
  int32_t first_child;
  int32_t child_count;
  int32_t parent;
  int32_t slot;
  int32_t leaf_children;
};

// One leaf block: up to 16 runs, lane l = one run (index -1 in an empty
// lane). d_end is the predecessor key r_end - q_end; score is written at
// activation, so an unscored lane is rejected by index >= B alone.
struct alignas(64) WideLeaf {
  int32_t q_begin[kWideLanes];
  int32_t q_end[kWideLanes];
  int32_t r_begin[kWideLanes];
  int32_t r_end[kWideLanes];
  int32_t w[kWideLanes];
  int32_t score[kWideLanes];
  int32_t index[kWideLanes];
  int32_t d_end[kWideLanes];
  int32_t parent;
  int32_t slot;
};

// The build's records: both split keys gathered once per wide level so the
// four median splits below it never dereference a run.
struct WideKey {
  int32_t q_end;
  int32_t d_end;
  int32_t id;
};
struct WideTask {
  int32_t lo;
  int32_t hi;
  int32_t depth; // binary depth of the range, which picks the split keys
  int32_t node;
};
// A child left on the search stack with the bound and max index its parent's
// pass computed for it.
struct WideVisit {
  int32_t bound;
  int32_t max_index;
  int32_t child;
  int32_t is_leaf;
};
#endif // FA_DENSE_EXACT_WIDE

// Per-thread scratch for dense_run_chain and chain_dense_colinear, kept at capacity across
// calls. Every vector is resized or assigned before it is read.
struct DenseChainScratch {
  ::fa::cpu::lr::FlatInt64Map<int32_t> open_runs;
  std::vector<int32_t> previous_r;
  std::vector<HotRun> hot;
  std::vector<int32_t> predecessors;
  std::vector<int32_t> marks;
  std::vector<std::pair<int32_t, int64_t>> ranked_ends;
  std::vector<int32_t> links;
  // The exact arm's tree and stacks; empty unless a pool takes that arm.
  std::vector<ExactNode> exact_tree;
  std::vector<int32_t> exact_ids;
  std::vector<int32_t> exact_leaf;
  std::vector<ExactVisit> exact_stack;
  std::vector<ExactBuild> exact_build;
  // End diagonal -> the largest run index activated on it: the exact arm's search seed.
  ::fa::cpu::lr::FlatInt64Map<int32_t> exact_diag_last;
#if FA_DENSE_EXACT_WIDE
  // The wide arm's tree, leaf blocks, build records and stack (exact_diag_last is shared),
  // and the last (gap, bw) whose drift was checked for monotonicity, with the verdict.
  std::vector<WideNode> wide_nodes;
  std::vector<WideLeaf> wide_leaves;
  std::vector<int32_t> wide_run_leaf;
  std::vector<WideKey> wide_keys;
  std::vector<int32_t> wide_ids;
  std::vector<WideTask> wide_tasks;
  std::vector<int32_t> wide_cuts;
  std::vector<WideVisit> wide_stack;
  float wide_gap = -1.0f;
  int32_t wide_bw = -1;
  bool wide_log_term = false;
#endif
};

DenseChainScratch& dense_chain_scratch() {
  static thread_local DenseChainScratch scratch;
  return scratch;
}

} // namespace

std::vector<DenseRun>
dense_collapse_runs_stream(const std::vector<Anchor>& pool, int32_t span,
                           std::vector<int32_t>& links) {
  links.assign(pool.size(), -1);
  std::vector<DenseRun> runs;
  if (pool.empty())
    return runs;
  // d -> index of that diagonal's open run. The pool is ascending in r, so each diagonal's
  // anchors arrive in ascending r and its cuts do not depend on other diagonals.
  DenseChainScratch& scratch = dense_chain_scratch();
  ::fa::cpu::lr::FlatInt64Map<int32_t>& open = scratch.open_runs;
  open.clear();
  // At most one slot per anchor; reserve twice that to keep probe chains short when
  // every run is a single anchor.
  open.reserve(pool.size() * 2);
  std::vector<int32_t>& previous_r = scratch.previous_r; // parallel to `runs`
  previous_r.clear();
  runs.reserve(pool.size() / 4 + 1);
  previous_r.reserve(pool.size() / 4 + 1);
  for (size_t at_pool = 0; at_pool < pool.size(); ++at_pool) {
    const Anchor& anchor = pool[at_pool];
    const int32_t index = static_cast<int32_t>(at_pool);
    const int64_t d = static_cast<int64_t>(anchor.r) - anchor.q;
    auto found = open.find(d);
    const bool inserted = found == open.end();
    int32_t& which = inserted ? (open[d] = -1) : found->second;
    if (!inserted && which >= 0) {
      const size_t at = static_cast<size_t>(which);
      const int64_t step = static_cast<int64_t>(anchor.r) - previous_r[at];
      if (step <= static_cast<int64_t>(span)) {
        // The k-mers overlap or touch, so the step is the newly matched
        // length.
        DenseRun& current = runs[at];
        current.w += static_cast<int32_t>(step);
        current.q_end = anchor.q + span;
        current.r_end = anchor.r + span;
        ++current.raw_count;
        links[static_cast<size_t>(current.last)] = index;
        current.last = index;
        previous_r[at] = anchor.r;
        continue;
      }
    }
    DenseRun current;
    current.d = d;
    current.q_begin = anchor.q;
    current.q_end = anchor.q + span;
    current.r_begin = anchor.r;
    current.r_end = anchor.r + span;
    current.w = span;
    current.raw_count = 1;
    current.first = current.last = index;
    which = static_cast<int32_t>(runs.size());
    runs.push_back(current);
    previous_r.push_back(anchor.r);
  }
  std::sort(runs.begin(), runs.end(), run_order_less);
  return runs;
}

// The portable exact arm (see dense_chain.h): a static binary k-d tree over the runs,
// aggregates over the scored ones, and a branch-and-bound descent per successor that
// prices the leaves it reaches with dense_step.
class ExactArm {
public:
  ExactArm(const std::vector<HotRun>& hot, const DenseChainParams& params,
           DenseChainScratch& scratch, int32_t count)
      : hot_(hot), params_(params), tree_(scratch.exact_tree),
        ids_(scratch.exact_ids), leaf_(scratch.exact_leaf),
        stack_(scratch.exact_stack), build_(scratch.exact_build),
        diag_last_(scratch.exact_diag_last), count_(count),
        reach_(std::min(params.max_dist_x, params.max_dist_y)) {}

  // One nth_element split per level, alternating the key between q_end and the end
  // diagonal. The index tie-break makes the tree a function of the run list alone.
  void build() {
    ids_.resize(static_cast<size_t>(count_));
    leaf_.resize(static_cast<size_t>(count_));
    for (int32_t i = 0; i < count_; ++i)
      ids_[static_cast<size_t>(i)] = i;
    tree_.clear();
    tree_.reserve(static_cast<size_t>(2 * count_));
    fresh_node(-1);
    build_.clear();
    build_.push_back(ExactBuild{0, count_, 0, 0});
    while (!build_.empty()) {
      const ExactBuild task = build_.back();
      build_.pop_back();
      if (task.hi - task.lo == 1) {
        const int32_t run = ids_[static_cast<size_t>(task.lo)];
        tree_[static_cast<size_t>(task.node)].run = run;
        leaf_[static_cast<size_t>(run)] = task.node;
        continue;
      }
      const int32_t mid = task.lo + (task.hi - task.lo) / 2;
      const bool by_diagonal = (task.depth & 1) != 0;
      const std::vector<HotRun>& hot = hot_;
      std::nth_element(ids_.begin() + task.lo, ids_.begin() + mid,
                       ids_.begin() + task.hi,
                       [&hot, by_diagonal](int32_t x, int32_t y) {
                         const int64_t kx =
                             by_diagonal ? diagonal(hot[static_cast<size_t>(x)])
                                         : hot[static_cast<size_t>(x)].q_end;
                         const int64_t ky =
                             by_diagonal ? diagonal(hot[static_cast<size_t>(y)])
                                         : hot[static_cast<size_t>(y)].q_end;
                         return kx != ky ? kx < ky : x < y;
                       });
      const int32_t left = fresh_node(task.node);
      const int32_t right = fresh_node(task.node);
      tree_[static_cast<size_t>(task.node)].left = left;
      tree_[static_cast<size_t>(task.node)].right = right;
      build_.push_back(ExactBuild{mid, task.hi, right, task.depth + 1});
      build_.push_back(ExactBuild{task.lo, mid, left, task.depth + 1});
    }
    // One slot per end diagonal, and a diagonal cannot outnumber the runs.
    diag_last_.clear();
    diag_last_.reserve(static_cast<size_t>(count_));
  }

  // Run i joins the active set, leaf to root, right after it is scored. An unscored node
  // keeps max_index -1 and is rejected whole.
  void activate(int32_t i) {
    const HotRun& r = hot_[static_cast<size_t>(i)];
    const int64_t d = diagonal(r);
    // Runs activate in ascending i, so this keeps the largest index on the diagonal.
    diag_last_[d] = i;
    for (int32_t at = leaf_[static_cast<size_t>(i)]; at >= 0;
         at = tree_[static_cast<size_t>(at)].parent) {
      ExactNode& n = tree_[static_cast<size_t>(at)];
      // The sweep activates in ascending i, so the newest is the largest.
      n.max_index = i;
      if (r.score > n.score)
        n.score = r.score;
      if (r.q_begin < n.min_q_begin)
        n.min_q_begin = r.q_begin;
      if (r.r_begin < n.min_r_begin)
        n.min_r_begin = r.r_begin;
      if (r.q_end < n.min_q_end)
        n.min_q_end = r.q_end;
      if (r.r_end < n.min_r_end)
        n.min_r_end = r.r_end;
      if (r.q_end > n.max_q_end)
        n.max_q_end = r.q_end;
      if (d < n.min_d)
        n.min_d = d;
      if (d > n.max_d)
        n.max_d = d;
    }
  }

  // `best_score` enters at B's own weight with `best_predecessor` -1 (the fresh start) and
  // leaves as the exact maximum, with the predecessor the linear scan would elect.
  void search(int32_t i, int32_t& best_score, int32_t& best_predecessor) {
    const HotRun& b = hot_[static_cast<size_t>(i)];
    // Seed the incumbent with the last run on B's own begin diagonal, where dd is 0 and the
    // drift vanishes. The election rule is a total order on (score, index), with the fresh
    // start highest, so evaluation order cannot change the result; a better incumbent
    // only prunes more.
    const auto seed =
        diag_last_.find(static_cast<int64_t>(b.r_begin) - b.q_begin);
    if (seed != diag_last_.end())
      elect(b, seed->second, best_score, best_predecessor);
    // A node that can only tie the best is entered only for a predecessor above the elected
    // one, and never while the best is the fresh start.
    const auto blocked = [&](int64_t bound, int32_t node) {
      return bound < static_cast<int64_t>(best_score) ||
             (bound == static_cast<int64_t>(best_score) &&
              (best_predecessor < 0 ||
               tree_[static_cast<size_t>(node)].max_index <= best_predecessor));
    };
    // Descend into the better child and stack the other. A child already blocked is never
    // stacked, since the incumbent only improves.
    stack_.clear();
    ++nodes_;
    ExactVisit visit{upper(0, b), 0};
    if (blocked(visit.bound, 0))
      return;
    for (;;) {
      const ExactNode& n = tree_[static_cast<size_t>(visit.node)];
      if (n.run >= 0) {
        elect(b, n.run, best_score, best_predecessor);
      } else {
        // The higher bound first; on a tie, the larger max index.
        int32_t first = n.left;
        int32_t second = n.right;
        ++nodes_;
        int64_t first_bound = upper(first, b);
        ++nodes_;
        int64_t second_bound = upper(second, b);
        if (first_bound < second_bound ||
            (first_bound == second_bound &&
             tree_[static_cast<size_t>(first)].max_index <
                 tree_[static_cast<size_t>(second)].max_index)) {
          std::swap(first, second);
          std::swap(first_bound, second_bound);
        }
        if (!blocked(second_bound, second))
          stack_.push_back(ExactVisit{second_bound, second});
        if (!blocked(first_bound, first)) {
          visit = ExactVisit{first_bound, first};
          continue;
        }
      }
      // Pop stacked siblings, re-testing each against the improved incumbent.
      for (;;) {
        if (stack_.empty())
          return;
        visit = stack_.back();
        stack_.pop_back();
        if (!blocked(visit.bound, visit.node))
          break;
      }
    }
  }

  int64_t nodes() const { return nodes_; }
  int64_t steps() const { return steps_; }

private:
  static constexpr int64_t kImpossible = INT64_MIN;

  // Prices one candidate with dense_step and elects it by the linear scan's rule: a higher
  // score wins, an equal score wins only from a higher index, and a tie never displaces
  // the fresh start.
  void elect(const HotRun& b, int32_t candidate, int32_t& best_score,
             int32_t& best_predecessor) {
    ++steps_;
    const HotRun& a = hot_[static_cast<size_t>(candidate)];
    const int32_t step = dense_step(b, a, params_);
    if (step == kNegInf)
      return;
    const int32_t score = step + a.score;
    if (score > best_score || (score == best_score && best_predecessor >= 0 &&
                               candidate > best_predecessor)) {
      best_score = score;
      best_predecessor = candidate;
    }
  }

  // dense_step's dd is |(B.r_begin - B.q_begin) - (A.r_end - A.q_end)|, so a
  // predecessor's key is its END diagonal and a successor's its BEGIN one.
  static int64_t diagonal(const HotRun& r) {
    return static_cast<int64_t>(r.r_end) - r.q_end;
  }

  int32_t fresh_node(int32_t parent) {
    ExactNode node;
    node.left = -1;
    node.right = -1;
    node.parent = parent;
    node.run = -1;
    node.max_index = -1;
    node.score = kNegInf;
    node.min_q_begin = INT32_MAX;
    node.min_r_begin = INT32_MAX;
    node.min_q_end = INT32_MAX;
    node.min_r_end = INT32_MAX;
    node.max_q_end = INT32_MIN;
    node.min_d = INT64_MAX;
    node.max_d = INT64_MIN;
    tree_.push_back(node);
    return static_cast<int32_t>(tree_.size()) - 1;
  }

  // Upper bound over a node's scored members. Each rejection is one of dense_step's tests
  // applied to the member extremum most likely to pass it.
  int64_t upper(int32_t node, const HotRun& b) const {
    const ExactNode& n = tree_[static_cast<size_t>(node)];
    if (n.max_index < 0)
      return kImpossible;
    // Forward starts and forward ends, against the members closest to passing.
    if (b.q_begin <= n.min_q_begin || b.r_begin <= n.min_r_begin ||
        b.q_end <= n.min_q_end || b.r_end <= n.min_r_end)
      return kImpossible;
    // dense_step's reach test compares dq with both max distances: a floor on q_end.
    if (static_cast<int64_t>(n.max_q_end) <
        static_cast<int64_t>(b.q_begin) + params_.span - reach_)
      return kImpossible;
    const int64_t d = static_cast<int64_t>(b.r_begin) - b.q_begin;
    const int64_t dd =
        d < n.min_d ? n.min_d - d : (d > n.max_d ? d - n.max_d : 0);
    if (dd > params_.bw)
      return kImpossible;
    // The gain at the smallest overlap the node's minima allow; dense_step's gain is
    // non-increasing in the overlap.
    const int64_t overlap = std::max<int64_t>(
        0, std::max<int64_t>(static_cast<int64_t>(n.min_q_end) - b.q_begin,
                             static_cast<int64_t>(n.min_r_end) - b.r_begin));
    const int64_t length = static_cast<int64_t>(b.q_end) - b.q_begin;
    int64_t gain;
    // No overlap, the common case: the gain is B's weight, without the division.
    if (overlap == 0 && length > 0) {
      gain = b.w;
    } else {
      const int64_t trimmed =
          overlap >= length
              ? 0
              : (static_cast<int64_t>(b.w) * (length - overlap)) / length;
      gain = std::min<int64_t>(trimmed,
                               static_cast<int64_t>(b.w) -
                                   std::min<int64_t>(overlap, params_.span));
    }
    if (gain < 0)
      gain = 0;
    // dense_step's own float multiply, so the bound cannot floor above the step; the
    // omitted skip and log terms are non-negative.
    const int64_t drift =
        static_cast<int64_t>(params_.chn_pen_gap * static_cast<float>(dd));
    const int64_t bound = static_cast<int64_t>(n.score) + gain;
    return bound - drift;
  }

  const std::vector<HotRun>& hot_;
  const DenseChainParams& params_;
  std::vector<ExactNode>& tree_;
  std::vector<int32_t>& ids_;
  std::vector<int32_t>& leaf_;
  std::vector<ExactVisit>& stack_;
  std::vector<ExactBuild>& build_;
  ::fa::cpu::lr::FlatInt64Map<int32_t>& diag_last_;
  const int32_t count_;
  const int64_t reach_;
  int64_t nodes_ = 0;
  int64_t steps_ = 0;
};

#if FA_DENSE_EXACT_WIDE
#if FA_DENSE_EXACT_WIDE_NEON
// NEON is in the AArch64 baseline, so the lane functions are compiled for the
// target the whole unit is compiled for and need no attribute.
#define FA_WIDE_TARGET
#else
#define FA_WIDE_TARGET __attribute__((target("avx512f")))
#endif

#if FA_DENSE_EXACT_WIDE_NEON
// NEON lanes: a 16-lane quantity is four int32x4_t, quarter q holding lanes 4q..4q+3.
// Masks stay in vector form and fold to one bit per lane at the end. Kept lanes are
// visited in ascending lane order, the order AVX-512's compress-store writes them in, so
// both lane sets follow the same search trajectory.

// Bit 4*l of the result is lane l's mask bit, so a lane's bit position is
// 4 * its lane index and the mask has at most one bit per nibble.
constexpr uint64_t kWideLaneBits = 0x1111111111111111ull;

// Two narrowing steps put the sixteen lanes' low bytes in one 16-byte vector
// (byte l = lane l, 0x00 or 0xff), and the shift-and-narrow folds each byte
// PAIR into a byte of two nibbles, landing lane l in nibble l of a 64-bit
// word; the AND leaves one bit per lane.
inline uint64_t wide_lane_bits(uint32x4_t m0, uint32x4_t m1, uint32x4_t m2,
                               uint32x4_t m3) {
  const uint16x8_t lo = vcombine_u16(vmovn_u32(m0), vmovn_u32(m1));
  const uint16x8_t hi = vcombine_u16(vmovn_u32(m2), vmovn_u32(m3));
  const uint8x16_t bytes = vcombine_u8(vmovn_u16(lo), vmovn_u16(hi));
  const uint8x8_t nibbles = vshrn_n_u16(vreinterpretq_u16_u8(bytes), 4);
  return vget_lane_u64(vreinterpret_u64_u8(nibbles), 0) & kWideLaneBits;
}

// mg_log2 over four lanes in the scalar operation order with no fused multiply-add, so the
// lanes equal the scalar floats. The exponent uses a logical shift, as mg_log2's unsigned
// `z.i >> 23`.
inline float32x4_t wide_log2_ps(float32x4_t x) {
  uint32x4_t bits = vreinterpretq_u32_f32(x);
  const float32x4_t exponent =
      vsubq_f32(vcvtq_f32_s32(vreinterpretq_s32_u32(
                    vandq_u32(vshrq_n_u32(bits, 23), vdupq_n_u32(255)))),
                vdupq_n_f32(128.0f));
  bits = vandq_u32(bits, vdupq_n_u32(~(255u << 23)));
  bits = vaddq_u32(bits, vdupq_n_u32(127u << 23));
  const float32x4_t mantissa = vreinterpretq_f32_u32(bits);
  float32x4_t poly = vmulq_f32(vdupq_n_f32(-0.34484843f), mantissa);
  poly = vaddq_f32(poly, vdupq_n_f32(2.02466578f));
  poly = vmulq_f32(poly, mantissa);
  poly = vsubq_f32(poly, vdupq_n_f32(0.67487759f));
  return vaddq_f32(exponent, poly);
}

// dense_step's drift at four lanes of dd >= 0 without the skip term (>= 0, so the bound
// stays above the step): (int)(gap * dd [+ 0.5 * log2(dd + 1)]), and 0 at dd == 0.
// vcvtq_s32_f32 truncates like dense_step's cast; gap * bw keeps the drift in int32 range.
inline int32x4_t wide_drift_lanes(int32x4_t dd, float gap, bool log_term) {
  float32x4_t drift = vmulq_f32(vdupq_n_f32(gap), vcvtq_f32_s32(dd));
  if (log_term) {
    const float32x4_t log_penalty =
        wide_log2_ps(vcvtq_f32_s32(vaddq_s32(dd, vdupq_n_s32(1))));
    drift = vaddq_f32(drift, vmulq_f32(vdupq_n_f32(0.5f), log_penalty));
  }
  return vbicq_s32(vcvtq_s32_f32(drift),
                   vreinterpretq_s32_u32(vceqzq_s32(dd)));
}

// dense_step's gain at four lanes of overlap >= 0: min(floor(w * (L - o) / L),
// w - min(o, span)), clamped at 0, with the quotient floored in double precision as in the
// AVX-512 lanes. The clamp to +-2^30 hands vcvtq_s64_f64 an in-range integer.
inline int32x4_t wide_gain_lanes(int32x4_t overlap, int32_t w, int64_t length,
                                 int32_t span) {
  const float64x2_t wd = vdupq_n_f64(static_cast<double>(w));
  const float64x2_t ld = vdupq_n_f64(static_cast<double>(length));
  const float64x2_t lo = vcvtq_f64_s64(vmovl_s32(vget_low_s32(overlap)));
  const float64x2_t hi = vcvtq_f64_s64(vmovl_s32(vget_high_s32(overlap)));
  const float64x2_t qlo =
      vrndmq_f64(vdivq_f64(vmulq_f64(wd, vsubq_f64(ld, lo)), ld));
  const float64x2_t qhi =
      vrndmq_f64(vdivq_f64(vmulq_f64(wd, vsubq_f64(ld, hi)), ld));
  const float64x2_t cap = vdupq_n_f64(1073741824.0);
  const float64x2_t floor_cap = vdupq_n_f64(-1073741824.0);
  const int32x2_t tlo =
      vmovn_s64(vcvtq_s64_f64(vminq_f64(cap, vmaxq_f64(floor_cap, qlo))));
  const int32x2_t thi =
      vmovn_s64(vcvtq_s64_f64(vminq_f64(cap, vmaxq_f64(floor_cap, qhi))));
  const int32x4_t trimmed = vcombine_s32(tlo, thi);
  const int32x4_t head_loss =
      vsubq_s32(vdupq_n_s32(w), vminq_s32(overlap, vdupq_n_s32(span)));
  return vmaxq_s32(vminq_s32(trimmed, head_loss), vdupq_n_s32(0));
}

#else

// mg_log2 over 16 lanes in the scalar operation order; with no FP contraction in this file
// the lanes equal the scalar floats.
FA_WIDE_TARGET inline __m512 wide_log2_ps(__m512 x) {
  __m512i bits = _mm512_castps_si512(x);
  const __m512 exponent =
      _mm512_sub_ps(_mm512_cvtepi32_ps(_mm512_and_si512(
                        _mm512_srli_epi32(bits, 23), _mm512_set1_epi32(255))),
                    _mm512_set1_ps(128.0f));
  bits = _mm512_and_si512(
      bits, _mm512_set1_epi32(static_cast<int32_t>(~(255u << 23))));
  bits = _mm512_add_epi32(bits, _mm512_set1_epi32(127 << 23));
  const __m512 mantissa = _mm512_castsi512_ps(bits);
  __m512 poly = _mm512_mul_ps(_mm512_set1_ps(-0.34484843f), mantissa);
  poly = _mm512_add_ps(poly, _mm512_set1_ps(2.02466578f));
  poly = _mm512_mul_ps(poly, mantissa);
  poly = _mm512_sub_ps(poly, _mm512_set1_ps(0.67487759f));
  return _mm512_add_ps(exponent, poly);
}

// dense_step's drift at 16 lanes of dd >= 0, the skip term left out (it is
// >= 0, so leaving it out keeps the bound above the step): (int)(gap * dd
// [+ 0.5 * log2(dd + 1)]), and 0 at dd == 0 as dense_step's branch has it.
FA_WIDE_TARGET inline __m512i wide_drift_lanes(__m512i dd, float gap,
                                               bool log_term) {
  __m512 drift = _mm512_mul_ps(_mm512_set1_ps(gap), _mm512_cvtepi32_ps(dd));
  if (log_term) {
    const __m512 log_penalty = wide_log2_ps(
        _mm512_cvtepi32_ps(_mm512_add_epi32(dd, _mm512_set1_epi32(1))));
    drift =
        _mm512_add_ps(drift, _mm512_mul_ps(_mm512_set1_ps(0.5f), log_penalty));
  }
  return _mm512_maskz_mov_epi32(
      _mm512_cmpneq_epi32_mask(dd, _mm512_setzero_si512()),
      _mm512_cvttps_epi32(drift));
}

// dense_step's gain at 16 lanes of overlap >= 0: min(floor(w * (L - o) / L),
// w - min(o, span)), clamped at 0. floor(fl(n / L)) equals n / L for
// n < 2^53 (a non-integer quotient sits at least 1/L from an integer and the
// rounding error is below 2^-53 * n / L), and n = w * (L - o) < 2^53 inside
// kWideLengthLimit. An overlap at or past L gives a non-positive quotient and
// the clamp, which is dense_step's `trimmed = 0` branch.
FA_WIDE_TARGET inline __m512i wide_gain_lanes(__m512i overlap, int32_t w,
                                              int64_t length, int32_t span) {
  const __m512d wd = _mm512_set1_pd(static_cast<double>(w));
  const __m512d ld = _mm512_set1_pd(static_cast<double>(length));
  const __m512d lo = _mm512_cvtepi32_pd(_mm512_castsi512_si256(overlap));
  const __m512d hi = _mm512_cvtepi32_pd(_mm512_extracti64x4_epi64(overlap, 1));
  const __m512d qlo = _mm512_roundscale_pd(
      _mm512_div_pd(_mm512_mul_pd(wd, _mm512_sub_pd(ld, lo)), ld),
      _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  const __m512d qhi = _mm512_roundscale_pd(
      _mm512_div_pd(_mm512_mul_pd(wd, _mm512_sub_pd(ld, hi)), ld),
      _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  // Clamped before the int conversion so no lane can overflow it.
  const __m512d cap = _mm512_set1_pd(1073741824.0);
  const __m512d floor_cap = _mm512_set1_pd(-1073741824.0);
  const __m256i tlo =
      _mm512_cvtpd_epi32(_mm512_min_pd(cap, _mm512_max_pd(floor_cap, qlo)));
  const __m256i thi =
      _mm512_cvtpd_epi32(_mm512_min_pd(cap, _mm512_max_pd(floor_cap, qhi)));
  const __m512i trimmed =
      _mm512_inserti64x4(_mm512_castsi256_si512(tlo), thi, 1);
  const __m512i head_loss = _mm512_sub_epi32(
      _mm512_set1_epi32(w), _mm512_min_epi32(overlap, _mm512_set1_epi32(span)));
  return _mm512_max_epi32(_mm512_min_epi32(trimmed, head_loss),
                          _mm512_setzero_si512());
}

#endif // FA_DENSE_EXACT_WIDE_NEON

class ExactWideArm {
public:
  ExactWideArm(const std::vector<HotRun>& hot, const DenseChainParams& params,
               DenseChainScratch& scratch, int32_t count)
      : hot_(hot), params_(params), scratch_(scratch),
        tree_(scratch.wide_nodes), leaves_(scratch.wide_leaves),
        run_leaf_(scratch.wide_run_leaf), keys_(scratch.wide_keys),
        ids_(scratch.wide_ids), tasks_(scratch.wide_tasks),
        cuts_(scratch.wide_cuts), stack_(scratch.wide_stack),
        diag_last_(scratch.exact_diag_last), count_(count),
        reach_(std::min(params.max_dist_x, params.max_dist_y)) {}

  static bool supported() {
#if FA_DENSE_EXACT_WIDE_NEON
    // ASIMD is mandatory on AArch64.
    return true;
#else
    static const bool ok = __builtin_cpu_supports("avx512f") != 0;
    return ok;
#endif
  }

  // The pool the lanes can represent: see kWideCoordinateLimit and kWideLengthLimit.
  bool in_domain() const {
    if (params_.bw >= kWideCoordinateLimit ||
        params_.span >= kWideCoordinateLimit ||
        reach_ >= kWideCoordinateLimit || params_.bw < 0 || params_.span < 0 ||
        reach_ < 0)
      return false;
    int64_t sum_w = 0;
    for (int32_t i = 0; i < count_; ++i) {
      const HotRun& h = hot_[static_cast<size_t>(i)];
      if (h.q_begin < 0 || h.r_begin < 0 || h.q_end >= kWideCoordinateLimit ||
          h.r_end >= kWideCoordinateLimit || h.w < 0 ||
          h.q_end - h.q_begin >= kWideLengthLimit || h.q_end <= h.q_begin)
        return false;
      sum_w += h.w;
    }
    return sum_w < kWideCoordinateLimit;
  }

  void build() {
    log_term_ = drift_monotone();
    ids_.resize(static_cast<size_t>(count_));
    for (int32_t i = 0; i < count_; ++i)
      ids_[static_cast<size_t>(i)] = i;
    tree_.clear();
    leaves_.clear();
    run_leaf_.assign(static_cast<size_t>(count_), -1);
    diag_last_.clear();
    diag_last_.reserve(static_cast<size_t>(count_));
    tree_.push_back(fresh_node(-1, -1));
    tasks_.clear();
    // The root takes total_levels mod 4 median levels (or 4) and every node below it four,
    // so a leaf block holds 8 to 16 runs instead of one or two when the pool is just above
    // a power of sixteen.
    int total_levels = 0;
    while ((static_cast<int64_t>(kWideLanes) << total_levels) < count_)
      ++total_levels;
    const int root_levels = total_levels % 4 == 0 ? 4 : total_levels % 4;
    tasks_.push_back(WideTask{0, count_, 0, 0});
    while (!tasks_.empty()) {
      const WideTask t = tasks_.back();
      tasks_.pop_back();
      const int32_t size = t.hi - t.lo;
      // Gather both keys once, then make this node's median levels over them: ExactArm's
      // k-d cut over the same depths.
      keys_.resize(static_cast<size_t>(size));
      for (int32_t at = 0; at < size; ++at) {
        const int32_t id = ids_[static_cast<size_t>(t.lo + at)];
        const HotRun& h = hot_[static_cast<size_t>(id)];
        keys_[static_cast<size_t>(at)] =
            WideKey{h.q_end, h.r_end - h.q_end, id};
      }
      cuts_.clear();
      const int levels = t.node == 0 ? root_levels : 4;
      partition(0, size, t.depth, levels);
      cuts_.push_back(size);
      for (int32_t at = 0; at < size; ++at)
        ids_[static_cast<size_t>(t.lo + at)] =
            keys_[static_cast<size_t>(at)].id;
      if (size <= (kWideLanes << levels)) {
        // `levels` median levels over at most 16 << levels runs leave ranges
        // of at most 16: one leaf block each.
        const int32_t first = static_cast<int32_t>(leaves_.size());
        int32_t child = 0;
        for (size_t k = 0; k + 1 < cuts_.size(); ++k) {
          const int32_t lo = cuts_[k], hi = cuts_[k + 1];
          if (hi <= lo)
            continue;
          WideLeaf leaf;
          for (int l = 0; l < kWideLanes; ++l) {
            leaf.q_begin[l] = 0;
            leaf.q_end[l] = 0;
            leaf.r_begin[l] = 0;
            leaf.r_end[l] = 0;
            leaf.w[l] = 0;
            leaf.score[l] = kNegInf;
            leaf.index[l] = -1;
            leaf.d_end[l] = 0;
          }
          leaf.parent = t.node;
          leaf.slot = child;
          for (int32_t at = lo; at < hi; ++at) {
            const int32_t run = ids_[static_cast<size_t>(t.lo + at)];
            const HotRun& h = hot_[static_cast<size_t>(run)];
            const int l = at - lo;
            leaf.q_begin[l] = h.q_begin;
            leaf.q_end[l] = h.q_end;
            leaf.r_begin[l] = h.r_begin;
            leaf.r_end[l] = h.r_end;
            leaf.w[l] = h.w;
            leaf.index[l] = run;
            leaf.d_end[l] = h.r_end - h.q_end;
            run_leaf_[static_cast<size_t>(run)] =
                static_cast<int32_t>(leaves_.size()) * kWideLanes + l;
          }
          leaves_.push_back(leaf);
          ++child;
        }
        WideNode& n = tree_[static_cast<size_t>(t.node)];
        n.leaf_children = 1;
        n.first_child = first;
        n.child_count = child;
      } else {
        const int32_t first = static_cast<int32_t>(tree_.size());
        int32_t child = 0;
        for (size_t k = 0; k + 1 < cuts_.size(); ++k) {
          const int32_t lo = cuts_[k], hi = cuts_[k + 1];
          if (hi <= lo)
            continue;
          tree_.push_back(fresh_node(t.node, child));
          ++child;
        }
        WideNode& n = tree_[static_cast<size_t>(t.node)];
        n.leaf_children = 0;
        n.first_child = first;
        n.child_count = child;
        // The children's ranges, in reverse so the left one is built first.
        int32_t c = child;
        for (size_t k = cuts_.size() - 1; k > 0; --k) {
          const int32_t lo = cuts_[k - 1], hi = cuts_[k];
          if (hi <= lo)
            continue;
          --c;
          tasks_.push_back(
              WideTask{t.lo + lo, t.lo + hi, t.depth + levels, first + c});
        }
      }
    }
  }

  // Run i joins the active set: its leaf lane takes its score, and its lane
  // in every ancestor node folds it in. At most one lane per level.
  void activate(int32_t i) {
    const HotRun& r = hot_[static_cast<size_t>(i)];
    const int32_t d = r.r_end - r.q_end;
    diag_last_[static_cast<int64_t>(d)] = i;
    const int32_t packed = run_leaf_[static_cast<size_t>(i)];
    WideLeaf& leaf = leaves_[static_cast<size_t>(packed / kWideLanes)];
    leaf.score[packed % kWideLanes] = r.score;
    int32_t node = leaf.parent;
    int32_t slot = leaf.slot;
    while (node >= 0) {
      WideNode& n = tree_[static_cast<size_t>(node)];
      n.max_index[slot] = i;
      if (r.score > n.score[slot])
        n.score[slot] = r.score;
      if (r.q_begin < n.min_q_begin[slot])
        n.min_q_begin[slot] = r.q_begin;
      if (r.r_begin < n.min_r_begin[slot])
        n.min_r_begin[slot] = r.r_begin;
      if (r.q_end < n.min_q_end[slot])
        n.min_q_end[slot] = r.q_end;
      if (r.r_end < n.min_r_end[slot])
        n.min_r_end[slot] = r.r_end;
      if (r.q_end > n.max_q_end[slot])
        n.max_q_end[slot] = r.q_end;
      if (d < n.min_d[slot])
        n.min_d[slot] = d;
      if (d > n.max_d[slot])
        n.max_d[slot] = d;
      slot = n.slot;
      node = n.parent;
    }
  }

  // As ExactArm::search, sixteen children per pass; a stacked child is re-tested against
  // the incumbent when popped.
  void search(int32_t i, int32_t& best_score, int32_t& best_predecessor) {
    const HotRun& b = hot_[static_cast<size_t>(i)];
    const auto seed =
        diag_last_.find(static_cast<int64_t>(b.r_begin) - b.q_begin);
    if (seed != diag_last_.end())
      elect(b, seed->second, best_score, best_predecessor);
    stack_.clear();
    expand(0, b, best_score, best_predecessor);
    while (!stack_.empty()) {
      const WideVisit v = stack_.back();
      stack_.pop_back();
      if (v.bound < best_score ||
          (v.bound == best_score &&
           (best_predecessor < 0 || v.max_index <= best_predecessor)))
        continue;
      if (v.is_leaf)
        scan_leaf(v.child, i, b, best_score, best_predecessor);
      else
        expand(v.child, b, best_score, best_predecessor);
    }
  }

  // Lane bounds computed (16 per pass) and dense_step calls, the wide arm's
  // reading of DenseRunResult::exact_nodes and exact_steps.
  int64_t bounds() const { return bounds_; }
  int64_t steps() const { return steps_; }
  bool log_term() const { return log_term_; }

private:
  // dense_step's drift at dd with no skip term, in its float order. The bound charges it at
  // a node's smallest dd, which is sound only if it is non-decreasing in dd; drift_monotone
  // checks that over [0, bw], once per (gap, bw) and thread.
  int32_t drift_at(int32_t dd) const {
    if (dd == 0)
      return 0;
    const float linear = params_.chn_pen_gap * static_cast<float>(dd) + 0.0f;
    const float log_penalty = mg_log2(static_cast<float>(dd + 1));
    return static_cast<int32_t>(linear + 0.5f * log_penalty);
  }

  bool drift_monotone() {
    if (scratch_.wide_gap == params_.chn_pen_gap &&
        scratch_.wide_bw == params_.bw)
      return scratch_.wide_log_term;
    bool monotone = params_.chn_pen_gap >= 0.0f;
    int32_t previous = 0;
    for (int32_t dd = 1; monotone && dd <= params_.bw; ++dd) {
      const int32_t current = drift_at(dd);
      if (current < previous)
        monotone = false;
      previous = current;
    }
    scratch_.wide_gap = params_.chn_pen_gap;
    scratch_.wide_bw = params_.bw;
    scratch_.wide_log_term = monotone;
    return monotone;
  }

  // Median splits over keys_[lo, hi), `levels` deep, the key alternating with the binary
  // depth as in ExactArm::build (q_end at even depths, the end diagonal at odd), ties broken
  // by id. Appends every range start to cuts_ in position order.
  void partition(int32_t lo, int32_t hi, int32_t depth, int levels) {
    if (levels == 0 || hi - lo <= 1) {
      cuts_.push_back(lo);
      return;
    }
    const int32_t mid = lo + (hi - lo) / 2;
    const bool by_diagonal = (depth & 1) != 0;
    std::nth_element(keys_.begin() + lo, keys_.begin() + mid,
                     keys_.begin() + hi,
                     [by_diagonal](const WideKey& x, const WideKey& y) {
                       const int32_t kx = by_diagonal ? x.d_end : x.q_end;
                       const int32_t ky = by_diagonal ? y.d_end : y.q_end;
                       return kx != ky ? kx < ky : x.id < y.id;
                     });
    partition(lo, mid, depth + 1, levels - 1);
    partition(mid, hi, depth + 1, levels - 1);
  }

  WideNode fresh_node(int32_t parent, int32_t slot) {
    WideNode n;
    for (int c = 0; c < kWideLanes; ++c) {
      n.max_index[c] = -1;
      n.score[c] = kNegInf;
      n.min_q_begin[c] = kWideCoordinateLimit;
      n.min_r_begin[c] = kWideCoordinateLimit;
      n.min_q_end[c] = kWideCoordinateLimit;
      n.min_r_end[c] = kWideCoordinateLimit;
      n.max_q_end[c] = -kWideCoordinateLimit;
      n.min_d[c] = kWideCoordinateLimit;
      n.max_d[c] = -kWideCoordinateLimit;
    }
    n.first_child = -1;
    n.child_count = 0;
    n.parent = parent;
    n.slot = slot;
    n.leaf_children = 0;
    return n;
  }

  // As ExactArm::elect.
  void elect(const HotRun& b, int32_t candidate, int32_t& best_score,
             int32_t& best_predecessor) {
    ++steps_;
    const HotRun& a = hot_[static_cast<size_t>(candidate)];
    const int32_t step = dense_step(b, a, params_);
    if (step == kNegInf)
      return;
    const int32_t score = step + a.score;
    if (score > best_score || (score == best_score && best_predecessor >= 0 &&
                               candidate > best_predecessor)) {
      best_score = score;
      best_predecessor = candidate;
    }
  }

  // ExactArm::upper over a node's sixteen children at once, with the drift charged at each
  // child's smallest dd (log term included when monotone). Unblocked children are pushed
  // in ascending (bound, max index) order, so the best is popped first.
#if FA_DENSE_EXACT_WIDE_NEON
  FA_WIDE_TARGET void expand(int32_t node, const HotRun& b, int32_t best_score,
                             int32_t best_predecessor) {
    bounds_ += kWideLanes;
    const WideNode& n = tree_[static_cast<size_t>(node)];
    const int32_t d_b = b.r_begin - b.q_begin;
    const int32_t reach_floor = static_cast<int32_t>(
        static_cast<int64_t>(b.q_begin) + params_.span - reach_);
    const int64_t length = static_cast<int64_t>(b.q_end) - b.q_begin;
    const int32x4_t vqb = vdupq_n_s32(b.q_begin);
    const int32x4_t vrb = vdupq_n_s32(b.r_begin);
    const int32x4_t vqe = vdupq_n_s32(b.q_end);
    const int32x4_t vre = vdupq_n_s32(b.r_end);
    const int32x4_t vd = vdupq_n_s32(d_b);
    const int32x4_t vreach = vdupq_n_s32(reach_floor);
    const int32x4_t vbw = vdupq_n_s32(params_.bw);
    const int32x4_t zero = vdupq_n_s32(0);
    uint32x4_t impossible[4];
    int32x4_t dd[4];
    int32x4_t overlap[4];
    uint32x4_t overlapping = vdupq_n_u32(0);
    for (int q = 0; q < 4; ++q) {
      const int at = 4 * q;
      const int32x4_t max_index = vld1q_s32(n.max_index + at);
      const int32x4_t min_q_begin = vld1q_s32(n.min_q_begin + at);
      const int32x4_t min_r_begin = vld1q_s32(n.min_r_begin + at);
      const int32x4_t min_q_end = vld1q_s32(n.min_q_end + at);
      const int32x4_t min_r_end = vld1q_s32(n.min_r_end + at);
      const int32x4_t max_q_end = vld1q_s32(n.max_q_end + at);
      const int32x4_t min_d = vld1q_s32(n.min_d + at);
      const int32x4_t max_d = vld1q_s32(n.max_d + at);
      // The AVX-512 predicates, with the same operand order.
      uint32x4_t imp = vcltq_s32(max_index, zero);
      imp = vorrq_u32(imp, vcleq_s32(vqb, min_q_begin));
      imp = vorrq_u32(imp, vcleq_s32(vrb, min_r_begin));
      imp = vorrq_u32(imp, vcleq_s32(vqe, min_q_end));
      imp = vorrq_u32(imp, vcleq_s32(vre, min_r_end));
      imp = vorrq_u32(imp, vcltq_s32(max_q_end, vreach));
      const int32x4_t dd_q = vmaxq_s32(
          vmaxq_s32(vsubq_s32(min_d, vd), vsubq_s32(vd, max_d)), zero);
      imp = vorrq_u32(imp, vcgtq_s32(dd_q, vbw));
      const int32x4_t overlap_q = vmaxq_s32(
          vmaxq_s32(vsubq_s32(min_q_end, vqb), vsubq_s32(min_r_end, vrb)),
          zero);
      // vtstq_s32(v, v) marks the nonzero lanes, as AVX-512's test-to-mask.
      overlapping = vorrq_u32(
          overlapping, vbicq_u32(vtstq_s32(overlap_q, overlap_q), imp));
      impossible[q] = imp;
      dd[q] = dd_q;
      overlap[q] = overlap_q;
    }
    // No admissible child overlaps B in the common case: the gain is then B's weight.
    const bool no_overlap = vmaxvq_u32(overlapping) == 0;
    const int32x4_t vw = vdupq_n_s32(b.w);
    const int32x4_t vbest = vdupq_n_s32(best_score);
    const int32x4_t vbestpred = vdupq_n_s32(best_predecessor);
    alignas(64) int32_t bound_lane[kWideLanes];
    uint32x4_t dropped[4];
    for (int q = 0; q < 4; ++q) {
      const int at = 4 * q;
      const int32x4_t gain =
          no_overlap ? vw
                     : wide_gain_lanes(overlap[q], b.w, length, params_.span);
      const int32x4_t drift =
          wide_drift_lanes(dd[q], params_.chn_pen_gap, log_term_);
      const int32x4_t bound =
          vsubq_s32(vaddq_s32(vld1q_s32(n.score + at), gain), drift);
      // Only the bound is read back; max index and child id come from the node.
      vst1q_s32(bound_lane + at, bound);
      uint32x4_t blk = vcltq_s32(bound, vbest);
      const uint32x4_t equal = vceqq_s32(bound, vbest);
      if (best_predecessor < 0)
        blk = vorrq_u32(blk, equal);
      else
        blk = vorrq_u32(
            blk, vandq_u32(equal, vcleq_s32(vld1q_s32(n.max_index + at),
                                            vbestpred)));
      // Fold both drop masks so the bit form is taken once.
      dropped[q] = vorrq_u32(impossible[q], blk);
    }
    // The valid children are the low child_count lanes; at 16 the shift would be by 64.
    const uint64_t valid =
        n.child_count >= kWideLanes
            ? kWideLaneBits
            : (kWideLaneBits & ((1ull << (4 * n.child_count)) - 1ull));
    uint64_t keep =
        valid & ~wide_lane_bits(dropped[0], dropped[1], dropped[2], dropped[3]);
    if (keep == 0)
      return;
    WideVisit order[kWideLanes];
    int kept = 0;
    // Ascending lane order, as AVX-512's compress-store.
    while (keep != 0) {
      const int lane = __builtin_ctzll(keep) >> 2;
      keep &= keep - 1;
      const WideVisit v{bound_lane[lane], n.max_index[lane],
                        n.first_child + lane, n.leaf_children};
      int p = kept;
      while (p > 0 && (order[p - 1].bound > v.bound ||
                       (order[p - 1].bound == v.bound &&
                        order[p - 1].max_index > v.max_index))) {
        order[p] = order[p - 1];
        --p;
      }
      order[p] = v;
      ++kept;
    }
    for (int c = 0; c < kept; ++c)
      stack_.push_back(order[c]);
  }
#else
  FA_WIDE_TARGET void expand(int32_t node, const HotRun& b, int32_t best_score,
                             int32_t best_predecessor) {
    bounds_ += kWideLanes;
    const WideNode& n = tree_[static_cast<size_t>(node)];
    const int32_t d_b = b.r_begin - b.q_begin;
    const int32_t reach_floor = static_cast<int32_t>(
        static_cast<int64_t>(b.q_begin) + params_.span - reach_);
    const int64_t length = static_cast<int64_t>(b.q_end) - b.q_begin;
    const __m512i max_index = _mm512_load_si512(n.max_index);
    const __m512i score = _mm512_load_si512(n.score);
    const __m512i min_q_begin = _mm512_load_si512(n.min_q_begin);
    const __m512i min_r_begin = _mm512_load_si512(n.min_r_begin);
    const __m512i min_q_end = _mm512_load_si512(n.min_q_end);
    const __m512i min_r_end = _mm512_load_si512(n.min_r_end);
    const __m512i max_q_end = _mm512_load_si512(n.max_q_end);
    const __m512i min_d = _mm512_load_si512(n.min_d);
    const __m512i max_d = _mm512_load_si512(n.max_d);
    const __m512i vqb = _mm512_set1_epi32(b.q_begin);
    const __m512i vrb = _mm512_set1_epi32(b.r_begin);
    const __m512i vqe = _mm512_set1_epi32(b.q_end);
    const __m512i vre = _mm512_set1_epi32(b.r_end);
    const __m512i vd = _mm512_set1_epi32(d_b);
    const __m512i zero = _mm512_setzero_si512();
    __mmask16 impossible = _mm512_cmplt_epi32_mask(max_index, zero);
    impossible |= _mm512_cmple_epi32_mask(vqb, min_q_begin);
    impossible |= _mm512_cmple_epi32_mask(vrb, min_r_begin);
    impossible |= _mm512_cmple_epi32_mask(vqe, min_q_end);
    impossible |= _mm512_cmple_epi32_mask(vre, min_r_end);
    impossible |=
        _mm512_cmplt_epi32_mask(max_q_end, _mm512_set1_epi32(reach_floor));
    const __m512i dd =
        _mm512_max_epi32(_mm512_max_epi32(_mm512_sub_epi32(min_d, vd),
                                          _mm512_sub_epi32(vd, max_d)),
                         zero);
    impossible |= _mm512_cmpgt_epi32_mask(dd, _mm512_set1_epi32(params_.bw));
    const __m512i overlap =
        _mm512_max_epi32(_mm512_max_epi32(_mm512_sub_epi32(min_q_end, vqb),
                                          _mm512_sub_epi32(min_r_end, vrb)),
                         zero);
    // No admissible child overlaps B in the common case: the gain is then B's weight.
    __m512i gain;
    if ((_mm512_test_epi32_mask(overlap, overlap) & ~impossible) == 0)
      gain = _mm512_set1_epi32(b.w);
    else
      gain = wide_gain_lanes(overlap, b.w, length, params_.span);
    const __m512i drift = wide_drift_lanes(dd, params_.chn_pen_gap, log_term_);
    const __m512i bound =
        _mm512_sub_epi32(_mm512_add_epi32(score, gain), drift);
    const __m512i vbest = _mm512_set1_epi32(best_score);
    __mmask16 blocked = _mm512_cmplt_epi32_mask(bound, vbest);
    const __mmask16 equal = _mm512_cmpeq_epi32_mask(bound, vbest);
    if (best_predecessor < 0)
      blocked |= equal;
    else
      blocked |= equal & _mm512_cmple_epi32_mask(
                             max_index, _mm512_set1_epi32(best_predecessor));
    const __mmask16 valid = static_cast<__mmask16>((1u << n.child_count) - 1u);
    const __mmask16 keep =
        static_cast<__mmask16>(valid & ~impossible & ~blocked);
    if (keep == 0)
      return;
    alignas(64) int32_t kept_bound[kWideLanes];
    alignas(64) int32_t kept_index[kWideLanes];
    alignas(64) int32_t kept_child[kWideLanes];
    const __m512i lane =
        _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    _mm512_mask_compressstoreu_epi32(kept_bound, keep, bound);
    _mm512_mask_compressstoreu_epi32(kept_index, keep, max_index);
    _mm512_mask_compressstoreu_epi32(
        kept_child, keep,
        _mm512_add_epi32(lane, _mm512_set1_epi32(n.first_child)));
    const int kept = __builtin_popcount(static_cast<unsigned>(keep));
    WideVisit order[kWideLanes];
    for (int c = 0; c < kept; ++c) {
      const WideVisit v{kept_bound[c], kept_index[c], kept_child[c],
                        n.leaf_children};
      int p = c;
      while (p > 0 && (order[p - 1].bound > v.bound ||
                       (order[p - 1].bound == v.bound &&
                        order[p - 1].max_index > v.max_index))) {
        order[p] = order[p - 1];
        --p;
      }
      order[p] = v;
    }
    for (int c = 0; c < kept; ++c)
      stack_.push_back(order[c]);
  }
#endif // FA_DENSE_EXACT_WIDE_NEON

  // The leaf pass: a lane is one run, so its overlap and dd are exact. Lanes that can
  // still win are priced by elect, each re-tested against the incumbent the previous
  // election may have raised.
#if FA_DENSE_EXACT_WIDE_NEON
  FA_WIDE_TARGET void scan_leaf(int32_t block, int32_t i, const HotRun& b,
                                int32_t& best_score,
                                int32_t& best_predecessor) {
    bounds_ += kWideLanes;
    const WideLeaf& L = leaves_[static_cast<size_t>(block)];
    const int32_t d_b = b.r_begin - b.q_begin;
    const int32_t reach_floor = static_cast<int32_t>(
        static_cast<int64_t>(b.q_begin) + params_.span - reach_);
    const int64_t length = static_cast<int64_t>(b.q_end) - b.q_begin;
    const int32x4_t vqb = vdupq_n_s32(b.q_begin);
    const int32x4_t vrb = vdupq_n_s32(b.r_begin);
    const int32x4_t vqe = vdupq_n_s32(b.q_end);
    const int32x4_t vre = vdupq_n_s32(b.r_end);
    const int32x4_t vd = vdupq_n_s32(d_b);
    const int32x4_t vi = vdupq_n_s32(i);
    const int32x4_t vreach = vdupq_n_s32(reach_floor);
    const int32x4_t vbw = vdupq_n_s32(params_.bw);
    const int32x4_t zero = vdupq_n_s32(0);
    uint32x4_t impossible[4];
    int32x4_t dd[4];
    int32x4_t overlap[4];
    uint32x4_t overlapping = vdupq_n_u32(0);
    for (int q = 0; q < 4; ++q) {
      const int at = 4 * q;
      const int32x4_t index = vld1q_s32(L.index + at);
      const int32x4_t q_begin = vld1q_s32(L.q_begin + at);
      const int32x4_t q_end = vld1q_s32(L.q_end + at);
      const int32x4_t r_begin = vld1q_s32(L.r_begin + at);
      const int32x4_t r_end = vld1q_s32(L.r_end + at);
      const int32x4_t d_end = vld1q_s32(L.d_end + at);
      // Empty lanes (index -1) and unscored ones (index >= B's) are rejected, then
      // dense_step's tests run on the lane's exact coordinates.
      uint32x4_t imp = vcltq_s32(index, zero);
      imp = vorrq_u32(imp, vcgeq_s32(index, vi));
      imp = vorrq_u32(imp, vcleq_s32(vqb, q_begin));
      imp = vorrq_u32(imp, vcleq_s32(vrb, r_begin));
      imp = vorrq_u32(imp, vcleq_s32(vqe, q_end));
      imp = vorrq_u32(imp, vcleq_s32(vre, r_end));
      imp = vorrq_u32(imp, vcltq_s32(q_end, vreach));
      const int32x4_t dd_q = vabsq_s32(vsubq_s32(vd, d_end));
      imp = vorrq_u32(imp, vcgtq_s32(dd_q, vbw));
      const int32x4_t overlap_q = vmaxq_s32(
          vmaxq_s32(vsubq_s32(q_end, vqb), vsubq_s32(r_end, vrb)), zero);
      overlapping = vorrq_u32(
          overlapping, vbicq_u32(vtstq_s32(overlap_q, overlap_q), imp));
      impossible[q] = imp;
      dd[q] = dd_q;
      overlap[q] = overlap_q;
    }
    const bool no_overlap = vmaxvq_u32(overlapping) == 0;
    const int32x4_t vw = vdupq_n_s32(b.w);
    const int32x4_t vbest = vdupq_n_s32(best_score);
    const int32x4_t vbestpred = vdupq_n_s32(best_predecessor);
    alignas(64) int32_t bound_lane[kWideLanes];
    uint32x4_t kept_mask[4];
    for (int q = 0; q < 4; ++q) {
      const int at = 4 * q;
      const int32x4_t gain =
          no_overlap ? vw
                     : wide_gain_lanes(overlap[q], b.w, length, params_.span);
      const int32x4_t drift =
          wide_drift_lanes(dd[q], params_.chn_pen_gap, log_term_);
      const int32x4_t bound =
          vsubq_s32(vaddq_s32(vld1q_s32(L.score + at), gain), drift);
      vst1q_s32(bound_lane + at, bound);
      uint32x4_t k = vcgtq_s32(bound, vbest);
      if (best_predecessor >= 0)
        k = vorrq_u32(k, vandq_u32(vceqq_s32(bound, vbest),
                                   vcgtq_s32(vld1q_s32(L.index + at),
                                             vbestpred)));
      // `& ~impossible` folded in, so the bit form is taken once.
      kept_mask[q] = vbicq_u32(k, impossible[q]);
    }
    uint64_t keep = wide_lane_bits(kept_mask[0], kept_mask[1], kept_mask[2],
                                   kept_mask[3]);
    if (keep == 0)
      return;
    // Ascending lane order again: elect raises the incumbent as it goes.
    while (keep != 0) {
      const int lane = __builtin_ctzll(keep) >> 2;
      keep &= keep - 1;
      if (bound_lane[lane] < best_score ||
          (bound_lane[lane] == best_score &&
           (best_predecessor < 0 || L.index[lane] <= best_predecessor)))
        continue;
      elect(b, L.index[lane], best_score, best_predecessor);
    }
  }
#else
  FA_WIDE_TARGET void scan_leaf(int32_t block, int32_t i, const HotRun& b,
                                int32_t& best_score,
                                int32_t& best_predecessor) {
    bounds_ += kWideLanes;
    const WideLeaf& L = leaves_[static_cast<size_t>(block)];
    const int32_t d_b = b.r_begin - b.q_begin;
    const int32_t reach_floor = static_cast<int32_t>(
        static_cast<int64_t>(b.q_begin) + params_.span - reach_);
    const int64_t length = static_cast<int64_t>(b.q_end) - b.q_begin;
    const __m512i index = _mm512_load_si512(L.index);
    const __m512i q_begin = _mm512_load_si512(L.q_begin);
    const __m512i q_end = _mm512_load_si512(L.q_end);
    const __m512i r_begin = _mm512_load_si512(L.r_begin);
    const __m512i r_end = _mm512_load_si512(L.r_end);
    const __m512i score = _mm512_load_si512(L.score);
    const __m512i d_end = _mm512_load_si512(L.d_end);
    const __m512i vqb = _mm512_set1_epi32(b.q_begin);
    const __m512i vrb = _mm512_set1_epi32(b.r_begin);
    const __m512i vqe = _mm512_set1_epi32(b.q_end);
    const __m512i vre = _mm512_set1_epi32(b.r_end);
    const __m512i vd = _mm512_set1_epi32(d_b);
    const __m512i zero = _mm512_setzero_si512();
    __mmask16 impossible = _mm512_cmplt_epi32_mask(index, zero);
    impossible |= _mm512_cmpge_epi32_mask(index, _mm512_set1_epi32(i));
    impossible |= _mm512_cmple_epi32_mask(vqb, q_begin);
    impossible |= _mm512_cmple_epi32_mask(vrb, r_begin);
    impossible |= _mm512_cmple_epi32_mask(vqe, q_end);
    impossible |= _mm512_cmple_epi32_mask(vre, r_end);
    impossible |=
        _mm512_cmplt_epi32_mask(q_end, _mm512_set1_epi32(reach_floor));
    const __m512i dd = _mm512_abs_epi32(_mm512_sub_epi32(vd, d_end));
    impossible |= _mm512_cmpgt_epi32_mask(dd, _mm512_set1_epi32(params_.bw));
    const __m512i overlap =
        _mm512_max_epi32(_mm512_max_epi32(_mm512_sub_epi32(q_end, vqb),
                                          _mm512_sub_epi32(r_end, vrb)),
                         zero);
    __m512i gain;
    if ((_mm512_test_epi32_mask(overlap, overlap) & ~impossible) == 0)
      gain = _mm512_set1_epi32(b.w);
    else
      gain = wide_gain_lanes(overlap, b.w, length, params_.span);
    const __m512i drift = wide_drift_lanes(dd, params_.chn_pen_gap, log_term_);
    const __m512i bound =
        _mm512_sub_epi32(_mm512_add_epi32(score, gain), drift);
    const __m512i vbest = _mm512_set1_epi32(best_score);
    __mmask16 keep = _mm512_cmpgt_epi32_mask(bound, vbest);
    if (best_predecessor >= 0)
      keep |=
          _mm512_cmpeq_epi32_mask(bound, vbest) &
          _mm512_cmpgt_epi32_mask(index, _mm512_set1_epi32(best_predecessor));
    keep = static_cast<__mmask16>(keep & ~impossible);
    if (keep == 0)
      return;
    alignas(64) int32_t kept_index[kWideLanes];
    alignas(64) int32_t kept_bound[kWideLanes];
    _mm512_mask_compressstoreu_epi32(kept_index, keep, index);
    _mm512_mask_compressstoreu_epi32(kept_bound, keep, bound);
    const int kept = __builtin_popcount(static_cast<unsigned>(keep));
    for (int c = 0; c < kept; ++c) {
      if (kept_bound[c] < best_score ||
          (kept_bound[c] == best_score &&
           (best_predecessor < 0 || kept_index[c] <= best_predecessor)))
        continue;
      elect(b, kept_index[c], best_score, best_predecessor);
    }
  }
#endif // FA_DENSE_EXACT_WIDE_NEON

  const std::vector<HotRun>& hot_;
  const DenseChainParams& params_;
  DenseChainScratch& scratch_;
  std::vector<WideNode>& tree_;
  std::vector<WideLeaf>& leaves_;
  std::vector<int32_t>& run_leaf_;
  std::vector<WideKey>& keys_;
  std::vector<int32_t>& ids_;
  std::vector<WideTask>& tasks_;
  std::vector<int32_t>& cuts_;
  std::vector<WideVisit>& stack_;
  ::fa::cpu::lr::FlatInt64Map<int32_t>& diag_last_;
  const int32_t count_;
  const int64_t reach_;
  bool log_term_ = false;
  int64_t bounds_ = 0;
  int64_t steps_ = 0;
};

#undef FA_WIDE_TARGET
#endif // FA_DENSE_EXACT_WIDE

// The run DP, as colinear_chain.cpp's chain_driver over runs, with two exact predecessor
// searches:
//   linear  a descending scan pruned by `upper`, the prefix max of scores. gain <= W_B and
//           drift >= 0, so the scan stops once upper[j] + W_B <= best; upper only falls as
//           j descends and best only grows, so the stop is exact.
//   exact   ExactWideArm or ExactArm, for pools of diag_min_runs runs or more: branch and
//           bound over every scored run, returning the linear scan's score and predecessor.
DenseRunResult dense_run_chain(const std::vector<DenseRun>& runs,
                               const DenseChainParams& params) {
  DenseRunResult result;
  const int64_t count = static_cast<int64_t>(runs.size());
  if (count == 0)
    return result;
  DenseChainScratch& scratch = dense_chain_scratch();
  // Coordinates copied in; score and upper are written by the sweep before they are read.
  std::vector<HotRun>& hot = scratch.hot;
  hot.resize(static_cast<size_t>(count));
  for (int64_t i = 0; i < count; ++i) {
    const DenseRun& run = runs[static_cast<size_t>(i)];
    HotRun& h = hot[static_cast<size_t>(i)];
    h.q_begin = run.q_begin;
    h.q_end = run.q_end;
    h.r_begin = run.r_begin;
    h.r_end = run.r_end;
    h.w = run.w;
  }
  std::vector<int32_t>& predecessors = scratch.predecessors;
  predecessors.resize(static_cast<size_t>(count));
  std::vector<int32_t>& marks = scratch.marks;
  marks.assign(static_cast<size_t>(count), 0);
  // Admissibility runs from A's end to B's start, so a long predecessor within
  // reach can start more than max_dist_x back. The look-back window adds the longest run's
  // r-extent so no admissible predecessor is dropped.
  int64_t max_run_r_len = 0;
  for (const DenseRun& run : runs)
    max_run_r_len = std::max<int64_t>(
        max_run_r_len, static_cast<int64_t>(run.r_end) - run.r_begin);
  const int64_t window =
      static_cast<int64_t>(params.max_dist_x) + max_run_r_len;
  // The break needs drift >= 0, which takes non-negative penalties (and dense_step's
  // clamp of the skip term on overlap).
  const bool prune_by_upper_bound =
      params.chn_pen_gap >= 0.0f && params.chn_pen_skip >= 0.0f;
  // The exact arm searches the whole scored prefix, so the linear scan does not run after it.
  const bool use_exact =
      params.diag_min_runs >= 0 &&
      count >= static_cast<int64_t>(params.diag_min_runs);

  // The arms hold references only, so an unused one allocates nothing.
  ExactArm exact(hot, params, scratch, static_cast<int32_t>(count));
  // exact == 1 takes the wide arm where the CPU and the pool allow it; 2 forces ExactArm.
#if FA_DENSE_EXACT_WIDE
  ExactWideArm wide(hot, params, scratch, static_cast<int32_t>(count));
  const bool use_wide = use_exact && params.exact == 1 &&
                        ExactWideArm::supported() && wide.in_domain();
#else
  const bool use_wide = false;
#endif
  if (use_wide) {
#if FA_DENSE_EXACT_WIDE
    wide.build();
#endif
  } else if (use_exact) {
    exact.build();
  }

  int64_t start = 0;
  int64_t max_index = -1;
  // The carried argmax of scores over [start, scan_seen - 1]; -1 before the first score.
  int64_t scan_arg = -1;
  int64_t scan_seen = 0;
  for (int64_t i = 0; i < count; ++i) {
    const DenseRun& bi = runs[static_cast<size_t>(i)];
    int32_t best_predecessor = -1;
    int32_t best_score = bi.w;
    while (start<i&& static_cast<int64_t>(bi.r_begin)> static_cast<int64_t>(
               runs[static_cast<size_t>(start)].r_begin) +
           window) {
      ++start;
    }
    int64_t j = i - 1;
    if (use_exact) {
      // The exact arm answers over the whole scored prefix, beyond the look-back window too.
#if FA_DENSE_EXACT_WIDE
      if (use_wide)
        wide.search(static_cast<int32_t>(i), best_score, best_predecessor);
      else
#endif
        exact.search(static_cast<int32_t>(i), best_score, best_predecessor);
      j = i;
    }
    for (; (!use_exact || j < i) && j >= start; --j) {
      const HotRun& a = hot[static_cast<size_t>(j)];
      if (prune_by_upper_bound && a.upper + bi.w <= best_score)
        break;
      const int32_t step = dense_step(bi, a, params);
      if (step == kNegInf)
        continue;
      const int32_t score = step + a.score;
      if (score > best_score) {
        best_score = score;
        best_predecessor = static_cast<int32_t>(j);
      }
    }
    // Where the scan stopped: start - 1 at the window's end, else the break index.
    const int64_t end_j = j;
    if (max_index < 0 || static_cast<int64_t>(bi.r_begin) -
                                 runs[static_cast<size_t>(max_index)].r_begin >
                             static_cast<int64_t>(params.max_dist_x)) {
      // chain_driver rescans [start, i-1] here. This window is wider than max_dist_x, so a
      // stale run can be re-elected and the rescan would repeat on every run, a quadratic
      // term. The argmax is carried instead: ties go to the largest index as in the
      // descending rescan, it stays valid until start passes it, and new scores are folded
      // in as the window extends.
      if (scan_arg < start) {
        int32_t maximum = kNegInf;
        scan_arg = -1;
        for (int64_t scan = i - 1; scan >= start; --scan) {
          if (maximum < hot[static_cast<size_t>(scan)].score) {
            maximum = hot[static_cast<size_t>(scan)].score;
            scan_arg = scan;
          }
        }
      } else {
        for (; scan_seen < i; ++scan_seen) {
          // `>=`: a later equal score wins, as the descending rescan keeps the largest index.
          if (scan_arg < 0 ||
              hot[static_cast<size_t>(scan_seen)].score >=
                  hot[static_cast<size_t>(scan_arg)].score)
            scan_arg = scan_seen;
        }
      }
      scan_seen = i;
      max_index = scan_arg;
    }
    if (max_index >= 0 && max_index < end_j) {
      const HotRun& a = hot[static_cast<size_t>(max_index)];
      const int32_t extra = dense_step(bi, a, params);
      if (extra != kNegInf && best_score < extra + a.score) {
        best_score = extra + a.score;
        best_predecessor = static_cast<int32_t>(max_index);
      }
    }
    hot[static_cast<size_t>(i)].score = best_score;
    predecessors[static_cast<size_t>(i)] = best_predecessor;
    if (use_wide) {
#if FA_DENSE_EXACT_WIDE
      wide.activate(static_cast<int32_t>(i));
#endif
    } else if (use_exact) {
      exact.activate(static_cast<int32_t>(i));
    }
    hot[static_cast<size_t>(i)].upper =
        (i > 0 && hot[static_cast<size_t>(i - 1)].upper > best_score)
            ? hot[static_cast<size_t>(i - 1)].upper
            : best_score;
    if (max_index < 0 ||
        (static_cast<int64_t>(bi.r_begin) -
                 runs[static_cast<size_t>(max_index)].r_begin <=
             static_cast<int64_t>(params.max_dist_x) &&
         hot[static_cast<size_t>(max_index)].score <
             hot[static_cast<size_t>(i)].score)) {
      max_index = i;
    }
  }
  std::vector<int32_t>& scores = result.scores;
  scores.resize(static_cast<size_t>(count));
  for (int64_t i = 0; i < count; ++i)
    scores[static_cast<size_t>(i)] = hot[static_cast<size_t>(i)].score;
  result.exact_runs = use_exact ? count : 0;
#if FA_DENSE_EXACT_WIDE
  result.exact_nodes = use_wide ? wide.bounds() : exact.nodes();
  result.exact_steps = use_wide ? wide.steps() : exact.steps();
#else
  result.exact_nodes = exact.nodes();
  result.exact_steps = exact.steps();
#endif
  result.exact_wide = use_wide;

  std::vector<std::pair<int32_t, int64_t>>& ranked_ends = scratch.ranked_ends;
  ranked_ends.clear();
  ranked_ends.reserve(static_cast<size_t>(count));
  for (int64_t i = 0; i < count; ++i) {
    if (scores[static_cast<size_t>(i)] >= params.min_sc)
      ranked_ends.emplace_back(scores[static_cast<size_t>(i)], i);
  }
  if (ranked_ends.empty())
    return result;
  std::sort(ranked_ends.begin(), ranked_ends.end(),
            [](const std::pair<int32_t, int64_t>& left,
               const std::pair<int32_t, int64_t>& right) {
              if (left.first != right.first)
                return left.first < right.first;
              return left.second < right.second;
            });
  std::fill(marks.begin(), marks.end(), 0);
  const int32_t max_drop = params.bw;
  for (int64_t k = static_cast<int64_t>(ranked_ends.size()) - 1; k >= 0; --k) {
    const int64_t chain_start = ranked_ends[static_cast<size_t>(k)].second;
    if (marks[static_cast<size_t>(chain_start)] != 0)
      continue;
    const int32_t score_end = ranked_ends[static_cast<size_t>(k)].first;
    const int64_t end_i = backtrack_end(max_drop, score_end, scores,
                                        predecessors, marks, chain_start);
    int64_t length = 0;
    for (int64_t walk = chain_start; walk != end_i;
         walk = predecessors[static_cast<size_t>(walk)])
      ++length;
    // The cut run. chain_driver excludes the node whose net-negative link the chain would
    // pay for, which forfeits one span for an anchor but a whole run here, never re-seeded.
    // When no kept chain owns the cut run, it starts this chain at its fresh-start weight.
    const bool adopt_cut_run =
        end_i >= 0 && marks[static_cast<size_t>(end_i)] == 0;
    if (adopt_cut_run)
      ++length;
    std::vector<int32_t> indices(static_cast<size_t>(length));
    size_t fill = static_cast<size_t>(length);
    int64_t i = chain_start;
    for (; i != end_i; i = predecessors[static_cast<size_t>(i)]) {
      indices[--fill] = static_cast<int32_t>(i);
      marks[static_cast<size_t>(i)] = 1;
    }
    int32_t score =
        i < 0 ? score_end : score_end - scores[static_cast<size_t>(i)];
    if (adopt_cut_run) {
      indices[--fill] = static_cast<int32_t>(end_i);
      marks[static_cast<size_t>(end_i)] = 1;
      score += runs[static_cast<size_t>(end_i)].w;
    }
    // min_cnt counts runs here, where chain_driver counts anchors.
    if (score >= params.min_sc &&
        static_cast<int32_t>(indices.size()) >= params.min_cnt) {
      DenseRunChain chain;
      chain.idx = std::move(indices);
      chain.score = score;
      result.chains.push_back(std::move(chain));
    }
  }
  return result;
}

ChainResult chain_dense_colinear(std::vector<Anchor> anchors,
                                 DenseChainParams params,
                                 DenseChainStats* stats) {
  ChainResult result;
  if (params.max_dist_x < params.bw)
    params.max_dist_x = params.bw;
  if (params.max_dist_y < params.bw)
    params.max_dist_y = params.bw;
  if (anchors.empty())
    return result;

  // Callers usually pass anchors strictly sorted by (r, q), the streaming builder's
  // precondition; sort only when an O(n) check fails.
  const auto rq_strictly_less = [](const Anchor& left, const Anchor& right) {
    if (left.r != right.r)
      return left.r < right.r;
    return left.q < right.q;
  };
  const bool rq_strictly_sorted =
      std::adjacent_find(
          anchors.begin(), anchors.end(),
          [&rq_strictly_less](const Anchor& left, const Anchor& right) {
            return !rq_strictly_less(left, right);
          }) == anchors.end();
  if (!rq_strictly_sorted)
    std::sort(anchors.begin(), anchors.end(), rq_strictly_less);
  result.anchors = std::move(anchors);

  std::vector<int32_t>& links = dense_chain_scratch().links;
  const std::vector<DenseRun> runs =
      dense_collapse_runs_stream(result.anchors, params.span, links);
  if (stats != nullptr) {
    stats->runs = static_cast<int64_t>(runs.size());
    stats->anchors = static_cast<int64_t>(result.anchors.size());
    // dense_run_chain's use_exact test.
    stats->used_exact_arm = params.diag_min_runs >= 0 &&
                            static_cast<int64_t>(runs.size()) >=
                                static_cast<int64_t>(params.diag_min_runs);
  }
  const DenseRunResult collapsed = dense_run_chain(runs, params);
  if (stats != nullptr) {
    stats->exact_runs = collapsed.exact_runs;
    stats->exact_nodes = collapsed.exact_nodes;
    stats->exact_steps = collapsed.exact_steps;
  }

  // Expand each run chain to its anchors. Within a run they ascend in r and q; across runs
  // a successor may overlap its predecessor, and its overlapped head, already discounted
  // by the trim, is dropped so the chain is strictly ascending in (q, r).
  result.chains.reserve(collapsed.chains.size());
  for (size_t which_chain = 0; which_chain < collapsed.chains.size();
       ++which_chain) {
    const DenseRunChain& chain = collapsed.chains[which_chain];
    std::vector<int32_t> idx;
    int32_t reserve = 0;
    for (const int32_t which : chain.idx)
      reserve += runs[static_cast<size_t>(which)].raw_count;
    idx.reserve(static_cast<size_t>(reserve));
    int32_t last_q = INT32_MIN;
    int32_t last_r = INT32_MIN;
    const auto emit_run = [&](const DenseRun& run) {
      for (int32_t at = run.first; at >= 0;
           at = links[static_cast<size_t>(at)]) {
        const Anchor& anchor = result.anchors[static_cast<size_t>(at)];
        if (anchor.q <= last_q || anchor.r <= last_r)
          continue;
        idx.push_back(at);
        last_q = anchor.q;
        last_r = anchor.r;
      }
    };
    for (size_t position = 0; position < chain.idx.size(); ++position)
      emit_run(runs[static_cast<size_t>(chain.idx[position])]);
    if (idx.empty())
      continue;
    result.chains.push_back(Chain{std::move(idx), chain.score});
  }
  return result;
}

} // namespace chaining
} // namespace cpu
} // namespace fa
