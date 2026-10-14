// Anchor, Chain and ChainResult, the chaining parameters, query intervals and the mask-level overlap test, and the chainers' declarations.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

/**** anchor ****/
// Anchor type shared by the DNA and RNA chainers.

namespace fa {
namespace cpu {
namespace chaining {

// As minimap2's MM_SEED_IGNORE, MM_SEED_TANDEM and MM_SEED_LONG_JOIN.
enum AnchorFlag : uint32_t {
  ANCHOR_IGNORE = 1u << 0,
  ANCHOR_TANDEM = 1u << 1,
  ANCHOR_LONG_JOIN = 1u << 2,
};

// `r` and `q` are reference/query coordinates. The recurrence uses pairwise
// differences, so uniform-span START and END coordinates are equivalent.
struct Anchor {
  int32_t r = 0;
  int32_t q = 0;
  int32_t span = 0;
  int32_t id = -1;
  uint32_t flags = 0;

  int32_t q_end() const { return q + span; }
  int32_t r_end() const { return r + span; }
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa

/**** result ****/

namespace fa {
namespace cpu {
namespace chaining {

struct Chain {
  std::vector<int32_t> idx;
  int32_t score = 0;
};

struct ChainResult {
  std::vector<Anchor> anchors;
  std::vector<Chain> chains;
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa

/**** colinear_params ****/

namespace fa {
namespace cpu {
namespace chaining {

// minimap2 single-segment DNA chaining defaults (map-ont).
struct ColinearChainParams {
  int32_t max_dist_x = 5000;
  int32_t max_dist_y = 5000;
  int32_t bw = 500;
  int32_t max_skip = 25;
  int32_t max_iter = 5000;
  int32_t min_cnt = 3;
  int32_t min_sc = 40;
  float chn_pen_gap = 0.168f;
  float chn_pen_skip = 0.0f;
};

// Spliced (cDNA) chaining, as minimap2's: forward reference gaps pay the cheaper of the
// linear and log penalties, spans are uncapped, max_dist_y is not floored to bw, and the
// chain Z-drop is off.
struct SplicedChainParams {
  int32_t max_dist_x = 5000;
  int32_t max_dist_y = 5000;
  int32_t bw = 500;
  int32_t max_skip = 25;
  int32_t max_iter = 5000;
  int32_t min_cnt = 3;
  int32_t min_sc = 40;
  float chn_pen_gap = 0.168f;
  float chn_pen_skip = 0.0f;
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa

/**** coordinates ****/
// Query intervals and strands for split mapping. An interval passed between stages is in
// forward-read coordinates; the strand frame is used only inside realization.

namespace fa {
namespace cpu {
namespace split {

// Strand of a placement: Forward is the read as given, Reverse its reverse complement.
enum class Strand : uint8_t { Forward = 0, Reverse = 1 };

// Conversions to and from the `bool is_reverse` form.
constexpr Strand strand_from_is_reverse(bool is_reverse) noexcept {
  return is_reverse ? Strand::Reverse : Strand::Forward;
}
constexpr bool is_reverse_strand(Strand s) noexcept {
  return s == Strand::Reverse;
}

// Half-open query interval [begin, end).
struct QueryInterval {
  int begin = 0;
  int end = 0;

  // 0 for an empty or inverted interval.
  constexpr int length() const noexcept { return end > begin ? end - begin : 0; }
  constexpr bool valid() const noexcept { return end > begin; }

  friend constexpr bool operator==(const QueryInterval &a,
                                   const QueryInterval &b) noexcept {
    return a.begin == b.begin && a.end == b.end;
  }
  friend constexpr bool operator!=(const QueryInterval &a,
                                   const QueryInterval &b) noexcept {
    return !(a == b);
  }
};

// [b, e) on a read of length L maps to [L - e, L - b). Its own inverse.
constexpr QueryInterval reverse_complement_interval(QueryInterval qi,
                                                    int read_len) noexcept {
  return QueryInterval{read_len - qi.end, read_len - qi.begin};
}

// Forward-read interval to the strand frame used for chaining and realization.
constexpr QueryInterval to_strand_frame(QueryInterval forward, int read_len,
                                        Strand s) noexcept {
  return s == Strand::Reverse ? reverse_complement_interval(forward, read_len)
                              : forward;
}

// Strand-frame interval back to forward-read coordinates.
constexpr QueryInterval to_forward(QueryInterval strand_frame, int read_len,
                                   Strand s) noexcept {
  return s == Strand::Reverse
             ? reverse_complement_interval(strand_frame, read_len)
             : strand_frame;
}

} // namespace split
} // namespace cpu
} // namespace fa

/**** query_geometry ****/
// Query-overlap test that decides whether two placements are rivals for the same read span or
// separate pieces of the read: minimap2's mask_level test (mm_set_parent).

namespace fa {
namespace cpu {
namespace split {

// Relationship between two query intervals for split-mapping selection.
enum class QueryRelation : uint8_t {
  DisjointSegment, // separately alignable piece of the read (supplementary)
  Rival            // alternative placement of the same span (secondary/f2)
};

// Overlap fraction of the shorter interval at or above which the pair are rivals. The
// threshold is rounded down.
struct QueryMaskPolicy {
  uint32_t numerator = 1;
  uint32_t denominator = 2;
};

inline constexpr QueryMaskPolicy kHalfFloor{1, 2};

// Rival when the overlap reaches the policy fraction of the shorter interval; no overlap or
// an empty interval is always DisjointSegment.
// Inline: called per rival chain from partition_chains (colinear_chain.cpp).
inline QueryRelation classify_query_relation(QueryInterval a, QueryInterval b,
                                             QueryMaskPolicy p = {}) noexcept {
  const int overlap =
      std::min(a.end, b.end) - std::max(a.begin, b.begin);
  if (overlap <= 0)
    return QueryRelation::DisjointSegment;
  const int shorter = std::min(a.length(), b.length());
  if (shorter <= 0)
    return QueryRelation::DisjointSegment;
  const int64_t scaled =
      static_cast<int64_t>(shorter) * static_cast<int64_t>(p.numerator);
  const int64_t den = static_cast<int64_t>(p.denominator);
  const int threshold = static_cast<int>(scaled / den);
  return overlap >= threshold ? QueryRelation::Rival
                              : QueryRelation::DisjointSegment;
}

} // namespace split
} // namespace cpu
} // namespace fa

/**** colinear_chain ****/
// Colinear chaining as minimap2's mg_lchain_dp, for DNA and spliced (cDNA) anchors.

namespace fa {
namespace cpu {
namespace chaining {

ChainResult chain_colinear(std::vector<Anchor> anchors,
                           ColinearChainParams params);

ChainResult chain_spliced_colinear(std::vector<Anchor> anchors,
                                   SplicedChainParams params);

int best_chain_index(const ChainResult& result);

}  // namespace chaining
}  // namespace cpu
}  // namespace fa

/**** partition ****/
// Picks the primary chain and counts its rivals (n_sub, f2) for MAPQ.

namespace fa {
namespace cpu {
namespace chaining {

struct ChainPartition {
  int primary = -1;
  int32_t f2 = 0;
  // Index in ChainResult::chains of the best-scoring rival (the chain f2 comes from);
  // -1 when no chain rivals the primary.
  int f2_index = -1;
  int n_sub = 0;
};

ChainPartition partition_chains(
    const ChainResult& result,
    ::fa::cpu::split::QueryMaskPolicy policy =
        ::fa::cpu::split::kHalfFloor);

}  // namespace chaining
}  // namespace cpu
}  // namespace fa

/**** dense_chain ****/
// Dense-pool colinear chaining: collapse an anchor pool to same-diagonal runs, run the
// minimap2 recurrence over the runs and expand the run chains back to anchors. With one
// anchor per run the recurrence is minimap2's anchor chain.

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
