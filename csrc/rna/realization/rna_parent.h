// minimap2's region bookkeeping (hit.c, align.c) over a small hit struct:
// mm_set_parent, mm_update_dp_max and mm_hit_sort's order. Pure; the caller builds the
// hits and reads the verdicts back by `owner`.
#pragma once

#include "rival_lifecycle.h" // RnaRankRecalAccounting

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// minimap2's MM_PARENT_UNSET.
inline constexpr int kRnaParentUnset = -1;

// One region (mm_reg1_t). Query coordinates are on the forward read.
struct RnaParentHit {
  int qs = 0;
  int qe = 0;
  int score = 0; // the region's own chain score
  int cnt = 0;
  int rid = -1;
  std::int64_t rs = 0;
  std::int64_t re = 0;
  // minimap2's r->p != 0: the region is aligned and dp_max / dp_max2 are set.
  bool has_dp = false;
  int dp_max = 0;
  int dp_max2 = 0;
  // mm_update_dp_max's input for an aligned region; accounting_ok false when it could
  // not be formed.
  RnaRankRecalAccounting accounting;
  bool accounting_ok = false;
  // minimap2's bookkeeping, written by rna_parent_set.
  int parent = kRnaParentUnset;
  int subsc = 0;
  int n_sub = 0;
  // The caller's keys: which hypothesis the region belongs to, its catalogue index (the
  // declared tie-break) and its segment within the hypothesis.
  int owner = -1;
  int catalogue = 0;
  int segment = 0;
};

// mm_set_parent with hard_mask_level 0 and no mask_len bound, without the alt-contig
// terms. Sets parent, subsc, n_sub and, on aligned pairs, dp_max2; r[0]'s n_sub is never
// reset (minimap2 hit.c:182).
void rna_parent_set(std::vector<RnaParentHit>& hits, float mask_level, int sub_diff);

// mm_update_dp_max with rank_min_len: reprices every aligned region's dp_max by
// mm_recal_max_dp when the read is at least min_len long, there are two aligned regions,
// the top spans at least frac of the read and the second is within frac of the top.
// Atomic: a region whose accounting is missing, or a top region without an identity,
// leaves every price raw. Returns whether the prices were rewritten.
bool rna_parent_update_dp_max(int qlen, int min_len, std::vector<RnaParentHit>& hits,
                              float frac, int a, int b);

// mm_hit_sort's order: dp_max for aligned regions, score otherwise, descending; ties to
// the lower catalogue index, then the lower segment. Stable.
void rna_parent_hit_sort(std::vector<RnaParentHit>& hits);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
