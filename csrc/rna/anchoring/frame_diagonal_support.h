// Which query frame a harvested candidate pool is diagonally coherent in, read off the
// pool before any chain runs. In the right frame an exon's anchors share a diagonal r - q.
// In the mirrored frame q' = length - span - q, and every anchor of a pool has the same
// span, so the mirror's diagonal is the anti-diagonal r + q shifted by a constant, which
// run and band clustering cannot see. One pool therefore answers for both frames. The
// fine stage elects each window's frame this way and chains both frames only when the
// verdict is ambiguous.
#pragma once

#include "skeleton_harvest.h"

#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Diagonal coherence of one pool in both query frames at one tolerance, over the
// final-pool anchors only. A run is a maximal group of sorted keys whose consecutive
// differences are within the tolerance: `mass` counts the anchors in runs of two or more,
// `best_run` is the largest run and `runs_ge3` counts runs of three or more. `band` is
// the most keys inside any window of width `tolerance`; unlike run mass it stays
// discriminating on a dense pool, where the wrong frame's keys chain into long runs.
struct FrameDiagonalSupport {
  int32_t tolerance = 0;
  uint32_t anchors = 0;
  uint32_t mass_nominated = 0;
  uint32_t best_run_nominated = 0;
  uint32_t runs_ge3_nominated = 0;
  uint32_t mass_mirror = 0;
  uint32_t best_run_mirror = 0;
  uint32_t runs_ge3_mirror = 0;
  uint32_t band_nominated = 0;
  uint32_t band_mirror = 0;
};

// Clusters the pool's final anchors in both key spaces. `tolerance` is in reference
// bases; a negative value is read as 0. A refused or empty pool reports zeros.
FrameDiagonalSupport compute_frame_diagonal_support(const CandidatePool& pool,
                                                    int32_t tolerance);

// kAmbiguous means the bands are too close to settle the frame cheaply; the window is
// then chained in both frames and settled by the chain scores.
enum class FrameElection { kNominated, kMirror, kAmbiguous };

// With margin = |band_n - band_m|, the window is decided when
//
//     margin >= min_margin  and  margin >= min_ratio * (band_n + band_m)
//
// and a decided window elects the mirror iff band_m > band_n, the tie rule of the
// chain-score comparison. The absolute bound keeps a tiny pool from deciding on noise;
// the relative bound keeps a dense pool from deciding on a small margin.
FrameElection elect_query_frame(const FrameDiagonalSupport& support);

// The same reading off the two strand-routed pools: `nominated` is the harvested lane's
// pool and `opposite` the other lane's, already mirrored to q' = length - span - q. The
// opposite pool's own diagonal key r - q' is therefore the mirror key, so both pools are
// read with r - q. Unlike the one-pool reading, the two masses come from different anchor
// sets, each lane capped on its own. A refused pool contributes zeros.
FrameDiagonalSupport
routed_frame_diagonal_support(const CandidatePool& nominated,
                              const CandidatePool& opposite, int32_t tolerance);

// Tolerance at which the fine stage computes the support, the same for rank-1 and rival
// windows.
constexpr int32_t kFrameElectionTolerance = 8;

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
