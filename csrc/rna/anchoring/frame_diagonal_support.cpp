#include "frame_diagonal_support.h"

#include <algorithm>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

namespace {

struct RunTotals {
  uint32_t mass = 0;
  uint32_t best = 0;
  uint32_t runs_ge3 = 0;
  uint32_t band = 0;
};

// Maximal runs of the sorted keys under `tolerance`, and the best band. Sorts `keys` in
// place.
RunTotals summarize_runs(std::vector<int64_t>& keys, int64_t tolerance) {
  RunTotals totals;
  std::sort(keys.begin(), keys.end());
  size_t begin = 0;
  while (begin < keys.size()) {
    size_t end = begin + 1;
    while (end < keys.size() && keys[end] - keys[end - 1] <= tolerance)
      ++end;
    const uint32_t run = static_cast<uint32_t>(end - begin);
    if (run >= 2)
      totals.mass += run;
    if (run >= 3)
      ++totals.runs_ge3;
    totals.best = std::max(totals.best, run);
    begin = end;
  }
  // Best band by two pointers over the sorted keys. A band's keys are pairwise within the
  // tolerance, so every band lies inside one run: band <= best.
  size_t low = 0;
  for (size_t high = 0; high < keys.size(); ++high) {
    while (keys[high] - keys[low] > tolerance)
      ++low;
    totals.band =
        std::max(totals.band, static_cast<uint32_t>(high - low + 1));
  }
  return totals;
}

// One pool's r - q keys over its final anchors, multiplicity kept. `keys` is caller
// scratch, left sorted.
struct LaneTotals {
  uint32_t anchors = 0;
  RunTotals runs;
};

LaneTotals summarize_pool_diagonal(const CandidatePool& pool, int64_t width,
                                   std::vector<int64_t>& keys) {
  keys.clear();
  keys.reserve(pool.anchors.size());
  for (const FineAnchor& anchor : pool.anchors) {
    if (!anchor.final_pool)
      continue;
    const int64_t reference = anchor.reference_begin;
    const int64_t query = anchor.query_begin;
    keys.push_back(reference - query);
  }
  LaneTotals totals;
  totals.anchors = static_cast<uint32_t>(keys.size());
  totals.runs = summarize_runs(keys, width);
  return totals;
}

} // namespace

FrameDiagonalSupport compute_frame_diagonal_support(const CandidatePool& pool,
                                                    int32_t tolerance) {
  FrameDiagonalSupport support;
  support.tolerance = tolerance;
  std::vector<int64_t> nominated;
  std::vector<int64_t> mirror;
  nominated.reserve(pool.anchors.size());
  mirror.reserve(pool.anchors.size());
  for (const FineAnchor& anchor : pool.anchors) {
    if (!anchor.final_pool)
      continue;
    const int64_t reference = anchor.reference_begin;
    const int64_t query = anchor.query_begin;
    nominated.push_back(reference - query);
    mirror.push_back(reference + query);
  }
  support.anchors = static_cast<uint32_t>(nominated.size());

  const int64_t width = tolerance > 0 ? tolerance : 0;
  const RunTotals in_frame = summarize_runs(nominated, width);
  support.mass_nominated = in_frame.mass;
  support.best_run_nominated = in_frame.best;
  support.runs_ge3_nominated = in_frame.runs_ge3;
  support.band_nominated = in_frame.band;
  const RunTotals mirrored = summarize_runs(mirror, width);
  support.mass_mirror = mirrored.mass;
  support.best_run_mirror = mirrored.best;
  support.runs_ge3_mirror = mirrored.runs_ge3;
  support.band_mirror = mirrored.band;
  return support;
}

FrameDiagonalSupport
routed_frame_diagonal_support(const CandidatePool& nominated,
                              const CandidatePool& opposite,
                              int32_t tolerance) {
  // Each pool is read in its own frame with the r - q key; the opposite pool's reading is
  // the mirror half (see the header). A refused pool contributes zeros.
  FrameDiagonalSupport support;
  support.tolerance = tolerance;
  const int64_t width = tolerance > 0 ? tolerance : 0;
  std::vector<int64_t> keys;
  if (!nominated.refused) {
    const LaneTotals in_lane = summarize_pool_diagonal(nominated, width, keys);
    support.anchors += in_lane.anchors;
    support.mass_nominated = in_lane.runs.mass;
    support.best_run_nominated = in_lane.runs.best;
    support.runs_ge3_nominated = in_lane.runs.runs_ge3;
    support.band_nominated = in_lane.runs.band;
  }
  if (!opposite.refused) {
    const LaneTotals other_lane =
        summarize_pool_diagonal(opposite, width, keys);
    support.anchors += other_lane.anchors;
    support.mass_mirror = other_lane.runs.mass;
    support.best_run_mirror = other_lane.runs.best;
    support.runs_ge3_mirror = other_lane.runs.runs_ge3;
    support.band_mirror = other_lane.runs.band;
  }
  return support;
}

namespace {

constexpr uint32_t kFrameElectionMinMargin = 8;
constexpr double kFrameElectionMinRatio = 0.5;

} // namespace

FrameElection elect_query_frame(const FrameDiagonalSupport& support) {
  const uint32_t nominated = support.band_nominated;
  const uint32_t mirror = support.band_mirror;
  const uint32_t margin =
      nominated > mirror ? nominated - mirror : mirror - nominated;
  if (margin < kFrameElectionMinMargin)
    return FrameElection::kAmbiguous;
  // Both bands count anchors of one pool, so the double comparison is exact.
  const double total =
      static_cast<double>(nominated) + static_cast<double>(mirror);
  if (static_cast<double>(margin) < kFrameElectionMinRatio * total)
    return FrameElection::kAmbiguous;
  return mirror > nominated ? FrameElection::kMirror
                            : FrameElection::kNominated;
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
