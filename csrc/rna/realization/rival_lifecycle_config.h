// Rival lifecycle options, in their own header so context.h can carry them without
// including the realizer.
#pragma once

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Defaults follow minimap2: pri_ratio and min_diff are its secondary retention band
// (min_diff = 2k = 30 at k = 15).

struct RnaRivalLifecycleConfig {
  // -p / --rival-min-diff: a competitor is admitted when its chain score is within
  // either band of the primary's (inclusive).
  double pri_ratio = 0.8;
  int min_diff = 30;
  // --realize-max: bound on rival realizations per read (minimap2's best_n). The scan
  // itself covers the whole catalogue.
  int realize_max = 5;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
