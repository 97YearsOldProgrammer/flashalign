// A stage of the mapping pipeline, which a preset's lane runs or leaves out
// (presets.h lane_runs_stage). An option that configures a stage the lane
// leaves out is refused.
#pragma once

namespace fa {
namespace cpu {
namespace options {

enum class LaneStage {
  None,            // run by every lane
  BaseLevelOutput, // a requested CIGAR: SAM, or PAF with one
  BaseAlignment,   // base-level alignment: DP scoring, bands, CIGAR operators
  SamOutput,       // the shape of SAM output
  Selection,       // primary and secondary selection
  Partition,       // the query partition
  OverlapPairs,    // which reads of an overlap pair print it
};

} // namespace options
} // namespace cpu
} // namespace fa
