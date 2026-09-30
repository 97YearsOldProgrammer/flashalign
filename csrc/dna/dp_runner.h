// The banded DP entry point shared by DNA realization and RNA exon fills.
#pragma once

#include "realization_controller.h"

#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {

// Dual-affine scoring and bands, built once by the caller.
struct DpScoringParams {
  int match = 2;     // minimap2 -A
  int mismatch = 4;  // minimap2 -B, positive
  int ambi = 1;      // minimap2 --score-N
  int gap_open1 = 4; // minimap2 -O
  int gap_extend1 = 2;  // minimap2 -E
  int gap_open2 = 24;   // minimap2 -O2
  int gap_extend2 = 1;  // minimap2 -E2
  int zdrop = 400; // minimap2 -z
  int tail_end_bonus =
      -1;       // minimap2 opt->end_bonus
  int bw = 500; // minimap2 -r, first value
  int bw_long = 20000; // minimap2 -r, second value
  // Local-inversion parameters, read only by DNA internal and supplementary
  // fills.
  int inversion_zdrop = 200;    // minimap2 zdrop_inv
  int inversion_max_gap = 5000; // the preset's max_gap
  int inversion_min_chain_score = 40;
  int inversion_min_dp_max = 80;
};

// The certified gaps of a verified region, region-relative. Passed only for
// a verified region; null for every other packet.
struct VerifiedRegionInputs {
  const ::fa::cpu::lr::RegionGap *gaps = nullptr;
  int gap_count = 0;
};

// Runs one realization packet through the controller. Slicing the packet and
// committing the resulting CIGAR are left to the caller. `inversion_probe_gate`,
// when set, gates the late probe of a fill whose probe is enabled.
realization::RealizationOutcome run_dna_long_realization(
    realization::RealizationRole role, const DpScoringParams &dp,
    realization::SequenceSlice query, realization::SequenceSlice target,
    int caller_requested_band, int interior_zdrop = -1, bool long_join = false,
    realization::ExecutionContext *execution_context = nullptr,
    bool inversion_probe_enabled = false,
    const VerifiedRegionInputs *verified_region = nullptr,
    const realization::InversionProbeGate *inversion_probe_gate = nullptr);

} // namespace lr
} // namespace cpu
} // namespace fa
