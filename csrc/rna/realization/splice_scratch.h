// Caller-owned per-worker realization scratch, like minimap2's mm_tbuf_t: grown on demand,
// never shrunk, and reused across packets, segments, hypotheses and reads. A null pointer
// keeps every buffer function-local.
#pragma once

#include "../../dp/splice_kernel.h"
#include "controller_metrics.h"
#include "splice_anchor.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// One ignore/long-join bit set the cluster filters added to one anchor.
struct SpliceAnchorFlagMutation {
  std::uint32_t index = 0;
  std::uint8_t added = 0;
};

// One segment's orientation-blind preparation: the terminal-anchor filter, the two
// cluster filters, the pivots and the terminal envelope. None of it reads the transcript
// orientation or the flags the cluster filters add, so the second hypothesis replays the
// first's.
struct SplicePreparedSegment {
  int selected_begin = 0;
  int selected_count = 0;
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;
  // The segment-result fields the preparation wrote before any refusal, so a replay
  // writes the same ones.
  bool wrote_filtered = false;
  bool wrote_envelope = false;
  int filtered_begin = 0;
  int filtered_count = 0;
  int first_reference_pivot = 0;
  int first_query_pivot = 0;
  int last_reference_pivot = 0;
  int last_query_pivot = 0;
  int envelope_reference_begin = 0;
  int envelope_query_begin = 0;
  int envelope_reference_end = 0;
  int envelope_query_end = 0;
  int left_neighbor = -1;
  int right_neighbor = -1;
  std::vector<SpliceAnchorFlagMutation> anchor_flags;
};

struct SpliceRealizationScratch {
  // ksw2 arena and the reused ez.cigar.
  dp::SpliceKernelWorkspace kernel;
  // The kernel's last result, reused so its CIGAR keeps its capacity.
  dp::SpliceKernelResult aligned;
  // One packet's known-junction mask.
  std::vector<std::uint8_t> junction;
  // The left packet's reversed query and reference flanks.
  std::vector<std::uint8_t> reversed_query;
  std::vector<std::uint8_t> reversed_reference;
  // The hypothesis' mutable copy of the frozen anchor view.
  std::vector<SpliceAnchor> working;
  // Segments prepared during this realization, keyed by selected range. Anchor geometry
  // is fixed for a realization, so entries stay valid for later segments and the second
  // hypothesis; realize_exact_anchor_path clears them for each new anchor path.
  std::vector<SplicePreparedSegment> prepared;
  // Flag bytes of the selected range before the cluster filters ran.
  std::vector<std::uint8_t> flags_before;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
