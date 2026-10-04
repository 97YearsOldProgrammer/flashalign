// Splice realization of one selected anchor segment. The controller works on a mutable
// copy of the anchor view and holds no global state; the DP kernel is dp/splice_kernel.h.
#pragma once

#include "../../dp/splice_kernel.h"
#include "controller_metrics.h"
#include "splice_anchor.h"
#include "splice_scratch.h"
#include "../annotation/known_junctions.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The transcript hypothesis a realization pass scores. None is minimap2's -un: no splice
// motif, so there is a single orientation-free pass.
enum class TranscriptOrientation : std::uint8_t {
  Forward = 1,
  Reverse = 2,
  None = 3,
};

enum class SplicePacketRole : std::uint8_t {
  Left = 1,
  Internal = 2,
  Right = 3,
};

struct SpliceControllerOptions {
  int k = 15;
  int match = 1;
  int mismatch = 2;
  int gap_open = 2;
  int gap_extend = 1;
  int long_gap_open = 32;
  int long_gap_extend = 0;
  int transition = 0;
  int ambiguous = 1;
  int junction_bonus = 9;
  int junction_penalty = 5;
  int zdrop = 200;
  int inversion_zdrop = 100;
  int end_bonus = -1;
  int maximum_gap = 2000;
  int maximum_reference_gap = 200000;
  int bandwidth = 200000;
  int long_bandwidth = 200000;
  int minimum_anchor_count = 3;
  int minimum_match_bases = 40;
  int minimum_dp_maximum = 80;
  int minimum_packet_length = 200;
  // Shortest reference gap called an intron, used only by the packet-shape accounting;
  // the kernel derives its own threshold from the long-gap scores.
  int minimum_intron = 20;
  int anchor_extension_length = 20;
  int anchor_extension_shift = 6;
  std::int64_t maximum_dp_cells = 0;
};

SpliceControllerOptions ordinary_long_read_splice_options() noexcept;
bool splice_controller_options_supported(
    const SpliceControllerOptions& options) noexcept;

struct SpliceControllerRequest {
  const std::vector<std::uint8_t>* query_forward = nullptr;
  const std::vector<std::uint8_t>* query_reverse = nullptr;
  const std::vector<std::uint8_t>* reference = nullptr;
  std::vector<SpliceAnchor>* working_anchors = nullptr;
  int selected_begin = 0;
  int selected_count = 0;
  int chain_score = 0;
  int initial_chain_score = 0;
  bool mapping_reverse = false;
  TranscriptOrientation transcript_orientation =
      TranscriptOrientation::Forward;
  SpliceControllerOptions options = ordinary_long_read_splice_options();
  KnownJunctionContigView known_junctions;
  // Optional per-worker scratch. Null keeps every buffer function-local.
  SpliceRealizationScratch* scratch = nullptr;
};

struct SpliceContinuation {
  int selected_begin = 0;
  int selected_count = 0;
  int chain_score = 0;
  int initial_chain_score = 0;

  bool present() const noexcept { return selected_count > 0; }
};

struct SpliceSegmentResult {
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;

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
  int left_neighbor_index = -1;
  int right_neighbor_index = -1;

  int query_begin = 0;
  int query_end = 0;
  int reference_begin = 0;
  int reference_end = 0;
  int region_anchor_count = 0;
  int chain_score = 0;
  int matches = 0;
  int block_length = 0;
  bool is_spliced = false;
  int dp_score = 0;
  int dp_maximum = 0;
  int dp_maximum_before_strand = 0;
  int ambiguous_bases = 0;
  std::vector<std::uint32_t> cigar;

  SpliceContinuation continuation;
  std::uint32_t split_flags = 0;
  RnaControllerMetrics metrics;
};

SpliceSegmentResult realize_splice_segment(
    const SpliceControllerRequest& request) noexcept;

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
