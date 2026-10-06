// The anchor path of one block and the pre-DP control minimap2's mm_align1
// applies to it: k-mer center adjustment, terminal DP windows, bad-seed
// filtering and the verified realization geometry. No DP runs here; the
// geometry planner reads bases only to test gaps against the certificate.
//
// Anchor coordinates are exact k-mer starts. minimap2 stores a minimizer end
// and mm_adjust_minier() subtracts k/2 to reach the center;
// adjusted_endpoint() applies the equivalent start-to-center shift.
#pragma once

#include "../chaining/anchor.h"
#include "region_gap.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace ordered_anchor {

enum class CoordinateConvention : std::uint8_t {
  ExactKmerStart,
};

struct OrderedAnchorPath {
  int reference_id = -1;
  bool reverse_complemented = false;
  int query_length = 0;
  int target_length = 0;
  int query_bound_begin = 0;
  int query_bound_end = 0;
  int target_bound_begin = 0;
  int target_bound_end = 0;
  int minimizer_k = 0;
  CoordinateConvention coordinate_convention =
      CoordinateConvention::ExactKmerStart;
  bool homopolymer_compressed = false;

  // Coordinates are in the oriented query and the forward reference.
  // [realization_begin, realization_end) is the subrange mm_fix_bad_ends
  // keeps; the anchors outside it still bound the terminal windows.
  std::vector<::fa::cpu::chaining::Anchor> selected;
  std::size_t realization_begin = 0;
  std::size_t realization_end = 0;
};

enum class ValidationCode : std::uint8_t {
  Valid,
  EmptySelectedPath,
  InvalidReference,
  InvalidBounds,
  UnsupportedCoordinateConvention,
  UnsupportedHomopolymerCompression,
  InvalidAnchorSpan,
  InconsistentMinimizerSpan,
  AnchorOutOfBounds,
  NonMonotoneSelectedPath,
  InvalidAnchorFlags,
  InvalidSelectedSubrange,
  InvalidControl,
};

struct ValidationResult {
  ValidationCode code = ValidationCode::Valid;
  std::size_t anchor_index = 0;

  explicit operator bool() const { return code == ValidationCode::Valid; }
};

ValidationResult validate(const OrderedAnchorPath& path);

struct AdjustedEndpoint {
  int query = 0;
  int target = 0;
};

// mm_adjust_minier for exact k-mer starts, without HPC:
// start + (k - 1) - floor(k / 2).
AdjustedEndpoint adjusted_endpoint(const OrderedAnchorPath& path,
                                   const ::fa::cpu::chaining::Anchor& anchor);

struct NormalizationControl {
  int bandwidth = 0;
  int min_chain_score = 0;
  int max_gap = 0;
  bool trim_left = true;
  bool trim_right = true;
};

struct NormalizationStats {
  std::size_t trimmed_front = 0;
  std::size_t trimmed_back = 0;
  std::size_t ignored = 0;
  std::size_t long_join = 0;
};

struct NormalizedPath {
  ValidationResult validation;
  OrderedAnchorPath path;
  NormalizationStats stats;

  explicit operator bool() const {
    return validation.code == ValidationCode::Valid;
  }
};

// Pre-DP normalization as minimap2's mm_fix_bad_ends, mm_filter_bad_seeds and
// mm_filter_bad_seeds_alt. End trimming is enabled per side; the bad-seed
// filters always run and flag anchors ANCHOR_IGNORE / ANCHOR_LONG_JOIN.
NormalizedPath normalize_for_realization(const OrderedAnchorPath& path,
                                         const NormalizationControl& control);

// The verified realization geometry (see kernel_packets.h). A block is laid
// out left to right as steps of three kinds.
enum class GeometryStepKind : std::uint8_t {
  // A verified region, emitted ungapped from its first anchor's center to its
  // last's, with no DP.
  Verified,
  // The fill between two verified regions with no retained anchor between
  // them.
  Seam,
  // One fill of a stretch: unverified material between verified regions that
  // are not neighbours, before the first or after the last. minimap2's loop
  // cuts a stretch into pieces.
  Piece,
};

struct GeometryRegion {
  std::size_t first = 0; // selected index of its first retained anchor
  std::size_t last = 0;  // selected index of its last retained anchor
  int anchors = 0;       // retained anchors in it
  bool verified = false;
  // Gaps that failed the certificate, counting all of them.
  int failing_gaps = 0;
};

// One executor call, in path order. Consecutive steps abut and tile
// [center(first retained anchor), center(last retained anchor)].
struct GeometryStep {
  GeometryStepKind kind = GeometryStepKind::Piece;
  // The right corner anchor, where a Z-drop split cuts (minimap2's `i`).
  std::size_t selected_index = 0;
  int query_begin = 0; // [center(left corner), center(right corner))
  int query_end = 0;
  int target_begin = 0;
  int target_end = 0;
  // Seam / Piece: the right corner carries ANCHOR_LONG_JOIN.
  bool long_join = false;
  // Seam / Piece: bw_long_eff(), or max(q_span, r_span) on a long join.
  // Verified: bw_long_eff(), recorded; no DP runs.
  int selected_band = -1;
  // Verified: its gaps, [gap_begin, gap_begin + gap_count) in
  // GeometryPlan::gaps, relative to the step start on both axes.
  std::size_t gap_begin = 0;
  std::size_t gap_count = 0;
  // Piece: its index in GeometryPlan::stretches and its ordinal within the
  // stretch. -1 elsewhere.
  int stretch = -1;
  int piece = -1;

  // Counters only. Verified: the region's retained anchors. Piece: the
  // retained anchors it crosses. Seam: 0.
  int anchor_count = 0;
  // Seam / Piece: max(next.q - prev.q_end, next.r - prev.r_end) between the
  // corner anchors, negative when they overlap. 0 on a verified region.
  int anchor_free = 0;
};

// One stretch of unverified material. Counters only.
struct GeometryStretch {
  std::size_t first_step = 0; // its first piece in GeometryPlan::steps
  int pieces = 0;
  int query_bases = 0;    // center(end corner) - center(start corner)
  int lone_anchors = 0;   // regions of one retained anchor inside it
  int failed_regions = 0; // regions of >= 2 anchors with a failing gap
  int failed_anchors = 0; // the retained anchors of those regions
  int failed_bases = 0;   // their query span, center to center
  int failing_gaps = 0;   // the gaps in them that failed the certificate
};

struct GeometryPlan {
  int initial_query_cursor = 0; // center(first retained anchor)
  int initial_target_cursor = 0;
  int final_query_cursor = 0; // center(last retained anchor)
  int final_target_cursor = 0;
  std::vector<GeometryRegion> regions;
  std::vector<GeometryStep> steps;
  std::vector<GeometryStretch> stretches;
  // The gaps of every verified region.
  std::vector<::fa::cpu::lr::RegionGap> gaps;
};

struct GeometryControl {
  int normal_long_band = 0; // DpMapOpt::bw_long_eff()
  // minimap2 opt->min_ksw_len: inside a stretch a piece ends at a retained
  // anchor once both spans since the previous corner reach it.
  int min_ksw_len = 0;
};

// Inputs of the ungapped certificate, borrowed for the call. The ungapped
// alignment of an L-base gap is optimal when its score S0 reaches
// (L - 1) * a - 2 * g(1): a gapped alignment of two equal-length windows keeps
// at most L - 1 match columns and pays two one-base gaps. L < 2 always passes.
// This allows up to 3 mismatches per gap with map-hifi scoring, 2 with map-ont.
struct GapCertificate {
  const std::uint8_t* query = nullptr;  // the whole oriented query (anchor.q)
  const std::uint8_t* target = nullptr; // the whole chromosome (anchor.r)
  // ksw2_simple_mat(a, b, sc_ambi): N scores sc_ambi, even against N.
  std::array<std::int8_t, 25> matrix{};
  int match = 0;        // a
  int one_base_gap = 0; // g(1) = min(q + e, q2 + e2)
};

// Plans the verified geometry of one normalized block:
//   1. interior anchors flagged ANCHOR_IGNORE or ANCHOR_TANDEM are skipped, as
//      in minimap2; the first and last anchors never are;
//   2. the retained anchors are cut into regions, and each region of two or
//      more anchors is tested gap by gap, a gap being the bases between
//      consecutive k-mers, [p.q_end(), n.q); bases of skipped anchors inside a
//      gap are tested like any other;
//   3. a verified region becomes a Verified step, neighbouring verified
//      regions are joined by a Seam, and everything else is a stretch cut into
//      Pieces by minimap2's loop.
// `path` must be valid normalize_for_realization output; nothing is
// re-checked.
GeometryPlan plan_verified_geometry(const OrderedAnchorPath& path,
                                    const GeometryControl& control,
                                    const GapCertificate& certificate);

struct TerminalWindowControl {
  int max_gap = 0;
  int match = 0;
  int gap_open = 0;
  int gap_extend = 0;
  // minimap2's cap at nearby seeds (mm_align1): where set, the left window
  // starts at or after (left_cap_query, left_cap_target) and the right one
  // ends at or before (right_cap_query, right_cap_target).
  bool left_cap = false;
  int left_cap_query = 0;
  int left_cap_target = 0;
  bool right_cap = false;
  int right_cap_query = 0;
  int right_cap_target = 0;
};

struct TerminalWindowPlan {
  ValidationResult validation;
  AdjustedEndpoint first_adjusted;
  AdjustedEndpoint last_adjusted;

  int query_start = 0;  // qs0
  int target_start = 0; // rs0
  int query_end = 0;    // qe0
  int target_end = 0;   // re0

  bool has_left_packet() const {
    return query_start < first_adjusted.query &&
           target_start < first_adjusted.target;
  }
  bool has_right_packet() const {
    return last_adjusted.query < query_end && last_adjusted.target < target_end;
  }

  explicit operator bool() const {
    return validation.code == ValidationCode::Valid;
  }
};

TerminalWindowPlan plan_terminal_windows(const OrderedAnchorPath& path,
                                         const TerminalWindowControl& control);

} // namespace ordered_anchor
} // namespace lr
} // namespace cpu
} // namespace fa
