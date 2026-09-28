#pragma once

#include "geometry_work.h"
#include "context.h"
#include "ordered_anchor_path.h"
#include "placement_family_adapter.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

struct DnaPlacementChainingResult;
struct DnaPlacementCandidateChain;

enum class DnaFamilyFailure : std::uint8_t {
  None,
  InvalidInput,
  PlacementChainRefused,
  NoSelectedFamily,
  NoAnchors,
  InvalidOrderedPath,
  InteriorRefused,
  TerminalRefused,
  FamilyValidation,
};

enum class DnaFamilyJoinKind : std::uint8_t {
  NotEligible,
  ZeroSeam,
  PureInsertion,
  PureDeletion,
  TwoAxis,
};

struct DnaFamilyJoinGeometry {
  int left_reference = -1;
  int right_reference = -1;
  bool left_reverse = false;
  bool right_reverse = false;
  int left_forward_begin = 0;
  int left_forward_end = 0;
  int right_forward_begin = 0;
  int right_forward_end = 0;
  int left_query_begin_cut = 0;
  int left_query_end_cut = 0;
  int right_query_begin_cut = 0;
  int right_query_end_cut = 0;
  int left_target_begin_cut = 0;
  int left_target_end_cut = 0;
  int right_target_begin_cut = 0;
  int right_target_end_cut = 0;
  int max_query_gap = 0;
  int max_reference_gap = 0;
  std::int64_t matrix_cell_cap = 0;
};

struct DnaFamilyJoinPlan {
  DnaFamilyJoinKind kind = DnaFamilyJoinKind::NotEligible;
  int query_gap = -1;
  int reference_gap = -1;
};

DnaFamilyJoinPlan plan_dna_family_join(
    const DnaFamilyJoinGeometry& geometry) noexcept;

// Whether the two anchors that meet at a family seam are the same anchor.
// Reverse-strand blocks run backwards in the oriented query, so their seam is
// left.begin / right.end rather than left.end / right.begin.
bool dna_family_seam_has_duplicate_anchor(
    const ordered_anchor::OrderedAnchorPath& left,
    const ordered_anchor::OrderedAnchorPath& right,
    bool reverse) noexcept;

// The anchors for a block whose oriented query interval [oriented_begin,
// oriented_end) holds none of the primary's anchors, in realization and in the
// map-only projection: those of the sibling path with the most anchors inside
// the interval (the first on ties), deduplicated; empty when fewer than two
// remain. `sibling` receives that path's index in sibling_paths, or -1.
std::vector<chaining::Anchor>
dna_selected_sibling_path(const DnaPlacementCandidateChain& evidence,
                          int oriented_begin, int oriented_end, int* sibling);

struct DnaFamilyRealizationRequest {
  const DnaPlacementFamily* family = nullptr;
  const DnaPlacementChainingResult* placement = nullptr;
  const std::vector<std::uint8_t>* forward_query = nullptr;
  const std::vector<std::uint8_t>* reverse_query = nullptr;
  // Forward-read query clip for the MAPQ's evidence-only rival realizations:
  // when clip_forward_end > clip_forward_begin every block is intersected
  // with [clip_forward_begin, clip_forward_end) before its anchors are
  // gathered. Zero on every path that emits records.
  int clip_forward_begin = 0;
  int clip_forward_end = 0;
  // Set only for the primary family: the primary record is extended past its
  // block edges toward the read ends, and block records it then covers are
  // demoted (see DnaFamilyRealizationOutcome::demoted).
  bool primary_family = false;
};

struct DnaFamilyRealizationOutcome {
  dna::Result output;
  // Primary score before seam settlement, used only to choose between
  // hypotheses. output.score is the score of the emitted CIGAR.
  int decision_score = 0;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  DnaFamilyFailure failure = DnaFamilyFailure::InvalidInput;
  int failed_block = -1;
  int block_count = 0;
  // Realized segments. Differs from block_count when a Z-drop split a block
  // into further segments, or a segment was cut at its start and realized
  // nothing.
  int segment_count = 0;
  // Interior fills truncated at a Z-drop, and how many of those continued in
  // a new segment.
  int zdrop_truncations = 0;
  int zdrop_splits = 0;
  // Local inversions (minimap2 mm_test_zdrop, mm_split_reg, mm_align1_inv).
  // Probes: reverse scans run over an interior fill's largest drop. Splits:
  // probes that found an inversion. Middles: inversion middles attempted.
  // Records: middles emitted as records.
  int inversion_probes = 0;
  int inversion_splits = 0;
  int inversion_middles = 0;
  int inversion_records = 0;
  int unit_count = 0;
  int bridge_attempts = 0;
  int bridge_accepts = 0;
  int pure_axis_bridges = 0;
  int zero_seams = 0;
  std::int64_t anchor_candidates = 0;
  int anchor_count = 0;
  int ksw2_attempts = 0;
  std::int64_t estimated_cells = 0;
  DnaGeometryWork geometry;

  // Extend-then-demote (request.primary_family). `demoted` holds the block
  // records the extended primary covers by at least half their query span,
  // removed from output.supplementary; they are emitted as secondaries with
  // MAPQ 0. `demoted_candidates` and `demoted_positions` are parallel to it,
  // the latter giving each record's index in output.supplementary had it
  // stayed. `pre_extension_primary` is the primary before the extension
  // (valid when `primary_extended`); MAPQ is scored on the family as it was
  // before this step.
  std::vector<AlignResult> demoted;
  std::vector<::fa::cpu::voting::CandidateId> demoted_candidates;
  std::vector<int> demoted_positions;
  bool primary_extended = false;
  AlignResult pre_extension_primary;


  bool accepted() const noexcept {
    return failure == DnaFamilyFailure::None && output.mapped();
  }
};

const char* dna_family_failure_name(DnaFamilyFailure failure) noexcept;

DnaFamilyRealizationOutcome realize_full_cigar_family(
    const DnaContext& context,
    const DnaFamilyRealizationRequest& request);

}  // namespace fa::cpu::lr
