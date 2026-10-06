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

enum class DnaFamilyFailure : std::uint8_t {
  None,
  InvalidInput,
  PlacementChainRefused,
  NoOwnerChain,
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

struct DnaFamilyRealizationRequest {
  const DnaPlacementFamily* family = nullptr;
  const DnaPlacementChainingResult* placement = nullptr;
  const std::vector<std::uint8_t>* forward_query = nullptr;
  const std::vector<std::uint8_t>* reverse_query = nullptr;
};

// The -c lane's realization of one read (region_realization.h) and its work.
struct DnaFamilyRealizationOutcome {
  dna::Result output;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  DnaFamilyFailure failure = DnaFamilyFailure::InvalidInput;
  int failed_block = -1;
  int block_count = 0;
  // Realized segments. Differs from block_count when a Z-drop split a region
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
  int ksw2_attempts = 0;
  std::int64_t estimated_cells = 0;
  DnaGeometryWork geometry;
};

const char* dna_family_failure_name(DnaFamilyFailure failure) noexcept;

}  // namespace fa::cpu::lr
