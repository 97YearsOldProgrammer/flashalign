// The realization steps of one anchor path (family_realization.cpp), which the
// -c lane's controller (region_realization.h) runs on each region.
#pragma once

#include "family_realization.h"
#include "ordered_anchor_path.h"
#include "../core/cigar.h"

#include <vector>

namespace fa::cpu::lr::family_internal {

namespace ordered = ::fa::cpu::lr::ordered_anchor;

struct BlockPlan {
  int chromosome = -1;
  bool reverse = false;
  ordered::OrderedAnchorPath path;
  ordered::GeometryPlan geometry;
  // The realized interior as packed runs; text is rendered once per record.
  ::fa::cpu::output::PackedCigar interior_cigar;
  int interior_score = 0;
  int query_begin_cut = 0;
  int query_end_cut = 0;
  int target_begin_cut = 0;
  int target_end_cut = 0;

  // Z-drop split state. As minimap2's mm_split_reg, a Z-dropped segment is
  // truncated at the maximum-scoring cell and its remaining anchors form a
  // new segment and record. Segments of one region share `split_group`.
  int split_group = -1;
  bool split_continuation = false;
  // Set on a continuation whose Z-drop the inversion probe accepted as a
  // local inversion (mm_test_zdrop code 2), as minimap2's r2->split_inv.
  bool split_inv = false;
  bool truncated_right = false;
  // A continuation's own anchors, taken from the parent segment's path.
  std::vector<chaining::Anchor> continuation_anchors;
};

// One interior Z-drop. `truncated`: the segment ends at the Z-drop cell.
// `has_continuation`: the rest of the path became a new segment.
struct BlockSplit {
  bool truncated = false;
  bool has_continuation = false;
  BlockPlan continuation;
};

// A realized segment: its oriented query and reference spans, score and
// CIGAR.
struct Unit {
  int chromosome = -1;
  bool reverse = false;
  int oriented_begin = 0;
  int oriented_end = 0;
  int target_begin = 0;
  int target_end = 0;
  int score = 0;
  ::fa::cpu::output::PackedCigar core_cigar;
};

void deduplicate_exact_anchors(std::vector<chaining::Anchor>& anchors);

// Plans the segment's verified geometry and runs every step; a Z-dropped step
// splits it (split.continuation holds the rest when enough anchors remain).
bool plan_and_realize_packets(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              BlockPlan& plan, DnaFamilyFailure& failure,
                              DnaFamilyRealizationOutcome& result,
                              BlockSplit& split);

Unit unit_from_block(const BlockPlan& block);

// The unit's left and right extensions toward the read ends, the windows
// capped where `caps` sets a cap (minimap2's cap at nearby seeds).
bool extend_unit_terminals(const DnaContext& context,
                           const DnaFamilyRealizationRequest& request,
                           const BlockPlan& block, Unit& unit,
                           DnaFamilyRealizationOutcome& result,
                           const ordered::TerminalWindowControl& caps);

bool validate_and_commit_unit(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              Unit& unit, AlignResult& output);

// minimap2's mm_align1_inv between a truncated segment and its continuation.
bool realize_inversion_middle(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              const Unit& incumbent, const Unit& continuation,
                              DnaFamilyRealizationOutcome& result,
                              AlignResult& output);

}  // namespace fa::cpu::lr::family_internal
