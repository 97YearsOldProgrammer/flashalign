#include "family_realization.h"

#include "chain_mapq.h" // kDnaChainMapqStudyBridgeOwner
#include "cigar_geometry.h"
#include "dp_runner.h"
#include "../chaining/dense_chain.h"
#include "inv_local_chain.h"
#include "ordered_anchor_path.h"
#include "placement_chaining.h"
#include "record_family.h"
#include "../core/cigar.h"
#include "../dp/control.h"     // dp_local_score (ksw_ll) for the inversion probe
#include "../dp/ksw2_align.h"  // ksw2_simple_mat
#include "../dp/params.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

namespace ordered = ::fa::cpu::lr::ordered_anchor;
namespace realization = ::fa::cpu::lr::realization;

struct BlockPlan {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  int chromosome = -1;
  bool reverse = false;
  int forward_begin = 0;
  int forward_end = 0;
  int oriented_begin = 0;
  int oriented_end = 0;
  ordered::OrderedAnchorPath path;
  ordered::GeometryPlan geometry;
  // The realized interior as packed runs; text is rendered once per record.
  ::fa::cpu::output::PackedCigar interior_cigar;
  int interior_score = 0;
  int query_begin_cut = 0;
  int query_end_cut = 0;
  int target_begin_cut = 0;
  int target_end_cut = 0;

  // Z-drop split state. As minimap2's mm_split_reg, a Z-dropped block is
  // truncated at the maximum-scoring cell and its remaining anchors form a
  // new segment and record. Segments of one block share `split_group` and
  // never bridge back together.
  int split_group = -1;
  bool split_continuation = false;
  // Set on a continuation whose Z-drop the inversion probe accepted as a
  // local inversion (mm_test_zdrop code 2), as minimap2's r2->split_inv.
  bool split_inv = false;
  bool truncated_right = false;
  // A continuation's own anchors, taken from the parent segment's path, and
  // its lower bounds.
  int bound_query_begin = 0;
  int bound_target_begin = 0;
  std::vector<chaining::Anchor> continuation_anchors;
};

// One interior Z-drop. `truncated`: the segment ends at the Z-drop cell.
// `has_continuation`: the rest of the path became a new segment.
struct BlockSplit {
  bool truncated = false;
  bool has_continuation = false;
  BlockPlan continuation;
};

struct Unit {
  std::vector<std::size_t> blocks;
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  int chromosome = -1;
  bool reverse = false;
  int forward_begin = 0;
  int forward_end = 0;
  int oriented_begin = 0;
  int oriented_end = 0;
  int target_begin = 0;
  int target_end = 0;
  int score = 0;
  ::fa::cpu::output::PackedCigar core_cigar;
  // Set when seam settling changed this unit's CIGAR, geometry or score; an
  // untouched unit reuses its pre-settle record.
  bool settled_dirty = false;
  // The extension on that oriented side stopped at a Z-drop short of its
  // window's end; extend_primary_and_demote does not extend past it.
  bool zdrop_oriented_left = false;
  bool zdrop_oriented_right = false;
};

DpScoringParams dp_scoring_from_opts(const ResolvedDnaOptions& opts) {
  DpScoringParams scoring;
  scoring.match = opts.cigar_dp_match;
  scoring.mismatch = opts.cigar_dp_mismatch;
  scoring.ambi = opts.cigar_dp_ambi;
  scoring.gap_open1 = opts.cigar_dp_gap_open1;
  scoring.gap_extend1 = opts.cigar_dp_gap_extend1;
  scoring.gap_open2 = opts.cigar_dp_gap_open2;
  scoring.gap_extend2 = opts.cigar_dp_gap_extend2;
  scoring.zdrop = opts.cigar_dp_tail_zdrop;
  scoring.tail_end_bonus = opts.cigar_dp_tail_end_bonus;
  scoring.bw = opts.cigar_dp_bw;
  scoring.bw_long = opts.cigar_dp_bw_long;
  scoring.inversion_zdrop = opts.cigar_dp_inversion_zdrop;
  scoring.inversion_max_gap = opts.cigar_dp_max_gap;
  scoring.inversion_min_chain_score = opts.min_chain_score;
  scoring.inversion_min_dp_max = opts.cigar_dp_min_dp_max;
  return scoring;
}

// The gap fills between anchors (a block's seams, stretch pieces and verified
// regions with the gap certificate that plans them, the bridge between
// blocks, and the inversion probe and middle a fill's Z-drop test starts) run
// under the fill row, -A -B -O -E -z --score-N. The read-end extensions, seam
// settling and EXTEND run under the preset's end row (cigar_dp_*), which also
// prices every path. The inversion gates of a fill read -S scaled to the fill
// row (options/resolve.cpp).
DpScoringParams fill_dp_scoring(const ResolvedDnaOptions& opts) {
  DpScoringParams scoring = dp_scoring_from_opts(opts);
  scoring.match = opts.fill_dp_match;
  scoring.mismatch = opts.fill_dp_mismatch;
  scoring.ambi = opts.fill_dp_ambi;
  scoring.gap_open1 = opts.fill_dp_gap_open1;
  scoring.gap_extend1 = opts.fill_dp_gap_extend1;
  scoring.gap_open2 = opts.fill_dp_gap_open2;
  scoring.gap_extend2 = opts.fill_dp_gap_extend2;
  scoring.zdrop = opts.fill_dp_tail_zdrop;
  scoring.inversion_zdrop = opts.fill_dp_inversion_zdrop;
  scoring.inversion_min_dp_max = opts.fill_dp_min_dp_max;
  return scoring;
}

// Restates a fill path's score under the end row, as the kernels score a
// path: the substitution matrix per aligned base and min(O1 + E1 * L,
// O2 + E2 * L) per gap run. `query` and `target` are the packet's slices; the
// path starts at their first bases. Nothing to do when the two rows agree.
void price_fill_path(const ResolvedDnaOptions& opts,
                     const std::uint8_t* query, int query_length,
                     const std::uint8_t* target, int target_length,
                     realization::RealizationOutcome& outcome) {
  const DpScoringParams fill = fill_dp_scoring(opts);
  const DpScoringParams end = dp_scoring_from_opts(opts);
  if (fill.match == end.match && fill.mismatch == end.mismatch &&
      fill.ambi == end.ambi && fill.gap_open1 == end.gap_open1 &&
      fill.gap_extend1 == end.gap_extend1 &&
      fill.gap_open2 == end.gap_open2 && fill.gap_extend2 == end.gap_extend2)
    return;
  const std::vector<std::uint32_t>& path = outcome.raw.packed_cigar;
  int query_used = 0;
  int target_used = 0;
  for (const std::uint32_t run : path) {
    const int op = static_cast<int>(run & 0xf);
    const int length = static_cast<int>(run >> 4);
    if (op == 0 || op == 1)
      query_used += length;
    if (op != 1)
      target_used += length;
  }
  if (query_used > query_length || target_used > target_length)
    return;
  std::array<std::int8_t, 25> matrix{};
  ::fa::cpu::ksw2_simple_mat(matrix.data(), end.match, end.mismatch, end.ambi);
  const auto gap = [&end](std::int64_t length) {
    return std::min(end.gap_open1 + end.gap_extend1 * length,
                    end.gap_open2 + end.gap_extend2 * length);
  };
  int q = 0;
  int t = 0;
  std::int64_t score = 0;
  for (const std::uint32_t run : path) {
    const int op = static_cast<int>(run & 0xf);
    const int length = static_cast<int>(run >> 4);
    if (op == 0) {
      for (int k = 0; k < length; ++k)
        score += matrix[static_cast<std::size_t>(
            std::min<int>(target[t + k], 4) * 5 +
            std::min<int>(query[q + k], 4))];
      q += length;
      t += length;
    } else if (op == 1) {
      score -= gap(length);
      q += length;
    } else {
      score -= gap(length);
      t += length;
    }
  }
  outcome.alignment.score = static_cast<int>(score);
}

void add_saturating(std::int64_t value, std::int64_t& total) {
  if (value <= 0 || total == std::numeric_limits<std::int64_t>::max())
    return;
  if (value > std::numeric_limits<std::int64_t>::max() - total)
    total = std::numeric_limits<std::int64_t>::max();
  else
    total += value;
}

void record_kernel_work(const realization::RealizationOutcome& outcome,
                        int query_span, int target_span,
                        DnaFamilyRealizationOutcome& result) {
  // A verified region runs no kernel; it counts bases, not DP.
  if (outcome.committed_executor == realization::ExecutorId::Verified) {
    ++result.geometry.verified_regions;
    result.geometry.verified_bases += outcome.trace.verified.bases;
    return;
  }
  for (std::size_t index = 0; index < outcome.trace.attempt_count; ++index) {
    const realization::KernelAttempt& attempt = outcome.trace.attempts[index];
    if (attempt.executor != realization::ExecutorId::Ksw2)
      continue;
    ++result.ksw2_attempts;
    const std::int64_t query_cells =
        static_cast<std::int64_t>(std::max(0, query_span)) + 1;
    const std::int64_t target_cells =
        static_cast<std::int64_t>(std::max(0, target_span)) + 1;
    const std::int64_t band_width =
        attempt.band < 0
            ? target_cells
            : std::min(target_cells,
                       static_cast<std::int64_t>(2) * attempt.band + 1);
    if (query_cells > 0 && band_width > 0 &&
        query_cells <= std::numeric_limits<std::int64_t>::max() / band_width) {
      add_saturating(query_cells * band_width, result.estimated_cells);
    } else {
      result.estimated_cells = std::numeric_limits<std::int64_t>::max();
    }
  }
}

// Counts mm_test_zdrop's inversion test: the reverse probe runs when the
// largest drop exceeds zdrop_inv, and code 2 means an inversion was found.
void record_inversion_probe(const realization::RealizationOutcome& outcome,
                            DnaFamilyRealizationOutcome& result) {
  if (outcome.trace.zdrop_test.reverse_probe_attempted)
    ++result.inversion_probes;
  if (outcome.trace.zdrop_test.code == 2)
    ++result.inversion_splits;
}

int query_begin_for_tile(const DnaPlacementFamily& family, int tile) {
  return ::fa::cpu::voting::query_tile_begin(tile, family.read_length,
                                             family.seed_length);
}

bool exact_anchor_equal(const chaining::Anchor& left,
                        const chaining::Anchor& right) {
  return left.q == right.q && left.r == right.r && left.span == right.span;
}

void deduplicate_exact_anchors(std::vector<chaining::Anchor>& anchors) {
  anchors.erase(std::unique(anchors.begin(), anchors.end(), exact_anchor_equal),
                anchors.end());
}

bool materialize_plans(const DnaContext& context,
                       const DnaFamilyRealizationRequest& request,
                       std::vector<BlockPlan>& plans) {
  const DnaPlacementFamily& family = *request.family;
  for (const auto& block : family.partition.selected.blocks) {
    if (block.candidate == ::fa::cpu::voting::kNullCandidate)
      continue;
    const DnaPlacementCandidate* candidate = family.find(block.candidate);
    if (candidate == nullptr)
      return false;
    BlockPlan plan;
    plan.candidate = block.candidate;
    plan.chromosome = candidate->peak.chr;
    plan.reverse = candidate->peak.is_rc;
    plan.forward_begin = query_begin_for_tile(family, block.query_tile_begin);
    plan.forward_end =
        block.query_tile_end == ::fa::cpu::voting::kQueryTileCount
            ? family.read_length
            : query_begin_for_tile(family, block.query_tile_end);
    // The evidence-only query clip. A single-candidate family's block spans
    // the whole read, so the clip is applied to the plan; everything
    // downstream reads these bounds. A block outside the window fails below.
    if (request.clip_forward_end > request.clip_forward_begin) {
      plan.forward_begin =
          std::max(plan.forward_begin, request.clip_forward_begin);
      plan.forward_end = std::min(plan.forward_end, request.clip_forward_end);
    }
    plan.oriented_begin = plan.reverse ? family.read_length - plan.forward_end
                                       : plan.forward_begin;
    plan.oriented_end = plan.reverse ? family.read_length - plan.forward_begin
                                     : plan.forward_end;
    if (plan.chromosome < 0 || plan.forward_begin < 0 ||
        plan.forward_end <= plan.forward_begin ||
        plan.forward_end > family.read_length || plan.oriented_begin < 0 ||
        plan.oriented_end <= plan.oriented_begin ||
        block.supporting_tiles < kDnaPostCommitRecordMinBlockTiles)
      return false;
    plans.push_back(std::move(plan));
  }
  return true;
}

// Splits a block at a Z-dropped seam or piece, as minimap2's mm_align1: the
// returned traceback, ending at the maximum-scoring cell, is kept as the
// segment's prefix and the anchors past that cell form the continuation. A
// fill dropped at its start (over the matrix cap) cuts at the step start.
// Returns false only on a malformed truncation.
bool plan_block_split(const DnaContext& context,
                      const DnaFamilyRealizationRequest& request,
                      BlockPlan& plan, const ordered::GeometryStep& step,
                      const realization::RealizationOutcome& outcome,
                      int query_cursor, int target_cursor, BlockSplit& split) {
  const int q_span = step.query_end - step.query_begin;
  const int r_span = step.target_end - step.target_begin;
  const int q_used =
      ::fa::cpu::output::packed_cigar_query_consumed(outcome.raw.packed_cigar);
  const int r_used = outcome.alignment.ref_consumed;
  // The truncated CIGAR must consume exactly max_query + 1 by
  // max_target + 1 (0 by 0 when the fill dropped at its start).
  if (outcome.alignment.ref_offset != 0 || q_used > q_span ||
      r_used > r_span || outcome.raw.max_query + 1 != q_used ||
      outcome.raw.max_target + 1 != r_used)
    return false;
  // A cut at the block start leaves the segment empty; the segment loop
  // drops it.
  const int cut_query = query_cursor + q_used;
  const int cut_target = target_cursor + r_used;

  ::fa::cpu::output::append_packed_cigar(plan.interior_cigar,
                                         outcome.raw.packed_cigar);
  plan.interior_score += outcome.alignment.score;
  plan.query_end_cut = cut_query;
  plan.target_end_cut = cut_target;
  plan.truncated_right = true;
  split.truncated = true;

  // Split after the last anchor ending at or before the cut, as minimap2.
  // Records here must be query-disjoint, so the continuation also skips
  // anchors reaching back over the cut.
  const std::size_t original_begin = plan.path.realization_begin;
  const std::size_t original_end = plan.path.realization_end;
  std::size_t start = original_begin + 1;
  for (std::size_t index = step.selected_index; index > original_begin;
       --index) {
    if (plan.path.selected[index - 1].r_end() <= cut_target) {
      start = index;
      break;
    }
  }
  while (start < original_end && (plan.path.selected[start].q < cut_query ||
                                  plan.path.selected[start].r < cut_target))
    ++start;
  if (start > original_end)
    start = original_end;
  // The segment's path now ends at the cut.
  plan.path.realization_end = start;

  const int remaining = static_cast<int>(original_end - start);
  // With fewer than opt->min_cnt anchors left there is no continuation, as
  // in minimap2.
  if (remaining < std::max(1, context.opts.cigar_dp_split_min_anchors))
    return true;

  const int read_length = request.family->read_length;
  BlockPlan& continuation = split.continuation;
  continuation.candidate = plan.candidate;
  continuation.chromosome = plan.chromosome;
  continuation.reverse = plan.reverse;
  continuation.split_group = plan.split_group;
  continuation.split_continuation = true;
  continuation.split_inv = outcome.trace.zdrop_test.code == 2;
  continuation.bound_query_begin = cut_query;
  continuation.bound_target_begin = cut_target;
  continuation.oriented_begin = cut_query;
  continuation.oriented_end = plan.oriented_end;
  // On the reverse strand the continuation is the lower forward interval.
  if (plan.reverse) {
    continuation.forward_begin = plan.forward_begin;
    continuation.forward_end = read_length - cut_query;
    plan.forward_begin = read_length - cut_query;
  } else {
    continuation.forward_begin = cut_query;
    continuation.forward_end = plan.forward_end;
    plan.forward_end = cut_query;
  }
  plan.oriented_end = cut_query;
  continuation.continuation_anchors.assign(
      plan.path.selected.begin() + static_cast<std::ptrdiff_t>(start),
      plan.path.selected.begin() + static_cast<std::ptrdiff_t>(original_end));
  split.has_continuation = true;
  return true;
}

std::vector<chaining::Anchor>
selected_sibling_path(const DnaPlacementCandidateChain& evidence,
                      int oriented_begin, int oriented_end,
                      int* sibling = nullptr) {
  std::vector<chaining::Anchor> best;
  std::vector<chaining::Anchor> inside;
  const std::vector<chaining::Anchor>* chosen = nullptr;
  if (sibling != nullptr)
    *sibling = -1;
  for (const std::vector<chaining::Anchor>& path : evidence.sibling_paths) {
    inside.clear();
    for (const chaining::Anchor& anchor : path) {
      if (anchor.q >= oriented_begin && anchor.q_end() <= oriented_end)
        inside.push_back(anchor);
    }
    if (inside.size() > best.size()) {
      best = inside;
      chosen = &path;
    }
  }
  if (best.size() < 2)
    return {};
  deduplicate_exact_anchors(best);
  if (best.size() < 2)
    return {};
  if (sibling != nullptr)
    *sibling = static_cast<int>(chosen - evidence.sibling_paths.data());
  return best;
}

// Plans the block's verified geometry and runs every step: a verified region
// on the Verified executor, a seam or piece as one gap fill.
bool plan_and_realize_packets(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              BlockPlan& plan, DnaFamilyFailure& failure,
                              DnaFamilyRealizationOutcome& result,
                              BlockSplit& split) {
  const auto& query =
      plan.reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(plan.chromosome)];
  const DpScoringParams scoring = fill_dp_scoring(context.opts);
  ::fa::cpu::DpMapOpt map_opt;
  map_opt.bw = context.opts.cigar_dp_bw;
  map_opt.bw_long = context.opts.cigar_dp_bw_long;
  ordered::GeometryControl control;
  control.normal_long_band = map_opt.bw_long_eff();
  control.min_ksw_len = context.opts.cigar_dp_min_ksw_len;
  // The same matrix the ksw2 executor uses.
  ordered::GapCertificate certificate;
  certificate.query = query.data();
  certificate.target = reference.data();
  ::fa::cpu::ksw2_simple_mat(certificate.matrix.data(), scoring.match,
                             scoring.mismatch, scoring.ambi);
  certificate.match = scoring.match;
  certificate.one_base_gap = std::min(scoring.gap_open1 + scoring.gap_extend1,
                                      scoring.gap_open2 + scoring.gap_extend2);
  plan.geometry =
      ordered::plan_verified_geometry(plan.path, control, certificate);
  plan.query_begin_cut = plan.geometry.initial_query_cursor;
  plan.query_end_cut = plan.geometry.final_query_cursor;
  plan.target_begin_cut = plan.geometry.initial_target_cursor;
  plan.target_end_cut = plan.geometry.final_target_cursor;

  int q_cursor = plan.query_begin_cut;
  int r_cursor = plan.target_begin_cut;
  // HiFi presets: a seam's or piece's late inversion probe runs only over a
  // local chain of the read's opposite-lane seeds in its drop window. A
  // bridge's probe, a veto, is not gated.
  const bool gate_probe = context.opts.inversion_probe_local_gate;
  DnaInvLocalChainSource local_source;
  realization::InversionProbeGate probe_gate;
  if (gate_probe) {
    // Placement keeps the read's seed list for the gate (map_read); checked
    // in release builds too.
    if (context.inversion_gate_seeds == nullptr) {
      std::fputs("flashalign: internal error: the inversion probe gate has "
                 "no seed list\n",
                 stderr);
      std::abort();
    }
    local_source.seeds = context.inversion_gate_seeds;
    local_source.index = context.ref.index;
    local_source.chromosome = plan.chromosome;
    local_source.block_reverse = plan.reverse;
    local_source.read_length = request.family->read_length;
    local_source.seed_length = request.family->seed_length;
    // The whole-query harvest's occurrence gate.
    local_source.occurrence_cap =
        context.opts.dna_pool_gate_occ > 0
            ? static_cast<std::uint32_t>(context.opts.dna_pool_gate_occ)
            : 0u;
    probe_gate.count = dna_inv_local_chain;
    probe_gate.source = &local_source;
    probe_gate.min_count = kDnaInvLocalMinAnchors;
  }
  for (const ordered::GeometryStep& step : plan.geometry.steps) {
    const int q_span = step.query_end - step.query_begin;
    const int r_span = step.target_end - step.target_begin;
    // A pure-axis step is materialized by the controller without an executor.
    const bool two_axis = q_span > 0 && r_span > 0;
    const bool verified = step.kind == ordered::GeometryStepKind::Verified;
    const bool gated = gate_probe && two_axis && !verified;
    if (gated) {
      local_source.query_offset = step.query_begin;
      local_source.target_offset = step.target_begin;
    }
    VerifiedRegionInputs region;
    region.gaps = step.gap_count > 0
                      ? plan.geometry.gaps.data() + step.gap_begin
                      : nullptr;
    region.gap_count = static_cast<int>(step.gap_count);
    const VerifiedRegionInputs* verified_region =
        verified && two_axis ? &region : nullptr;
    if (step.long_join && two_axis)
      ++result.geometry.join_calls;
    // Seams and pieces are counted here, verified regions in
    // record_kernel_work; a stretch is counted at its first piece.
    const ordered::GeometryStretch* first_piece_of = nullptr;
    if (step.kind == ordered::GeometryStepKind::Seam) {
      ++result.geometry.seams;
    } else if (step.kind == ordered::GeometryStepKind::Piece) {
      ++result.geometry.stretch_pieces;
      if (step.piece == 0) {
        first_piece_of =
            &plan.geometry.stretches[static_cast<std::size_t>(step.stretch)];
        ++result.geometry.stretches;
        result.geometry.stretch_bases += first_piece_of->query_bases;
        result.geometry.failed_regions += first_piece_of->failed_regions;
        result.geometry.failed_region_anchors +=
            first_piece_of->failed_anchors;
        result.geometry.failed_region_bases += first_piece_of->failed_bases;
        result.geometry.failing_gaps += first_piece_of->failing_gaps;
      }
    }
    realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaInternalFill, scoring,
        realization::make_query_slice(
            q_span > 0 ? query.data() + step.query_begin : nullptr, q_span,
            step.query_begin, plan.path.query_length, plan.reverse),
        realization::make_target_slice(
            r_span > 0 ? reference.data() + step.target_begin : nullptr, r_span,
            step.target_begin),
        two_axis ? step.selected_band : -1, scoring.zdrop, step.long_join, nullptr,
        /*inversion_probe_enabled=*/true, verified_region, gated ? &probe_gate : nullptr);
    price_fill_path(context.opts, query.data() + step.query_begin, q_span,
                    reference.data() + step.target_begin, r_span, outcome);
    record_kernel_work(outcome, q_span, r_span, result);
    record_inversion_probe(outcome, result);
    const bool successful =
        outcome.kind == realization::OutcomeKind::Aligned ||
        outcome.kind == realization::OutcomeKind::PureInsertion ||
        outcome.kind == realization::OutcomeKind::PureDeletion;
    if (!successful) {
      if (outcome.kind == realization::OutcomeKind::SplitRequested) {
        if (plan_block_split(context, request, plan, step, outcome, q_cursor,
                             r_cursor, split))
          return true;
      }
      failure = DnaFamilyFailure::InteriorRefused;
      return false;
    }
    // A successful step consumes exactly its spans.
    ::fa::cpu::output::append_packed_cigar(plan.interior_cigar,
                                           outcome.raw.packed_cigar);
    plan.interior_score += outcome.alignment.score;
    q_cursor = step.query_end;
    r_cursor = step.target_end;
  }
  // The steps tile the block, so the cursor must end at the end cut.
  if (q_cursor != plan.query_end_cut || r_cursor != plan.target_end_cut) {
    failure = DnaFamilyFailure::InvalidOrderedPath;
    return false;
  }
  return true;
}

bool acquire_and_realize_block(const DnaContext& context,
                               const DnaFamilyRealizationRequest& request,
                               BlockPlan& plan, bool trim_left, bool trim_right,
                               DnaFamilyFailure& failure,
                               DnaFamilyRealizationOutcome& result,
                               BlockSplit& split) {
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(plan.chromosome)];
  const DnaPlacementCandidateChain* evidence =
      request.placement == nullptr ? nullptr
                                   : request.placement->find(plan.candidate);
  if (evidence == nullptr || !evidence->exact ||
      evidence->status != DnaPlacementChainStatus::Accepted) {
    failure = DnaFamilyFailure::PlacementChainRefused;
    return false;
  }
  std::vector<chaining::Anchor> anchors;
  if (plan.split_continuation) {
    // Already deduplicated and counted by the parent segment.
    anchors = std::move(plan.continuation_anchors);
  } else {
    anchors.reserve(evidence->primary.size());
    for (const chaining::Anchor& anchor : evidence->primary) {
      if (anchor.q >= plan.oriented_begin &&
          anchor.q_end() <= plan.oriented_end)
        anchors.push_back(anchor);
    }
    deduplicate_exact_anchors(anchors);
    if (anchors.empty()) {
      anchors = selected_sibling_path(*evidence, plan.oriented_begin,
                                      plan.oriented_end);
    }
    result.anchor_candidates += evidence->interval_hits;
    result.anchor_count += static_cast<int>(anchors.size());
  }
  if (anchors.empty()) {
    failure = DnaFamilyFailure::NoAnchors;
    return false;
  }

  ordered::OrderedAnchorPath raw;
  raw.reference_id = plan.chromosome;
  raw.reverse_complemented = plan.reverse;
  raw.query_length = request.family->read_length;
  raw.target_length = static_cast<int>(reference.size());
  // A continuation starts at the Z-drop cut, so its left terminal window
  // stays out of the previous segment.
  raw.query_bound_begin =
      plan.split_continuation ? plan.bound_query_begin : plan.oriented_begin;
  raw.query_bound_end = plan.oriented_end;
  raw.target_bound_begin =
      plan.split_continuation ? plan.bound_target_begin : 0;
  raw.target_bound_end = static_cast<int>(reference.size());
  raw.minimizer_k = context.ref.index->k();
  raw.selected = std::move(anchors);
  raw.realization_end = raw.selected.size();

  // As mm_fix_bad_ends(r, a, bw, min_chain_score * 2, ...).
  ordered::NormalizationControl normalization{
      std::max(1, context.opts.cigar_dp_bw), context.opts.min_chain_score,
      std::max(1, context.opts.cigar_dp_max_gap)};
  normalization.trim_left = trim_left;
  normalization.trim_right = trim_right;
  const ordered::NormalizedPath normalized =
      ordered::normalize_for_realization(raw, normalization);
  if (!normalized) {
    failure = DnaFamilyFailure::InvalidOrderedPath;
    return false;
  }
  plan.path = normalized.path;
  return plan_and_realize_packets(context, request, plan, failure, result,
                                  split);
}

bool seam_has_duplicate_anchor(const BlockPlan& left, const BlockPlan& right) {
  return dna_family_seam_has_duplicate_anchor(left.path, right.path,
                                              left.reverse);
}

bool bridge_geometry(const DnaContext& context, const BlockPlan& left,
                     const BlockPlan& right, int& q_gap, int& r_gap) {
  if (seam_has_duplicate_anchor(left, right))
    return false;
  ::fa::cpu::DpMapOpt opt;
  const DnaFamilyJoinPlan plan = plan_dna_family_join(DnaFamilyJoinGeometry{
      left.chromosome, right.chromosome, left.reverse, right.reverse,
      left.forward_begin, left.forward_end, right.forward_begin,
      right.forward_end, left.query_begin_cut, left.query_end_cut,
      right.query_begin_cut, right.query_end_cut, left.target_begin_cut,
      left.target_end_cut, right.target_begin_cut, right.target_end_cut,
      context.opts.cigar_dp_max_gap, context.opts.cigar_dp_max_gap,
      opt.max_sw_mat});
  q_gap = plan.query_gap;
  r_gap = plan.reference_gap;
  return plan.kind != DnaFamilyJoinKind::NotEligible;
}

bool attempt_bridge(const DnaContext& context,
                    const DnaFamilyRealizationRequest& request,
                    const BlockPlan& left, const BlockPlan& right,
                    ::fa::cpu::output::PackedCigar& cigar, int& score,
                    DnaFamilyRealizationOutcome& result) {
  int q_gap = 0;
  int r_gap = 0;
  if (!bridge_geometry(context, left, right, q_gap, r_gap))
    return false;
  ++result.bridge_attempts;
  if (q_gap == 0 && r_gap == 0) {
    ++result.zero_seams;
    ++result.bridge_accepts;
    return true;
  }

  const bool reverse = left.reverse;
  const int q_begin = reverse ? right.query_end_cut : left.query_end_cut;
  const int r_begin = reverse ? right.target_end_cut : left.target_end_cut;
  const auto& query = reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(left.chromosome)];
  const DpScoringParams scoring = fill_dp_scoring(context.opts);
  const bool two_axis = q_gap > 0 && r_gap > 0;
  realization::RealizationOutcome outcome = run_dna_long_realization(
      realization::RealizationRole::DnaInternalFill, scoring,
      realization::make_query_slice(
          q_gap > 0 ? query.data() + q_begin : nullptr, q_gap, q_begin,
          request.family->read_length, reverse),
      realization::make_target_slice(
          r_gap > 0 ? reference.data() + r_begin : nullptr, r_gap, r_begin),
      two_axis ? std::max(q_gap, r_gap) : -1, scoring.zdrop,
      /*long_join=*/two_axis, nullptr, /*inversion_probe_enabled=*/true);
  price_fill_path(context.opts, query.data() + q_begin, q_gap,
                  reference.data() + r_begin, r_gap, outcome);
  record_kernel_work(outcome, q_gap, r_gap, result);
  record_inversion_probe(outcome, result);
  if (outcome.trace.zdrop_test.code == 2)
    return false;
  const bool expected =
      two_axis    ? outcome.kind == realization::OutcomeKind::Aligned
      : q_gap > 0 ? outcome.kind == realization::OutcomeKind::PureInsertion
                  : outcome.kind == realization::OutcomeKind::PureDeletion;
  if (!expected ||
      ::fa::cpu::output::packed_cigar_query_consumed(
          outcome.raw.packed_cigar) != q_gap ||
      outcome.alignment.ref_offset != 0 ||
      outcome.alignment.ref_consumed != r_gap)
    return false;
  cigar = outcome.raw.packed_cigar;
  score = outcome.alignment.score;
  ++result.bridge_accepts;
  if (!two_axis)
    ++result.pure_axis_bridges;
  return true;
}

Unit unit_from_block(const BlockPlan& block, std::size_t index) {
  Unit unit;
  unit.blocks.push_back(index);
  unit.candidate = block.candidate;
  unit.chromosome = block.chromosome;
  unit.reverse = block.reverse;
  unit.forward_begin = block.forward_begin;
  unit.forward_end = block.forward_end;
  unit.oriented_begin = block.query_begin_cut;
  unit.oriented_end = block.query_end_cut;
  unit.target_begin = block.target_begin_cut;
  unit.target_end = block.target_end_cut;
  unit.score = block.interior_score;
  unit.core_cigar = block.interior_cigar;
  return unit;
}

void accept_join(Unit& unit, const BlockPlan& current,
                 std::size_t current_index,
                 const ::fa::cpu::output::PackedCigar& bridge,
                 int bridge_score) {
  if (!unit.reverse) {
    ::fa::cpu::output::append_packed_cigar(unit.core_cigar, bridge);
    ::fa::cpu::output::append_packed_cigar(unit.core_cigar,
                                           current.interior_cigar);
    unit.oriented_end = current.query_end_cut;
    unit.target_end = current.target_end_cut;
  } else {
    ::fa::cpu::output::PackedCigar joined = current.interior_cigar;
    ::fa::cpu::output::append_packed_cigar(joined, bridge);
    ::fa::cpu::output::prepend_packed_cigar(unit.core_cigar, std::move(joined));
    unit.oriented_begin = current.query_begin_cut;
    unit.target_begin = current.target_begin_cut;
  }
  unit.blocks.push_back(current_index);
  unit.forward_end = current.forward_end;
  unit.score += bridge_score + current.interior_score;
}

// An extension the Z-drop stopped before its query window ran out. A Zdropped
// extension whose best cell is the window's last query base (the drop fell in
// the reference overhang) did reach its edge.
bool zdrop_stopped_short(const realization::RealizationOutcome& outcome,
                         int query_used, int query_span) {
  return outcome.kind == realization::OutcomeKind::Zdropped &&
         query_used < query_span;
}

bool extend_unit_terminals(const DnaContext& context,
                           const DnaFamilyRealizationRequest& request,
                           const std::vector<BlockPlan>& blocks, Unit& unit,
                           DnaFamilyRealizationOutcome& result) {
  const std::size_t reference_left_index =
      unit.reverse ? unit.blocks.back() : unit.blocks.front();
  const std::size_t reference_right_index =
      unit.reverse ? unit.blocks.front() : unit.blocks.back();
  const BlockPlan& left_block = blocks[reference_left_index];
  const BlockPlan& right_block = blocks[reference_right_index];
  const auto& query =
      unit.reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(unit.chromosome)];
  const DpScoringParams scoring = dp_scoring_from_opts(context.opts);
  const ordered::TerminalWindowControl terminal_control{
      std::max(1, context.opts.cigar_dp_max_gap), scoring.match,
      scoring.gap_open1, scoring.gap_extend1};
  const ordered::TerminalWindowPlan left_terminals =
      ordered::plan_terminal_windows(left_block.path, terminal_control);
  const ordered::TerminalWindowPlan right_terminals =
      ordered::plan_terminal_windows(right_block.path, terminal_control);
  if (!left_terminals || !right_terminals)
    return false;

  // As minimap2 (r->split_inv ? zdrop_inv : zdrop): the left terminal of a
  // continuation is its Z-drop seam, and after an inversion split it extends
  // at the tighter inversion Z-drop so it does not run through the inversion.
  DpScoringParams seam_scoring = scoring;
  const bool seam_is_inversion =
      left_block.split_continuation && left_block.split_inv;
  if (seam_is_inversion)
    seam_scoring.zdrop = scoring.inversion_zdrop;

  if (left_terminals.has_left_packet()) {
    const int q_span =
        left_terminals.first_adjusted.query - left_terminals.query_start;
    const int r_span =
        left_terminals.first_adjusted.target - left_terminals.target_start;
    const realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaLeftExtension, seam_scoring,
        realization::make_query_slice(query.data() + left_terminals.query_start,
                                      q_span, left_terminals.query_start,
                                      request.family->read_length,
                                      unit.reverse),
        realization::make_target_slice(reference.data() +
                                           left_terminals.target_start,
                                       r_span, left_terminals.target_start),
        context.opts.cigar_dp_bw);
    record_kernel_work(outcome, q_span, r_span, result);
    const int q_used = ::fa::cpu::output::packed_cigar_query_consumed(
        outcome.raw.packed_cigar);
    if (zdrop_stopped_short(outcome, q_used, q_span))
      unit.zdrop_oriented_left = true;
    if ((outcome.kind == realization::OutcomeKind::Aligned ||
         outcome.kind == realization::OutcomeKind::Zdropped) &&
        !outcome.raw.packed_cigar.empty() && q_used > 0 && q_used <= q_span &&
        outcome.alignment.ref_offset == 0 &&
        outcome.alignment.ref_consumed >= 0 &&
        outcome.alignment.ref_consumed <= r_span) {
      ::fa::cpu::output::prepend_packed_cigar(unit.core_cigar,
                                              outcome.raw.packed_cigar);
      unit.oriented_begin = left_terminals.first_adjusted.query - q_used;
      unit.target_begin =
          left_terminals.first_adjusted.target - outcome.alignment.ref_consumed;
      unit.score += outcome.alignment.score;
    }
  }

  // As minimap2, a Z-dropped block gets no right extension: it ends at the
  // truncation cell.
  if (right_terminals.has_right_packet() && !right_block.truncated_right) {
    const int q_span =
        right_terminals.query_end - right_terminals.last_adjusted.query;
    const int r_span =
        right_terminals.target_end - right_terminals.last_adjusted.target;
    const realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaRightExtension, scoring,
        realization::make_query_slice(
            query.data() + right_terminals.last_adjusted.query, q_span,
            right_terminals.last_adjusted.query, request.family->read_length,
            unit.reverse),
        realization::make_target_slice(
            reference.data() + right_terminals.last_adjusted.target, r_span,
            right_terminals.last_adjusted.target),
        context.opts.cigar_dp_bw);
    record_kernel_work(outcome, q_span, r_span, result);
    const int q_used = ::fa::cpu::output::packed_cigar_query_consumed(
        outcome.raw.packed_cigar);
    if (zdrop_stopped_short(outcome, q_used, q_span))
      unit.zdrop_oriented_right = true;
    if ((outcome.kind == realization::OutcomeKind::Aligned ||
         outcome.kind == realization::OutcomeKind::Zdropped) &&
        !outcome.raw.packed_cigar.empty() && q_used > 0 && q_used <= q_span &&
        outcome.alignment.ref_offset == 0 &&
        outcome.alignment.ref_consumed >= 0 &&
        outcome.alignment.ref_consumed <= r_span) {
      ::fa::cpu::output::append_packed_cigar(unit.core_cigar,
                                             outcome.raw.packed_cigar);
      unit.oriented_end = right_terminals.last_adjusted.query + q_used;
      unit.target_end =
          right_terminals.last_adjusted.target + outcome.alignment.ref_consumed;
      unit.score += outcome.alignment.score;
    }
  }
  return true;
}

int unit_forward_begin(const Unit& unit, int read_length) {
  return unit.reverse ? read_length - unit.oriented_end : unit.oriented_begin;
}

int unit_forward_end(const Unit& unit, int read_length) {
  return unit.reverse ? read_length - unit.oriented_begin : unit.oriented_end;
}

// The reference budget of a query span, sized as the terminal windows are:
// its length plus the gap its best match reward can pay for, capped at
// max_gap.
int seam_reference_budget(const DpScoringParams& scoring, int query_span,
                          int max_gap) {
  std::int64_t span = query_span;
  const std::int64_t reward = span * static_cast<std::int64_t>(scoring.match);
  if (reward > scoring.gap_open1 && scoring.gap_extend1 > 0)
    span += (reward - scoring.gap_open1) / scoring.gap_extend1;
  span = std::min<std::int64_t>(span, std::max(1, max_gap));
  return static_cast<int>(
      std::min<std::int64_t>(span, std::numeric_limits<int>::max()));
}

// Extends one end of a unit by an extension DP over at most
// min(budget, max_gap) query bases and seam_reference_budget reference bases.
// Returns the query bases added (0 when refused, leaving the unit unchanged)
// and marks the side when a Z-drop stopped the DP short. Used by seam
// settling and by the primary's extension toward the read ends.
int extend_unit_side(const DnaContext& context,
                     const DnaFamilyRealizationRequest& request,
                     const std::vector<BlockPlan>& blocks, Unit& unit,
                     bool forward_right, int budget,
                     DnaFamilyRealizationOutcome& result) {
  if (budget <= 0)
    return 0;
  const bool oriented_right = unit.reverse ? !forward_right : forward_right;
  const BlockPlan& facing =
      blocks[forward_right ? unit.blocks.back() : unit.blocks.front()];
  // Never extend past a Z-drop truncation or back across a split.
  if (oriented_right ? facing.truncated_right : facing.split_continuation)
    return 0;

  const auto& query =
      unit.reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(unit.chromosome)];
  const DpScoringParams scoring = dp_scoring_from_opts(context.opts);
  const int max_gap = std::max(1, context.opts.cigar_dp_max_gap);
  const int capped_budget = std::min(budget, max_gap);
  const int reference_budget =
      seam_reference_budget(scoring, capped_budget, max_gap);

  if (oriented_right) {
    const int query_span = std::min(capped_budget, request.family->read_length -
                                                       unit.oriented_end);
    const int reference_span = std::min(
        reference_budget, static_cast<int>(reference.size()) - unit.target_end);
    if (query_span <= 0 || reference_span <= 0)
      return 0;
    const realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaRightExtension, scoring,
        realization::make_query_slice(
            query.data() + unit.oriented_end, query_span, unit.oriented_end,
            request.family->read_length, unit.reverse),
        realization::make_target_slice(reference.data() + unit.target_end,
                                       reference_span, unit.target_end),
        context.opts.cigar_dp_bw);
    record_kernel_work(outcome, query_span, reference_span, result);
    const int query_used = ::fa::cpu::output::packed_cigar_query_consumed(
        outcome.raw.packed_cigar);
    if (zdrop_stopped_short(outcome, query_used, query_span))
      unit.zdrop_oriented_right = true;
    if ((outcome.kind != realization::OutcomeKind::Aligned &&
         outcome.kind != realization::OutcomeKind::Zdropped) ||
        outcome.raw.packed_cigar.empty() || query_used <= 0 ||
        query_used > query_span || outcome.alignment.ref_offset != 0 ||
        outcome.alignment.ref_consumed < 0 ||
        outcome.alignment.ref_consumed > reference_span) {
      return 0;
    }
    ::fa::cpu::output::append_packed_cigar(unit.core_cigar,
                                           outcome.raw.packed_cigar);
    unit.oriented_end += query_used;
    unit.target_end += outcome.alignment.ref_consumed;
    unit.score += outcome.alignment.score;
    unit.settled_dirty = true;
    return query_used;
  }

  const int query_span = std::min(capped_budget, unit.oriented_begin);
  const int reference_span = std::min(reference_budget, unit.target_begin);
  if (query_span <= 0 || reference_span <= 0)
    return 0;
  const int query_begin = unit.oriented_begin - query_span;
  const int reference_begin = unit.target_begin - reference_span;
  const realization::RealizationOutcome outcome = run_dna_long_realization(
      realization::RealizationRole::DnaLeftExtension, scoring,
      realization::make_query_slice(query.data() + query_begin, query_span,
                                    query_begin, request.family->read_length,
                                    unit.reverse),
      realization::make_target_slice(reference.data() + reference_begin,
                                     reference_span, reference_begin),
      context.opts.cigar_dp_bw);
  record_kernel_work(outcome, query_span, reference_span, result);
  const int query_used = ::fa::cpu::output::packed_cigar_query_consumed(
      outcome.raw.packed_cigar);
  if (zdrop_stopped_short(outcome, query_used, query_span))
    unit.zdrop_oriented_left = true;
  if ((outcome.kind != realization::OutcomeKind::Aligned &&
       outcome.kind != realization::OutcomeKind::Zdropped) ||
      outcome.raw.packed_cigar.empty() || query_used <= 0 ||
      query_used > query_span || outcome.alignment.ref_offset != 0 ||
      outcome.alignment.ref_consumed < 0 ||
      outcome.alignment.ref_consumed > reference_span) {
    return 0;
  }
  ::fa::cpu::output::prepend_packed_cigar(unit.core_cigar,
                                          outcome.raw.packed_cigar);
  unit.oriented_begin -= query_used;
  unit.target_begin -= outcome.alignment.ref_consumed;
  unit.score += outcome.alignment.score;
  unit.settled_dirty = true;
  return query_used;
}

// Extend one end into an unclaimed query seam (budget = the seam's gap).
int settle_seam_side(const DnaContext& context,
                     const DnaFamilyRealizationRequest& request,
                     const std::vector<BlockPlan>& blocks, Unit& unit,
                     bool forward_right, int budget,
                     DnaFamilyRealizationOutcome& result) {
  return extend_unit_side(context, request, blocks, unit, forward_right,
                          budget, result);
}

// Settles the query gaps between units in forward-read order. Each side
// extends against its own reference, so junctions across contigs or strands
// work too. The left record extends first and the right gets the remainder,
// so records never overlap in query.
void settle_family_junction_seams(const DnaContext& context,
                                  const DnaFamilyRealizationRequest& request,
                                  const std::vector<BlockPlan>& blocks,
                                  std::vector<Unit>& units,
                                  DnaFamilyRealizationOutcome& result) {
  const int read_length = request.family->read_length;
  for (std::size_t index = 1; index < units.size(); ++index) {
    Unit& left = units[index - 1];
    Unit& right = units[index];
    const int gap = unit_forward_begin(right, read_length) -
                    unit_forward_end(left, read_length);
    if (gap <= 0)
      continue;
    const int left_settled = settle_seam_side(
        context, request, blocks, left, /*forward_right=*/true, gap, result);
    settle_seam_side(context, request, blocks, right,
                     /*forward_right=*/false, gap - left_settled, result);
  }
  // Refresh the forward interval from the oriented one.
  for (Unit& unit : units) {
    unit.forward_begin = unit_forward_begin(unit, read_length);
    unit.forward_end = unit_forward_end(unit, read_length);
  }
}

bool validate_and_commit_unit(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              Unit& unit, AlignResult& output) {
  if (unit.core_cigar.empty() || unit.oriented_begin < 0 ||
      unit.oriented_end <= unit.oriented_begin || unit.target_begin < 0 ||
      unit.target_end <= unit.target_begin)
    return false;
  const int forward_begin =
      unit.reverse ? request.family->read_length - unit.oriented_end
                   : unit.oriented_begin;
  const int forward_end =
      unit.reverse ? request.family->read_length - unit.oriented_begin
                   : unit.oriented_end;
  // The one text render of the record's CIGAR.
  const std::string clipped =
      ::fa::cpu::output::clipped_reference_oriented_cigar(
          ::fa::cpu::output::packed_cigar_text(unit.core_cigar),
          request.family->read_length, forward_begin, forward_end,
          unit.reverse);
  if (clipped.empty())
    return false;

  output = {};
  output.read_len = request.family->read_length;
  output.chromosome =
      (*context.ref.names)[static_cast<std::size_t>(unit.chromosome)];
  output.pos = unit.target_begin;
  output.is_reverse = unit.reverse;
  output.origin = AlignmentOrigin::DnaCigarFamily;
  output.score = unit.score;
  output.cigar = clipped;
  const auto& query =
      unit.reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(unit.chromosome)];
  return commit_fixed_cigar_geometry(output, query, reference,
                                     context.opts.cigar_replay_request) &&
         output.query_start == forward_begin &&
         output.query_end == forward_end &&
         output.target_end == unit.target_end;
}

// minimap2 extends the primary toward the read ends until a Z-drop stops it
// and only then assigns roles on the aligned spans (mm_set_parent), so a copy
// the extended primary overlaps becomes a secondary. Here each block's
// extension stops at its query edge, so this step does the same afterwards,
// for the primary family only:
//  1. each side of the primary that borders another block record is extended
//     toward the read end with extend_unit_side, unless it already Z-dropped
//     or is a split edge (minimap2's cap at nearby seeds is not applied);
//  2. the extended unit is revalidated as one record; on failure the family
//     is left as it was;
//  3. every other block record the extended primary covers by at least half
//     of its own span is moved to result.demoted, to be emitted as a
//     secondary with MAPQ 0. The rest stay supplementaries, overlap allowed.
// Map-only output has no DP, so there such a read keeps its split records.
// `primary_unit` is the unit that became the primary record.
void extend_primary_and_demote(const DnaContext& context,
                               const DnaFamilyRealizationRequest& request,
                               const std::vector<BlockPlan>& segments,
                               const std::vector<Unit>& units,
                               std::size_t primary_unit,
                               DnaRecordFamily& assembled,
                               DnaFamilyRealizationOutcome& result) {
  if (primary_unit >= units.size())
    return;
  const int read_length = request.family->read_length;
  Unit extended = units[primary_unit];
  int added = 0;
  for (const bool forward_right : {false, true}) {
    // Units are in forward-read order.
    const bool borders =
        forward_right ? primary_unit + 1 < units.size() : primary_unit > 0;
    if (!borders)
      continue;
    const bool oriented_right =
        extended.reverse ? !forward_right : forward_right;
    const BlockPlan& facing =
        segments[forward_right ? extended.blocks.back()
                               : extended.blocks.front()];
    const bool split_edge =
        oriented_right ? facing.truncated_right : facing.split_continuation;
    const bool zdropped = oriented_right ? extended.zdrop_oriented_right
                                         : extended.zdrop_oriented_left;
    if (split_edge || zdropped)
      continue;
    const int budget =
        forward_right ? read_length - unit_forward_end(extended, read_length)
                      : unit_forward_begin(extended, read_length);
    added += extend_unit_side(context, request, segments, extended,
                              forward_right, budget, result);
  }
  if (added == 0)
    return;
  AlignResult record;
  if (!validate_and_commit_unit(context, request, extended, record))
    return;
  const int covered_begin = record.query_start;
  const int covered_end = record.query_end;
  std::vector<DnaSegmentRecord> kept;
  kept.reserve(assembled.supplementary.size());
  for (std::size_t position = 0; position < assembled.supplementary.size();
       ++position) {
    DnaSegmentRecord& other = assembled.supplementary[position];
    const int span = other.alignment.query_end - other.alignment.query_start;
    const int overlap = std::min(covered_end, other.alignment.query_end) -
                        std::max(covered_begin, other.alignment.query_start);
    if (other.candidate != ::fa::cpu::voting::kNullCandidate && span > 0 &&
        overlap > 0 &&
        2 * static_cast<std::int64_t>(overlap) >=
            static_cast<std::int64_t>(span)) {
      result.demoted.push_back(std::move(other.alignment));
      result.demoted_candidates.push_back(other.candidate);
      result.demoted_positions.push_back(static_cast<int>(position));
    } else {
      kept.push_back(std::move(other));
    }
  }
  assembled.supplementary = std::move(kept);
  result.pre_extension_primary = assembled.primary;
  result.primary_extended = true;
  assembled.primary = std::move(record);
}

// The inversion middles join the supplementaries, which are then sorted by
// (query_start, query_end, chromosome, pos). Recomputes where each demoted
// record would sit in that sort had it stayed, since MAPQ is scored on the
// family before demotion.
void place_demoted_among_middles(const std::vector<DnaSegmentRecord>& kept,
                                 const std::vector<AlignResult>& middles,
                                 DnaFamilyRealizationOutcome& result) {
  struct Entry {
    const AlignResult* record;
    int demoted; // index into result.demoted, or -1
  };
  std::vector<Entry> entries;
  entries.reserve(kept.size() + result.demoted.size() + middles.size());
  std::size_t next_kept = 0;
  std::size_t next_demoted = 0;
  const std::size_t blocks = kept.size() + result.demoted.size();
  for (std::size_t position = 0; position < blocks; ++position) {
    if (next_demoted < result.demoted.size() &&
        result.demoted_positions[next_demoted] ==
            static_cast<int>(position)) {
      entries.push_back(
          {&result.demoted[next_demoted], static_cast<int>(next_demoted)});
      ++next_demoted;
    } else if (next_kept < kept.size()) {
      entries.push_back({&kept[next_kept].alignment, -1});
      ++next_kept;
    }
  }
  for (const AlignResult& middle : middles)
    entries.push_back({&middle, -1});
  std::stable_sort(entries.begin(), entries.end(),
                   [](const Entry& left, const Entry& right) {
                     return std::tie(left.record->query_start,
                                     left.record->query_end,
                                     left.record->chromosome,
                                     left.record->pos) <
                            std::tie(right.record->query_start,
                                     right.record->query_end,
                                     right.record->chromosome,
                                     right.record->pos);
                   });
  for (std::size_t index = 0; index < entries.size(); ++index)
    if (entries[index].demoted >= 0)
      result.demoted_positions[static_cast<std::size_t>(
          entries[index].demoted)] = static_cast<int>(index);
}

// minimap2's mm_align1_inv. `incumbent` is minimap2's r1 (truncated at the
// Z-drop cell) and `continuation` its r2 (split off with split_inv):
//
//   ql = r1->rev? r1->qs - r2->qe : r2->qs - r1->qe   (forward-read coords)
//   tl = r2->rs - r1->re
//   qseq = r1->rev? &qseq0[0][r2->qe] : &qseq0[1][qlen - r2->qs]
//
// The query gap in the orientation opposite the halves is aligned against the
// reference gap: both slices are reversed for one ksw_ll local score, whose
// end cell gives the start offsets, and the middle is then an extension from
// there at band bw * 1.5 with no end bonus. minimap2's split and parent flags
// are replaced by the equivalent split-group adjacency, and the ksw_ll
// endpoints are also range-checked.
bool realize_inversion_middle(const DnaContext& context,
                              const DnaFamilyRealizationRequest& request,
                              const Unit& incumbent, const Unit& continuation,
                              DnaFamilyRealizationOutcome& result,
                              AlignResult& output) {
  const int read_length = request.family->read_length;
  const bool rev = incumbent.reverse;
  const int incumbent_begin = unit_forward_begin(incumbent, read_length);
  const int incumbent_end = unit_forward_end(incumbent, read_length);
  const int continuation_begin = unit_forward_begin(continuation, read_length);
  const int continuation_end = unit_forward_end(continuation, read_length);
  const int ql = rev ? incumbent_begin - continuation_end
                     : continuation_begin - incumbent_end;
  const int tl = continuation.target_begin - incumbent.target_end;

  const DpScoringParams scoring = fill_dp_scoring(context.opts);
  // As mm_align1_inv: min_chain_score <= ql, tl <= max_gap.
  if (ql < scoring.inversion_min_chain_score || ql > scoring.inversion_max_gap)
    return false;
  if (tl < scoring.inversion_min_chain_score || tl > scoring.inversion_max_gap)
    return false;

  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(incumbent.chromosome)];
  const int target_origin = incumbent.target_end;
  if (target_origin < 0 ||
      target_origin + tl > static_cast<int>(reference.size()))
    return false;
  // The middle is on the opposite strand, so its CIGAR is written against the
  // other oriented read.
  const bool middle_reverse = !rev;
  const auto& oriented_query =
      middle_reverse ? *request.reverse_query : *request.forward_query;
  const int query_origin =
      rev ? continuation_end : read_length - continuation_begin;
  if (query_origin < 0 || query_origin + ql > read_length)
    return false;

  ++result.inversion_middles;

  // Reverse both slices, run ksw_ll, and turn the end cell into start
  // offsets in the unreversed slices.
  std::vector<std::uint8_t> query_reversed(static_cast<std::size_t>(ql));
  for (int index = 0; index < ql; ++index) {
    query_reversed[static_cast<std::size_t>(index)] =
        oriented_query[static_cast<std::size_t>(query_origin + ql - 1 - index)];
  }
  std::vector<std::uint8_t> target_reversed(static_cast<std::size_t>(tl));
  for (int index = 0; index < tl; ++index) {
    target_reversed[static_cast<std::size_t>(index)] =
        reference[static_cast<std::size_t>(target_origin + tl - 1 - index)];
  }
  std::array<std::int8_t, 25> matrix{};
  ::fa::cpu::ksw2_simple_mat(matrix.data(), scoring.match, scoring.mismatch,
                             scoring.ambi);
  const ::fa::cpu::dp::DpLocalScoreResult probe = ::fa::cpu::dp::dp_local_score(
      query_reversed.data(), ql, target_reversed.data(), tl, matrix.data(),
      scoring.gap_open1, scoring.gap_extend1);
  if (!probe.attempted || probe.query_end < 0 || probe.query_end >= ql ||
      probe.target_end < 0 || probe.target_end >= tl ||
      probe.score < scoring.inversion_min_dp_max)
    return false;
  const int q_off = ql - (probe.query_end + 1);
  const int t_off = tl - (probe.target_end + 1);
  const int query_span = ql - q_off;
  const int target_span = tl - t_off;
  if (query_span <= 0 || target_span <= 0)
    return false;

  const int oriented_begin = query_origin + q_off;
  const int target_begin = target_origin + t_off;
  realization::RealizationOutcome outcome = run_dna_long_realization(
      realization::RealizationRole::DnaLocalInversionMiddle, scoring,
      realization::make_query_slice(oriented_query.data() + oriented_begin,
                                    query_span, oriented_begin, read_length,
                                    middle_reverse),
      realization::make_target_slice(reference.data() + target_begin,
                                     target_span, target_begin),
      static_cast<int>(scoring.bw * 1.5));
  price_fill_path(context.opts, oriented_query.data() + oriented_begin,
                  query_span, reference.data() + target_begin, target_span,
                  outcome);
  record_kernel_work(outcome, query_span, target_span, result);
  // As mm_align1_inv, any traceback is accepted; a Z-dropped extension still
  // ends at its maximum cell.
  if (outcome.kind != realization::OutcomeKind::Aligned &&
      outcome.kind != realization::OutcomeKind::Zdropped)
    return false;
  const int query_used = ::fa::cpu::output::packed_cigar_query_consumed(
      outcome.raw.packed_cigar);
  const int target_used = outcome.alignment.ref_consumed;
  if (outcome.raw.packed_cigar.empty() || outcome.alignment.ref_offset != 0 ||
      query_used <= 0 || query_used > query_span || target_used <= 0 ||
      target_used > target_span || outcome.raw.max_query + 1 != query_used ||
      outcome.raw.max_target + 1 != target_used)
    return false;

  Unit middle;
  middle.candidate = incumbent.candidate;
  middle.chromosome = incumbent.chromosome;
  middle.reverse = middle_reverse;
  middle.oriented_begin = oriented_begin;
  middle.oriented_end = oriented_begin + query_used;
  middle.target_begin = target_begin;
  middle.target_end = target_begin + target_used;
  // The max-cell score, as minimap2's dp_score = ez->max.
  middle.score = outcome.alignment.score;
  middle.core_cigar = outcome.raw.packed_cigar;
  middle.forward_begin = unit_forward_begin(middle, read_length);
  middle.forward_end = unit_forward_end(middle, read_length);
  if (!validate_and_commit_unit(context, request, middle, output))
    return false;
  // Marks the record for the inversion MAPQ rule (min of the flanks), tp:A:I
  // and the supplementary flag.
  output.origin = AlignmentOrigin::DnaLocalInversion;
  ++result.inversion_records;
  return true;
}

// Collects the inversion middles of a realized family. An inversion split
// always lies between two adjacent units (split seams are never bridged); the
// incumbent is the left unit on the forward strand and the right one on the
// reverse strand.
std::vector<AlignResult>
collect_inversion_middles(const DnaContext& context,
                          const DnaFamilyRealizationRequest& request,
                          const std::vector<BlockPlan>& segments,
                          const std::vector<Unit>& units,
                          DnaFamilyRealizationOutcome& result) {
  std::vector<AlignResult> middles;
  for (std::size_t index = 0; index + 1 < units.size(); ++index) {
    const Unit& left = units[index];
    const Unit& right = units[index + 1];
    if (left.chromosome != right.chromosome || left.reverse != right.reverse ||
        left.blocks.empty() || right.blocks.empty())
      continue;
    const BlockPlan& left_segment = segments[left.blocks.back()];
    const BlockPlan& right_segment = segments[right.blocks.front()];
    if (left_segment.split_group != right_segment.split_group)
      continue;
    const BlockPlan& continuation_segment =
        left.reverse ? left_segment : right_segment;
    const BlockPlan& incumbent_segment =
        left.reverse ? right_segment : left_segment;
    if (!continuation_segment.split_continuation ||
        !continuation_segment.split_inv || !incumbent_segment.truncated_right)
      continue;
    AlignResult inverted;
    if (realize_inversion_middle(context, request,
                                 left.reverse ? right : left,
                                 left.reverse ? left : right, result, inverted))
      middles.push_back(std::move(inverted));
  }
  return middles;
}

} // namespace

DnaFamilyJoinPlan
plan_dna_family_join(const DnaFamilyJoinGeometry& geometry) noexcept {
  DnaFamilyJoinPlan plan;
  if (geometry.left_reference < 0 ||
      geometry.left_reference != geometry.right_reference ||
      geometry.left_reverse != geometry.right_reverse ||
      geometry.left_forward_begin < 0 ||
      geometry.left_forward_end <= geometry.left_forward_begin ||
      geometry.right_forward_end <= geometry.right_forward_begin ||
      geometry.left_forward_end > geometry.right_forward_begin ||
      geometry.max_query_gap < 0 || geometry.max_reference_gap < 0)
    return plan;
  if (!geometry.left_reverse) {
    plan.query_gap =
        geometry.right_query_begin_cut - geometry.left_query_end_cut;
    plan.reference_gap =
        geometry.right_target_begin_cut - geometry.left_target_end_cut;
  } else {
    plan.query_gap =
        geometry.left_query_begin_cut - geometry.right_query_end_cut;
    plan.reference_gap =
        geometry.left_target_begin_cut - geometry.right_target_end_cut;
  }
  if (plan.query_gap < 0 || plan.reference_gap < 0 ||
      plan.query_gap > geometry.max_query_gap ||
      plan.reference_gap > geometry.max_reference_gap) {
    plan.query_gap = -1;
    plan.reference_gap = -1;
    return plan;
  }
  if (plan.query_gap > 0 && plan.reference_gap > 0) {
    if (geometry.matrix_cell_cap <= 0 ||
        static_cast<std::int64_t>(plan.query_gap) >
            geometry.matrix_cell_cap /
                static_cast<std::int64_t>(plan.reference_gap)) {
      plan.query_gap = -1;
      plan.reference_gap = -1;
      return plan;
    }
    plan.kind = DnaFamilyJoinKind::TwoAxis;
    return plan;
  }
  if (plan.query_gap > 0)
    plan.kind = DnaFamilyJoinKind::PureInsertion;
  else if (plan.reference_gap > 0)
    plan.kind = DnaFamilyJoinKind::PureDeletion;
  else
    plan.kind = DnaFamilyJoinKind::ZeroSeam;
  return plan;
}

bool dna_family_seam_has_duplicate_anchor(
    const ordered::OrderedAnchorPath& left,
    const ordered::OrderedAnchorPath& right, bool reverse) noexcept {
  if (left.selected.empty() || right.selected.empty() ||
      left.realization_begin >= left.realization_end ||
      right.realization_begin >= right.realization_end ||
      left.realization_end > left.selected.size() ||
      right.realization_end > right.selected.size())
    return false;
  const std::size_t left_index =
      reverse ? left.realization_begin : left.realization_end - 1;
  const std::size_t right_index =
      reverse ? right.realization_end - 1 : right.realization_begin;
  return exact_anchor_equal(left.selected[left_index],
                            right.selected[right_index]);
}

std::vector<chaining::Anchor>
dna_selected_sibling_path(const DnaPlacementCandidateChain& evidence,
                          int oriented_begin, int oriented_end,
                          int* sibling) {
  return selected_sibling_path(evidence, oriented_begin, oriented_end,
                               sibling);
}

const char* dna_family_failure_name(DnaFamilyFailure failure) noexcept {
  switch (failure) {
  case DnaFamilyFailure::None:
    return "none";
  case DnaFamilyFailure::InvalidInput:
    return "invalid_input";
  case DnaFamilyFailure::PlacementChainRefused:
    return "placement_chain_refused";
  case DnaFamilyFailure::NoSelectedFamily:
    return "no_selected_family";
  case DnaFamilyFailure::NoAnchors:
    return "no_anchors";
  case DnaFamilyFailure::InvalidOrderedPath:
    return "invalid_ordered_path";
  case DnaFamilyFailure::InteriorRefused:
    return "interior_refused";
  case DnaFamilyFailure::TerminalRefused:
    return "terminal_refused";
  case DnaFamilyFailure::FamilyValidation:
    return "family_validation";
  }
  return "invalid_failure";
}

DnaFamilyRealizationOutcome
realize_full_cigar_family(const DnaContext& context,
                          const DnaFamilyRealizationRequest& request) {
  DnaFamilyRealizationOutcome result;
  if (request.placement == nullptr || !request.placement->accepted) {
    result.failure = DnaFamilyFailure::PlacementChainRefused;
    return result;
  }
  if (request.family == nullptr || request.forward_query == nullptr ||
      request.reverse_query == nullptr || context.ref.index == nullptr ||
      context.ref.names == nullptr || context.ref.encoded == nullptr ||
      !request.family->valid || request.family->read_length <= 0 ||
      static_cast<int>(request.forward_query->size()) !=
          request.family->read_length ||
      static_cast<int>(request.reverse_query->size()) !=
          request.family->read_length) {
    return result;
  }

  std::vector<BlockPlan> blocks;
  if (!materialize_plans(context, request, blocks)) {
    result.failure = DnaFamilyFailure::InvalidInput;
    return result;
  }
  result.block_count = static_cast<int>(blocks.size());
  if (blocks.empty()) {
    result.failure = DnaFamilyFailure::NoSelectedFamily;
    return result;
  }
  // A block realizes into one or more segments: a Z-dropped interior fill
  // truncates the segment and the rest continues in a new one (minimap2's
  // mm_split_reg).
  std::vector<BlockPlan> segments;
  segments.reserve(blocks.size());
  for (std::size_t index = 0; index < blocks.size(); ++index) {
    const std::size_t group_begin = segments.size();
    BlockPlan current = std::move(blocks[index]);
    current.split_group = static_cast<int>(index);
    // The bad-end trim falls on the read's two ends, never on a seam. The
    // flags name ends of the oriented path, so a reverse block swaps them.
    const bool read_start = index == 0;
    const bool read_end = index + 1 == blocks.size();
    const bool trim_left = current.reverse ? read_end : read_start;
    const bool trim_right = current.reverse ? read_start : read_end;
    for (;;) {
      DnaFamilyFailure failure = DnaFamilyFailure::None;
      BlockSplit split;
      // A continuation's left end is the Z-drop junction; never trim it.
      if (!acquire_and_realize_block(context, request, current,
                                     trim_left && !current.split_continuation,
                                     trim_right, failure, result, split)) {
        result.failure = failure;
        result.failed_block = static_cast<int>(index);
        return result;
      }
      const bool truncated = split.truncated;
      const bool continues = split.has_continuation;
      // A segment cut at its start realized nothing and is dropped, as
      // minimap2's mm_filter_regs drops such a region.
      if (!(truncated && current.query_end_cut == current.query_begin_cut))
        segments.push_back(std::move(current));
      if (truncated)
        ++result.zdrop_truncations;
      if (!continues)
        break;
      ++result.zdrop_splits;
      current = std::move(split.continuation);
    }
    // Keep segments in forward-read order; reverse blocks produce them
    // backwards.
    if (segments.size() > group_begin && segments[group_begin].reverse)
      std::reverse(segments.begin() + static_cast<std::ptrdiff_t>(group_begin),
                   segments.end());
  }
  result.segment_count = static_cast<int>(segments.size());
  // Every block was cut at its start with too few anchors left.
  if (segments.empty()) {
    result.failure = DnaFamilyFailure::NoAnchors;
    return result;
  }

  std::vector<Unit> units;
  units.push_back(unit_from_block(segments.front(), 0));
  for (std::size_t index = 1; index < segments.size(); ++index) {
    const BlockPlan& previous = segments[index - 1];
    const BlockPlan& current = segments[index];
    ::fa::cpu::output::PackedCigar bridge;
    int bridge_score = 0;
    // Two segments of one block are the sides of a Z-drop and stay separate
    // records.
    const bool split_seam = previous.split_group == current.split_group;
    if (!split_seam && attempt_bridge(context, request, previous, current,
                                      bridge, bridge_score, result)) {
      accept_join(units.back(), current, index, bridge, bridge_score);
      // With kDnaChainMapqStudyBridgeOwner set, the joined unit takes the
      // candidate with the higher chain score; otherwise the leftmost
      // block's candidate keeps it.
      if ((dna_chain_mapq_study_bits(context.opts.chain_mapq_hifi_margin) &
           kDnaChainMapqStudyBridgeOwner) != 0 &&
          request.placement != nullptr) {
        Unit& joined = units.back();
        const DnaPlacementCandidateChain* incoming =
            request.placement->find(current.candidate);
        const DnaPlacementCandidateChain* owner =
            request.placement->find(joined.candidate);
        const int incoming_score =
            incoming != nullptr ? incoming->chain_score : 0;
        const int owner_score = owner != nullptr ? owner->chain_score : 0;
        if (incoming_score > owner_score)
          joined.candidate = current.candidate;
      }
    } else {
      units.push_back(unit_from_block(current, index));
    }
  }
  result.unit_count = static_cast<int>(units.size());

  std::vector<DnaSegmentRecord> records;
  records.reserve(units.size());
  std::vector<AlignResult> pre_settle_alignments;
  pre_settle_alignments.reserve(units.size());
  const auto commit_unit = [&](Unit& unit) -> bool {
    AlignResult record;
    if (!validate_and_commit_unit(context, request, unit, record))
      return false;
    pre_settle_alignments.push_back(record);
    records.push_back({std::move(record), unit.candidate});
    return true;
  };
  // The pre-settle records decide disjointness, roles and the primary's
  // decision score.
  for (Unit& unit : units) {
    if (!extend_unit_terminals(context, request, segments, unit, result)) {
      result.failure = DnaFamilyFailure::TerminalRefused;
      return result;
    }
    if (!commit_unit(unit)) {
      result.failure = DnaFamilyFailure::FamilyValidation;
      return result;
    }
  }

  // Seam settlement. Every changed unit is revalidated first; one failure
  // keeps the whole pre-settle family.
  settle_family_junction_seams(context, request, segments, units, result);
  std::vector<AlignResult> settled_alignments;
  settled_alignments.reserve(units.size());
  bool settled_valid = true;
  for (std::size_t index = 0; index < units.size(); ++index) {
    Unit& unit = units[index];
    // An unchanged unit would commit to the same record.
    if (!unit.settled_dirty) {
      settled_alignments.push_back(pre_settle_alignments[index]);
      continue;
    }
    AlignResult settled;
    if (!validate_and_commit_unit(context, request, unit, settled)) {
      settled_valid = false;
      break;
    }
    settled_alignments.push_back(std::move(settled));
  }
  // As minimap2, inversion middles are realized once both halves are
  // aligned. Settling never moves a split edge, so the settled units give the
  // same seam geometry as the committed records.
  const std::vector<AlignResult> inversion_middles =
      collect_inversion_middles(context, request, segments, units, result);
  DnaRecordFamily assembled = assemble_dna_record_family(std::move(records));
  if (!assembled.valid) {
    result.failure = DnaFamilyFailure::FamilyValidation;
    return result;
  }
  result.decision_score = assembled.primary.score;

  // Substitute the settled records into the role slots. Each slot must match
  // exactly one pre-settle record and each record must be used once;
  // otherwise the pre-settle family is kept whole. `primary_unit` is the unit
  // that became the primary, or units.size() when nothing was substituted.
  std::size_t primary_unit = units.size();
  if (settled_valid &&
      settled_alignments.size() == pre_settle_alignments.size()) {
    const auto same_slot = [](const AlignResult& left,
                              const AlignResult& right) {
      return std::tie(left.query_start, left.query_end, left.chromosome,
                      left.pos, left.target_end, left.is_reverse, left.score) ==
             std::tie(right.query_start, right.query_end, right.chromosome,
                      right.pos, right.target_end, right.is_reverse,
                      right.score);
    };
    std::vector<AlignResult*> role_slots;
    role_slots.reserve(1 + assembled.supplementary.size());
    role_slots.push_back(&assembled.primary);
    for (DnaSegmentRecord& supplementary : assembled.supplementary)
      role_slots.push_back(&supplementary.alignment);

    std::vector<std::size_t> substitution;
    substitution.reserve(role_slots.size());
    std::vector<bool> used(pre_settle_alignments.size(), false);
    bool bijective = role_slots.size() == pre_settle_alignments.size();
    for (const AlignResult* slot : role_slots) {
      std::size_t match = pre_settle_alignments.size();
      int matches = 0;
      for (std::size_t index = 0; index < pre_settle_alignments.size();
           ++index) {
        if (same_slot(*slot, pre_settle_alignments[index])) {
          match = index;
          ++matches;
        }
      }
      if (matches != 1) {
        bijective = false;
        break;
      }
      if (used[match]) {
        bijective = false;
        break;
      }
      used[match] = true;
      substitution.push_back(match);
    }
    if (bijective && std::all_of(used.begin(), used.end(),
                                 [](bool value) { return value; })) {
      for (std::size_t slot = 0; slot < role_slots.size(); ++slot)
        *role_slots[slot] = std::move(settled_alignments[substitution[slot]]);
      primary_unit = substitution.front();
    }
  }
  // Before the inversion middles join, so they are never demoted.
  if (request.primary_family && units.size() > 1)
    extend_primary_and_demote(context, request, segments, units, primary_unit,
                              assembled, result);
  // An inversion middle is always a supplementary. It lies inside the query
  // gap between its halves, so the family stays query-disjoint. It owns no
  // block and takes its MAPQ from its flanks.
  if (!inversion_middles.empty()) {
    if (!result.demoted.empty())
      place_demoted_among_middles(assembled.supplementary, inversion_middles,
                                  result);
    for (const AlignResult& middle : inversion_middles)
      assembled.supplementary.push_back(
          {middle, ::fa::cpu::voting::kNullCandidate});
    std::sort(assembled.supplementary.begin(), assembled.supplementary.end(),
              [](const DnaSegmentRecord& left, const DnaSegmentRecord& right) {
                return std::tie(left.alignment.query_start,
                                left.alignment.query_end,
                                left.alignment.chromosome,
                                left.alignment.pos) <
                       std::tie(right.alignment.query_start,
                                right.alignment.query_end,
                                right.alignment.chromosome,
                                right.alignment.pos);
              });
  }
  static_cast<AlignResult&>(result.output) = std::move(assembled.primary);
  result.primary_candidate = assembled.primary_candidate;
  result.output.supplementary.clear();
  result.output.supplementary_candidates.clear();
  for (DnaSegmentRecord& supplementary : assembled.supplementary) {
    result.output.supplementary.push_back(std::move(supplementary.alignment));
    result.output.supplementary_candidates.push_back(supplementary.candidate);
  }
  result.failure = DnaFamilyFailure::None;
  return result;
}

} // namespace fa::cpu::lr
