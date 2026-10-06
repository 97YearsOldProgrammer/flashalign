#include "family_realization.h"
#include "family_realization_internal.h"

#include "cigar_geometry.h"
#include "dp_runner.h"
#include "../chaining/dense_chain.h"
#include "inv_local_chain.h"
#include "ordered_anchor_path.h"
#include "placement_chaining.h"
#include "../core/cigar.h"
#include "../dp/control.h"     // dp_local_score (ksw_ll) for the inversion probe
#include "../dp/ksw2_align.h"  // ksw2_simple_mat
#include "../dp/params.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

namespace ordered = ::fa::cpu::lr::ordered_anchor;
namespace realization = ::fa::cpu::lr::realization;

using family_internal::BlockPlan;
using family_internal::BlockSplit;
using family_internal::Unit;

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
// regions with the gap certificate that plans them, and the inversion probe
// and middle a fill's Z-drop test starts) run under the fill row, -A -B -O -E
// -z --score-N. The read-end extensions run under the preset's end row
// (cigar_dp_*), which also prices every path. The inversion gates of a fill read -S scaled to the fill
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

bool exact_anchor_equal(const chaining::Anchor& left,
                        const chaining::Anchor& right) {
  return left.q == right.q && left.r == right.r && left.span == right.span;
}

// Splits a block at a Z-dropped seam or piece, as minimap2's mm_align1: the
// returned traceback, ending at the maximum-scoring cell, is kept as the
// segment's prefix and the anchors past that cell form the continuation. A
// fill dropped at its start (over the matrix cap) cuts at the step start.
// Returns false only on a malformed truncation.
bool plan_block_split(const DnaContext& context, BlockPlan& plan,
                      const ordered::GeometryStep& step,
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

  // Split after the last anchor ending at or before the cut, as minimap2;
  // the continuation also skips anchors reaching back over the cut.
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

  BlockPlan& continuation = split.continuation;
  continuation.chromosome = plan.chromosome;
  continuation.reverse = plan.reverse;
  continuation.split_group = plan.split_group;
  continuation.split_continuation = true;
  continuation.split_inv = outcome.trace.zdrop_test.code == 2;
  // Every anchor past the cut, the ones the bad-end trim left out included,
  // as minimap2's mm_split_reg.
  continuation.continuation_anchors.assign(
      plan.path.selected.begin() + static_cast<std::ptrdiff_t>(start),
      plan.path.selected.end());
  split.has_continuation = true;
  return true;
}

int unit_forward_begin(const Unit& unit, int read_length) {
  return unit.reverse ? read_length - unit.oriented_end : unit.oriented_begin;
}

int unit_forward_end(const Unit& unit, int read_length) {
  return unit.reverse ? read_length - unit.oriented_begin : unit.oriented_end;
}

} // namespace

namespace family_internal {

void deduplicate_exact_anchors(std::vector<chaining::Anchor>& anchors) {
  anchors.erase(std::unique(anchors.begin(), anchors.end(), exact_anchor_equal),
                anchors.end());
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
  // local chain of the read's opposite-lane seeds in its drop window.
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
        if (plan_block_split(context, plan, step, outcome, q_cursor, r_cursor,
                             split))
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

Unit unit_from_block(const BlockPlan& block) {
  Unit unit;
  unit.chromosome = block.chromosome;
  unit.reverse = block.reverse;
  unit.oriented_begin = block.query_begin_cut;
  unit.oriented_end = block.query_end_cut;
  unit.target_begin = block.target_begin_cut;
  unit.target_end = block.target_end_cut;
  unit.score = block.interior_score;
  unit.core_cigar = block.interior_cigar;
  return unit;
}

bool extend_unit_terminals(const DnaContext& context,
                           const DnaFamilyRealizationRequest& request,
                           const BlockPlan& block, Unit& unit,
                           DnaFamilyRealizationOutcome& result,
                           const ordered::TerminalWindowControl& caps) {
  const auto& query =
      unit.reverse ? *request.reverse_query : *request.forward_query;
  const auto& reference =
      (*context.ref.encoded)[static_cast<std::size_t>(unit.chromosome)];
  const DpScoringParams scoring = dp_scoring_from_opts(context.opts);
  ordered::TerminalWindowControl terminal_control{
      std::max(1, context.opts.cigar_dp_max_gap), scoring.match,
      scoring.gap_open1, scoring.gap_extend1};
  terminal_control.left_cap = caps.left_cap;
  terminal_control.left_cap_query = caps.left_cap_query;
  terminal_control.left_cap_target = caps.left_cap_target;
  terminal_control.right_cap = caps.right_cap;
  terminal_control.right_cap_query = caps.right_cap_query;
  terminal_control.right_cap_target = caps.right_cap_target;
  const ordered::TerminalWindowPlan terminals =
      ordered::plan_terminal_windows(block.path, terminal_control);
  if (!terminals)
    return false;

  // As minimap2 (r->split_inv ? zdrop_inv : zdrop): the left terminal of a
  // continuation is its Z-drop seam, and after an inversion split it extends
  // at the tighter inversion Z-drop so it does not run through the inversion.
  DpScoringParams seam_scoring = scoring;
  const bool seam_is_inversion =
      block.split_continuation && block.split_inv;
  if (seam_is_inversion)
    seam_scoring.zdrop = scoring.inversion_zdrop;

  if (terminals.has_left_packet()) {
    const int q_span =
        terminals.first_adjusted.query - terminals.query_start;
    const int r_span =
        terminals.first_adjusted.target - terminals.target_start;
    const realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaLeftExtension, seam_scoring,
        realization::make_query_slice(query.data() + terminals.query_start,
                                      q_span, terminals.query_start,
                                      request.family->read_length,
                                      unit.reverse),
        realization::make_target_slice(reference.data() +
                                           terminals.target_start,
                                       r_span, terminals.target_start),
        context.opts.cigar_dp_bw);
    record_kernel_work(outcome, q_span, r_span, result);
    const int q_used = ::fa::cpu::output::packed_cigar_query_consumed(
        outcome.raw.packed_cigar);
    if ((outcome.kind == realization::OutcomeKind::Aligned ||
         outcome.kind == realization::OutcomeKind::Zdropped) &&
        !outcome.raw.packed_cigar.empty() && q_used > 0 && q_used <= q_span &&
        outcome.alignment.ref_offset == 0 &&
        outcome.alignment.ref_consumed >= 0 &&
        outcome.alignment.ref_consumed <= r_span) {
      ::fa::cpu::output::prepend_packed_cigar(unit.core_cigar,
                                              outcome.raw.packed_cigar);
      unit.oriented_begin = terminals.first_adjusted.query - q_used;
      unit.target_begin =
          terminals.first_adjusted.target - outcome.alignment.ref_consumed;
      unit.score += outcome.alignment.score;
    }
  }

  // As minimap2, a Z-dropped block gets no right extension: it ends at the
  // truncation cell.
  if (terminals.has_right_packet() && !block.truncated_right) {
    const int q_span =
        terminals.query_end - terminals.last_adjusted.query;
    const int r_span =
        terminals.target_end - terminals.last_adjusted.target;
    const realization::RealizationOutcome outcome = run_dna_long_realization(
        realization::RealizationRole::DnaRightExtension, scoring,
        realization::make_query_slice(
            query.data() + terminals.last_adjusted.query, q_span,
            terminals.last_adjusted.query, request.family->read_length,
            unit.reverse),
        realization::make_target_slice(
            reference.data() + terminals.last_adjusted.target, r_span,
            terminals.last_adjusted.target),
        context.opts.cigar_dp_bw);
    record_kernel_work(outcome, q_span, r_span, result);
    const int q_used = ::fa::cpu::output::packed_cigar_query_consumed(
        outcome.raw.packed_cigar);
    if ((outcome.kind == realization::OutcomeKind::Aligned ||
         outcome.kind == realization::OutcomeKind::Zdropped) &&
        !outcome.raw.packed_cigar.empty() && q_used > 0 && q_used <= q_span &&
        outcome.alignment.ref_offset == 0 &&
        outcome.alignment.ref_consumed >= 0 &&
        outcome.alignment.ref_consumed <= r_span) {
      ::fa::cpu::output::append_packed_cigar(unit.core_cigar,
                                             outcome.raw.packed_cigar);
      unit.oriented_end = terminals.last_adjusted.query + q_used;
      unit.target_end =
          terminals.last_adjusted.target + outcome.alignment.ref_consumed;
      unit.score += outcome.alignment.score;
    }
  }
  return true;
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
// there at band bw * 1.5 with no end bonus. The ksw_ll endpoints are also
// range-checked.
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
  middle.chromosome = incumbent.chromosome;
  middle.reverse = middle_reverse;
  middle.oriented_begin = oriented_begin;
  middle.oriented_end = oriented_begin + query_used;
  middle.target_begin = target_begin;
  middle.target_end = target_begin + target_used;
  // The max-cell score, as minimap2's dp_score = ez->max.
  middle.score = outcome.alignment.score;
  middle.core_cigar = outcome.raw.packed_cigar;
  if (!validate_and_commit_unit(context, request, middle, output))
    return false;
  // Marks the record for the inversion MAPQ rule (min of the flanks), tp:A:I
  // and the supplementary flag.
  output.origin = AlignmentOrigin::DnaLocalInversion;
  ++result.inversion_records;
  return true;
}

} // namespace family_internal

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

const char* dna_family_failure_name(DnaFamilyFailure failure) noexcept {
  switch (failure) {
  case DnaFamilyFailure::None:
    return "none";
  case DnaFamilyFailure::InvalidInput:
    return "invalid_input";
  case DnaFamilyFailure::PlacementChainRefused:
    return "placement_chain_refused";
  case DnaFamilyFailure::NoOwnerChain:
    return "no_owner_chain";
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

} // namespace fa::cpu::lr
