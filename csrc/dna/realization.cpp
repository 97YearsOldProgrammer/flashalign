// CIGAR realization of the family: regions from the final ownership, trimmed, split on Z-drop, their gaps filled between anchors, seams, the inversion probe, the committed geometry.
#include "dna.h"

#include "../chain/chain.h"
#include "../core/core.h"
#include "../dp/dp.h" // dp_local_score (ksw_ll) for the inversion probe; ksw2_simple_mat
#include "../index/format.h"
#include "../seed/seed.h"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/**** geometry_work ****/

namespace fa::cpu::lr {

void DnaGeometryWork::add(const DnaGeometryWork& other) {
  verified_regions += other.verified_regions;
  verified_bases += other.verified_bases;
  seams += other.seams;
  stretches += other.stretches;
  stretch_pieces += other.stretch_pieces;
  stretch_bases += other.stretch_bases;
  failed_regions += other.failed_regions;
  failed_region_anchors += other.failed_region_anchors;
  failed_region_bases += other.failed_region_bases;
  failing_gaps += other.failing_gaps;
  join_calls += other.join_calls;
}

} // namespace fa::cpu::lr

/**** cigar_geometry ****/

namespace fa::cpu::lr {

// Replays a fixed CIGAR and fills the record's coordinates and accounting;
// false if the CIGAR does not replay. `request` asks for cs:Z / MD:Z / =X
// output (minimap2 --cs / --MD / --eqx), built during the same replay; with
// --eqx the CIGAR is replaced by its =/X spelling.
static bool commit_fixed_cigar_geometry(
    AlignResult& result, const std::vector<std::uint8_t>& query,
    const std::vector<std::uint8_t>& reference,
    const ::fa::cpu::output::CigarReplayRequest& request = {}) {
  if (result.cigar.empty() || result.read_len < 0 ||
      static_cast<int>(query.size()) < result.read_len || result.pos < 0) {
    return false;
  }
  const ::fa::cpu::output::CigarReplay replay = ::fa::cpu::output::replay_cigar(
      result.cigar, query.data(), result.read_len, reference.data(),
      static_cast<int>(reference.size()), result.pos, request);
  const ::fa::cpu::output::ForwardQueryInterval query_interval =
      ::fa::cpu::output::forward_query_interval_from_cigar(
          result.cigar, result.read_len, result.is_reverse);
  if (!replay.valid || replay.query_consumed != result.read_len ||
      replay.target_start != result.pos || !query_interval.valid) {
    return false;
  }

  result.query_start = query_interval.q_lo;
  result.query_end = query_interval.q_hi;
  result.target_end = replay.target_end;
  result.matches = replay.matches;
  result.mismatches = replay.mismatches;
  result.insertions = replay.insertions;
  result.deletions = replay.deletions;
  result.ambiguities = replay.ambiguities;
  result.edit_distance = replay.edit_distance;
  result.block_len = replay.block_len;
  result.cs = replay.cs;
  result.md = replay.md;
  if (request.eqx)
    result.cigar = replay.eqx_cigar;
  result.alignment_accounting_valid = true;
  return true;
}

}  // namespace fa::cpu::lr

/**** family_realization_internal ****/
// The realization steps of one anchor path (realization.cpp), which the
// -c lane's controller (realization.cpp) runs on each region.

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
  // A segment of a secondary-role region.
  bool rival = false;
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

/**** inv_local_chain ****/

namespace fa::cpu::lr {

namespace {

// Anchors in the longest chain of the opposite lane's seeds lying wholly inside
// the fill-local window [query_begin, query_end) x [target_begin, target_end).
// An anchor is (q', r): q' the seed's start in its lane, where the block's
// oriented window [a, b) is [L - b, L - a), and r its posting's start on the
// contig. Chained anchors have q' and r strictly increasing. `source` is a
// DnaInvLocalChainSource; the signature is dp/dp.h DpInversionProbeGate.
int dna_inv_local_chain(const void* opaque, int query_begin, int query_end,
                        int target_begin, int target_end) {
  const DnaInvLocalChainSource& source =
      *static_cast<const DnaInvLocalChainSource*>(opaque);
  const int k = source.seed_length;
  // The window in the opposite lane's frame, the reverse complement of the
  // block's.
  const int lane_begin =
      source.read_length - (source.query_offset + query_end);
  const int lane_end =
      source.read_length - (source.query_offset + query_begin);
  const int reference_begin = source.target_offset + target_begin;
  const int reference_end = source.target_offset + target_end;
  if (lane_end - lane_begin < k || reference_end - reference_begin < k)
    return 0;
  const bool lane_reverse = !source.block_reverse;
  const std::vector<RetainedSeedRef>& lane =
      lane_reverse ? source.seeds->fine_reverse()
                   : source.seeds->fine_forward();
  // Both lanes are in ascending read_pos.
  auto seed = std::lower_bound(
      lane.begin(), lane.end(), lane_begin,
      [](const RetainedSeedRef& ref, int position) {
        return ref.seed.read_pos < position;
      });
  std::vector<std::pair<int, int>> anchors;  // (q', r)
  for (; seed != lane.end() && seed->seed.read_pos <= lane_end - k; ++seed) {
    // Postings starting in [reference_begin, reference_end - k].
    const KmerPostingIntervalView interval = source.seeds->slice(
        *source.index, seed->entry, source.chromosome,
        static_cast<std::uint32_t>(reference_begin),
        static_cast<std::uint32_t>(reference_end - k + 1));
    if (!interval.found() || interval.count == 0 ||
        (source.occurrence_cap != 0 &&
         interval.global_count > source.occurrence_cap))
      continue;
    for (std::uint32_t posting = 0; posting < interval.count; ++posting) {
      const PackedRefPos packed = interval.positions.packed_at(posting);
      if (!packed_ref_orientation_compatible(seed->seed.z, packed,
                                             lane_reverse))
        continue;
      anchors.emplace_back(seed->seed.read_pos,
                           static_cast<int>(packed_ref_local(packed)));
    }
  }
  if (anchors.empty()) return 0;
  std::sort(anchors.begin(), anchors.end());
  // Longest chain, exact over the predecessors within kDnaInvLocalGap in q'.
  std::vector<int> length(anchors.size(), 1);
  int best = 1;
  std::size_t low = 0;
  for (std::size_t j = 0; j < anchors.size(); ++j) {
    const int qj = anchors[j].first;
    const int rj = anchors[j].second;
    while (anchors[low].first < qj - kDnaInvLocalGap) ++low;
    const int dj = rj - qj;
    for (std::size_t i = low; i < j; ++i) {
      const int qi = anchors[i].first;
      const int ri = anchors[i].second;
      if (qi >= qj || ri >= rj || rj - ri > kDnaInvLocalGap ||
          std::abs((ri - qi) - dj) > kDnaInvLocalDiagonal)
        continue;
      length[j] = std::max(length[j], length[i] + 1);
    }
    best = std::max(best, length[j]);
  }
  return best;
}

}  // namespace

}  // namespace fa::cpu::lr

/**** family_realization ****/

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
// row (options/options.cpp).
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
  // A rival region ends at the cut.
  if (plan.rival && context.opts.rival_split_ends)
    return true;

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
  // Every anchor past the cut, the ones the bad-end trim left out included:
  // the window after the cut, as mm_split_reg's (as + n, cnt - n).
  continuation.path.selected = {plan.path.selected.data + start,
                                plan.path.selected.size() - start};
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

/**** region_realization ****/

namespace fa::cpu::lr {
namespace {

namespace ordered = ::fa::cpu::lr::ordered_anchor;
using family_internal::BlockPlan;
using family_internal::BlockSplit;
using family_internal::Unit;

// minimap2's mask_level, mask_len, rank_min_len and rank_frac (options.c).
// Its best_n (-N) is context.opts.secondary_max.
constexpr float kMaskLevel = 0.5f;
constexpr int kMaskLen = INT_MAX;
constexpr int kRankMinLen = 500;
constexpr float kRankFrac = 0.9f;

// One chain of the selection: its anchors in chain order, deduplicated, and
// its forward-query and reference spans (mm_reg_set_coor).
struct Chain {
  const DnaSelectionItem* item = nullptr;
  std::vector<chaining::Anchor> anchors;
  int q_begin = 0;
  int q_end = 0;
  int ref_begin = 0;
  int ref_end = 0;
};

// A region mm_select_sub kept: an owner with its chain-level subsc and n_sub,
// or an attached chain.
struct Region {
  int item = 0;
  bool owner = false;
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct Record {
  AlignResult alignment;
  int region = 0;
  int contig = -1;
  bool inversion = false;
  int score = 0;
  int cnt = 0;
  int dp_max = 0;
  int dp_max2 = 0;
  // The record dp_max2 came from, -1 for none.
  int dp_max2_from = -1;
  // The sweep's dp_max before the rank rescale (minimap2's dp_max0), its own
  // and that of the record dp_max2 came from.
  int dp_max0 = 0;
  int dp_max2_0 = 0;
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
  int parent = -1;
  // Its rank in an exact dp_max tie (step 4).
  std::uint64_t tie = 0;
  // Its entry of the outcome's rivals and the entry dp_max2 came from; -1 for
  // none.
  int rival = -1;
  int dp_max2_rival = -1;
};

bool identical_chains(const Chain& a, const Chain& b) {
  return a.q_begin == b.q_begin && a.q_end == b.q_end &&
         a.item->contig == b.item->contig && a.ref_begin == b.ref_begin &&
         a.ref_end == b.ref_end;
}

// The read's seeds the whole-query harvest admits (key in the index, not over
// the pool gate), minimap2's mini_pos, as sorted forward-read starts, both
// lanes.
std::vector<int> admitted_seed_starts(const DnaContext& context,
                                      const RetainedSeedDensity& seeds,
                                      int read_length, int seed_length) {
  const internal::PoolGate gate = internal::pool_gate(context);
  std::vector<int> starts;
  starts.reserve(seeds.fine_forward().size() + seeds.fine_reverse().size());
  for (const RetainedSeedRef& seed : seeds.fine_forward())
    if (!gate.drops(seeds.occurrence(seed.entry), seed.rescued))
      starts.push_back(seed.seed.read_pos);
  for (const RetainedSeedRef& seed : seeds.fine_reverse())
    if (!gate.drops(seeds.occurrence(seed.entry), seed.rescued))
      starts.push_back(read_length - seed.seed.read_pos - seed_length);
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  return starts;
}

// mm_est_err (esterr.c 30-63) on one chain: from its first anchor in
// forward-read order, the walk over the read's admitted seeds matching its
// anchors in order, one seed more at each end with room on both axes, and
// 1 - (matched / total)^(1 / k); -1 when the first anchor is not a seed.
// Every seed spans k, so minimap2's avg_k is k. The right-end test reads
// q_begin, as minimap2's reads qs.
float chain_divergence(const Chain& chain, const std::vector<int>& seeds,
                       int read_length, int seed_length, int contig_length,
                       bool reverse) {
  const int cnt = static_cast<int>(chain.anchors.size());
  const auto forward_start = [&](int k) {
    const chaining::Anchor& anchor =
        reverse ? chain.anchors[static_cast<std::size_t>(cnt - 1 - k)]
                : chain.anchors[static_cast<std::size_t>(k)];
    return reverse ? read_length - anchor.q - anchor.span : anchor.q;
  };
  if (cnt == 0 || seeds.empty()) return -1.0f;
  const auto first =
      std::lower_bound(seeds.begin(), seeds.end(), forward_start(0));
  if (first == seeds.end() || *first != forward_start(0)) return -1.0f;
  const int n = static_cast<int>(seeds.size());
  const int st = static_cast<int>(first - seeds.begin());
  int en = st;
  int n_match = 1;
  for (int k = 1, j = st + 1; j < n && k < cnt; ++j) {
    if (forward_start(k) == seeds[static_cast<std::size_t>(j)]) {
      ++k;
      en = j;
      ++n_match;
    }
  }
  int n_tot = en - st + 1;
  const float avg_k = static_cast<float>(seed_length);
  if (chain.q_begin > avg_k && chain.ref_begin > avg_k) ++n_tot;
  if (read_length - chain.q_begin > avg_k &&
      contig_length - chain.ref_end > avg_k)
    ++n_tot;
  return n_match >= n_tot
             ? 0.0f
             : static_cast<float>(
                   1.0 - std::pow(static_cast<double>(n_match) / n_tot,
                                  1.0 / avg_k));
}

// A chain the strand rule kept: its and its owner's divergence and
// mm_filter_strand_retained's verdict.
struct StrandVerdict {
  int item = 0;
  int owner_item = 0;
  float div = 0.0f;
  float owner_div = 0.0f;
  bool kept = true;
};

// hit.c 255-281 with check_strand: every owner, and an attached chain within
// pri_ratio or 2k of its owner, not identical to it, or on the other strand
// above 0.8 * max_gap, at most best_n per read. The roles are in the parent
// walk's order, minimap2's chain order. With `density`, a chain the strand
// rule keeps then meets mm_filter_strand_retained (hit.c 283-299) on the read's
// admitted seeds, its best_n slot spent either way (map.c 333-336).
std::vector<Region> retain_regions(const DnaContext& context,
                                   const DnaKeptSelection& selection,
                                   const std::vector<Chain>& chains,
                                   int seed_length, int read_length,
                                   const RetainedSeedDensity* density) {
  std::vector<Region> regions;
  const float pri_ratio = static_cast<float>(context.opts.pri_ratio);
  const int min_diff = 2 * seed_length;
  const int min_strand_sc =
      static_cast<int>(context.opts.cigar_dp_max_gap * 0.8);
  int n_2nd = 0;
  const auto valid = [&](const DnaOwnershipRole& role) {
    return role.item >= 0 &&
           static_cast<std::size_t>(role.item) < chains.size() &&
           role.owner_item >= 0 &&
           static_cast<std::size_t>(role.owner_item) < chains.size() &&
           !chains[static_cast<std::size_t>(role.item)].anchors.empty();
  };
  // The read's admitted seeds, built at the first chain the filter tests.
  std::vector<int> seeds;
  bool seeds_built = false;
  for (const DnaOwnershipRole& role : selection.roles) {
    if (!valid(role)) continue;
    const Chain& chain = chains[static_cast<std::size_t>(role.item)];
    if (role.owner_item == role.item) {
      regions.push_back({role.item, true, role.subsc, role.n_sub,
                         role.pool_subsc, role.pool_n_sub});
      continue;
    }
    const std::size_t owner = static_cast<std::size_t>(role.owner_item);
    const Chain& parent = chains[owner];
    const int score = chain.item->score;
    const int parent_score = parent.item->score;
    if ((static_cast<float>(score) >=
             static_cast<float>(parent_score) * pri_ratio ||
         score + min_diff >= parent_score) &&
        n_2nd < context.opts.secondary_max) {
      if (!identical_chains(chain, parent)) {
        regions.push_back({role.item, false});
        ++n_2nd;
      }
    } else if (n_2nd < context.opts.secondary_max && score > min_strand_sc &&
               chain.item->reverse != parent.item->reverse) {
      ++n_2nd;
      if (density != nullptr) {
        if (!seeds_built) {
          seeds = admitted_seed_starts(context, *density, read_length,
                                       seed_length);
          seeds_built = true;
        }
        StrandVerdict verdict;
        verdict.item = role.item;
        verdict.owner_item = role.owner_item;
        verdict.div = chain_divergence(
            chain, seeds, read_length, seed_length,
            static_cast<int>(context.ref.contig_length(chain.item->contig)),
            chain.item->reverse);
        verdict.owner_div = chain_divergence(
            parent, seeds, read_length, seed_length,
            static_cast<int>(context.ref.contig_length(parent.item->contig)),
            parent.item->reverse);
        verdict.kept = verdict.div < verdict.owner_div * 5.0f ||
                       verdict.div < 0.01f;
        if (!verdict.kept) continue;
      }
      regions.push_back({role.item, false});
    }
  }
  return regions;
}

// The kept selection's chains, parallel to its items: each item's anchors
// deduplicated, with its forward-query and reference spans; an item with no
// anchors or a contig at or past `contigs` stays empty. Then retain_regions
// over them.
std::vector<Region> retain_selection(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    int read_length, int seed_length, std::size_t contigs,
    std::vector<Chain>& chains) {
  const DnaKeptSelection& selection = placement.kept_selection;
  chains.assign(selection.items.size(), Chain{});
  for (std::size_t i = 0; i < selection.items.size(); ++i) {
    Chain& chain = chains[i];
    chain.item = &selection.items[i];
    const std::vector<chaining::Anchor>* anchors =
        placement.selection_anchors(selection.items[i]);
    if (anchors == nullptr || anchors->empty() || chain.item->contig < 0 ||
        static_cast<std::size_t>(chain.item->contig) >= contigs)
      continue;
    chain.anchors = *anchors;
    family_internal::deduplicate_exact_anchors(chain.anchors);
    const chaining::Anchor& first = chain.anchors.front();
    const chaining::Anchor& last = chain.anchors.back();
    chain.q_begin = chain.item->reverse ? read_length - last.q - last.span
                                        : first.q;
    chain.q_end = chain.item->reverse ? read_length - first.q
                                      : last.q + last.span;
    chain.ref_begin = first.r;
    chain.ref_end = last.r + last.span;
  }
  return retain_regions(context, selection, chains, seed_length, read_length,
                        placement.strand_filter_seeds);
}

// minimap2's cap at nearby seeds (align.c 706-768). minimap2 lays the
// regions' anchors out per contig and strand in order of their first anchor's
// reference position (lchain.c compact_a); a segment holding [first, end) of
// its region's anchors looks back from its first anchor and forward from its
// last, and the (min_cnt + 1)-th anchor lying beyond it on both axes caps the
// window at its distance.
class NearbySeeds {
 public:
  NearbySeeds(const std::vector<Region>& regions,
              const std::vector<Chain>& chains)
      : regions_(regions), chains_(chains), order_(regions.size()),
        position_(regions.size()) {
    for (std::size_t r = 0; r < regions.size(); ++r) order_[r] = r;
    std::sort(order_.begin(), order_.end(),
              [this](std::size_t left, std::size_t right) {
                return key(left) < key(right);
              });
    for (std::size_t p = 0; p < order_.size(); ++p) position_[order_[p]] = p;
  }

  ordered::TerminalWindowControl caps(std::size_t region, std::size_t first,
                                      int min_cnt) const {
    ordered::TerminalWindowControl out;
    const std::vector<chaining::Anchor>& own = anchors(region);
    const chaining::Anchor& head = own[first];
    const chaining::Anchor& tail = own.back();
    const int rs0 = head.r;
    const int qs0 = head.q;
    const int re0 = tail.r + tail.span;
    const int qe0 = tail.q + tail.span;
    int count = 0;
    const auto before = [&](const chaining::Anchor& anchor) {
      if (anchor.r >= rs0 || anchor.q >= qs0 || ++count <= min_cnt)
        return false;
      const int l = std::max(rs0 - anchor.r, qs0 - anchor.q);
      out.left_cap = true;
      out.left_cap_query = qs0 - l;
      out.left_cap_target = std::max(0, rs0 - l);
      return true;
    };
    bool found = false;
    for (std::size_t i = first; !found && i-- > 0;) found = before(own[i]);
    for (std::size_t p = position_[region]; !found && p-- > 0;) {
      if (!same_group(order_[p], region)) break;
      const std::vector<chaining::Anchor>& other = anchors(order_[p]);
      for (std::size_t i = other.size(); !found && i-- > 0;)
        found = before(other[i]);
    }
    count = 0;
    const auto after = [&](const chaining::Anchor& anchor) {
      const int x = anchor.r + anchor.span;
      const int y = anchor.q + anchor.span;
      if (x <= re0 || y <= qe0 || ++count <= min_cnt) return false;
      const int l = std::max(x - re0, y - qe0);
      out.right_cap = true;
      out.right_cap_query = qe0 + l;
      out.right_cap_target = re0 + l;
      return true;
    };
    found = false;
    for (std::size_t p = position_[region] + 1;
         !found && p < order_.size(); ++p) {
      if (!same_group(order_[p], region)) break;
      for (const chaining::Anchor& anchor : anchors(order_[p]))
        if ((found = after(anchor))) break;
    }
    return out;
  }

 private:
  const std::vector<chaining::Anchor>& anchors(std::size_t region) const {
    return chains_[static_cast<std::size_t>(regions_[region].item)].anchors;
  }
  std::tuple<int, bool, int, std::size_t> key(std::size_t region) const {
    const DnaSelectionItem& item =
        *chains_[static_cast<std::size_t>(regions_[region].item)].item;
    return {item.contig, item.reverse, anchors(region).front().r, region};
  }
  bool same_group(std::size_t a, std::size_t b) const {
    const DnaSelectionItem& x =
        *chains_[static_cast<std::size_t>(regions_[a].item)].item;
    const DnaSelectionItem& y =
        *chains_[static_cast<std::size_t>(regions_[b].item)].item;
    return x.contig == y.contig && x.reverse == y.reverse;
  }

  const std::vector<Region>& regions_;
  const std::vector<Chain>& chains_;
  std::vector<std::size_t> order_;
  std::vector<std::size_t> position_;
};

// mm_filter_regs (hit.c 301-320); max_clip_ratio is 1, which never fires.
// Returns the records dropped.
int filter_records(const DnaContext& context, std::vector<Record>& records) {
  const int min_cnt = context.opts.cigar_dp_split_min_anchors;
  const std::size_t before = records.size();
  records.erase(
      std::remove_if(records.begin(), records.end(),
                     [&](const Record& record) {
                       return (!record.inversion && record.cnt < min_cnt) ||
                              record.alignment.matches <
                                  context.opts.min_chain_score ||
                              record.dp_max < context.opts.cigar_dp_min_dp_max;
                     }),
      records.end());
  return static_cast<int>(before - records.size());
}

// mm_update_dp_max (align.c 1022-1046).
void update_dp_max(const DnaContext& context,
                   const DnaFamilyRealizationRequest& request, int read_length,
                   std::vector<Record>& records) {
  if (records.size() < 2) return;
  int max = -1;
  int max2 = -1;
  int max_i = -1;
  for (std::size_t i = 0; i < records.size(); ++i) {
    const int dp = records[i].dp_max;
    if (dp > max) {
      max2 = max;
      max = dp;
      max_i = static_cast<int>(i);
    } else if (dp > max2) {
      max2 = dp;
    }
  }
  if (max_i < 0 || max < 0 || max2 < 0) return;
  const AlignResult& best = records[static_cast<std::size_t>(max_i)].alignment;
  if (best.query_end - best.query_start <
      static_cast<double>(read_length) * static_cast<double>(kRankFrac))
    return;
  if (max2 < static_cast<double>(max) * static_cast<double>(kRankFrac)) return;
  double identity = -1.0;
  dna_record_event_identity(context, best, *request.forward_query,
                            *request.reverse_query, identity);
  const double b2 = dna_rank_b2(identity, context.opts.cigar_dp_match,
                                context.opts.cigar_dp_mismatch);
  for (Record& record : records)
    record.dp_max = std::max(
        0, dna_record_recal_dp_max(context, record.alignment,
                                   *request.forward_query,
                                   *request.reverse_query, b2));
}

// mm_set_parent (hit.c 125-186) on the aligned forward-read spans: parent,
// subsc, n_sub and dp_max2. A record's subsc carries over from its chain, and
// the first record keeps its chain's n_sub, as in minimap2.
void set_parents(const DnaContext& context,
                 const std::vector<Chain>& chains,
                 const std::vector<Region>& regions,
                 std::vector<Record>& records) {
  const int sub_diff =
      2 * context.opts.cigar_dp_match + context.opts.cigar_dp_mismatch;
  const auto pool_of = [&](const Record& record) {
    const DnaSelectionItem& item =
        *chains[static_cast<std::size_t>(
                    regions[static_cast<std::size_t>(record.region)].item)]
             .item;
    return !record.inversion ? item.candidate : -1;
  };
  if (records.empty()) return;
  std::vector<std::size_t> w{0};
  std::vector<std::pair<int, int>> cover;
  records[0].parent = 0;
  for (std::size_t i = 1; i < records.size(); ++i) {
    Record& ri = records[i];
    const int si = ri.alignment.query_start;
    const int ei = ri.alignment.query_end;
    int uncov_len = 0;
    cover.clear();
    for (const std::size_t j : w) {
      int sj = records[j].alignment.query_start;
      int ej = records[j].alignment.query_end;
      if (ej <= si || sj >= ei) continue;
      if (sj < si) sj = si;
      if (ej > ei) ej = ei;
      cover.emplace_back(sj, ej);
    }
    bool secondary = false;
    if (!cover.empty()) {
      std::sort(cover.begin(), cover.end());
      int x = si;
      for (const auto& [begin, end] : cover) {
        if (begin > x) uncov_len += begin - x;
        x = end > x ? end : x;
      }
      if (ei > x) uncov_len += ei - x;
      for (const std::size_t j : w) {
        Record& rp = records[j];
        const int sj = rp.alignment.query_start;
        const int ej = rp.alignment.query_end;
        if (ej <= si || sj >= ei) continue;
        const int min = std::min(ej - sj, ei - si);
        const int max = std::max(ej - sj, ei - si);
        const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                               : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
        if (static_cast<float>(ol) / static_cast<float>(min) -
                    static_cast<float>(uncov_len) / static_cast<float>(max) <=
                kMaskLevel ||
            uncov_len > kMaskLen)
          continue;
        ri.parent = rp.parent;
        const double sci = ri.score;
        rp.subsc = std::max(rp.subsc, sci);
        bool cnt_sub = ri.cnt >= rp.cnt;
        if (rp.contig != ri.contig || rp.alignment.pos != ri.alignment.pos ||
            rp.alignment.target_end != ri.alignment.target_end || ol != min) {
          if (ri.dp_max > rp.dp_max2) {
            rp.dp_max2 = ri.dp_max;
            rp.dp_max2_from = static_cast<int>(i);
            rp.dp_max2_0 = ri.dp_max0;
          }
          if (rp.dp_max - ri.dp_max <= sub_diff) cnt_sub = true;
        }
        if (cnt_sub) ++rp.n_sub;
        const int pool = pool_of(rp);
        if (pool >= 0 && pool == pool_of(ri)) {
          rp.pool_subsc = std::max(rp.pool_subsc, sci);
          if (cnt_sub) ++rp.pool_n_sub;
        }
        secondary = true;
        break;
      }
    }
    if (!secondary) {
      w.push_back(i);
      ri.parent = static_cast<int>(i);
      ri.n_sub = 0;
      ri.pool_n_sub = 0;
    }
  }
}

// mm_select_sub (hit.c 255-281) without check_strand, on chain scores, then
// mm_sync_regs.
void select_records(const DnaContext& context, int seed_length,
                    std::vector<Record>& records) {
  const float pri_ratio = static_cast<float>(context.opts.pri_ratio);
  const int min_diff = 2 * seed_length;
  std::vector<int> index(records.size(), -1);
  std::vector<bool> keep(records.size(), false);
  int n_2nd = 0;
  for (std::size_t i = 0; i < records.size(); ++i) {
    const Record& ri = records[i];
    const Record& rp = records[static_cast<std::size_t>(ri.parent)];
    if (ri.parent == static_cast<int>(i) || ri.inversion) {
      keep[i] = true;
    } else if ((static_cast<float>(ri.score) >=
                    static_cast<float>(rp.score) * pri_ratio ||
                ri.score + min_diff >= rp.score) &&
               n_2nd < context.opts.secondary_max) {
      const AlignResult& a = ri.alignment;
      const AlignResult& b = rp.alignment;
      if (!(a.query_start == b.query_start && a.query_end == b.query_end &&
            ri.contig == rp.contig && a.pos == b.pos &&
            a.target_end == b.target_end)) {
        keep[i] = true;
        ++n_2nd;
      }
    }
  }
  std::vector<Record> kept;
  for (std::size_t i = 0; i < records.size(); ++i) {
    if (!keep[i]) continue;
    index[i] = static_cast<int>(kept.size());
    kept.push_back(std::move(records[i]));
  }
  for (Record& record : kept) {
    record.parent = index[static_cast<std::size_t>(record.parent)];
    if (record.dp_max2_from >= 0)
      record.dp_max2_from = index[static_cast<std::size_t>(record.dp_max2_from)];
  }
  records = std::move(kept);
}

}  // namespace

std::vector<int> dna_retained_secondary_items(
    const DnaContext& context, const DnaPlacementChainingResult& placement) {
  std::vector<Chain> chains;
  const std::vector<Region> regions = retain_selection(
      context, placement, placement.family.read_length,
      placement.family.seed_length,
      static_cast<std::size_t>(context.ref.contig_count()), chains);
  std::vector<int> items;
  for (const Region& region : regions)
    if (!region.owner) items.push_back(region.item);
  return items;
}

DnaRegionOutcome realize_dna_regions(
    const DnaContext& context, const DnaFamilyRealizationRequest& request) {
  DnaRegionOutcome out;
  DnaFamilyRealizationOutcome& result = out.family;
  if (request.placement == nullptr || !request.placement->accepted) {
    result.failure = request.placement != nullptr &&
                             request.placement->no_owner_chain
                         ? DnaFamilyFailure::NoOwnerChain
                         : DnaFamilyFailure::PlacementChainRefused;
    return out;
  }
  if (request.family == nullptr || request.forward_query == nullptr ||
      request.reverse_query == nullptr || context.ref.index == nullptr ||
      context.ref.names == nullptr || context.ref.encoded == nullptr ||
      !request.family->valid || request.family->read_length <= 0 ||
      static_cast<int>(request.forward_query->size()) !=
          request.family->read_length ||
      static_cast<int>(request.reverse_query->size()) !=
          request.family->read_length)
    return out;
  const DnaPlacementChainingResult& placement = *request.placement;
  const int read_length = request.family->read_length;
  const int seed_length = request.family->seed_length;

  // 1. Retain.
  std::vector<Chain> chains;
  const std::vector<Region> regions =
      retain_selection(context, placement, read_length, seed_length,
                       context.ref.encoded->size(), chains);
  result.block_count = static_cast<int>(regions.size());
  if (regions.empty()) {
    result.failure = DnaFamilyFailure::NoSelectedFamily;
    return out;
  }

  // 2. Realize each region from its whole anchor vector. A refused segment
  // is dropped; the segments before it stand.
  const NearbySeeds nearby(regions, chains);
  const int min_cnt = context.opts.cigar_dp_split_min_anchors;
  const ordered::NormalizationControl normalization{
      std::max(1, context.opts.cigar_dp_bw), context.opts.min_chain_score,
      std::max(1, context.opts.cigar_dp_max_gap)};
  DnaFamilyFailure first_failure = DnaFamilyFailure::None;
  std::vector<Record> records;
  for (std::size_t r = 0; r < regions.size(); ++r) {
    const Region& region = regions[r];
    const Chain& chain = chains[static_cast<std::size_t>(region.item)];
    const DnaSelectionItem& item = *chain.item;
    const std::vector<std::uint8_t>& reference =
        (*context.ref.encoded)[static_cast<std::size_t>(item.contig)];
    const auto refuse = [&](DnaFamilyFailure failure) {
      if (first_failure == DnaFamilyFailure::None) {
        first_failure = failure;
        result.failed_block = static_cast<int>(r);
      }
    };
    BlockPlan current;
    current.chromosome = item.contig;
    current.reverse = item.reverse;
    current.split_group = static_cast<int>(r);
    current.rival = !region.owner;
    // The region's one anchor array, flagged in place as minimap2 flags a[];
    // each segment's path is a window of it.
    std::vector<chaining::Anchor> anchors = chain.anchors;
    current.path.selected = {anchors.data(), anchors.size()};
    int score = item.score;
    int cnt = static_cast<int>(anchors.size());
    bool previous_committed = false;
    Unit previous_unit;
    for (;;) {
      ordered::OrderedAnchorPath& path = current.path;
      path.reference_id = item.contig;
      path.reverse_complemented = item.reverse;
      path.query_length = read_length;
      path.target_length = static_cast<int>(reference.size());
      path.query_bound_begin = 0;
      path.query_bound_end = read_length;
      path.target_bound_begin = 0;
      path.target_bound_end = static_cast<int>(reference.size());
      path.minimizer_k = context.ref.index->k();
      path.realization_begin = 0;
      path.realization_end = path.selected.size();
      if (!ordered::normalize_for_realization(path, normalization)) {
        refuse(DnaFamilyFailure::InvalidOrderedPath);
        break;
      }
      BlockSplit split;
      DnaFamilyFailure failure = DnaFamilyFailure::None;
      if (!family_internal::plan_and_realize_packets(context, request, current,
                                                     failure, result, split)) {
        refuse(failure);
        break;
      }
      ++result.segment_count;
      // mm_split_reg (hit.c 105-123): the rest takes its share of the score.
      int segment_score = score;
      int segment_cnt = cnt;
      int next_score = 0;
      int next_cnt = 0;
      if (split.has_continuation) {
        next_cnt = static_cast<int>(split.continuation.path.selected.size());
        next_score = static_cast<int>(
            static_cast<float>(score) *
                (static_cast<float>(next_cnt) / static_cast<float>(cnt)) +
            .499);
        segment_cnt = cnt - next_cnt;
        segment_score = score - next_score;
      }
      if (split.truncated) ++result.zdrop_truncations;
      bool committed = false;
      Unit unit;
      // A segment cut at its start realized nothing (mm_filter_regs).
      if (!(split.truncated &&
            current.query_end_cut == current.query_begin_cut)) {
        unit = family_internal::unit_from_block(current);
        const ordered::TerminalWindowControl caps = nearby.caps(
            r, chain.anchors.size() - current.path.selected.size(), min_cnt);
        AlignResult alignment;
        if (!family_internal::extend_unit_terminals(context, request, current,
                                                    unit, result, caps)) {
          refuse(DnaFamilyFailure::TerminalRefused);
        } else if (!family_internal::validate_and_commit_unit(
                       context, request, unit, alignment)) {
          refuse(DnaFamilyFailure::FamilyValidation);
        } else {
          // mm_align1_inv between two pieces of a chain-level primary
          // (align.c 1102-1108, 925-926).
          if (region.owner && current.split_continuation &&
              current.split_inv && previous_committed) {
            Record middle;
            if (family_internal::realize_inversion_middle(
                    context, request, previous_unit, unit, result,
                    middle.alignment)) {
              middle.region = static_cast<int>(r);
              middle.contig = item.contig;
              middle.inversion = true;
              records.push_back(std::move(middle));
            }
          }
          Record record;
          record.alignment = std::move(alignment);
          record.alignment.chain_anchors = segment_cnt;
          record.alignment.chain_score = segment_score;
          record.region = static_cast<int>(r);
          record.contig = item.contig;
          record.score = segment_score;
          record.cnt = segment_cnt;
          record.subsc = region.subsc;
          record.n_sub = region.n_sub;
          record.pool_subsc = region.pool_subsc;
          record.pool_n_sub = region.pool_n_sub;
          records.push_back(std::move(record));
          committed = true;
        }
      }
      previous_committed = committed;
      previous_unit = std::move(unit);
      if (!split.has_continuation) break;
      ++result.zdrop_splits;
      current = std::move(split.continuation);
      score = next_score;
      cnt = next_cnt;
    }
  }
  result.unit_count = static_cast<int>(records.size());

  // 3. Filter, and for a read of rank_min_len or more rescale and filter
  // again (align.c 1113-1117).
  for (Record& record : records) {
    record.dp_max = dna_record_dp_max_segment(
        context, record.alignment, *request.forward_query,
        *request.reverse_query);
    record.dp_max0 = record.dp_max;
  }
  filter_records(context, records);
  if (read_length >= kRankMinLen) {
    update_dp_max(context, request, read_length, records);
    filter_records(context, records);
  }
  if (records.empty()) {
    out.floored = result.unit_count > 0;
    result.failure = out.floored ? DnaFamilyFailure::None
                     : first_failure != DnaFamilyFailure::None
                         ? first_failure
                         : DnaFamilyFailure::NoAnchors;
    return out;
  }

  // 4. Re-select on the aligned spans, in AS order, dp_max breaking its ties
  // (mm_hit_sort orders by dp_max alone). In an exact tie the read-seeded hash
  // of the record's contig, strand and reference start decides (seed/seed.h),
  // then contig and start.
  for (Record& record : records)
    record.tie = tie_locus_hash(context.vote_tie_seed, record.contig,
                                record.alignment.is_reverse,
                                record.alignment.pos);
  std::stable_sort(records.begin(), records.end(),
                   [](const Record& left, const Record& right) {
                     if (left.alignment.score != right.alignment.score)
                       return left.alignment.score > right.alignment.score;
                     if (left.dp_max != right.dp_max)
                       return left.dp_max > right.dp_max;
                     return std::tie(left.tie, left.contig,
                                     left.alignment.pos) <
                            std::tie(right.tie, right.contig,
                                     right.alignment.pos);
                   });
  set_parents(context, chains, regions, records);
  // Every secondary is a MAPQ rival, the ones select_records drops included,
  // so a dp_max2 source outlives the drop.
  for (std::size_t i = 0; i < records.size(); ++i) {
    Record& record = records[i];
    if (record.parent == static_cast<int>(i) || record.inversion) continue;
    const Chain& chain = chains[static_cast<std::size_t>(
        regions[static_cast<std::size_t>(record.region)].item)];
    DnaRealizedRival rival;
    rival.candidate = chain.item->candidate;
    rival.path = chain.item->path;
    rival.region = record.region;
    rival.score0 = chain.item->score;
    rival.anchors0 = static_cast<int>(chain.anchors.size());
    rival.chain_q_begin = chain.q_begin;
    rival.chain_q_end = chain.q_end;
    rival.chain_ref_begin = chain.ref_begin;
    rival.chain_ref_end = chain.ref_end;
    record.rival = static_cast<int>(out.rivals.size());
    out.rivals.push_back(rival);
  }
  for (Record& record : records)
    if (record.dp_max2_from >= 0)
      record.dp_max2_rival =
          records[static_cast<std::size_t>(record.dp_max2_from)].rival;
  select_records(context, seed_length, records);

  // The first primary is the SAM primary (mm_set_sam_pri), the other
  // primaries its supplementaries, each secondary a hypothesis of its own.
  dna::Result& output = result.output;
  bool head = true;
  for (std::size_t i = 0; i < records.size(); ++i) {
    Record& record = records[i];
    const Chain& chain = chains[static_cast<std::size_t>(
        regions[static_cast<std::size_t>(record.region)].item)];
    const DnaSelectionItem& item = *chain.item;
    if (record.parent != static_cast<int>(i)) {
      record.alignment.mapq = 0;
      output.secondary.push_back(std::move(record.alignment));
      continue;
    }
    DnaRegionPrice price;
    price.candidate = record.inversion ? -1 : item.candidate;
    price.inversion = record.inversion;
    price.score = record.score;
    price.cnt = record.cnt;
    price.score0 = item.score;
    price.anchors0 = static_cast<int>(chain.anchors.size());
    price.chain_q_begin = chain.q_begin;
    price.chain_q_end = chain.q_end;
    price.chain_ref_begin = chain.ref_begin;
    price.chain_ref_end = chain.ref_end;
    price.dp_max = record.dp_max;
    price.dp_max2 = record.dp_max2;
    price.dp_max0 = record.dp_max0;
    price.dp_max2_0 = record.dp_max2_0;
    if (record.dp_max2_from >= 0) {
      const Record& from = records[static_cast<std::size_t>(record.dp_max2_from)];
      const DnaSelectionItem& source =
          *chains[static_cast<std::size_t>(
                      regions[static_cast<std::size_t>(from.region)].item)]
               .item;
      if (!from.inversion)
        price.dp_max2_candidate = source.candidate;
    }
    price.subsc = record.subsc;
    price.n_sub = record.n_sub;
    price.pool_subsc = record.pool_subsc;
    price.pool_n_sub = record.pool_n_sub;
    price.dp_max2_rival = record.dp_max2_rival;
    const ::fa::cpu::voting::CandidateId candidate =
        price.candidate >= 0
            ? static_cast<::fa::cpu::voting::CandidateId>(price.candidate)
            : ::fa::cpu::voting::kNullCandidate;
    if (head) {
      static_cast<AlignResult&>(output) = std::move(record.alignment);
      result.primary_candidate = candidate;
      head = false;
    } else {
      output.supplementary.push_back(std::move(record.alignment));
      output.supplementary_candidates.push_back(candidate);
    }
    out.prices.push_back(price);
  }
  output.primary_part = -1;
  result.failure = DnaFamilyFailure::None;
  return out;
}

}  // namespace fa::cpu::lr
