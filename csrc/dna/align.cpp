// The DP driver: one packet per gap, an executor per kind (ksw2 or verified base comparison), and minimap2's post-DP terms of a record.
#include "dna.h"

#include "../core/core.h"
#include "../dp/dp.h" // ::fa::cpu::ksw2_simple_mat + KSW flags; DpMapOpt defaults/effective bands

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>

/**** kernel_packets ****/
// The DNA realization geometry between a block's first and last anchors. No
// code here: the planner is plan_verified_geometry() (dna.h),
// the executors are in align.cpp and the gap fill is
// dp_fill_gap().
//
// The anchors are the dense chain's path inside the block after
// mm_fix_bad_ends and filter_bad_seeds, in chain order. An interior anchor
// flagged ANCHOR_IGNORE or ANCHOR_TANDEM is skipped as in minimap2; the first
// and last anchors never are. Every corner is a k-mer center, as minimap2's
// mm_adjust_minier.
//
// A region is a maximal run of consecutive retained anchors on one diagonal;
// a long join starts a new one. The gap between two neighbouring anchors of a
// region is the bases between their k-mers (none when they overlap). A region
// of two or more anchors is verified when every gap passes the ungapped
// certificate
//
//   S0 >= (L - 1) * a - 2 * g(1),   g(l) = min(q + e * l, q2 + e2 * l)
//
// with L the gap length and S0 its ungapped score: a path that leaves the
// diagonal keeps at most L - 1 match columns and pays two one-base gaps. A
// verified region is emitted ungapped, with no DP, and cannot Z-drop. It is
// optimal among paths through every anchor center of the region, not among all
// paths between its ends.
//
// A seam joins two neighbouring verified regions with one gap fill from the
// left region's last center to the right region's first. It is minimap2's
// fill: pinned-global DP at bw_long_eff(), or max(q_span, r_span) across a
// long join, with the approximate pass, the Z-drop test and one exact retry.
//
// A stretch is everything else: a region with any failing gap, a lone anchor,
// the material between verified regions that are not neighbours. minimap2's
// loop cuts it into pieces, each ending at the stretch end, at a long join or
// at the first retained anchor where both spans since the previous corner reach
// min_ksw_len. Each piece is one fill, like a seam.
//
// ANCHOR_TANDEM marks every anchor of a seed whose query neighbour has the
// same key, minimap2's MM_SEED_TANDEM. In a tandem array such a seed may sit on
// the wrong copy and a corner would pin the path there, so flagged interior
// anchors are skipped and one fill crosses the run.
//
// A Z-drop in a fill splits the block at the fill's maximum cell: the prefix
// is kept and the anchors past the cut form a new segment while at least
// cigar_dp_split_min_anchors remain. A fill over the matrix cap
// (DpMapOpt::max_sw_mat) is returned Z-dropped at its start, so the block
// splits there.

/**** realization_controller ****/

namespace fa {
namespace cpu {
namespace lr {
namespace realization {

namespace {

EndpointObjective objective_for_role(RealizationRole role);

RawAlignmentResult normalize_raw(const ::fa::cpu::dp::DpAlnResult& raw) {
  RawAlignmentResult out;
  out.has_cigar = raw.has_cigar;
  out.zdropped = raw.zdropped;
  out.reach_end = raw.reach_end;
  out.max_score = raw.max;
  out.max_query = raw.max_q;
  out.max_target = raw.max_t;
  out.query_end_score = raw.mqe;
  out.query_end_target = raw.mqe_t;
  out.target_end_score = raw.mte;
  out.target_end_query = raw.mte_q;
  out.corner_score = raw.corner_score;
  out.packed_cigar = raw.cigar;
  return out;
}

RawAlignmentResult normalize_raw(::fa::cpu::dp::DpAlnResult&& raw) {
  std::vector<std::uint32_t> packed_cigar = std::move(raw.cigar);
  RawAlignmentResult out =
      normalize_raw(static_cast<const ::fa::cpu::dp::DpAlnResult&>(raw));
  out.packed_cigar = std::move(packed_cigar);
  return out;
}

KernelId normalize_kernel(::fa::cpu::dp::DpKernelId kernel) {
  switch (kernel) {
  case ::fa::cpu::dp::DpKernelId::Extz2:
    return KernelId::Extz2;
  case ::fa::cpu::dp::DpKernelId::Extd2:
    return KernelId::Extd2;
  }
  return KernelId::None;
}

ZdropTestEvidence
normalize_zdrop_test(const ::fa::cpu::dp::DpZdropResult& test) {
  ZdropTestEvidence out;
  out.tested = test.tested;
  out.code = test.code;
  out.max_drop = test.max_drop;
  out.target_begin = test.target_begin;
  out.target_end = test.target_end;
  out.query_begin = test.query_begin;
  out.query_end = test.query_end;
  switch (test.inversion_status) {
  case ::fa::cpu::dp::DpInversionProbeStatus::Disabled:
    out.inversion_status = InversionProbeStatus::Disabled;
    break;
  case ::fa::cpu::dp::DpInversionProbeStatus::DropBelowInversionThreshold:
    out.inversion_status = InversionProbeStatus::DropBelowInversionThreshold;
    break;
  case ::fa::cpu::dp::DpInversionProbeStatus::InvalidBounds:
    out.inversion_status = InversionProbeStatus::InvalidBounds;
    break;
  case ::fa::cpu::dp::DpInversionProbeStatus::WeakReverseScore:
    out.inversion_status = InversionProbeStatus::WeakReverseScore;
    break;
  case ::fa::cpu::dp::DpInversionProbeStatus::Accepted:
    out.inversion_status = InversionProbeStatus::Accepted;
    break;
  case ::fa::cpu::dp::DpInversionProbeStatus::Gated:
    out.inversion_status = InversionProbeStatus::Gated;
    break;
  }
  out.reverse_probe_attempted = test.reverse_probe.attempted;
  out.reverse_score = test.reverse_probe.score;
  out.reverse_query_end = test.reverse_probe.query_end;
  out.reverse_target_end = test.reverse_probe.target_end;
  return out;
}

struct AttemptCollector {
  ExecutorId executor = ExecutorId::None;
  ExecutionTrace* trace = nullptr;
};

void collect_attempt(void* user_data, ::fa::cpu::dp::DpKernelId kernel,
                     int flags, int band, int zdrop, bool approximate,
                     const ::fa::cpu::dp::DpAlnResult& raw) {
  auto* collector = static_cast<AttemptCollector*>(user_data);
  if (collector == nullptr || collector->trace == nullptr)
    return;
  KernelAttempt attempt;
  attempt.executor = collector->executor;
  attempt.kernel = normalize_kernel(kernel);
  attempt.flags = static_cast<std::uint32_t>(flags);
  attempt.band = band;
  attempt.zdrop = zdrop;
  attempt.approximate = approximate;
  attempt.exact = !approximate;
  attempt.raw = normalize_raw(raw);
  collector->trace->push(std::move(attempt));
}

::fa::cpu::DpMapOpt map_opt(const RealizationRequest& request) {
  ::fa::cpu::DpMapOpt opt;
  opt.a = request.scoring.match;
  opt.b = request.scoring.mismatch;
  opt.sc_ambi = request.scoring.ambiguity;
  opt.q = request.scoring.gap1.open;
  opt.e = request.scoring.gap1.extend;
  opt.q2 = request.scoring.gap2.open;
  opt.e2 = request.scoring.gap2.extend;
  opt.bw = request.configured_band;
  opt.bw_long = request.configured_long_band;
  opt.zdrop = request.zdrop;
  opt.zdrop_inv = request.inversion_zdrop;
  opt.end_bonus = request.end_bonus;
  opt.max_sw_mat = request.matrix_cell_cap;
  return opt;
}

std::array<std::int8_t, 25> simple_matrix(const RealizationRequest& request) {
  std::array<std::int8_t, 25> matrix{};
  ::fa::cpu::ksw2_simple_mat(matrix.data(), request.scoring.match,
                             request.scoring.mismatch,
                             request.scoring.ambiguity);
  return matrix;
}

bool role_matches_objective(const RealizationRequest& request) {
  return request.objective == objective_for_role(request.role);
}

bool valid_role(RealizationRole role) {
  switch (role) {
  case RealizationRole::DnaLeftExtension:
  case RealizationRole::DnaInternalFill:
  case RealizationRole::DnaRightExtension:
  case RealizationRole::DnaSupplementaryFill:
  case RealizationRole::DnaLocalInversionMiddle:
  case RealizationRole::RnaExonFill:
    return true;
  }
  return false;
}

bool ksw2_supports(const RealizationRequest& request, SupportReason& reason) {
  if (request.verified) {
    // A verified region is never re-banded into a DP as a fallback.
    reason = SupportReason::UnsupportedPacketClass;
    return false;
  }
  if (!role_matches_objective(request)) {
    reason = SupportReason::UnsupportedObjective;
    return false;
  }
  const bool expected_reversal =
      request.objective == EndpointObjective::OpenEndedLeft;
  if (request.query.reversed_for_executor != expected_reversal ||
      request.target.reversed_for_executor != expected_reversal) {
    reason = SupportReason::UnsupportedOrientation;
    return false;
  }
  if (request.scoring.model == ScoreModel::Splice) {
    reason = SupportReason::UnsupportedScoreModel;
    return false;
  }
  const bool single_affine =
      request.scoring.gap1.open == request.scoring.gap2.open &&
      request.scoring.gap1.extend == request.scoring.gap2.extend;
  if ((request.scoring.model == ScoreModel::SingleAffine) != single_affine) {
    reason = SupportReason::UnsupportedScoreModel;
    return false;
  }
  if (request.scoring.substitution_matrix != simple_matrix(request)) {
    reason = SupportReason::UnsupportedSubstitutionMatrix;
    return false;
  }

  const ::fa::cpu::DpMapOpt opt = map_opt(request);
  // The caller's band must match the one derived here.
  const int expected_band =
      request.role == RealizationRole::DnaLocalInversionMiddle
          ? static_cast<int>(request.configured_band * 1.5)
      : request.objective == EndpointObjective::PinnedGlobal
          ? (request.long_join
                 ? std::max(request.query.length, request.target.length)
                 : opt.bw_long_eff())
          : opt.bw_eff();
  if (request.selected_band != expected_band) {
    reason = SupportReason::UnsupportedBand;
    return false;
  }

  std::uint32_t expected_flags = 0;
  if (request.objective == EndpointObjective::PinnedGlobal) {
    expected_flags = KSW_EZ_APPROX_MAX;
  } else {
    expected_flags = KSW_EZ_EXTZ_ONLY;
    if (request.objective == EndpointObjective::OpenEndedLeft)
      expected_flags |= KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR;
  }
  if (request.ksw_flags != expected_flags) {
    reason = SupportReason::UnsupportedFlags;
    return false;
  }

  reason = SupportReason::Supported;
  return true;
}

ExecutorResult ksw2_execute(const RealizationRequest& request,
                            ExecutionContext& context, ExecutionTrace& trace) {
  (void)context;
  const ::fa::cpu::DpMapOpt opt = map_opt(request);
  AttemptCollector collector{ExecutorId::Ksw2, &trace};
  ::fa::cpu::dp::DpAttemptObserver observer;
  observer.callback = collect_attempt;
  observer.user_data = &collector;

  ::fa::cpu::dp::DpAlnResult raw;
  if (request.objective == EndpointObjective::PinnedGlobal) {
    ::fa::cpu::dp::DpInversionProbeControl inversion_control;
    inversion_control.enabled = request.inversion_probe_enabled;
    inversion_control.max_gap = request.inversion_max_gap;
    inversion_control.min_chain_score = request.inversion_min_chain_score;
    inversion_control.min_dp_max = request.inversion_min_dp_max;
    inversion_control.gate = request.inversion_probe_gate.count;
    inversion_control.gate_source = request.inversion_probe_gate.source;
    inversion_control.gate_min_count = request.inversion_probe_gate.min_count;
    raw = ::fa::cpu::dp::dp_fill_gap(
        opt, request.query.length, request.query.data, request.target.length,
        request.target.data, request.long_join, &observer,
        request.inversion_probe_enabled ? &inversion_control : nullptr);
  } else if (request.role == RealizationRole::DnaLocalInversionMiddle) {
    raw = ::fa::cpu::dp::dp_align_inversion_middle(
        opt, request.query.length, request.query.data, request.target.length,
        request.target.data, &observer);
  } else {
    raw = ::fa::cpu::dp::dp_extend(
        opt, request.query.length, request.query.data, request.target.length,
        request.target.data,
        request.objective == EndpointObjective::OpenEndedLeft, &observer);
  }
  trace.zdrop_test = normalize_zdrop_test(raw.zdrop_test);
  trace.exact_retry = raw.exact_retry;
  trace.exact_retry_zdrop = raw.exact_retry_zdrop;

  ExecutorResult out;
  out.alignment = ::fa::cpu::dp::dp_to_gotoh(
      raw, request.query.data, request.query.length, request.target.data,
      request.target.length,
      request.objective == EndpointObjective::PinnedGlobal,
      /*render_text=*/false);
  out.raw = normalize_raw(std::move(raw));
  return out;
}

const ExecutorEntry* find_executor(const ExecutorTable& table, ExecutorId id) {
  for (std::size_t i = 0; i < table.size; ++i) {
    if (table.entries[i].id == id)
      return &table.entries[i];
  }
  return nullptr;
}

bool valid_orientation(SequenceOrientation orientation) {
  switch (orientation) {
  case SequenceOrientation::Forward:
  case SequenceOrientation::ReverseComplement:
    return true;
  }
  return false;
}

bool valid_slice(const SequenceSlice& slice) {
  if (slice.length < 0 || slice.original_begin < 0 ||
      slice.original_end < slice.original_begin ||
      !valid_orientation(slice.orientation))
    return false;
  if (slice.original_end - slice.original_begin != slice.length)
    return false;
  return slice.length == 0 || slice.data != nullptr;
}

bool valid_gap_scoring(const ScoringContract& scoring) {
  if (scoring.model == ScoreModel::Splice || scoring.gap1.open < 0 ||
      scoring.gap1.extend < 0 || scoring.gap2.open < 0 ||
      scoring.gap2.extend < 0)
    return false;
  const bool single_affine = scoring.gap1.open == scoring.gap2.open &&
                             scoring.gap1.extend == scoring.gap2.extend;
  return (scoring.model == ScoreModel::SingleAffine) == single_affine;
}

bool materialize_pure_axis(const RealizationRequest& request,
                           RealizationOutcome& outcome) {
  const bool insertion = request.query.length > 0 && request.target.length == 0;
  const int length = insertion ? request.query.length : request.target.length;
  if (length <= 0 || static_cast<std::uint64_t>(length) >
                         (std::numeric_limits<std::uint32_t>::max() >> 4))
    return false;
  const std::int64_t penalty1 =
      static_cast<std::int64_t>(request.scoring.gap1.open) +
      static_cast<std::int64_t>(request.scoring.gap1.extend) * length;
  const std::int64_t penalty2 =
      static_cast<std::int64_t>(request.scoring.gap2.open) +
      static_cast<std::int64_t>(request.scoring.gap2.extend) * length;
  const std::int64_t penalty = std::min(penalty1, penalty2);
  if (penalty > std::numeric_limits<int>::max())
    return false;

  const std::uint32_t operation = insertion ? 1u : 2u;
  outcome.alignment.ref_offset = 0;
  outcome.alignment.ref_consumed = insertion ? 0 : length;
  outcome.alignment.matches = 0;
  outcome.alignment.score = -static_cast<int>(penalty);
  outcome.raw.has_cigar = true;
  outcome.raw.corner_score = outcome.alignment.score;
  outcome.raw.packed_cigar.push_back((static_cast<std::uint32_t>(length) << 4) |
                                     operation);
  return true;
}

bool matrix_cap_exceeded(const RealizationRequest& request) {
  if (request.matrix_cell_cap <= 0 || request.query.length <= 0 ||
      request.target.length <= 0)
    return false;
  return static_cast<std::int64_t>(request.query.length) >
         request.matrix_cell_cap /
             static_cast<std::int64_t>(request.target.length);
}

void notify(ExecutionContext& context, const RealizationRequest& request,
            const RealizationOutcome& outcome) {
  if (context.observer.callback != nullptr)
    context.observer.callback(context.observer.user_data, request, outcome);
}

RealizationOutcome finish(ExecutionContext& context,
                          const RealizationRequest& request,
                          RealizationOutcome outcome) {
  notify(context, request, outcome);
  return outcome;
}

bool supports_entry(const ExecutorEntry* entry,
                    const RealizationRequest& request, SupportReason& reason) {
  if (entry == nullptr || entry->supports == nullptr ||
      entry->execute == nullptr) {
    reason = SupportReason::ExecutorUnavailable;
    return false;
  }
  return entry->supports(request, reason);
}

bool packed_consumption(const RawAlignmentResult& raw, int& query_consumed,
                        int& target_consumed) {
  query_consumed = 0;
  target_consumed = 0;
  for (const std::uint32_t packed : raw.packed_cigar) {
    const int length = static_cast<int>(packed >> 4);
    if (length <= 0)
      return false;
    switch (packed & UINT32_C(0xf)) {
    case 0: // M
      query_consumed += length;
      target_consumed += length;
      break;
    case 1: // I
      query_consumed += length;
      break;
    case 2: // D
      target_consumed += length;
      break;
    default:
      return false;
    }
  }
  return true;
}

// The Verified executor. A verified region lies on one diagonal and each of
// its gaps passed the ungapped certificate (GapCertificate), so it is scored
// base by base and emitted as one M run; the commit derives mismatches, cs and
// MD by replay. No DP, no Z-drop. minimap2 writes the certificate's ceiling
// with len - 2; len - 1 is the exact bound, and with map-hifi and map-ont
// scoring both make the same decisions.

int matrix_score(const RealizationRequest& request, std::uint8_t query_base,
                 std::uint8_t target_base) {
  const std::size_t q = query_base > 4 ? 4 : query_base;
  const std::size_t t = target_base > 4 ? 4 : target_base;
  return request.scoring.substitution_matrix[q * 5 + t];
}

void push_cigar_op(std::vector<std::uint32_t>& packed, std::uint32_t op,
                   int length) {
  if (length <= 0)
    return;
  if (!packed.empty() && (packed.back() & UINT32_C(0xf)) == op) {
    packed.back() += static_cast<std::uint32_t>(length) << 4;
    return;
  }
  packed.push_back((static_cast<std::uint32_t>(length) << 4) | op);
}

bool verified_supports(const RealizationRequest& request,
                       SupportReason& reason) {
  if (!request.verified) {
    reason = SupportReason::UnsupportedPacketClass;
    return false;
  }
  if (request.role != RealizationRole::DnaInternalFill ||
      request.objective != EndpointObjective::PinnedGlobal) {
    reason = SupportReason::UnsupportedObjective;
    return false;
  }
  if (request.query.reversed_for_executor ||
      request.target.reversed_for_executor) {
    reason = SupportReason::UnsupportedOrientation;
    return false;
  }
  if (request.scoring.model == ScoreModel::Splice) {
    reason = SupportReason::UnsupportedScoreModel;
    return false;
  }
  // Must be the matrix the certificate scored the gaps with.
  if (request.scoring.substitution_matrix != simple_matrix(request)) {
    reason = SupportReason::UnsupportedSubstitutionMatrix;
    return false;
  }
  if (request.ksw_flags != KSW_EZ_APPROX_MAX) {
    reason = SupportReason::UnsupportedFlags;
    return false;
  }
  // A verified region stays on one diagonal.
  if (request.query.length != request.target.length ||
      request.query.length <= 0 || request.long_join) {
    reason = SupportReason::UnsupportedPacketClass;
    return false;
  }
  if (request.gap_count < 0 ||
      (request.gap_count > 0 && request.gaps == nullptr)) {
    reason = SupportReason::InvalidPacket;
    return false;
  }
  for (int index = 0; index < request.gap_count; ++index) {
    const ::fa::cpu::lr::RegionGap& gap = request.gaps[index];
    if (gap.begin < 0 || gap.end <= gap.begin ||
        gap.end > request.query.length ||
        (index > 0 && gap.begin < request.gaps[index - 1].end)) {
      reason = SupportReason::InvalidPacket;
      return false;
    }
  }
  reason = SupportReason::Supported;
  return true;
}

ExecutorResult verified_execute(const RealizationRequest& request,
                                ExecutionContext& context,
                                ExecutionTrace& trace) {
  (void)context;
  const int length = request.query.length;
  const std::uint8_t* const query = request.query.data;
  const std::uint8_t* const target = request.target.data;

  VerifiedWork work;
  work.bases = length;

  const int match = request.scoring.match;

  std::vector<std::uint32_t> packed;
  std::int64_t score = 0;
  int matches = 0;

  const auto emit_ungapped = [&](int from, int to) {
    for (int i = from; i < to; ++i) {
      score += matrix_score(request, query[i], target[i]);
      if (query[i] == target[i])
        ++matches;
    }
    push_cigar_op(packed, 0u, to - from);
  };

  // Bases under the exact-match anchors: when memcmp confirms they are equal
  // and none is ambiguous, each scores `match`. Anything else takes the
  // per-base loop, which gives the same result.
  const auto emit_certified = [&](int from, int to) {
    const int covered = to - from;
    if (covered <= 0)
      return;
    if (std::memcmp(query + from, target + from,
                    static_cast<std::size_t>(covered)) == 0) {
      std::uint8_t acc = 0;
      for (int i = from; i < to; ++i)
        acc |= query[i];
      if (acc < 4) {
        score += static_cast<std::int64_t>(match) * covered;
        matches += covered;
        push_cigar_op(packed, 0u, covered);
        return;
      }
    }
    emit_ungapped(from, to);
  };

  // Certified gaps are optimal ungapped, so they are scored, not aligned.
  int cursor = 0;
  for (int index = 0; index < request.gap_count; ++index) {
    const ::fa::cpu::lr::RegionGap gap = request.gaps[index];
    emit_certified(cursor, gap.begin);
    emit_ungapped(gap.begin, gap.end);
    ++work.gaps;
    cursor = gap.end;
  }
  emit_certified(cursor, length);

  ExecutorResult out;
  out.raw.has_cigar = !packed.empty();
  out.raw.zdropped = false;
  out.raw.reach_end = true;
  out.raw.packed_cigar = std::move(packed);
  const int total = static_cast<int>(std::min<std::int64_t>(
      std::numeric_limits<int>::max(),
      std::max<std::int64_t>(std::numeric_limits<int>::min(), score)));
  out.raw.max_score = total;
  out.raw.max_query = length - 1;
  out.raw.max_target = length - 1;
  out.raw.query_end_score = total;
  out.raw.query_end_target = length - 1;
  out.raw.target_end_score = total;
  out.raw.target_end_query = length - 1;
  out.raw.corner_score = total;
  out.alignment.ref_offset = 0;
  out.alignment.ref_consumed = length;
  out.alignment.matches = matches;
  out.alignment.score = total;
  out.alignment.zdropped = false;
  out.alignment.max_q = -1;
  out.alignment.max_t = -1;
  trace.verified = work;
  return out;
}

void execute_supported(const ExecutorEntry& entry,
                       const RealizationRequest& request,
                       ExecutionContext& context, RealizationOutcome& outcome) {
  ExecutorResult result = entry.execute(request, context, outcome.trace);
  outcome.raw = std::move(result.raw);
  outcome.alignment = std::move(result.alignment);
  outcome.committed_executor = entry.id;
  outcome.support_reason = SupportReason::Supported;
  if (outcome.raw.zdropped) {
    outcome.kind = request.role == RealizationRole::DnaInternalFill ||
                           request.role == RealizationRole::DnaSupplementaryFill
                       ? OutcomeKind::SplitRequested
                       : OutcomeKind::Zdropped;
  } else {
    // A result that did not drop must carry a traceback, and a pinned-global
    // one must consume its whole rectangle; anything else is InvalidResult.
    int query_consumed = 0;
    int target_consumed = 0;
    const bool valid_traceback =
        outcome.raw.has_cigar && !outcome.raw.packed_cigar.empty() &&
        packed_consumption(outcome.raw, query_consumed, target_consumed);
    const bool complete_global =
        request.objective != EndpointObjective::PinnedGlobal ||
        (query_consumed == request.query.length &&
         target_consumed == request.target.length &&
         outcome.alignment.ref_offset == 0 &&
         outcome.alignment.ref_consumed == request.target.length);
    outcome.kind = valid_traceback && complete_global
                       ? OutcomeKind::Aligned
                       : OutcomeKind::InvalidResult;
  }
  // Fold op codes past N into M, after validation saw the raw codes.
  for (std::uint32_t& run : outcome.raw.packed_cigar)
    if ((run & UINT32_C(0xf)) > 3)
      run &= ~UINT32_C(0xf);
}

} // namespace

void ExecutionTrace::push(KernelAttempt attempt) {
  if (attempt_count >= attempts.size()) {
    overflow = true;
    return;
  }
  attempts[attempt_count++] = std::move(attempt);
}

namespace {

EndpointObjective objective_for_role(RealizationRole role) {
  switch (role) {
  case RealizationRole::DnaLeftExtension:
    return EndpointObjective::OpenEndedLeft;
  case RealizationRole::DnaRightExtension:
    return EndpointObjective::OpenEndedRight;
  case RealizationRole::DnaLocalInversionMiddle:
    return EndpointObjective::OpenEndedRight;
  case RealizationRole::DnaInternalFill:
  case RealizationRole::DnaSupplementaryFill:
  case RealizationRole::RnaExonFill:
    return EndpointObjective::PinnedGlobal;
  }
  return EndpointObjective::PinnedGlobal;
}

} // namespace

const char* role_name(RealizationRole role) {
  switch (role) {
  case RealizationRole::DnaLeftExtension:
    return "DNA_LEFT_EXTENSION";
  case RealizationRole::DnaInternalFill:
    return "DNA_INTERNAL_FILL";
  case RealizationRole::DnaRightExtension:
    return "DNA_RIGHT_EXTENSION";
  case RealizationRole::DnaSupplementaryFill:
    return "DNA_SUPPLEMENTARY_FILL";
  case RealizationRole::DnaLocalInversionMiddle:
    return "DNA_LOCAL_INVERSION_MIDDLE";
  case RealizationRole::RnaExonFill:
    return "RNA_EXON_FILL";
  }
  return "INVALID_ROLE";
}

SequenceSlice make_query_slice(const std::uint8_t* data, int length,
                               int strand_begin, int read_length,
                               bool reverse_complemented) {
  SequenceSlice slice;
  slice.data = data;
  slice.length = length;
  slice.orientation = reverse_complemented
                          ? SequenceOrientation::ReverseComplement
                          : SequenceOrientation::Forward;
  if (reverse_complemented) {
    slice.original_begin = read_length - strand_begin - length;
    slice.original_end = read_length - strand_begin;
  } else {
    slice.original_begin = strand_begin;
    slice.original_end = strand_begin + length;
  }
  return slice;
}

SequenceSlice make_target_slice(const std::uint8_t* data, int length,
                                int reference_begin) {
  SequenceSlice slice;
  slice.data = data;
  slice.length = length;
  slice.original_begin = reference_begin;
  slice.original_end = reference_begin + length;
  return slice;
}

namespace {

const ControllerPolicy& production_policy() {
  static const ControllerPolicy policy{};
  return policy;
}

const ExecutorTable& production_executors() {
  static const std::array<ExecutorEntry, 2> entries = {
      ExecutorEntry{ExecutorId::Ksw2, "KSW2", ksw2_supports, ksw2_execute},
      ExecutorEntry{ExecutorId::Verified, "VERIFIED", verified_supports,
                    verified_execute}};
  static const ExecutorTable table{entries.data(), entries.size()};
  return table;
}

} // namespace

RealizationOutcome run(const RealizationRequest& request,
                       const ControllerPolicy& policy,
                       const ExecutorTable& executors,
                       ExecutionContext& context) {
  RealizationOutcome outcome;
  // A verified region always goes to the Verified executor; everything else
  // to the policy's executor.
  const ExecutorId requested_executor =
      request.verified ? ExecutorId::Verified : policy.selected_executor;
  outcome.requested_executor = requested_executor;
  outcome.selected_executor = requested_executor;
  outcome.selection_reason = requested_executor == ExecutorId::Ksw2
                                 ? SelectionReason::ProductionDefault
                                 : SelectionReason::RequestedExecutor;

  if (!valid_role(request.role) || !valid_slice(request.query) ||
      !valid_slice(request.target) || !role_matches_objective(request) ||
      !valid_gap_scoring(request.scoring)) {
    outcome.kind = OutcomeKind::InvalidPacket;
    outcome.decision = ControllerDecision::RejectInvalid;
    outcome.support_reason = SupportReason::InvalidPacket;
    return finish(context, request, std::move(outcome));
  }

  if (request.query.length == 0 && request.target.length == 0) {
    outcome.kind = OutcomeKind::InvalidPacket;
    outcome.decision = ControllerDecision::RejectInvalid;
    outcome.support_reason = SupportReason::InvalidPacket;
    return finish(context, request, std::move(outcome));
  }
  if (request.target.length == 0) {
    if (!materialize_pure_axis(request, outcome)) {
      outcome.kind = OutcomeKind::InvalidPacket;
      outcome.decision = ControllerDecision::RejectInvalid;
      outcome.support_reason = SupportReason::InvalidPacket;
      return finish(context, request, std::move(outcome));
    }
    outcome.kind = OutcomeKind::PureInsertion;
    outcome.decision = ControllerDecision::PureInsertion;
    outcome.support_reason = SupportReason::Supported;
    return finish(context, request, std::move(outcome));
  }
  if (request.query.length == 0) {
    if (!materialize_pure_axis(request, outcome)) {
      outcome.kind = OutcomeKind::InvalidPacket;
      outcome.decision = ControllerDecision::RejectInvalid;
      outcome.support_reason = SupportReason::InvalidPacket;
      return finish(context, request, std::move(outcome));
    }
    outcome.kind = OutcomeKind::PureDeletion;
    outcome.decision = ControllerDecision::PureDeletion;
    outcome.support_reason = SupportReason::Supported;
    return finish(context, request, std::move(outcome));
  }

  const bool inversion_role =
      request.role == RealizationRole::DnaInternalFill ||
      request.role == RealizationRole::DnaSupplementaryFill;
  const bool invalid_inversion_control =
      request.inversion_probe_enabled &&
      (!inversion_role || request.inversion_zdrop < 0 ||
       request.inversion_max_gap <= 0 ||
       request.inversion_min_chain_score <= 0 ||
       request.inversion_min_dp_max < 0);
  if (request.caller_requested_band < 0 || request.configured_band < 0 ||
      request.configured_long_band < 0 || request.selected_band < 0 ||
      invalid_inversion_control) {
    outcome.kind = OutcomeKind::InvalidPacket;
    outcome.decision = ControllerDecision::RejectInvalid;
    outcome.support_reason = SupportReason::InvalidPacket;
    return finish(context, request, std::move(outcome));
  }

  // Check executor support before the matrix cap, so a malformed packet is
  // not reported as a cap refusal.
  SupportReason reason = SupportReason::ExecutorUnavailable;
  const ExecutorEntry* selected = find_executor(executors, requested_executor);
  if (!supports_entry(selected, request, reason)) {
    outcome.fallback_reason = reason;
    const ExecutorEntry* fallback = find_executor(executors, ExecutorId::Ksw2);
    SupportReason fallback_support = SupportReason::ExecutorUnavailable;
    if (requested_executor == ExecutorId::Ksw2 ||
        !supports_entry(fallback, request, fallback_support)) {
      outcome.kind = OutcomeKind::UnsupportedExecutor;
      outcome.decision = ControllerDecision::UnsupportedExecutor;
      outcome.support_reason =
          requested_executor == ExecutorId::Ksw2 ? reason : fallback_support;
      return finish(context, request, std::move(outcome));
    }
    selected = fallback;
    outcome.selection_reason = SelectionReason::UnsupportedExecutorFallback;
    outcome.selected_executor = ExecutorId::Ksw2;
  }

  // The matrix cap. Fills are not refused here: the kernel returns an
  // over-cap fill Z-dropped at its start (as minimap2's mm_align_pair) and
  // the block splits there. Extensions and the inversion middle are refused,
  // which their callers read as no extension and no middle.
  const bool fill_role = request.role == RealizationRole::DnaInternalFill ||
                         request.role == RealizationRole::DnaSupplementaryFill;
  if (!fill_role && matrix_cap_exceeded(request)) {
    outcome.kind = OutcomeKind::RefusedMatrixCap;
    outcome.decision = ControllerDecision::RefuseMatrixCap;
    outcome.support_reason = SupportReason::Supported;
    // The result ksw_reset_extz() would leave.
    outcome.raw.zdropped = true;
    outcome.raw.query_end_score = KSW_NEG_INF;
    outcome.raw.target_end_score = KSW_NEG_INF;
    outcome.raw.corner_score = KSW_NEG_INF;
    return finish(context, request, std::move(outcome));
  }

  execute_supported(*selected, request, context, outcome);

  if (outcome.kind == OutcomeKind::InvalidResult) {
    outcome.decision = ControllerDecision::RejectInvalidResult;
  } else {
    outcome.decision = outcome.committed_executor == ExecutorId::Ksw2
                           ? ControllerDecision::FireKsw2
                           : ControllerDecision::FireExecutor;
  }

  return finish(context, request, std::move(outcome));
}

} // namespace realization
} // namespace lr
} // namespace cpu
} // namespace fa

/**** dp_runner ****/

namespace fa {
namespace cpu {
namespace lr {

// As in minimap2: flank extensions use the band bw_eff() and gap fills
// bw_long_eff() (max(q_span, r_span) across a long join); a gap fill runs the
// KSW_EZ_APPROX_MAX pass first; a fill whose matrix exceeds opt.max_sw_mat is
// not run but returned Z-dropped.
realization::RealizationOutcome run_dna_long_realization(
    realization::RealizationRole role, const DpScoringParams &dp,
    realization::SequenceSlice query, realization::SequenceSlice target,
    int caller_requested_band, int interior_zdrop, bool long_join, realization::ExecutionContext *execution_context,
    bool inversion_probe_enabled, const VerifiedRegionInputs *verified_region,
    const realization::InversionProbeGate *inversion_probe_gate) {
  realization::RealizationRequest request;
  request.role = role;
  request.objective = realization::objective_for_role(role);
  request.query = query;
  request.target = target;
  const bool left =
      request.objective == realization::EndpointObjective::OpenEndedLeft;
  request.query.reversed_for_executor = left;
  request.target.reversed_for_executor = left;

  request.scoring.match = dp.match;
  request.scoring.mismatch = dp.mismatch;
  request.scoring.ambiguity = dp.ambi;
  ::fa::cpu::ksw2_simple_mat(request.scoring.substitution_matrix.data(),
                             dp.match, dp.mismatch, dp.ambi);
  request.scoring.gap1 = {dp.gap_open1, dp.gap_extend1};
  request.scoring.gap2 = {dp.gap_open2, dp.gap_extend2};
  request.scoring.model =
      dp.gap_open1 == dp.gap_open2 && dp.gap_extend1 == dp.gap_extend2
          ? realization::ScoreModel::SingleAffine
          : realization::ScoreModel::DualAffine;

  ::fa::cpu::DpMapOpt opt;
  opt.bw = dp.bw;
  opt.bw_long = dp.bw_long;
  request.caller_requested_band = caller_requested_band;
  request.configured_band = dp.bw;
  request.configured_long_band = dp.bw_long;
  request.selected_band =
      role == realization::RealizationRole::DnaLocalInversionMiddle
          ? static_cast<int>(dp.bw * 1.5)
          : request.objective == realization::EndpointObjective::PinnedGlobal
          ? (long_join ? std::max(query.length, target.length)
                       : opt.bw_long_eff())
          : opt.bw_eff();
  request.zdrop =
      request.objective == realization::EndpointObjective::PinnedGlobal
          ? interior_zdrop
          : dp.zdrop;
  request.inversion_zdrop = dp.inversion_zdrop;
  const bool inversion_role =
      role == realization::RealizationRole::DnaInternalFill ||
      role == realization::RealizationRole::DnaSupplementaryFill;
  request.inversion_probe_enabled =
      inversion_role && interior_zdrop >= 0 && inversion_probe_enabled;
  request.inversion_max_gap = dp.inversion_max_gap;
  request.inversion_min_chain_score = dp.inversion_min_chain_score;
  request.inversion_min_dp_max = dp.inversion_min_dp_max;
  if (inversion_probe_gate != nullptr)
    request.inversion_probe_gate = *inversion_probe_gate;
  request.end_bonus =
      request.objective == realization::EndpointObjective::PinnedGlobal ||
              role == realization::RealizationRole::DnaLocalInversionMiddle
          ? -1
          : dp.tail_end_bonus;
  request.matrix_cell_cap = opt.max_sw_mat;
  request.long_join = long_join;
  request.ksw_flags =
      request.objective == realization::EndpointObjective::PinnedGlobal
          ? KSW_EZ_APPROX_MAX
          : KSW_EZ_EXTZ_ONLY;
  if (left)
    request.ksw_flags |= KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR;
  if (verified_region != nullptr) {
    request.verified = true;
    request.gaps = verified_region->gaps;
    request.gap_count = verified_region->gap_count;
  }

  realization::ExecutionContext local_context;
  realization::ExecutionContext &context =
      execution_context == nullptr ? local_context : *execution_context;
  return realization::run(request, realization::production_policy(),
                          realization::production_executors(), context);
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** postdp_scoring ****/

namespace fa::cpu::lr {

namespace {

// The terms minimap2's mm_update_extra derives from one alignment. Bases are
// encoded 0..3; above 3 is ambiguous.
struct DpMaxSweep {
  int dp_max = 0; // max-scoring sub-segment, log-gap cost, (int)(max + .499)
  int mlen = 0;   // matched bases, ambiguous excluded
  int blen = 0;   // aligned block, ambiguous excluded
  int n_ambi = 0;
  int n_gap = 0;  // gap bases
  int n_gapo = 0; // gap opens
  // Sum over gaps of log2(1 + len): minibwa's gap cost is b2 * this.
  double log_gap_sum = 0.0;
};

// mm_update_extra's sweep with log gaps: add the substitution score per
// aligned base, subtract q + e * log2(1 + len) per gap, reset at 0 and keep
// the maximum. `query` is the oriented query from base 0 and `reference`
// starts at the record's pos.
DpMaxSweep dp_max_sweep(const std::vector<std::pair<int, char>>& ops,
                        const std::uint8_t* query,
                        const std::uint8_t* reference, int match, int mismatch,
                        int gap_open, int gap_extend) noexcept {
  DpMaxSweep out;
  if (query == nullptr || reference == nullptr)
    return out;
  // Signs as in ksw_gen_simple_mat.
  const double a = std::abs(match);
  const double b = -std::abs(mismatch);
  const double ambi = -1.0;
  const double q = static_cast<double>(gap_open);
  const double e = static_cast<double>(gap_extend);
  double s = 0.0, max = 0.0;
  std::size_t qoff = 0, toff = 0;
  for (const auto& op : ops) {
    const int len = op.first;
    if (len <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X': {
      int n_ambi = 0, n_diff = 0;
      for (int l = 0; l < len; ++l) {
        const int cq = query[qoff + static_cast<std::size_t>(l)];
        const int ct = reference[toff + static_cast<std::size_t>(l)];
        if (ct > 3 || cq > 3) {
          ++n_ambi;
          s += ambi;
        } else if (ct != cq) {
          ++n_diff;
          s += b;
        } else {
          s += a;
        }
        if (s < 0)
          s = 0;
        else
          max = max > s ? max : s;
      }
      out.blen += len - n_ambi;
      out.mlen += len - (n_ambi + n_diff);
      out.n_ambi += n_ambi;
      qoff += static_cast<std::size_t>(len);
      toff += static_cast<std::size_t>(len);
      break;
    }
    case 'I': {
      int n_ambi = 0;
      for (int l = 0; l < len; ++l)
        if (query[qoff + static_cast<std::size_t>(l)] > 3)
          ++n_ambi;
      out.blen += len - n_ambi;
      out.n_ambi += n_ambi;
      const double lg = std::log2(1.0 + static_cast<double>(len));
      s -= q + e * lg;
      if (s < 0)
        s = 0;
      out.n_gap += len;
      out.n_gapo += 1;
      out.log_gap_sum += lg;
      qoff += static_cast<std::size_t>(len);
      break;
    }
    case 'D': {
      int n_ambi = 0;
      for (int l = 0; l < len; ++l)
        if (reference[toff + static_cast<std::size_t>(l)] > 3)
          ++n_ambi;
      out.blen += len - n_ambi;
      out.n_ambi += n_ambi;
      const double lg = std::log2(1.0 + static_cast<double>(len));
      s -= q + e * lg;
      if (s < 0)
        s = 0;
      out.n_gap += len;
      out.n_gapo += 1;
      out.log_gap_sum += lg;
      toff += static_cast<std::size_t>(len);
      break;
    }
    case 'N':
      toff += static_cast<std::size_t>(len);
      break;
    case 'S':
      qoff += static_cast<std::size_t>(len);
      break;
    default: // H, P: consume nothing
      break;
    }
  }
  out.dp_max = truncate_to_int(max + 0.499);
  return out;
}

// minimap2's mm_event_identity: mlen / (blen + n_ambi - n_gap + n_gapo), or
// 0 on an empty denominator.
double event_identity(const DpMaxSweep& sweep) noexcept {
  const int denominator =
      sweep.blen + sweep.n_ambi - sweep.n_gap + sweep.n_gapo;
  if (denominator <= 0)
    return 0.0;
  return static_cast<double>(sweep.mlen) / static_cast<double>(denominator);
}

// b2 as in minimap2's mm_update_dp_max: 0.5 / max(0.02, 1 - identity), with
// its b2 * a < b guard.
double rescore_b2(double best_identity, int match, int mismatch) noexcept {
  double div = 1.0 - best_identity;
  if (div < 0.02)
    div = 0.02;
  double b2 = 0.5 / div;
  if (b2 * static_cast<double>(match) < static_cast<double>(mismatch))
    b2 = static_cast<double>(match) / static_cast<double>(mismatch);
  return b2;
}

// The record's contig index by name, or -1 when unknown.
int record_contig_index(const DnaContext& context, const std::string& name) {
  if (context.ref.names == nullptr)
    return -1;
  const auto& names = *context.ref.names;
  for (std::size_t index = 0; index < names.size(); ++index)
    if (names[index] == name)
      return static_cast<int>(index);
  return -1;
}

struct RecordPrice {
  DpMaxSweep sweep;
  // False when the record has no CIGAR or its contig is unknown, so that
  // dna_record_dp_max_segment can report -1 instead of 0.
  bool swept = false;
};

RecordPrice price_record(const DnaContext& context, const AlignResult& record,
                         const std::vector<std::uint8_t>& fwd,
                         const std::vector<std::uint8_t>& rc) {
  RecordPrice price;
  const int contig = record_contig_index(context, record.chromosome);
  if (context.ref.encoded != nullptr && contig >= 0 &&
      static_cast<std::size_t>(contig) < context.ref.encoded->size() &&
      record.pos >= 0 && !record.cigar.empty()) {
    const std::vector<std::uint8_t>& reference =
        (*context.ref.encoded)[static_cast<std::size_t>(contig)];
    const std::vector<std::uint8_t>& query = record.is_reverse ? rc : fwd;
    price.sweep = dp_max_sweep(::fa::cpu::output::parse_cigar_ops(record.cigar),
                               query.data(), reference.data() + record.pos,
                               context.opts.cigar_dp_match,
                               context.opts.cigar_dp_mismatch,
                               context.opts.cigar_dp_gap_open1,
                               context.opts.cigar_dp_gap_extend1);
    price.swept = true;
  }
  return price;
}

} // namespace

int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc) {
  if (!record.mapped())
    return -1;
  const RecordPrice price = price_record(context, record, fwd, rc);
  return price.swept ? price.sweep.dp_max : -1;
}

bool dna_record_event_identity(const DnaContext& context,
                               const AlignResult& record,
                               const std::vector<std::uint8_t>& fwd,
                               const std::vector<std::uint8_t>& rc,
                               double& identity) {
  if (!record.mapped())
    return false;
  const RecordPrice price = price_record(context, record, fwd, rc);
  if (!price.swept)
    return false;
  identity = event_identity(price.sweep);
  return true;
}

int dna_record_recal_dp_max(const DnaContext& context,
                            const AlignResult& record,
                            const std::vector<std::uint8_t>& fwd,
                            const std::vector<std::uint8_t>& rc, double b2) {
  if (!record.mapped())
    return -1;
  const RecordPrice price = price_record(context, record, fwd, rc);
  if (!price.swept)
    return -1;
  const DpMaxSweep& sweep = price.sweep;
  const double gap_cost =
      b2 * static_cast<double>(sweep.n_gapo) + sweep.log_gap_sum;
  const int n_mis = sweep.blen + sweep.n_ambi - sweep.mlen - sweep.n_gap;
  const int recal = truncate_to_int(
      static_cast<double>(std::abs(context.opts.cigar_dp_match)) *
          (static_cast<double>(sweep.mlen) - b2 * n_mis - gap_cost) +
      .499);
  return std::max(0, recal);
}

double dna_rank_b2(double identity, int match, int mismatch) noexcept {
  return rescore_b2(identity, match, mismatch);
}

} // namespace fa::cpu::lr
