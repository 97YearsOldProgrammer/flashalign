// The realization controller: validates one DP packet, picks an executor
// (ksw2, or base comparison for a verified region) and returns a packet-local
// result. Slicing, concatenation, clipping and commit are the caller's.
#pragma once

#include "../dp/result.h"    // ::fa::cpu::GotohCigarResult
#include "region_gap.h"      // ::fa::cpu::lr::RegionGap

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace realization {

enum class RealizationRole : std::uint8_t {
  DnaLeftExtension,
  DnaInternalFill,
  DnaRightExtension,
  DnaSupplementaryFill,
  DnaLocalInversionMiddle,
  RnaExonFill,
};

enum class EndpointObjective : std::uint8_t {
  OpenEndedLeft,
  PinnedGlobal,
  OpenEndedRight,
};

enum class SequenceOrientation : std::uint8_t {
  Forward,
  ReverseComplement,
};

struct SequenceSlice {
  const std::uint8_t* data = nullptr;
  int length = 0;
  int original_begin = 0;
  int original_end = 0;
  SequenceOrientation orientation = SequenceOrientation::Forward;
  bool reversed_for_executor = false;
};

struct AffineGapRow {
  int open = 0;
  int extend = 0;
};

enum class ScoreModel : std::uint8_t {
  SingleAffine,
  DualAffine,
  Splice,
};

struct ScoringContract {
  int match = 0;
  int mismatch = 0;
  int ambiguity = 0;
  std::array<std::int8_t, 25> substitution_matrix{};
  AffineGapRow gap1;
  AffineGapRow gap2;
  ScoreModel model = ScoreModel::DualAffine;
};

struct RealizationRequest {
  RealizationRole role = RealizationRole::DnaInternalFill;
  EndpointObjective objective = EndpointObjective::PinnedGlobal;
  SequenceSlice query;
  SequenceSlice target;
  ScoringContract scoring;

  int caller_requested_band = -1;
  int configured_band = -1;
  int configured_long_band = -1;
  int selected_band = -1;
  int zdrop = -1;
  int inversion_zdrop = -1;
  bool inversion_probe_enabled = false;
  int inversion_max_gap = 0;
  int inversion_min_chain_score = 0;
  int inversion_min_dp_max = 0;
  int end_bonus = -1;
  std::int64_t matrix_cell_cap = 0;
  bool long_join = false;
  std::uint32_t ksw_flags = 0;

  // A verified region: one diagonal, every base under an exact-match anchor
  // or in one of the certified gaps below. Runs on ExecutorId::Verified.
  bool verified = false;
  // Region-relative certified gaps, borrowed for the call. Empty when the
  // region's anchors overlap end to end.
  const ::fa::cpu::lr::RegionGap* gaps = nullptr;
  int gap_count = 0;
};

enum class ControllerDecision : std::uint8_t {
  FireKsw2,
  FireExecutor,
  PureInsertion,
  PureDeletion,
  RefuseMatrixCap,
  RejectInvalid,
  RejectInvalidResult,
  UnsupportedExecutor,
};

enum class OutcomeKind : std::uint8_t {
  Aligned,
  Zdropped,
  SplitRequested,
  PureInsertion,
  PureDeletion,
  RefusedMatrixCap,
  InvalidPacket,
  InvalidResult,
  UnsupportedExecutor,
};

enum class ExecutorId : std::uint8_t {
  None,
  Ksw2,
  // Base comparison over a verified region; no DP.
  Verified,
};

// The ksw2 kernels dp/control.h can run (DpKernelId), plus None.
enum class KernelId : std::uint8_t {
  None,
  Extz2,
  Extd2,
};

enum class SupportReason : std::uint8_t {
  Supported,
  InvalidPacket,
  UnsupportedObjective,
  UnsupportedOrientation,
  UnsupportedScoreModel,
  UnsupportedSubstitutionMatrix,
  UnsupportedFlags,
  UnsupportedBand,
  ExecutorUnavailable,
  // A verified region handed to ksw2, or anything else handed to Verified.
  UnsupportedPacketClass,
};

enum class SelectionReason : std::uint8_t {
  ProductionDefault,
  RequestedExecutor,
  UnsupportedExecutorFallback,
};

struct RawAlignmentResult {
  bool has_cigar = false;
  bool zdropped = false;
  bool reach_end = false;
  int max_score = 0;
  int max_query = -1;
  int max_target = -1;
  int query_end_score = 0;
  int query_end_target = -1;
  int target_end_score = 0;
  int target_end_query = -1;
  int corner_score = 0;
  std::vector<std::uint32_t> packed_cigar;
};

struct KernelAttempt {
  ExecutorId executor = ExecutorId::None;
  KernelId kernel = KernelId::None;
  std::uint32_t flags = 0;
  int band = -1;
  int zdrop = -1;
  bool approximate = false;
  bool exact = false;
  // raw.packed_cigar is empty even when raw.has_cigar: an attempt keeps only
  // metadata and scores.
  RawAlignmentResult raw;
};

enum class InversionProbeStatus : std::uint8_t {
  Disabled,
  DropBelowInversionThreshold,
  InvalidBounds,
  WeakReverseScore,
  Accepted,
};

struct ZdropTestEvidence {
  bool tested = false;
  int code = 0;
  int max_drop = 0;
  int target_begin = -1;
  int target_end = -1;
  int query_begin = -1;
  int query_end = -1;
  InversionProbeStatus inversion_status = InversionProbeStatus::Disabled;
  bool reverse_probe_attempted = false;
  int reverse_score = 0;
  int reverse_query_end = -1;
  int reverse_target_end = -1;
};

// Work of one verified region: bases scored and gaps emitted. It runs no
// kernel, so it records no attempt.
struct VerifiedWork {
  std::int64_t bases = 0;
  int gaps = 0;
};

struct ExecutionTrace {
  static constexpr std::size_t kMaxAttempts = 4;
  std::array<KernelAttempt, kMaxAttempts> attempts{};
  std::size_t attempt_count = 0;
  bool overflow = false;
  ZdropTestEvidence zdrop_test;
  bool exact_retry = false;
  int exact_retry_zdrop = -1;
  VerifiedWork verified;

  void push(KernelAttempt attempt);
};

// The traceback is `raw.packed_cigar` (BAM-packed, ops past N folded to M);
// `alignment.cigar` stays empty and `alignment` supplies only the accounting
// (ref_offset, ref_consumed, matches, score, zdropped, max_q, max_t).
struct RealizationOutcome {
  OutcomeKind kind = OutcomeKind::InvalidPacket;
  ControllerDecision decision = ControllerDecision::RejectInvalid;
  ExecutorId requested_executor = ExecutorId::None;
  ExecutorId selected_executor = ExecutorId::None;
  ExecutorId committed_executor = ExecutorId::None;
  SelectionReason selection_reason = SelectionReason::ProductionDefault;
  SupportReason support_reason = SupportReason::InvalidPacket;
  SupportReason fallback_reason = SupportReason::Supported;
  RawAlignmentResult raw;
  ::fa::cpu::GotohCigarResult alignment{};
  ExecutionTrace trace;
};

struct ControllerObserver {
  using Callback = void (*)(void* user_data, const RealizationRequest& request,
                            const RealizationOutcome& outcome);

  Callback callback = nullptr;
  void* user_data = nullptr;
};

struct ExecutionContext {
  ControllerObserver observer;
};

struct ExecutorResult {
  RawAlignmentResult raw;
  ::fa::cpu::GotohCigarResult alignment{};
};

using SupportsFn = bool (*)(const RealizationRequest& request,
                            SupportReason& reason);
using ExecuteFn = ExecutorResult (*)(const RealizationRequest& request,
                                     ExecutionContext& context,
                                     ExecutionTrace& trace);

struct ExecutorEntry {
  ExecutorId id = ExecutorId::None;
  const char* name = "none";
  SupportsFn supports = nullptr;
  ExecuteFn execute = nullptr;
};

struct ExecutorTable {
  const ExecutorEntry* entries = nullptr;
  std::size_t size = 0;
};

struct ControllerPolicy {
  ExecutorId selected_executor = ExecutorId::Ksw2;
};

EndpointObjective objective_for_role(RealizationRole role);
const char* role_name(RealizationRole role);

SequenceSlice make_query_slice(const std::uint8_t* data, int length,
                               int strand_begin, int read_length,
                               bool reverse_complemented);
SequenceSlice make_target_slice(const std::uint8_t* data, int length,
                                int reference_begin);

const ControllerPolicy& production_policy();
const ExecutorTable& production_executors();

RealizationOutcome run(const RealizationRequest& request,
                       const ControllerPolicy& policy,
                       const ExecutorTable& executors,
                       ExecutionContext& context);

} // namespace realization
} // namespace lr
} // namespace cpu
} // namespace fa
