// Saturating work counters for the RNA splice controller.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class RnaControllerRefusal : std::uint8_t {
  None = 0,
  InvalidRequest,
  UnsupportedProfile,
  UnrepresentableDomain,
  InvalidBundle,
  EmptyAnchorView,
  InvalidAnchor,
  InvalidSelectedRange,
  InvalidSequenceSlice,
  KernelFailure,
  EmptyCigar,
  CigarReconciliation,
  InvalidContinuation,
  ContinuationLimit,
  IncompleteHypothesis,
  RegionFilter,
  AllocationFailure,
  InternalInvariant,
  Count
};

constexpr std::size_t rna_controller_refusal_count() noexcept {
  return static_cast<std::size_t>(RnaControllerRefusal::Count);
}

constexpr const char* rna_controller_refusal_name(
    RnaControllerRefusal reason) noexcept {
  switch (reason) {
    case RnaControllerRefusal::None: return "none";
    case RnaControllerRefusal::InvalidRequest: return "invalid_request";
    case RnaControllerRefusal::UnsupportedProfile: return "unsupported_profile";
    case RnaControllerRefusal::UnrepresentableDomain: return "unrepresentable_domain";
    case RnaControllerRefusal::InvalidBundle: return "invalid_bundle";
    case RnaControllerRefusal::EmptyAnchorView: return "empty_anchor_view";
    case RnaControllerRefusal::InvalidAnchor: return "invalid_anchor";
    case RnaControllerRefusal::InvalidSelectedRange: return "invalid_selected_range";
    case RnaControllerRefusal::InvalidSequenceSlice: return "invalid_sequence_slice";
    case RnaControllerRefusal::KernelFailure: return "kernel_failure";
    case RnaControllerRefusal::EmptyCigar: return "empty_cigar";
    case RnaControllerRefusal::CigarReconciliation: return "cigar_reconciliation";
    case RnaControllerRefusal::InvalidContinuation: return "invalid_continuation";
    case RnaControllerRefusal::ContinuationLimit: return "continuation_limit";
    case RnaControllerRefusal::IncompleteHypothesis: return "incomplete_hypothesis";
    case RnaControllerRefusal::RegionFilter: return "region_filter";
    case RnaControllerRefusal::AllocationFailure: return "allocation_failure";
    case RnaControllerRefusal::InternalInvariant: return "internal_invariant";
    case RnaControllerRefusal::Count: return "count";
  }
  return "internal_invariant";
}

inline std::uint64_t saturating_add(std::uint64_t left,
                                    std::uint64_t right) noexcept {
  const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
  return right > maximum - left ? maximum : left + right;
}

inline std::uint64_t saturating_product(std::uint64_t left,
                                        std::uint64_t right) noexcept {
  if (left == 0 || right == 0) return 0;
  const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
  return left > maximum / right ? maximum : left * right;
}

struct RnaControllerMetrics {
  std::uint64_t top_level_controller_calls = 0;
  std::uint64_t transcript_hypothesis_passes = 0;
  std::uint64_t left_dp_calls = 0;
  std::uint64_t internal_dp_calls = 0;
  std::uint64_t right_dp_calls = 0;
  std::uint64_t approximate_passes = 0;
  std::uint64_t exact_retries = 0;
  std::uint64_t zdrop_signals = 0;
  std::uint64_t split_continuations = 0;
  std::uint64_t annotation_packets_masked = 0;
  std::uint64_t annotation_endpoint_bytes_set = 0;
  std::uint64_t estimated_dp_cells = 0;
  std::uint64_t query_bases_presented = 0;
  std::uint64_t reference_bases_presented = 0;
  std::uint64_t maximum_reference_slice = 0;
  // Internal packets by shape: intron-crossing when the reference run exceeds the query
  // run by at least one minimum intron, colinear otherwise. Counted once per logical
  // packet (an exact retry is not recounted), so the calls do not sum to internal_dp_calls.
  std::uint64_t internal_intron_crossing_calls = 0;
  std::uint64_t internal_colinear_calls = 0;
  std::uint64_t internal_intron_crossing_cells = 0;
  std::uint64_t internal_colinear_cells = 0;
  // The largest internal packet, and the packets whose traceback exceeds the DP arena and
  // is heap-allocated instead.
  std::uint64_t maximum_internal_packet_cells = 0;
  std::uint64_t internal_packets_over_arena = 0;
  std::array<std::uint64_t, rna_controller_refusal_count()> refusals{};

  void observe_refusal(RnaControllerRefusal reason) noexcept {
    const std::size_t index = static_cast<std::size_t>(reason);
    if (reason != RnaControllerRefusal::None && index < refusals.size())
      refusals[index] = saturating_add(refusals[index], 1);
  }

  void observe_dp(std::uint64_t query_bases, std::uint64_t reference_bases,
                  int role) noexcept {
    if (role == 1) left_dp_calls = saturating_add(left_dp_calls, 1);
    else if (role == 2) internal_dp_calls = saturating_add(internal_dp_calls, 1);
    else if (role == 3) right_dp_calls = saturating_add(right_dp_calls, 1);
    query_bases_presented = saturating_add(query_bases_presented, query_bases);
    reference_bases_presented = saturating_add(
        reference_bases_presented, reference_bases);
    estimated_dp_cells = saturating_add(
        estimated_dp_cells, saturating_product(query_bases, reference_bases));
    if (reference_bases > maximum_reference_slice)
      maximum_reference_slice = reference_bases;
  }

  // Charged once per internal packet before the kernel runs; independent of observe_dp.
  void observe_internal_packet(std::uint64_t query_bases,
                               std::uint64_t reference_bases,
                               bool intron_crossing,
                               bool traceback_over_arena) noexcept {
    const std::uint64_t cells = saturating_product(query_bases, reference_bases);
    if (intron_crossing) {
      internal_intron_crossing_calls =
          saturating_add(internal_intron_crossing_calls, 1);
      internal_intron_crossing_cells =
          saturating_add(internal_intron_crossing_cells, cells);
    } else {
      internal_colinear_calls = saturating_add(internal_colinear_calls, 1);
      internal_colinear_cells = saturating_add(internal_colinear_cells, cells);
    }
    if (cells > maximum_internal_packet_cells)
      maximum_internal_packet_cells = cells;
    if (traceback_over_arena)
      internal_packets_over_arena =
          saturating_add(internal_packets_over_arena, 1);
  }

  void merge(const RnaControllerMetrics& other) noexcept {
#define FA_RNA_CONTROLLER_SUM(field) \
    field = saturating_add(field, other.field)
    FA_RNA_CONTROLLER_SUM(top_level_controller_calls);
    FA_RNA_CONTROLLER_SUM(transcript_hypothesis_passes);
    FA_RNA_CONTROLLER_SUM(left_dp_calls);
    FA_RNA_CONTROLLER_SUM(internal_dp_calls);
    FA_RNA_CONTROLLER_SUM(right_dp_calls);
    FA_RNA_CONTROLLER_SUM(approximate_passes);
    FA_RNA_CONTROLLER_SUM(exact_retries);
    FA_RNA_CONTROLLER_SUM(zdrop_signals);
    FA_RNA_CONTROLLER_SUM(split_continuations);
    FA_RNA_CONTROLLER_SUM(annotation_packets_masked);
    FA_RNA_CONTROLLER_SUM(annotation_endpoint_bytes_set);
    FA_RNA_CONTROLLER_SUM(estimated_dp_cells);
    FA_RNA_CONTROLLER_SUM(query_bases_presented);
    FA_RNA_CONTROLLER_SUM(reference_bases_presented);
    FA_RNA_CONTROLLER_SUM(internal_intron_crossing_calls);
    FA_RNA_CONTROLLER_SUM(internal_colinear_calls);
    FA_RNA_CONTROLLER_SUM(internal_intron_crossing_cells);
    FA_RNA_CONTROLLER_SUM(internal_colinear_cells);
    FA_RNA_CONTROLLER_SUM(internal_packets_over_arena);
#undef FA_RNA_CONTROLLER_SUM
    if (other.maximum_reference_slice > maximum_reference_slice)
      maximum_reference_slice = other.maximum_reference_slice;
    if (other.maximum_internal_packet_cells > maximum_internal_packet_cells)
      maximum_internal_packet_cells = other.maximum_internal_packet_cells;
    for (std::size_t i = 0; i < refusals.size(); ++i)
      refusals[i] = saturating_add(refusals[i], other.refusals[i]);
  }
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
