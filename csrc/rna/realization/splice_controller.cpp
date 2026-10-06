// Splice realization of one anchor segment, following minimap2 2.30's long-read splice
// path (align.c): terminal-anchor and gap-cluster filters, packet boundaries, Z-drop
// replay and CIGAR normalization.

#include "splice_controller.h"

#include "../../core/checked_range.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <utility>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

constexpr int kCigarMatch = 0;
constexpr int kCigarInsertion = 1;
constexpr int kCigarDeletion = 2;
constexpr int kCigarReferenceSkip = 3;
constexpr std::uint32_t kMaximumCigarRun = UINT32_MAX >> 4;
// The DP arena's per-chunk ceiling (KSW_ARENA_CAPACITY_LIMIT); a larger packet is
// heap-allocated instead.
constexpr std::uint64_t kArenaCapacityLimit = std::uint64_t{32} << 20;

// The traceback matrix ksw_exts2 allocates for one packet, reproduced so the controller
// can price a packet without running it:
// n_col = (min(qlen,tlen) + 15) / 16 + 1 and
// mem2  = ((qlen + tlen - 1) * n_col + 1) * 16.
// Score-only passes allocate none of it, but every internal packet is traced back.
std::uint64_t traceback_matrix_bytes(int query_length,
                                     int reference_length) noexcept {
  const std::uint64_t query = static_cast<std::uint64_t>(query_length);
  const std::uint64_t reference = static_cast<std::uint64_t>(reference_length);
  const std::uint64_t rows = query < reference ? query : reference;
  const std::uint64_t columns = (rows + 15) / 16 + 1;
  return saturating_product(
      saturating_add(
          saturating_product(query + reference - 1, columns), 1),
      16);
}

bool has_flag(const SpliceAnchor& anchor, SpliceAnchorFlag flag) noexcept {
  return (anchor.flags & static_cast<std::uint8_t>(flag)) != 0;
}

void add_flag(SpliceAnchor& anchor, SpliceAnchorFlag flag) noexcept {
  anchor.flags |= static_cast<std::uint8_t>(flag);
}

int query_end(const SpliceAnchor& anchor) noexcept {
  return anchor.query_begin + anchor.span - 1;
}

int reference_end(const SpliceAnchor& anchor) noexcept {
  return anchor.reference_begin + anchor.span - 1;
}

bool valid_range(int begin, int end, std::size_t length) noexcept {
  return mapping::checked_half_open_range(
      begin, end, static_cast<std::uint64_t>(length));
}

bool valid_span(int begin, int count, std::size_t length) noexcept {
  return count >= 0 && mapping::checked_count_range(
      begin, static_cast<std::uint64_t>(count),
      static_cast<std::uint64_t>(length));
}

bool valid_nt4(const std::uint8_t* sequence, int length) noexcept {
  if (!sequence || length <= 0) return false;
  for (int index = 0; index < length; ++index)
    if (sequence[index] > 4) return false;
  return true;
}

RnaControllerRefusal kernel_refusal(dp::SpliceKernelFailure failure) noexcept {
  if (failure == dp::SpliceKernelFailure::AllocationFailure)
    return RnaControllerRefusal::AllocationFailure;
  if (failure == dp::SpliceKernelFailure::InvalidRequest)
    return RnaControllerRefusal::InvalidRequest;
  return RnaControllerRefusal::KernelFailure;
}

bool valid_selected_range(int begin, int count, std::size_t size) noexcept {
  return count > 0 && mapping::checked_count_range(
      begin, static_cast<std::uint64_t>(count),
      static_cast<std::uint64_t>(size));
}

bool options_are_supported(const SpliceControllerOptions& value) noexcept {
  const SpliceControllerOptions expected =
      ordinary_long_read_splice_options();
  const std::int64_t first_gap =
      static_cast<std::int64_t>(value.gap_open) + value.gap_extend;
  const std::int64_t second_gap =
      static_cast<std::int64_t>(value.long_gap_open) + value.long_gap_extend;
  const bool valid_dual_gap =
      value.gap_extend > value.long_gap_extend &&
      first_gap < second_gap;
  const bool safe_kernel_row =
      value.match > 0 && value.match <= 127 &&
      value.mismatch > 0 && value.mismatch <= 127 &&
      value.ambiguous >= 0 && value.ambiguous < value.mismatch &&
      value.gap_open > 0 && value.gap_extend > 0 &&
      value.long_gap_open >= 0 && value.long_gap_extend == 0 &&
      valid_dual_gap && first_gap + second_gap <= 127 &&
      value.junction_bonus >= 0 && value.junction_bonus <= 127 &&
      value.junction_penalty >= 0 && value.junction_penalty <= 127;
  return value.k > 0 && safe_kernel_row &&
      value.transition == expected.transition &&
      value.zdrop >= 0 && value.inversion_zdrop >= 0 &&
      value.zdrop >= value.inversion_zdrop && value.end_bonus >= -1 &&
      value.maximum_gap > 0 && value.maximum_reference_gap > 0 &&
      value.minimum_anchor_count == expected.minimum_anchor_count &&
      value.minimum_match_bases == expected.minimum_match_bases &&
      value.minimum_dp_maximum >= 0 &&
      value.minimum_packet_length == expected.minimum_packet_length &&
      value.minimum_intron > 0 &&
      value.anchor_extension_length == expected.anchor_extension_length &&
      value.anchor_extension_shift == expected.anchor_extension_shift &&
      value.maximum_dp_cells == expected.maximum_dp_cells;
}

dp::SpliceKernelScoring kernel_scoring(
    const SpliceControllerOptions& options) noexcept {
  dp::SpliceKernelScoring scoring;
  scoring.match = options.match;
  scoring.mismatch = options.mismatch;
  scoring.gap_open = options.gap_open;
  scoring.gap_extend = options.gap_extend;
  scoring.long_gap_open = options.long_gap_open;
  scoring.ambiguous = options.ambiguous;
  scoring.junction_bonus = options.junction_bonus;
  scoring.junction_penalty = options.junction_penalty;
  return scoring;
}

void refuse(SpliceSegmentResult& result,
            RnaControllerRefusal reason) noexcept {
  result.refused = true;
  result.refusal = reason;
  result.cigar.clear();
  result.continuation = {};
  result.metrics.observe_refusal(reason);
}

bool validate_request(const SpliceControllerRequest& request,
                      SpliceSegmentResult& result) noexcept {
  if (!request.query_forward || !request.query_reverse ||
      !request.reference || !request.working_anchors) {
    refuse(result, RnaControllerRefusal::InvalidRequest);
    return false;
  }
  if (!options_are_supported(request.options)) {
    refuse(result, RnaControllerRefusal::UnsupportedProfile);
    return false;
  }
  const std::size_t query_length = request.query_forward->size();
  const std::size_t reference_length = request.reference->size();
  if (query_length == 0 || reference_length == 0 ||
      request.query_reverse->size() != query_length) {
    refuse(result, RnaControllerRefusal::InvalidRequest);
    return false;
  }
  const std::size_t int_max = static_cast<std::size_t>(
      std::numeric_limits<int>::max());
  if (query_length > int_max || reference_length > int_max ||
      request.working_anchors->size() > int_max) {
    refuse(result, RnaControllerRefusal::UnrepresentableDomain);
    return false;
  }
  if (!valid_selected_range(request.selected_begin, request.selected_count,
                            request.working_anchors->size())) {
    refuse(result, RnaControllerRefusal::InvalidSelectedRange);
    return false;
  }
  // The sequences are nt4 by construction; validate pointers, lengths, the domain and
  // every slice below instead of rescanning a chromosome per request.
  for (const SpliceAnchor& anchor : *request.working_anchors) {
    const std::int64_t query_limit =
        static_cast<std::int64_t>(anchor.query_begin) + anchor.span;
    const std::int64_t reference_limit =
        static_cast<std::int64_t>(anchor.reference_begin) + anchor.span;
    if (anchor.span != request.options.k || anchor.query_begin < 0 ||
        anchor.reference_begin < 0 ||
        query_limit > static_cast<std::int64_t>(query_length) ||
        reference_limit > static_cast<std::int64_t>(reference_length)) {
      refuse(result, RnaControllerRefusal::InvalidAnchor);
      return false;
    }
  }
  const std::vector<SpliceAnchor>& anchors = *request.working_anchors;
  const int selected_end = request.selected_begin + request.selected_count;
  for (int index = request.selected_begin + 1; index < selected_end; ++index) {
    if (query_end(anchors[index]) <= query_end(anchors[index - 1]) ||
        reference_end(anchors[index]) <= reference_end(anchors[index - 1])) {
      refuse(result, RnaControllerRefusal::InvalidAnchor);
      return false;
    }
  }
  return true;
}

bool anchor_extension_score(
    const SpliceControllerRequest& request, const SpliceAnchor& anchor,
    const std::array<std::int8_t, 25>& matrix, int& score,
    RnaControllerRefusal& refusal) noexcept {
  const int query_length =
      static_cast<int>(request.query_forward->size());
  const int reference_length = static_cast<int>(request.reference->size());
  int reference_begin = anchor.reference_begin;
  int reference_limit = reference_end(anchor) + 1;
  int query_begin = anchor.query_begin;
  int query_limit = query_end(anchor) + 1;
  const int extension = request.options.anchor_extension_length;
  reference_begin = reference_begin - extension > 0
      ? reference_begin - extension : 0;
  query_begin = query_begin - extension > 0 ? query_begin - extension : 0;
  reference_limit = reference_limit + extension < reference_length
      ? reference_limit + extension : reference_length;
  query_limit = query_limit + extension < query_length
      ? query_limit + extension : query_length;
  if (!valid_range(query_begin, query_limit,
                   request.query_forward->size()) ||
      !valid_range(reference_begin, reference_limit,
                   request.reference->size())) {
    refusal = RnaControllerRefusal::InvalidSequenceSlice;
    return false;
  }
  const std::vector<std::uint8_t>& query = request.mapping_reverse
      ? *request.query_reverse : *request.query_forward;
  if (!valid_nt4(query.data() + query_begin, query_limit - query_begin) ||
      !valid_nt4(request.reference->data() + reference_begin,
                 reference_limit - reference_begin)) {
    refusal = RnaControllerRefusal::InvalidRequest;
    return false;
  }
  const dp::LocalAlignmentScore aligned = dp::local_alignment_score(
      query.data() + query_begin, query_limit - query_begin,
      request.reference->data() + reference_begin,
      reference_limit - reference_begin, matrix,
      request.options.gap_open, request.options.gap_extend);
  if (!aligned.succeeded) {
    refusal = kernel_refusal(aligned.failure);
    return false;
  }
  score = aligned.score;
  return true;
}

bool filter_weak_terminal_anchors(
    const SpliceControllerRequest& request,
    const std::array<std::int8_t, 25>& matrix,
    int& filtered_begin, int& filtered_count,
    RnaControllerRefusal& refusal) noexcept {
  filtered_begin = request.selected_begin;
  filtered_count = request.selected_count;
  if (request.selected_count < 3) return true;
  const std::vector<SpliceAnchor>& anchors = *request.working_anchors;
  const int original_end = request.selected_begin + request.selected_count;

  const SpliceAnchor& first = anchors[request.selected_begin];
  const SpliceAnchor& second = anchors[request.selected_begin + 1];
  const double left_log_gap = std::log(static_cast<double>(
      reference_end(second) - reference_end(first)));
  if (static_cast<double>(first.span) <
      left_log_gap + request.options.anchor_extension_shift) {
    int score = 0;
    if (!anchor_extension_score(request, first, matrix, score, refusal))
      return false;
    if (static_cast<double>(score) / matrix[0] <
        left_log_gap + request.options.anchor_extension_shift) {
      ++filtered_begin;
      --filtered_count;
    }
  }

  const SpliceAnchor& last = anchors[original_end - 1];
  const SpliceAnchor& before_last = anchors[original_end - 2];
  const double right_log_gap = std::log(static_cast<double>(
      reference_end(last) - reference_end(before_last)));
  if (static_cast<double>(last.span) <
      right_log_gap + request.options.anchor_extension_shift) {
    int score = 0;
    if (!anchor_extension_score(request, last, matrix, score, refusal))
      return false;
    if (static_cast<double>(score) / matrix[0] <
        right_log_gap + request.options.anchor_extension_shift)
      --filtered_count;
  }
  if (filtered_count <= 0) {
    refusal = RnaControllerRefusal::EmptyAnchorView;
    return false;
  }
  return true;
}

std::vector<int> collect_long_gaps(
    const std::vector<SpliceAnchor>& anchors, int begin, int count,
    int minimum_gap) {
  std::vector<int> positions;
  for (int offset = 1; offset < count; ++offset) {
    const SpliceAnchor& current = anchors[begin + offset];
    const SpliceAnchor& previous = anchors[begin + offset - 1];
    const std::int64_t gap =
        static_cast<std::int64_t>(query_end(current)) - query_end(previous) -
        (static_cast<std::int64_t>(reference_end(current)) -
         reference_end(previous));
    if (gap < -minimum_gap || gap > minimum_gap)
      positions.push_back(offset);
  }
  if (positions.size() <= 1) positions.clear();
  return positions;
}

void filter_balanced_gap_clusters(std::vector<SpliceAnchor>& anchors,
                                  int begin, int count,
                                  int minimum_gap, int difference_threshold,
                                  int maximum_extension_length,
                                  int maximum_extension_count) {
  const std::vector<int> gaps =
      collect_long_gaps(anchors, begin, count, minimum_gap);
  if (gaps.empty()) return;
  int best = 0;
  int best_start = -1;
  int best_end = -1;
  for (int gap_index = 0;; ++gap_index) {
    if (gap_index == static_cast<int>(gaps.size()) ||
        gap_index >= best_end) {
      if (best_end > 0) {
        for (int offset = gaps[best_start]; offset < gaps[best_end]; ++offset)
          add_flag(anchors[begin + offset], kSpliceAnchorIgnore);
      }
      best = 0;
      best_start = best_end = -1;
      if (gap_index == static_cast<int>(gaps.size())) break;
    }
    const int offset = gaps[gap_index];
    const SpliceAnchor& current = anchors[begin + offset];
    const SpliceAnchor& previous = anchors[begin + offset - 1];
    std::int64_t gap =
        static_cast<std::int64_t>(query_end(current)) - query_end(previous) -
        (static_cast<std::int64_t>(reference_end(current)) -
         reference_end(previous));
    std::int64_t inserted = gap > 0 ? gap : 0;
    std::int64_t deleted = gap < 0 ? -gap : 0;
    const int query_start = query_end(previous);
    const int reference_start = reference_end(previous);
    int maximum_difference = 0;
    int maximum_difference_index = -1;
    for (int lookahead = gap_index + 1;
         lookahead < static_cast<int>(gaps.size()) &&
         lookahead <= gap_index + maximum_extension_count; ++lookahead) {
      const int next_offset = gaps[lookahead];
      const SpliceAnchor& next = anchors[begin + next_offset];
      if (query_end(next) - query_start > maximum_extension_length ||
          reference_end(next) - reference_start > maximum_extension_length)
        break;
      const SpliceAnchor& next_previous = anchors[begin + next_offset - 1];
      gap = static_cast<std::int64_t>(query_end(next)) -
          query_end(next_previous) -
          (static_cast<std::int64_t>(reference_end(next)) -
           reference_end(next_previous));
      if (gap > 0) inserted += gap;
      else deleted -= gap;
      const std::int64_t difference =
          inserted + deleted - std::llabs(inserted - deleted);
      if (difference > maximum_difference) {
        maximum_difference = static_cast<int>(difference);
        maximum_difference_index = lookahead;
      }
    }
    if (maximum_difference > difference_threshold &&
        maximum_difference > best) {
      best = maximum_difference;
      best_start = gap_index;
      best_end = maximum_difference_index;
    }
  }
}

void filter_alternating_gap_clusters(std::vector<SpliceAnchor>& anchors,
                                     int begin, int count,
                                     int minimum_gap,
                                     int maximum_extension) {
  const std::vector<int> gaps =
      collect_long_gaps(anchors, begin, count, minimum_gap);
  if (gaps.empty()) return;
  for (int gap_index = 0; gap_index < static_cast<int>(gaps.size());) {
    const int offset = gaps[gap_index];
    const SpliceAnchor& current = anchors[begin + offset];
    const SpliceAnchor& previous = anchors[begin + offset - 1];
    std::int64_t gap_one =
        static_cast<std::int64_t>(query_end(current)) - query_end(previous) -
        (static_cast<std::int64_t>(reference_end(current)) -
         reference_end(previous));
    gap_one = std::llabs(gap_one);
    int previous_reference_end = reference_end(current);
    int previous_query_end = query_end(current);
    int lookahead = gap_index + 1;
    for (; lookahead < static_cast<int>(gaps.size()); ++lookahead) {
      const int next_offset = gaps[lookahead];
      const SpliceAnchor& next = anchors[begin + next_offset];
      if (query_end(next) - previous_query_end > maximum_extension ||
          reference_end(next) - previous_reference_end > maximum_extension)
        break;
      const SpliceAnchor& next_previous = anchors[begin + next_offset - 1];
      std::int64_t gap_two =
          static_cast<std::int64_t>(query_end(next)) -
          query_end(next_previous) -
          (static_cast<std::int64_t>(reference_end(next)) -
           reference_end(next_previous));
      const std::int64_t reference_limit =
          static_cast<std::int64_t>(reference_end(next_previous)) +
          next_previous.span;
      const std::int64_t query_limit =
          static_cast<std::int64_t>(query_end(next_previous)) +
          next_previous.span;
      const std::int64_t shared = std::min(
          reference_limit - previous_reference_end,
          query_limit - previous_query_end);
      gap_two = std::llabs(gap_two);
      if (shared > gap_one + gap_two) break;
      previous_reference_end = reference_end(next);
      previous_query_end = query_end(next);
      gap_one = gap_two;
    }
    if (lookahead > gap_index + 1) {
      const int end = gaps[lookahead - 1];
      for (int candidate = gaps[gap_index]; candidate < end; ++candidate)
        add_flag(anchors[begin + candidate], kSpliceAnchorIgnore);
      add_flag(anchors[begin + end], kSpliceAnchorLongJoin);
    }
    gap_index = lookahead;
  }
}

struct Envelope {
  int reference_begin = 0;
  int query_begin = 0;
  int reference_end = 0;
  int query_end = 0;
  int left_neighbor = -1;
  int right_neighbor = -1;
};

Envelope make_terminal_envelope(const SpliceControllerRequest& request,
                                int first_reference_pivot,
                                int first_query_pivot,
                                int last_reference_pivot,
                                int last_query_pivot) noexcept {
  const std::vector<SpliceAnchor>& anchors = *request.working_anchors;
  const int query_length =
      static_cast<int>(request.query_forward->size());
  const int reference_length = static_cast<int>(request.reference->size());
  const SpliceAnchor& original_first = anchors[request.selected_begin];
  const SpliceAnchor& original_last =
      anchors[request.selected_begin + request.selected_count - 1];
  Envelope envelope;
  envelope.reference_begin = original_first.reference_begin;
  envelope.query_begin = original_first.query_begin;
  int reference_constraint = 0;
  int query_constraint = 0;
  int qualifying = 0;
  for (int index = request.selected_begin - 1; index >= 0; --index) {
    const int candidate_reference = anchors[index].reference_begin;
    const int candidate_query = anchors[index].query_begin;
    if (candidate_reference < envelope.reference_begin &&
        candidate_query < envelope.query_begin &&
        ++qualifying > request.options.minimum_anchor_count) {
      const int distance = std::max(
          envelope.reference_begin - candidate_reference,
          envelope.query_begin - candidate_query);
      reference_constraint = envelope.reference_begin - distance;
      query_constraint = envelope.query_begin - distance;
      if (reference_constraint < 0) reference_constraint = 0;
      envelope.left_neighbor = index;
      break;
    }
  }
  if (first_query_pivot > 0 && first_reference_pivot > 0) {
    int length = std::min(first_query_pivot,
                          request.options.maximum_gap);
    query_constraint = std::max(query_constraint,
                                first_query_pivot - length);
    envelope.query_begin = std::min(envelope.query_begin,
                                    query_constraint);
    if (length * request.options.match > request.options.gap_open) {
      length += (length * request.options.match -
                 request.options.gap_open) / request.options.gap_extend;
    }
    length = std::min(length, request.options.maximum_gap);
    length = std::min(length, first_reference_pivot);
    reference_constraint = std::max(reference_constraint,
                                    first_reference_pivot - length);
    envelope.reference_begin = std::min(envelope.reference_begin,
                                        reference_constraint);
    envelope.reference_begin = std::min(envelope.reference_begin,
                                        first_reference_pivot);
  } else {
    envelope.reference_begin = first_reference_pivot;
    envelope.query_begin = first_query_pivot;
  }

  envelope.reference_end = reference_end(original_last) + 1;
  envelope.query_end = query_end(original_last) + 1;
  int reference_limit = reference_length;
  int query_limit = query_length;
  qualifying = 0;
  const int after_selected = request.selected_begin + request.selected_count;
  for (int index = after_selected;
       index < static_cast<int>(anchors.size()); ++index) {
    const int candidate_reference = reference_end(anchors[index]) + 1;
    const int candidate_query = query_end(anchors[index]) + 1;
    if (candidate_reference > envelope.reference_end &&
        candidate_query > envelope.query_end &&
        ++qualifying > request.options.minimum_anchor_count) {
      const int distance = std::max(
          candidate_reference - envelope.reference_end,
          candidate_query - envelope.query_end);
      reference_limit = static_cast<int>(std::min<std::int64_t>(
          reference_length,
          static_cast<std::int64_t>(envelope.reference_end) + distance));
      query_limit = static_cast<int>(std::min<std::int64_t>(
          query_length,
          static_cast<std::int64_t>(envelope.query_end) + distance));
      envelope.right_neighbor = index;
      break;
    }
  }
  if (last_query_pivot < query_length &&
      last_reference_pivot < reference_length) {
    int length = std::min(query_length - last_query_pivot,
                          request.options.maximum_gap);
    query_limit = std::min(query_limit, last_query_pivot + length);
    envelope.query_end = std::max(envelope.query_end, query_limit);
    if (length * request.options.match > request.options.gap_open) {
      length += (length * request.options.match -
                 request.options.gap_open) / request.options.gap_extend;
    }
    length = std::min(length, request.options.maximum_gap);
    length = std::min(length, reference_length - last_reference_pivot);
    reference_limit = std::min(reference_limit,
                               last_reference_pivot + length);
    envelope.reference_end = std::max(envelope.reference_end,
                                      reference_limit);
  } else {
    envelope.reference_end = last_reference_pivot;
    envelope.query_end = last_query_pivot;
  }
  return envelope;
}

// The orientation-blind half of a segment, computed once per selected range and
// replayed by the second transcript hypothesis; the caller applies `prepared` to its
// result on both passes. Not noexcept: its vectors may allocate, and
// realize_splice_segment turns bad_alloc into a refusal.
void prepare_segment(const SpliceControllerRequest& request,
                     const std::array<std::int8_t, 25>& matrix,
                     std::vector<std::uint8_t>& flags_before,
                     SplicePreparedSegment& prepared) {
  std::vector<SpliceAnchor>& anchors = *request.working_anchors;
  const std::vector<std::uint8_t>& mapping_query = request.mapping_reverse
      ? *request.query_reverse : *request.query_forward;
  prepared.selected_begin = request.selected_begin;
  prepared.selected_count = request.selected_count;
  prepared.refused = false;
  prepared.refusal = RnaControllerRefusal::None;
  prepared.wrote_filtered = false;
  prepared.wrote_envelope = false;
  prepared.anchor_flags.clear();

  int filtered_begin = 0;
  int filtered_count = 0;
  RnaControllerRefusal refusal = RnaControllerRefusal::None;
  if (!filter_weak_terminal_anchors(request, matrix, filtered_begin,
                                    filtered_count, refusal)) {
    prepared.refused = true;
    prepared.refusal = refusal;
    return;
  }
  flags_before.resize(static_cast<std::size_t>(filtered_count));
  for (int offset = 0; offset < filtered_count; ++offset)
    flags_before[static_cast<std::size_t>(offset)] =
        anchors[filtered_begin + offset].flags;
  filter_balanced_gap_clusters(
      anchors, filtered_begin, filtered_count, 10, 40,
      request.options.maximum_gap >> 1, 10);
  filter_alternating_gap_clusters(
      anchors, filtered_begin, filtered_count, 30,
      request.options.maximum_gap >> 1);
  for (int offset = 0; offset < filtered_count; ++offset) {
    const std::uint8_t before = flags_before[static_cast<std::size_t>(offset)];
    const std::uint8_t after = anchors[filtered_begin + offset].flags;
    if (after != before)
      prepared.anchor_flags.push_back(
          {static_cast<std::uint32_t>(filtered_begin + offset),
           static_cast<std::uint8_t>(after & ~before)});
  }

  const SpliceAnchor& first = anchors[filtered_begin];
  const SpliceAnchor& last = anchors[filtered_begin + filtered_count - 1];
  const int reference_cursor =
      reference_end(first) - (request.options.k >> 1);
  const int query_cursor = query_end(first) - (request.options.k >> 1);
  const int last_reference_pivot =
      reference_end(last) - (request.options.k >> 1);
  const int last_query_pivot = query_end(last) - (request.options.k >> 1);
  if (!valid_range(query_cursor, last_query_pivot, mapping_query.size()) ||
      !valid_range(reference_cursor, last_reference_pivot,
                   request.reference->size())) {
    prepared.refused = true;
    prepared.refusal = RnaControllerRefusal::InvalidAnchor;
    return;
  }
  prepared.wrote_filtered = true;
  prepared.filtered_begin = filtered_begin;
  prepared.filtered_count = filtered_count;
  prepared.first_reference_pivot = reference_cursor;
  prepared.first_query_pivot = query_cursor;
  prepared.last_reference_pivot = last_reference_pivot;
  prepared.last_query_pivot = last_query_pivot;

  const Envelope envelope = make_terminal_envelope(
      request, reference_cursor, query_cursor, last_reference_pivot,
      last_query_pivot);
  prepared.wrote_envelope = true;
  prepared.envelope_reference_begin = envelope.reference_begin;
  prepared.envelope_query_begin = envelope.query_begin;
  prepared.envelope_reference_end = envelope.reference_end;
  prepared.envelope_query_end = envelope.query_end;
  prepared.left_neighbor = envelope.left_neighbor;
  prepared.right_neighbor = envelope.right_neighbor;
  if (!valid_range(envelope.query_begin, envelope.query_end,
                   mapping_query.size()) ||
      !valid_range(envelope.reference_begin, envelope.reference_end,
                   request.reference->size()) ||
      envelope.reference_end <= envelope.reference_begin) {
    prepared.refused = true;
    prepared.refusal = RnaControllerRefusal::InvalidSequenceSlice;
  }
}

void append_cigar(std::vector<std::uint32_t>& destination,
                  const std::vector<std::uint32_t>& source) {
  if (source.empty()) return;
  const std::uint32_t destination_length = destination.empty()
      ? 0 : destination.back() >> 4;
  const std::uint32_t source_length = source.front() >> 4;
  std::size_t source_begin = 0;
  if (!destination.empty() &&
      (destination.back() & 0xfu) == (source.front() & 0xfu) &&
      source_length <= kMaximumCigarRun - destination_length) {
    destination.back() += source_length << 4;
    source_begin = 1;
  }
  for (std::size_t index = source_begin; index < source.size(); ++index)
    destination.push_back(source[index]);
}

bool replay_zdrop(const SpliceControllerOptions& options,
                  const std::array<std::int8_t, 25>& matrix,
                  const std::uint8_t* query, int query_length,
                  const std::uint8_t* reference, int reference_length,
                  const std::vector<std::uint32_t>& cigar,
                  bool& signalled) noexcept {
  std::int64_t score = 0;
  std::int64_t maximum = std::numeric_limits<std::int32_t>::min();
  int maximum_reference = -1;
  int maximum_query = -1;
  int reference_offset = 0;
  int query_offset = 0;
  std::int64_t maximum_drop = 0;
  for (std::uint32_t encoded : cigar) {
    const int operation = static_cast<int>(encoded & 0xfu);
    const int length = static_cast<int>(encoded >> 4);
    if (length <= 0) return false;
    if (operation == kCigarMatch) {
      if (!valid_span(query_offset, length,
                      static_cast<std::size_t>(query_length)) ||
          !valid_span(reference_offset, length,
                      static_cast<std::size_t>(reference_length)))
        return false;
      for (int base = 0; base < length; ++base) {
        score += matrix[reference[reference_offset + base] * 5 +
                        query[query_offset + base]];
        const int i = reference_offset + base;
        const int j = query_offset + base;
        if (score < maximum) {
          const int reference_distance = i - maximum_reference;
          const int query_distance = j - maximum_query;
          const int diagonal = std::abs(reference_distance - query_distance);
          const std::int64_t drop = maximum - score -
              static_cast<std::int64_t>(diagonal) * options.gap_extend;
          if (drop > maximum_drop) maximum_drop = drop;
        } else {
          maximum = score;
          maximum_reference = i;
          maximum_query = j;
        }
      }
      reference_offset += length;
      query_offset += length;
    } else if (operation == kCigarInsertion ||
               operation == kCigarDeletion ||
               operation == kCigarReferenceSkip) {
      score -= options.gap_open +
               static_cast<std::int64_t>(options.gap_extend) * length;
      if (operation == kCigarInsertion) {
        if (!valid_span(query_offset, length,
                        static_cast<std::size_t>(query_length))) return false;
        query_offset += length;
      } else {
        if (!valid_span(reference_offset, length,
                        static_cast<std::size_t>(reference_length))) return false;
        reference_offset += length;
      }
      if (score < maximum) {
        const int reference_distance =
            reference_offset - maximum_reference;
        const int query_distance = query_offset - maximum_query;
        const int diagonal = std::abs(reference_distance - query_distance);
        const std::int64_t drop = maximum - score -
            static_cast<std::int64_t>(diagonal) * options.gap_extend;
        if (drop > maximum_drop) maximum_drop = drop;
      } else {
        maximum = score;
        maximum_reference = reference_offset;
        maximum_query = query_offset;
      }
    } else {
      return false;
    }
  }
  if (query_offset != query_length || reference_offset != reference_length)
    return false;
  signalled = maximum_drop > options.zdrop;
  return true;
}

// minimap2's mg_log2 approximation. The input is 1 + a positive gap length, so >= 2.
float approximate_log2(float value) noexcept {
  union Bits {
    float floating;
    std::uint32_t integer;
  } bits{value};
  float result = static_cast<float>((bits.integer >> 23) & 255u) - 128.0f;
  bits.integer &= ~(255u << 23);
  bits.integer += 127u << 23;
  result += (-0.34484843f * bits.floating + 2.02466578f) *
      bits.floating - 0.67487759f;
  return result;
}

bool consumes_reference(int operation) noexcept {
  return operation == kCigarMatch || operation == kCigarDeletion ||
         operation == kCigarReferenceSkip;
}

bool consumes_query(int operation) noexcept {
  return operation == kCigarMatch || operation == kCigarInsertion;
}

// The reference span [begin, begin + span) lies wholly inside the anchor. An empty
// terminal exon (a record that would begin or end with `N`) is inside every anchor.
bool exon_inside_anchor(int begin, int span,
                        const SpliceAnchor& anchor) noexcept {
  if (span <= 0) return true;
  return begin >= anchor.reference_begin &&
         begin + span - 1 <= reference_end(anchor);
}

// A junction at a record end is evidence only when the exon outside it is. Packets are
// cut at the middle of the terminal seed (anchor_end - k/2, as minimap2's
// mm_adjust_minier), so the first or last internal packet can open a junction inside
// the very seed the outer exon rests on, leaving a pivot-cut fragment with an
// unsupported junction beside it. minimap2's record ends are extension DPs read back at
// their score maximum, which never lies past a junction. So clip the record past a
// terminal junction whose outer exon lies wholly inside the terminal anchor, and settle
// the new end on an aligned base: a record never begins or ends with N.
bool clip_unsupported_terminal_junctions(
    std::vector<std::uint32_t>& cigar, const SpliceAnchor& first_anchor,
    const SpliceAnchor& last_anchor, int& query_begin, int& query_end,
    int& reference_begin, int& reference_limit) noexcept {
  bool clipped_front = false;
  for (std::size_t guard = 0; guard < cigar.size(); ++guard) {
    std::size_t junction = 0;
    bool found = false;
    int exon_reference = 0;
    int exon_query = 0;
    for (std::size_t index = 0; index < cigar.size(); ++index) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int length = static_cast<int>(cigar[index] >> 4);
      if (operation == kCigarReferenceSkip) {
        junction = index;
        found = true;
        break;
      }
      if (consumes_reference(operation)) exon_reference += length;
      if (consumes_query(operation)) exon_query += length;
    }
    if (!found ||
        !exon_inside_anchor(reference_begin, exon_reference, first_anchor))
      break;
    reference_begin += exon_reference + static_cast<int>(cigar[junction] >> 4);
    query_begin += exon_query;
    cigar.erase(cigar.begin(),
                cigar.begin() + static_cast<std::ptrdiff_t>(junction) + 1);
    clipped_front = true;
  }
  // The end a clip leaves behind settles on an aligned base: a record that
  // starts or stops on a gap places reference or query bases against nothing.
  while (clipped_front && !cigar.empty() &&
         static_cast<int>(cigar.front() & 0xfu) != kCigarMatch) {
    const int operation = static_cast<int>(cigar.front() & 0xfu);
    const int length = static_cast<int>(cigar.front() >> 4);
    if (consumes_reference(operation)) reference_begin += length;
    if (consumes_query(operation)) query_begin += length;
    cigar.erase(cigar.begin());
  }

  bool clipped_back = false;
  for (std::size_t guard = 0; guard < cigar.size(); ++guard) {
    std::size_t junction = 0;
    bool found = false;
    int exon_reference = 0;
    int exon_query = 0;
    for (std::size_t index = cigar.size(); index-- > 0;) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int length = static_cast<int>(cigar[index] >> 4);
      if (operation == kCigarReferenceSkip) {
        junction = index;
        found = true;
        break;
      }
      if (consumes_reference(operation)) exon_reference += length;
      if (consumes_query(operation)) exon_query += length;
    }
    if (!found ||
        !exon_inside_anchor(reference_limit - exon_reference, exon_reference,
                            last_anchor))
      break;
    reference_limit -= exon_reference + static_cast<int>(cigar[junction] >> 4);
    query_end -= exon_query;
    cigar.erase(cigar.begin() + static_cast<std::ptrdiff_t>(junction),
                cigar.end());
    clipped_back = true;
  }
  while (clipped_back && !cigar.empty() &&
         static_cast<int>(cigar.back() & 0xfu) != kCigarMatch) {
    const int operation = static_cast<int>(cigar.back() & 0xfu);
    const int length = static_cast<int>(cigar.back() >> 4);
    if (consumes_reference(operation)) reference_limit -= length;
    if (consumes_query(operation)) query_end -= length;
    cigar.pop_back();
  }
  return !cigar.empty();
}

bool normalize_cigar(std::vector<std::uint32_t>& cigar,
                     const std::uint8_t* query,
                     const std::uint8_t* reference,
                     int query_span, int reference_span,
                     bool mapping_reverse,
                     int& projected_query_begin,
                     int& projected_query_end,
                     int& reference_begin,
                     int& query_shift,
                     int& reference_shift) {
  query_shift = reference_shift = 0;
  if (cigar.empty()) return false;
  // As in minimap2, a single-operation CIGAR skips normalization, including the
  // leading-gap rule.
  if (cigar.size() <= 1) return true;
  {
    int reference_offset = 0;
    int query_offset = 0;
    bool shrink = false;
    for (std::size_t index = 0; index < cigar.size(); ++index) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int length = static_cast<int>(cigar[index] >> 4);
      if (length == 0) shrink = true;
      if (operation == kCigarMatch) {
        if (!valid_span(query_offset, length,
                        static_cast<std::size_t>(query_span)) ||
            !valid_span(reference_offset, length,
                        static_cast<std::size_t>(reference_span)))
          return false;
        reference_offset += length;
        query_offset += length;
      } else if (operation == kCigarInsertion ||
                 operation == kCigarDeletion) {
        if (index > 0 && index + 1 < cigar.size() &&
            (cigar[index - 1] & 0xfu) == kCigarMatch &&
            (cigar[index + 1] & 0xfu) == kCigarMatch) {
          const int previous_length =
              static_cast<int>(cigar[index - 1] >> 4);
          int shift = 0;
          if (operation == kCigarInsertion) {
            if (!valid_span(query_offset, length,
                            static_cast<std::size_t>(query_span)))
              return false;
            for (; shift < previous_length; ++shift) {
              if (query[query_offset - 1 - shift] !=
                  query[query_offset + length - 1 - shift])
                break;
            }
          } else {
            if (!valid_span(reference_offset, length,
                            static_cast<std::size_t>(reference_span)))
              return false;
            for (; shift < previous_length; ++shift) {
              if (reference[reference_offset - 1 - shift] !=
                  reference[reference_offset + length - 1 - shift])
                break;
            }
          }
          const std::uint32_t next_length = cigar[index + 1] >> 4;
          if (shift > 0 && static_cast<std::uint32_t>(shift) <=
                  kMaximumCigarRun - next_length) {
            cigar[index - 1] -= static_cast<std::uint32_t>(shift) << 4;
            cigar[index + 1] += static_cast<std::uint32_t>(shift) << 4;
            query_offset -= shift;
            reference_offset -= shift;
          }
          if (shift == previous_length && static_cast<std::uint32_t>(shift) <=
                  kMaximumCigarRun - next_length) shrink = true;
        }
        if (operation == kCigarInsertion) {
          if (!valid_span(query_offset, length,
                          static_cast<std::size_t>(query_span))) return false;
          query_offset += length;
        } else {
          if (!valid_span(reference_offset, length,
                          static_cast<std::size_t>(reference_span))) return false;
          reference_offset += length;
        }
      } else if (operation == kCigarReferenceSkip) {
        if (!valid_span(reference_offset, length,
                        static_cast<std::size_t>(reference_span))) return false;
        reference_offset += length;
      } else {
        return false;
      }
      if (query_offset > query_span || reference_offset > reference_span)
        return false;
    }
    if (query_offset != query_span || reference_offset != reference_span)
      return false;

    for (std::size_t index = 0; index + 2 < cigar.size();) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int next_operation = static_cast<int>(cigar[index + 1] & 0xfu);
      if (operation > 0 && operation + next_operation == 3) {
        std::array<std::uint32_t, 3> totals{};
        std::size_t limit = index;
        for (; limit < cigar.size(); ++limit) {
          const int gap_operation = static_cast<int>(cigar[limit] & 0xfu);
          const std::uint32_t length = cigar[limit] >> 4;
          if (gap_operation == kCigarInsertion ||
              gap_operation == kCigarDeletion || length == 0) {
            if (gap_operation > kCigarDeletion) return false;
            totals[static_cast<std::size_t>(gap_operation)] += length;
          } else {
            break;
          }
        }
        if (totals[kCigarInsertion] > 0 && totals[kCigarDeletion] > 0 &&
            totals[kCigarInsertion] <= kMaximumCigarRun &&
            totals[kCigarDeletion] <= kMaximumCigarRun &&
            limit - index > 2) {
          cigar[index] = totals[kCigarInsertion] << 4 | kCigarInsertion;
          cigar[index + 1] = totals[kCigarDeletion] << 4 | kCigarDeletion;
          for (std::size_t zero = index + 2; zero < limit; ++zero)
            cigar[zero] &= 0xfu;
          shrink = true;
        }
        index = limit + 1;
      } else {
        ++index;
      }
    }
    if (shrink) {
      cigar.erase(std::remove_if(cigar.begin(), cigar.end(),
                                 [](std::uint32_t encoded) {
                                   return (encoded >> 4) == 0;
                                 }),
                  cigar.end());
      std::vector<std::uint32_t> merged;
      merged.reserve(cigar.size());
      for (std::uint32_t encoded : cigar) {
        const std::uint32_t encoded_length = encoded >> 4;
        const std::uint32_t merged_length = merged.empty()
            ? 0 : merged.back() >> 4;
        if (!merged.empty() &&
            (merged.back() & 0xfu) == (encoded & 0xfu) &&
            encoded_length <= kMaximumCigarRun - merged_length)
          merged.back() += encoded_length << 4;
        else
          merged.push_back(encoded);
      }
      cigar.swap(merged);
    }
  }
  if (cigar.empty()) return false;
  const int leading_operation = static_cast<int>(cigar.front() & 0xfu);
  if (leading_operation == kCigarInsertion ||
      leading_operation == kCigarDeletion) {
    const int length = static_cast<int>(cigar.front() >> 4);
    if (leading_operation == kCigarInsertion) {
      if (mapping_reverse) projected_query_end -= length;
      else projected_query_begin += length;
      query_shift = length;
    } else {
      reference_begin += length;
      reference_shift = length;
    }
    cigar.erase(cigar.begin());
  }
  return !cigar.empty();
}

bool update_alignment_summary(
    SpliceSegmentResult& result,
    const std::vector<std::uint8_t>& mapping_query,
    const std::vector<std::uint8_t>& reference,
    int oriented_query_begin, int oriented_query_end,
    const std::array<std::int8_t, 25>& matrix,
    const SpliceControllerOptions& options,
    bool mapping_reverse) {
  if (!valid_range(oriented_query_begin, oriented_query_end,
                   mapping_query.size()) ||
      !valid_range(result.reference_begin, result.reference_end,
                   reference.size()))
    return false;
  int query_shift = 0;
  int reference_shift = 0;
  if (!normalize_cigar(
          result.cigar, mapping_query.data() + oriented_query_begin,
          reference.data() + result.reference_begin,
          oriented_query_end - oriented_query_begin,
          result.reference_end - result.reference_begin,
          mapping_reverse, result.query_begin, result.query_end,
          result.reference_begin, query_shift, reference_shift))
    return false;

  const std::uint8_t* query =
      mapping_query.data() + oriented_query_begin + query_shift;
  const std::uint8_t* target =
      reference.data() + result.reference_begin;
  const int query_span = result.query_end - result.query_begin;
  const int reference_span = result.reference_end - result.reference_begin;
  int query_offset = 0;
  int reference_offset = 0;
  double score = 0.0;
  double maximum = 0.0;
  result.matches = 0;
  result.block_length = 0;
  result.ambiguous_bases = 0;
  result.is_spliced = false;
  for (std::uint32_t encoded : result.cigar) {
    const int operation = static_cast<int>(encoded & 0xfu);
    const int length = static_cast<int>(encoded >> 4);
    if (length <= 0) return false;
    if (operation == kCigarMatch) {
      if (!valid_span(query_offset, length,
                      static_cast<std::size_t>(query_span)) ||
          !valid_span(reference_offset, length,
                      static_cast<std::size_t>(reference_span)))
        return false;
      int ambiguous = 0;
      int differences = 0;
      for (int base = 0; base < length; ++base) {
        const int query_base = query[query_offset + base];
        const int reference_base = target[reference_offset + base];
        if (reference_base > 3 || query_base > 3) ++ambiguous;
        else if (reference_base != query_base) ++differences;
        score += matrix[reference_base * 5 + query_base];
        if (score < 0.0) score = 0.0;
        else maximum = std::max(maximum, score);
      }
      result.block_length += length - ambiguous;
      result.matches += length - ambiguous - differences;
      result.ambiguous_bases += ambiguous;
      reference_offset += length;
      query_offset += length;
    } else if (operation == kCigarInsertion) {
      if (!valid_span(query_offset, length,
                      static_cast<std::size_t>(query_span)))
        return false;
      int ambiguous = 0;
      for (int base = 0; base < length; ++base)
        if (query[query_offset + base] > 3) ++ambiguous;
      result.block_length += length - ambiguous;
      result.ambiguous_bases += ambiguous;
      score -= options.gap_open + static_cast<double>(options.gap_extend) *
          approximate_log2(static_cast<float>(1 + length));
      if (score < 0.0) score = 0.0;
      query_offset += length;
    } else if (operation == kCigarDeletion) {
      if (!valid_span(reference_offset, length,
                      static_cast<std::size_t>(reference_span)))
        return false;
      int ambiguous = 0;
      for (int base = 0; base < length; ++base)
        if (target[reference_offset + base] > 3) ++ambiguous;
      result.block_length += length - ambiguous;
      result.ambiguous_bases += ambiguous;
      score -= options.gap_open + static_cast<double>(options.gap_extend) *
          approximate_log2(static_cast<float>(1 + length));
      if (score < 0.0) score = 0.0;
      reference_offset += length;
    } else if (operation == kCigarReferenceSkip) {
      if (!valid_span(reference_offset, length,
                      static_cast<std::size_t>(reference_span)))
        return false;
      result.is_spliced = true;
      reference_offset += length;
    } else {
      return false;
    }
  }
  if (query_offset != query_span || reference_offset != reference_span)
    return false;
  result.dp_maximum = result.dp_maximum_before_strand =
      static_cast<int>(maximum + .499);
  return true;
}

struct KernelRunner {
  KernelRunner(const SpliceControllerRequest& input,
               SpliceSegmentResult& output,
               dp::SpliceKernelWorkspace& kernel,
               const std::array<std::int8_t, 25>& score_matrix)
      : request(input), result(output), workspace(kernel),
        matrix(score_matrix) {}

  bool run(SplicePacketRole role,
           int query_begin, int query_end,
           int reference_begin, int reference_end,
           const std::uint8_t* query,
           const std::uint8_t* reference,
           const std::uint8_t* junction,
           int flags, int zdrop, int end_bonus,
           bool approximate, bool exact_retry,
           dp::SpliceKernelResult& aligned) {
    const int query_length = query_end - query_begin;
    const int reference_length = reference_end - reference_begin;
    if (query_length <= 0 || reference_length <= 0 || !query || !reference) {
      refuse(result, RnaControllerRefusal::InvalidSequenceSlice);
      return false;
    }
    if (!valid_nt4(query, query_length) ||
        !valid_nt4(reference, reference_length)) {
      refuse(result, RnaControllerRefusal::InvalidRequest);
      return false;
    }
    dp::SpliceKernelRequest kernel_request;
    kernel_request.query = query;
    kernel_request.query_length = query_length;
    kernel_request.reference = reference;
    kernel_request.reference_length = reference_length;
    kernel_request.junction = junction;
    kernel_request.score_matrix = matrix.data();
    kernel_request.zdrop = zdrop;
    kernel_request.end_bonus = end_bonus;
    kernel_request.flags = flags;
    kernel_request.scoring = kernel_scoring(request.options);
    result.metrics.observe_dp(
        static_cast<std::uint64_t>(query_length),
        static_cast<std::uint64_t>(reference_length),
        static_cast<int>(role));
    if (approximate) ++result.metrics.approximate_passes;
    if (exact_retry) ++result.metrics.exact_retries;
    workspace.align(kernel_request, aligned);
    if (!aligned.succeeded) {
      refuse(result, kernel_refusal(aligned.failure));
      return false;
    }
    return true;
  }

  const SpliceControllerRequest& request;
  SpliceSegmentResult& result;
  dp::SpliceKernelWorkspace& workspace;
  const std::array<std::int8_t, 25>& matrix;
};

}  // namespace

SpliceControllerOptions ordinary_long_read_splice_options() noexcept {
  return {};
}

bool splice_controller_options_supported(
    const SpliceControllerOptions& options) noexcept {
  return options_are_supported(options);
}

SpliceSegmentResult realize_splice_segment(
    const SpliceControllerRequest& request) noexcept {
  SpliceSegmentResult result;
  // One per worker thread when the caller owns scratch, otherwise one per segment.
  std::optional<SpliceRealizationScratch> owned_scratch;
  if (!request.scratch) owned_scratch.emplace();
  SpliceRealizationScratch& scratch =
      request.scratch ? *request.scratch : *owned_scratch;
  try {
    if (!validate_request(request, result)) return result;
    std::vector<SpliceAnchor>& anchors = *request.working_anchors;
    const std::vector<std::uint8_t>& mapping_query = request.mapping_reverse
        ? *request.query_reverse : *request.query_forward;
    const int query_length = static_cast<int>(mapping_query.size());
    const int reference_length = static_cast<int>(request.reference->size());
    const auto scoring = kernel_scoring(request.options);
    const auto matrix = dp::make_splice_score_matrix(
        scoring, request.options.transition);

    // Prepare the segment once per selected range: both hypotheses see the
    // same anchors, so only base_flags below is orientation-dependent.
    const SplicePreparedSegment* found = nullptr;
    for (const SplicePreparedSegment& entry : scratch.prepared) {
      if (entry.selected_begin == request.selected_begin &&
          entry.selected_count == request.selected_count) {
        found = &entry;
        break;
      }
    }
    if (!found) {
      SplicePreparedSegment computed;
      prepare_segment(request, matrix, scratch.flags_before, computed);
      scratch.prepared.push_back(std::move(computed));
      found = &scratch.prepared.back();
    } else {
      // Replay the cluster filters' bits onto this hypothesis' own copy.
      for (const SpliceAnchorFlagMutation& mutation : found->anchor_flags)
        anchors[mutation.index].flags |= mutation.added;
    }
    const SplicePreparedSegment& prepared = *found;
    if (prepared.wrote_filtered) {
      result.filtered_begin = prepared.filtered_begin;
      result.filtered_count = prepared.filtered_count;
      result.first_reference_pivot = prepared.first_reference_pivot;
      result.first_query_pivot = prepared.first_query_pivot;
      result.last_reference_pivot = prepared.last_reference_pivot;
      result.last_query_pivot = prepared.last_query_pivot;
    }
    if (prepared.wrote_envelope) {
      result.envelope_reference_begin = prepared.envelope_reference_begin;
      result.envelope_query_begin = prepared.envelope_query_begin;
      result.envelope_reference_end = prepared.envelope_reference_end;
      result.envelope_query_end = prepared.envelope_query_end;
      result.left_neighbor_index = prepared.left_neighbor;
      result.right_neighbor_index = prepared.right_neighbor;
    }
    if (prepared.refused) {
      refuse(result, prepared.refusal);
      return result;
    }
    const int filtered_begin = prepared.filtered_begin;
    const int filtered_count = prepared.filtered_count;
    int reference_cursor = prepared.first_reference_pivot;
    int query_cursor = prepared.first_query_pivot;
    const int last_reference_pivot = prepared.last_reference_pivot;
    const int last_query_pivot = prepared.last_query_pivot;
    Envelope envelope;
    envelope.reference_begin = prepared.envelope_reference_begin;
    envelope.query_begin = prepared.envelope_query_begin;
    envelope.reference_end = prepared.envelope_reference_end;
    envelope.query_end = prepared.envelope_query_end;
    envelope.left_neighbor = prepared.left_neighbor;
    envelope.right_neighbor = prepared.right_neighbor;

    // kSpliceFlank is always set, as in minimap2's splice presets, but only ksw2's
    // ungraded splice branch reads it, and the graded model (kSpliceComplex, set in
    // dp/splice_kernel.cpp) switches that branch off; motif credit comes from the graded
    // table in ksw2_exts2_sse.c. The direction flag set below from the transcript
    // orientation is what the graded table reads.
    int base_flags = dp::kSpliceFlank;
    if (request.transcript_orientation == TranscriptOrientation::Forward)
      base_flags |= request.mapping_reverse
          ? dp::kSpliceReverse : dp::kSpliceForward;
    else if (request.transcript_orientation == TranscriptOrientation::Reverse)
      base_flags |= request.mapping_reverse
          ? dp::kSpliceForward : dp::kSpliceReverse;
    else if (request.transcript_orientation == TranscriptOrientation::None) {
      // minimap2 -un: neither direction flag, so no donor/acceptor signal is scored;
      // kSpliceFlank stays set and is inert.
    } else {
      refuse(result, RnaControllerRefusal::InvalidRequest);
      return result;
    }

    KernelRunner kernel(request, result, scratch.kernel, matrix);
    // The junction mask lives in scratch, so the annotation store stays immutable. With
    // no annotation it would be all zeroes, which the kernel reads as no annotated
    // junction, so it is built only when there is one.
    std::vector<std::uint8_t>& junction_scratch = scratch.junction;
    const bool materialize_junctions = !request.known_junctions.empty();
    int oriented_query_begin = query_cursor;
    int oriented_reference_begin = reference_cursor;
    int oriented_query_end = query_cursor;
    int oriented_reference_end = reference_cursor;

    if (query_cursor > 0 && reference_cursor > 0) {
      const int query_packet_length =
          query_cursor - envelope.query_begin;
      const int reference_packet_length =
          reference_cursor - envelope.reference_begin;
      if (query_packet_length <= 0 || reference_packet_length <= 0) {
        refuse(result, RnaControllerRefusal::InvalidSequenceSlice);
        return result;
      }
      std::vector<std::uint8_t>& reversed_query = scratch.reversed_query;
      std::vector<std::uint8_t>& reversed_reference =
          scratch.reversed_reference;
      reversed_query.assign(mapping_query.begin() + envelope.query_begin,
                            mapping_query.begin() + query_cursor);
      reversed_reference.assign(
          request.reference->begin() + envelope.reference_begin,
          request.reference->begin() + reference_cursor);
      std::reverse(reversed_query.begin(), reversed_query.end());
      std::reverse(reversed_reference.begin(), reversed_reference.end());
      const std::uint8_t* junction = nullptr;
      if (materialize_junctions) {
        const std::uint64_t annotated_endpoint_bytes =
            materialize_known_junction_mask(
                request.known_junctions, envelope.reference_begin,
                reference_cursor, junction_scratch);
        if (annotated_endpoint_bytes) {
          ++result.metrics.annotation_packets_masked;
          result.metrics.annotation_endpoint_bytes_set +=
              annotated_endpoint_bytes;
        }
        // Reference and mask undergo the same left-packet reversal.
        std::reverse(junction_scratch.begin(), junction_scratch.end());
        junction = junction_scratch.data();
      }
      dp::SpliceKernelResult& aligned = scratch.aligned;
      if (!kernel.run(
              SplicePacketRole::Left,
              envelope.query_begin, query_cursor,
              envelope.reference_begin, reference_cursor,
              reversed_query.data(), reversed_reference.data(), junction,
              base_flags | dp::kSpliceExtensionOnly |
                  dp::kSpliceRightAlign | dp::kSpliceReverseCigar,
              request.options.zdrop, request.options.end_bonus,
              false, false, aligned))
        return result;
      append_cigar(result.cigar, aligned.cigar);
      if (!aligned.cigar.empty()) result.dp_score += aligned.maximum;
      const int reference_extension = aligned.reached_end
          ? aligned.query_end_reference + 1
          : aligned.maximum_reference + 1;
      const int query_extension = aligned.reached_end
          ? query_packet_length : aligned.maximum_query + 1;
      oriented_reference_begin = reference_cursor - reference_extension;
      oriented_query_begin = query_cursor - query_extension;
      if (oriented_reference_begin < 0 || oriented_query_begin < 0) {
        refuse(result, RnaControllerRefusal::KernelFailure);
        return result;
      }
    }

    bool dropped = false;
    for (int offset = 1; offset < filtered_count; ++offset) {
      SpliceAnchor& anchor = anchors[filtered_begin + offset];
      if ((has_flag(anchor, kSpliceAnchorIgnore) ||
           has_flag(anchor, kSpliceAnchorTandem)) &&
          offset != filtered_count - 1)
        continue;
      const int packet_reference_end =
          reference_end(anchor) - (request.options.k >> 1);
      const int packet_query_end =
          query_end(anchor) - (request.options.k >> 1);
      oriented_reference_end = packet_reference_end;
      oriented_query_end = packet_query_end;
      if (offset == filtered_count - 1 ||
          has_flag(anchor, kSpliceAnchorLongJoin) ||
          (packet_query_end - query_cursor >=
               request.options.minimum_packet_length &&
           packet_reference_end - reference_cursor >=
               request.options.minimum_packet_length)) {
        if (!valid_range(query_cursor, packet_query_end,
                         mapping_query.size()) ||
            !valid_range(reference_cursor, packet_reference_end,
                         request.reference->size()) ||
            packet_query_end == query_cursor ||
            packet_reference_end == reference_cursor) {
          refuse(result, RnaControllerRefusal::InvalidSequenceSlice);
          return result;
        }
        const int query_packet_length = packet_query_end - query_cursor;
        const int reference_packet_length =
            packet_reference_end - reference_cursor;
        // The packet's shape, priced before the kernel sees it. A reference run
        // that outstrips the query run by at least one minimum intron is where
        // a splice has to be found; anything shorter is colinear work.
        result.metrics.observe_internal_packet(
            static_cast<std::uint64_t>(query_packet_length),
            static_cast<std::uint64_t>(reference_packet_length),
            static_cast<std::int64_t>(reference_packet_length) -
                    query_packet_length >=
                request.options.minimum_intron,
            traceback_matrix_bytes(query_packet_length,
                                   reference_packet_length) >
                kArenaCapacityLimit);
        const std::uint8_t* junction = nullptr;
        if (materialize_junctions) {
          const std::uint64_t annotated_endpoint_bytes =
              materialize_known_junction_mask(
                  request.known_junctions, reference_cursor,
                  packet_reference_end, junction_scratch);
          if (annotated_endpoint_bytes) {
            ++result.metrics.annotation_packets_masked;
            result.metrics.annotation_endpoint_bytes_set +=
                annotated_endpoint_bytes;
          }
          junction = junction_scratch.data();
        }
        dp::SpliceKernelResult& aligned = scratch.aligned;
        if (!kernel.run(
                SplicePacketRole::Internal,
                query_cursor, packet_query_end,
                reference_cursor, packet_reference_end,
                mapping_query.data() + query_cursor,
                request.reference->data() + reference_cursor, junction,
                base_flags | dp::kSpliceApproximateMaximum,
                request.options.zdrop, -1, true, false, aligned))
          return result;

        bool zdrop_signal = false;
        if (!replay_zdrop(
                request.options, matrix,
                mapping_query.data() + query_cursor, query_packet_length,
                request.reference->data() + reference_cursor,
                reference_packet_length, aligned.cigar, zdrop_signal)) {
          refuse(result, RnaControllerRefusal::CigarReconciliation);
          return result;
        }
        if (zdrop_signal) {
          ++result.metrics.zdrop_signals;
          if (!kernel.run(
                  SplicePacketRole::Internal,
                  query_cursor, packet_query_end,
                  reference_cursor, packet_reference_end,
                  mapping_query.data() + query_cursor,
                  request.reference->data() + reference_cursor, junction,
                  base_flags,
                  request.options.zdrop, -1, false, true, aligned))
            return result;
        }
        append_cigar(result.cigar, aligned.cigar);
        if (aligned.zdropped) {
          int split_offset = offset - 1;
          for (; split_offset >= 0; --split_offset) {
            if (reference_end(anchors[filtered_begin + split_offset]) <=
                reference_cursor + aligned.maximum_reference)
              break;
          }
          if (split_offset < 0) split_offset = 0;
          dropped = true;
          result.dp_score += aligned.maximum;
          oriented_reference_end = reference_cursor +
              aligned.maximum_reference + 1;
          oriented_query_end = query_cursor + aligned.maximum_query + 1;
          if (filtered_count - (split_offset + 1) >=
              request.options.minimum_anchor_count) {
            const int first_continuation =
                filtered_begin + split_offset + 1;
            const int retained_count =
                first_continuation - request.selected_begin;
            if (retained_count <= 0 ||
                retained_count >= request.selected_count) {
              refuse(result, RnaControllerRefusal::InvalidContinuation);
              return result;
            }
            result.continuation.selected_begin = first_continuation;
            result.continuation.selected_count =
                request.selected_count - retained_count;
            const float ratio =
                static_cast<float>(result.continuation.selected_count) /
                static_cast<float>(request.selected_count);
            const float scaled = static_cast<float>(request.chain_score) * ratio;
            result.continuation.chain_score =
                static_cast<int>(static_cast<double>(scaled) + .499);
            result.continuation.initial_chain_score =
                request.initial_chain_score;
            result.region_anchor_count = retained_count;
            result.chain_score = request.chain_score -
                result.continuation.chain_score;
            result.split_flags = 1u | (2u << 2);
            ++result.metrics.split_continuations;
          }
          break;
        }
        result.dp_score += aligned.score;
        reference_cursor = packet_reference_end;
        query_cursor = packet_query_end;
      }
    }

    if (!dropped && last_query_pivot < envelope.query_end &&
        last_reference_pivot < envelope.reference_end) {
      const int query_packet_length =
          envelope.query_end - last_query_pivot;
      const int reference_packet_length =
          envelope.reference_end - last_reference_pivot;
      if (query_packet_length <= 0 || reference_packet_length <= 0) {
        refuse(result, RnaControllerRefusal::InvalidSequenceSlice);
        return result;
      }
      const std::uint8_t* junction = nullptr;
      if (materialize_junctions) {
        const std::uint64_t annotated_endpoint_bytes =
            materialize_known_junction_mask(
                request.known_junctions, last_reference_pivot,
                envelope.reference_end, junction_scratch);
        if (annotated_endpoint_bytes) {
          ++result.metrics.annotation_packets_masked;
          result.metrics.annotation_endpoint_bytes_set +=
              annotated_endpoint_bytes;
        }
        junction = junction_scratch.data();
      }
      dp::SpliceKernelResult& aligned = scratch.aligned;
      if (!kernel.run(
              SplicePacketRole::Right,
              last_query_pivot, envelope.query_end,
              last_reference_pivot, envelope.reference_end,
              mapping_query.data() + last_query_pivot,
              request.reference->data() + last_reference_pivot, junction,
              base_flags | dp::kSpliceExtensionOnly,
              request.options.zdrop, request.options.end_bonus,
              false, false, aligned))
        return result;
      append_cigar(result.cigar, aligned.cigar);
      if (!aligned.cigar.empty()) result.dp_score += aligned.maximum;
      const int reference_extension = aligned.reached_end
          ? aligned.query_end_reference + 1
          : aligned.maximum_reference + 1;
      const int query_extension = aligned.reached_end
          ? query_packet_length : aligned.maximum_query + 1;
      oriented_reference_end = last_reference_pivot + reference_extension;
      oriented_query_end = last_query_pivot + query_extension;
    }

    // Terminal settlement: a junction at a record end must have evidence
    // outside itself (clip_unsupported_terminal_junctions above). This runs on
    // the assembled CIGAR, after every packet and before the bounds are
    // published, so it moves the record's two ends and nothing else.
    if (!clip_unsupported_terminal_junctions(
            result.cigar, anchors[filtered_begin],
            anchors[filtered_begin + filtered_count - 1],
            oriented_query_begin, oriented_query_end,
            oriented_reference_begin, oriented_reference_end)) {
      refuse(result, RnaControllerRefusal::EmptyCigar);
      return result;
    }

    if (oriented_query_end > query_length ||
        oriented_reference_end > reference_length ||
        oriented_query_begin < 0 || oriented_reference_begin < 0 ||
        oriented_query_end < oriented_query_begin ||
        oriented_reference_end < oriented_reference_begin) {
      refuse(result, RnaControllerRefusal::InvalidSequenceSlice);
      return result;
    }
    result.reference_begin = oriented_reference_begin;
    result.reference_end = oriented_reference_end;
    if (!request.mapping_reverse) {
      result.query_begin = oriented_query_begin;
      result.query_end = oriented_query_end;
    } else {
      result.query_begin = query_length - oriented_query_end;
      result.query_end = query_length - oriented_query_begin;
    }
    if (result.region_anchor_count == 0) {
      result.region_anchor_count = request.selected_count;
      result.chain_score = request.chain_score;
    }
    if (result.cigar.empty()) {
      refuse(result, RnaControllerRefusal::EmptyCigar);
      return result;
    }
    if (!update_alignment_summary(
            result, mapping_query, *request.reference,
            oriented_query_begin, oriented_query_end,
            matrix, request.options, request.mapping_reverse)) {
      refuse(result, RnaControllerRefusal::CigarReconciliation);
      return result;
    }
    return result;
  } catch (const std::bad_alloc&) {
    refuse(result, RnaControllerRefusal::AllocationFailure);
  } catch (...) {
    refuse(result, RnaControllerRefusal::InternalInvariant);
  }
  return result;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
