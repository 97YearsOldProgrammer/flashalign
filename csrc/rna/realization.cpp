// A fixed anchor path realized under one or both transcript hypotheses: the anchor view, the segments, CIGAR, cs and MD; known junctions from a BED file and their per-contig mask.
#include "rna.h"

#include "../core/core.h"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <zlib.h>

/**** splice_anchor_view ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

constexpr uint64_t kFnvOffset = 1469598103934665603ull;
constexpr uint64_t kFnvPrime = 1099511628211ull;

uint64_t mix_u64(uint64_t hash, uint64_t value) {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (byte * 8)) & 0xffu;
    hash *= kFnvPrime;
  }
  return hash;
}

uint64_t mix_vector(uint64_t hash, const std::vector<uint32_t>& values) {
  hash = mix_u64(hash, values.size());
  for (uint32_t value : values) hash = mix_u64(hash, value);
  return hash;
}

uint64_t hash_view(const SpliceAnchorView& view) {
  uint64_t hash = kFnvOffset;
  hash = mix_u64(hash, static_cast<uint32_t>(view.selected_begin));
  hash = mix_u64(hash, static_cast<uint32_t>(view.selected_count));
  for (size_t i = 0; i < view.anchors.size(); ++i) {
    const SpliceAnchor& anchor = view.anchors[i];
    const int64_t query_end = static_cast<int64_t>(anchor.query_begin) +
                              anchor.span - 1;
    const int64_t reference_end =
        static_cast<int64_t>(anchor.reference_begin) + anchor.span - 1;
    const uint64_t conceptual_x =
        (view.mapping_reverse ? (1ull << 63) : 0ull) |
        static_cast<uint32_t>(reference_end);
    uint64_t conceptual_y = static_cast<uint64_t>(anchor.span) << 32 |
                            static_cast<uint32_t>(query_end);
    if ((anchor.flags & kSpliceAnchorLongJoin) != 0)
      conceptual_y |= 1ull << 40;
    if ((anchor.flags & kSpliceAnchorIgnore) != 0)
      conceptual_y |= 1ull << 41;
    if ((anchor.flags & kSpliceAnchorTandem) != 0)
      conceptual_y |= 1ull << 42;
    hash = mix_u64(hash, conceptual_x);
    hash = mix_u64(hash, conceptual_y);
    hash = mix_u64(hash, anchor.stable_id);
    hash = mix_u64(hash, view.view_to_raw_index[i]);
  }
  return mix_vector(hash, view.excluded_interior_raw_indices);
}

// Builds the realization view of a selected anchor path; a non-empty `error` means the
// bundle was rejected.
SpliceAnchorView make_splice_anchor_view(
    const ExactAnchorPath& bundle, int index_k, std::string* error) {
  SpliceAnchorView view;
  AnchoringRefusal validation_refusal = AnchoringRefusal::None;
  if (!validate_exact_anchor_path(bundle, &validation_refusal)) {
    if (error) *error = anchoring_refusal_name(validation_refusal);
    return view;
  }
  if (index_k <= 0 || index_k > 255) {
    if (error) *error = "SPLICE_VIEW_INDEX_K";
    return view;
  }
  int raw_count = 0;
  int selected_count = 0;
  if (!mapping::checked_size_to_int(bundle.selected_pool_anchors.size(), raw_count) ||
      !mapping::checked_size_to_int(
          bundle.selected_raw_indices.size(), selected_count)) {
    if (error) *error = "SPLICE_VIEW_COUNT_DOMAIN";
    return view;
  }
  view.anchors.reserve(static_cast<std::size_t>(raw_count));
  view.view_to_raw_index.reserve(static_cast<std::size_t>(raw_count));

  const uint32_t first = bundle.selected_raw_indices.front();
  const uint32_t last = bundle.selected_raw_indices.back();
  std::vector<uint8_t> selected(bundle.selected_pool_anchors.size(), 0);
  for (uint32_t index : bundle.selected_raw_indices) selected[index] = 1;

  view.mapping_reverse = bundle.mapping_reverse;
  auto append = [&](uint32_t raw_index) {
    const SelectedPoolAnchor& raw = bundle.selected_pool_anchors[raw_index];
    if (raw.span != index_k) return false;
    const int64_t query_end = static_cast<int64_t>(raw.query_begin) +
                              raw.span - 1;
    const int64_t reference_end =
        static_cast<int64_t>(raw.reference_begin) + raw.span - 1;
    if (query_end < 0 || query_end > UINT32_MAX || reference_end < 0 ||
        reference_end > UINT32_MAX)
      return false;

    SpliceAnchor anchor;
    anchor.stable_id = raw.stable_id;
    anchor.query_begin = raw.query_begin;
    anchor.reference_begin = raw.reference_begin;
    anchor.span = raw.span;
    if (raw.long_join_flag) anchor.flags |= kSpliceAnchorLongJoin;
    if (raw.ignore_flag) anchor.flags |= kSpliceAnchorIgnore;
    if (raw.tandem_flag) anchor.flags |= kSpliceAnchorTandem;
    view.anchors.push_back(anchor);
    view.view_to_raw_index.push_back(raw_index);
    return true;
  };

  for (uint32_t index = 0; index < first; ++index) {
    if (!append(index)) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  if (!mapping::checked_size_to_int(view.anchors.size(), view.selected_begin)) {
    if (error) *error = "SPLICE_VIEW_COUNT_DOMAIN";
    return {};
  }
  for (uint32_t index : bundle.selected_raw_indices) {
    if (!append(index)) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  view.selected_count = selected_count;
  for (uint32_t index = first; index <= last; ++index) {
    if (!selected[index]) view.excluded_interior_raw_indices.push_back(index);
  }
  for (uint64_t index = static_cast<uint64_t>(last) + 1;
       index < bundle.selected_pool_anchors.size(); ++index) {
    if (!append(static_cast<uint32_t>(index))) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  view.content_hash = hash_view(view);
  if (error) error->clear();
  return view;
}

}  // namespace

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** splice_realizer ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

SpliceRegionReject
classify_splice_region(const SpliceRegionEvidence& evidence, int query_length,
                       const SpliceControllerOptions& options) noexcept {
  // Validity comes first: a quality reject may leave the region in place, a malformed
  // range may not.
  if (evidence.query_begin < 0 || evidence.query_end < evidence.query_begin ||
      evidence.query_end > query_length)
    return SpliceRegionReject::QueryRange;
  if (evidence.region_anchor_count < options.minimum_anchor_count)
    return SpliceRegionReject::AnchorCount;
  if (evidence.matches < options.minimum_match_bases)
    return SpliceRegionReject::MatchBases;
  if (evidence.dp_maximum < options.minimum_dp_maximum)
    return SpliceRegionReject::DpMaximum;
  return SpliceRegionReject::None;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

bool bind_splice_realization_reference(
    RnaSpliceRealizationRequest& request,
    const std::vector<std::vector<std::uint8_t>>& references,
    const std::vector<std::string>& reference_names,
    const KnownJunctionStore* known_junctions, int reference_id) noexcept {
  if (reference_id < 0 || references.size() != reference_names.size() ||
      static_cast<std::size_t>(reference_id) >= references.size())
    return false;

  const std::size_t id = static_cast<std::size_t>(reference_id);
  request.reference = &references[id];
  request.reference_name = reference_names[id];
  request.known_junctions = known_junctions
                                ? known_junctions->contig(id)
                                : KnownJunctionContigView{};
  return true;
}

namespace {

void refuse(RnaSpliceRealizationResult& result,
            RnaControllerRefusal reason) noexcept {
  result.refused = true;
  result.refusal = reason;
  result.segments.clear();
  result.metrics.observe_refusal(reason);
}

std::string numeric_cigar_body(
    const std::vector<std::uint32_t>& cigar) {
  static constexpr char kOperations[] = "MIDNSHP=XB";
  std::string text;
  for (std::uint32_t encoded : cigar) {
    const std::uint32_t operation = encoded & 0xfu;
    const std::uint32_t length = encoded >> 4;
    if (length == 0 || operation >= sizeof(kOperations) - 1) return {};
    ::fa::cpu::output::append_cigar_run(
        text, static_cast<int>(length), kOperations[operation]);
  }
  return text;
}

bool numeric_cigar_reconciles(const SpliceSegmentResult& raw) noexcept {
  std::int64_t query = 0;
  std::int64_t reference = 0;
  for (std::uint32_t encoded : raw.cigar) {
    const int operation = static_cast<int>(encoded & 0xfu);
    const std::int64_t length = encoded >> 4;
    // An operation numeric_cigar_body has no letter for cannot be projected.
    if (length <= 0 || operation >= 10) return false;
    if (operation == 0 || operation == 1 ||
        operation == 7 || operation == 8)
      query += length;
    if (operation == 0 || operation == 2 || operation == 3 ||
        operation == 7 || operation == 8)
      reference += length;
  }
  return query == raw.query_end - raw.query_begin &&
      reference == raw.reference_end - raw.reference_begin;
}

// Everything about a realized segment except its public CIGAR string, which is built
// only for the winning hypothesis. The checks here reject everything that would make that
// later build fail.
bool project_segment(const RnaSpliceRealizationRequest& request,
                     const ExactAnchorPath& bundle,
                     const SpliceSegmentResult& raw,
                     AlignResult& segment) {
  if (raw.cigar.empty() || !numeric_cigar_reconciles(raw)) return false;
  if (bundle.query_length < 0 || raw.query_begin < 0 ||
      raw.query_end < raw.query_begin ||
      raw.query_end > bundle.query_length)
    return false;
  segment.chromosome = request.reference_name;
  segment.pos = raw.reference_begin;
  segment.read_len = bundle.query_length;
  segment.is_reverse = bundle.mapping_reverse;
  segment.query_start = raw.query_begin;
  segment.query_end = raw.query_end;
  segment.target_end = raw.reference_end;
  segment.score = raw.dp_score;
  // ms:i is minimap2's dp_max0: the segment's running maximum before the transcript-strand
  // adjustment, which only the evidence copy receives (apply_transcript_strand).
  segment.dp_max_segment = raw.dp_maximum_before_strand;
  segment.matches = raw.matches;
  segment.block_len = raw.block_length;
  segment.ambiguities = raw.ambiguous_bases;
  // NM:i is minimap2's blen - mlen + n_ambi: the controller leaves ambiguous columns out
  // of block_length and matches, and NM charges each once.
  segment.edit_distance =
      std::max(0, raw.block_length - raw.matches) + raw.ambiguous_bases;
  segment.alignment_accounting_valid = true;
  for (std::uint32_t encoded : raw.cigar) {
    const int operation = static_cast<int>(encoded & 0xfu);
    const int length = static_cast<int>(encoded >> 4);
    if (operation == 1) segment.insertions += length;
    else if (operation == 2) segment.deletions += length;
  }
  // Unambiguous substitutions only; the rank recalibration adds the ambiguities back
  // itself (hit.cpp n_mis).
  segment.mismatches = std::max(0, raw.block_length - raw.matches -
                                       segment.insertions - segment.deletions);
  return true;
}

// The public reference-oriented CIGAR of one segment, built for the winner only.
bool materialize_segment_cigar(const ExactAnchorPath& bundle,
                               const std::vector<std::uint32_t>& cigar,
                               AlignResult& segment) {
  const std::string body = numeric_cigar_body(cigar);
  if (body.empty()) return false;
  segment.cigar = ::fa::cpu::output::clipped_reference_oriented_cigar(
      body, bundle.query_length, segment.query_start, segment.query_end,
      bundle.mapping_reverse);
  return !segment.cigar.empty();
}

// cs:Z / MD:Z strings and =/X CIGAR (minimap2 --cs / --MD / --eqx) of one elected
// segment, by core/core.cpp's replay as in DNA (`~` intron token included). The replay's
// accounting must agree with the controller's; a disagreement refuses the realization
// rather than emit a difference string that contradicts NM:i.
bool attach_difference_strings(const RnaSpliceRealizationRequest& request,
                               const ExactAnchorPath& bundle,
                               AlignResult& segment) {
  if (!request.cigar_replay_request.any()) return true;
  const std::vector<std::uint8_t>& query =
      bundle.mapping_reverse ? *request.query_reverse : *request.query_forward;
  const ::fa::cpu::output::CigarReplay replay =
      ::fa::cpu::output::replay_cigar(
          segment.cigar, query.data(), bundle.query_length,
          request.reference->data(),
          static_cast<int>(request.reference->size()), segment.pos,
          request.cigar_replay_request);
  if (!replay.valid || replay.target_start != segment.pos ||
      replay.target_end != segment.target_end ||
      replay.matches != segment.matches ||
      replay.block_len != segment.block_len ||
      replay.ambiguities != segment.ambiguities ||
      replay.edit_distance != segment.edit_distance)
    return false;
  segment.cs = replay.cs;
  segment.md = replay.md;
  // --eqx: the same walk's =/X spelling of the CIGAR.
  if (request.cigar_replay_request.eqx)
    segment.cigar = replay.eqx_cigar;
  return true;
}

RnaSpliceHypothesisResult run_hypothesis(
    const RnaSpliceRealizationRequest& request,
    const ExactAnchorPath& bundle,
    const SpliceAnchorView& frozen_view,
    TranscriptOrientation orientation) {
  RnaSpliceHypothesisResult hypothesis;
  hypothesis.orientation = orientation;
  hypothesis.metrics.transcript_hypothesis_passes = 1;
  // The controller mutates its anchors, so each hypothesis works on its own copy of the
  // frozen view, held in the caller's scratch when there is one.
  std::vector<SpliceAnchor> owned_working;
  std::vector<SpliceAnchor>& working =
      request.scratch ? request.scratch->working : owned_working;
  working = frozen_view.anchors;
  int selected_begin = frozen_view.selected_begin;
  int selected_count = frozen_view.selected_count;
  int chain_score = bundle.anchor_path_score;
  int initial_chain_score = bundle.anchor_path_score;
  constexpr int kMaximumSegments = 1024;

  for (int segment_id = 0; selected_count > 0; ++segment_id) {
    if (segment_id >= kMaximumSegments) {
      hypothesis.refused = true;
      hypothesis.refusal = RnaControllerRefusal::ContinuationLimit;
      hypothesis.metrics.observe_refusal(hypothesis.refusal);
      break;
    }
    SpliceControllerRequest segment_request;
    segment_request.query_forward = request.query_forward;
    segment_request.query_reverse = request.query_reverse;
    segment_request.reference = request.reference;
    segment_request.working_anchors = &working;
    segment_request.selected_begin = selected_begin;
    segment_request.selected_count = selected_count;
    segment_request.chain_score = chain_score;
    segment_request.initial_chain_score = initial_chain_score;
    segment_request.mapping_reverse = bundle.mapping_reverse;
    segment_request.transcript_orientation = orientation;
    segment_request.options = request.options;
    segment_request.known_junctions = request.known_junctions;
    segment_request.scratch = request.scratch;
    SpliceSegmentResult raw = realize_splice_segment(segment_request);
    hypothesis.metrics.merge(raw.metrics);
    if (raw.refused || raw.cigar.empty()) {
      hypothesis.refused = true;
      hypothesis.refusal = raw.refused
          ? raw.refusal : RnaControllerRefusal::EmptyCigar;
      if (!raw.refused) hypothesis.metrics.observe_refusal(hypothesis.refusal);
      break;
    }

    AlignResult alignment;
    if (!project_segment(request, bundle, raw, alignment)) {
      hypothesis.refused = true;
      hypothesis.refusal = RnaControllerRefusal::CigarReconciliation;
      hypothesis.metrics.observe_refusal(hypothesis.refusal);
      break;
    }
    if (segment_id == 0) {
      hypothesis.primary_dp_score = raw.dp_score;
      hypothesis.primary_dp_maximum = raw.dp_maximum;
    }
    hypothesis.segments.push_back(std::move(alignment));
    hypothesis.segment_cigars.push_back(std::move(raw.cigar));
    hypothesis.region_evidence.push_back({
        raw.query_begin, raw.query_end, raw.region_anchor_count,
        raw.matches, raw.dp_maximum, raw.is_spliced, raw.chain_score});

    if (!raw.continuation.present()) break;
    if (raw.continuation.selected_begin <= selected_begin ||
        raw.continuation.selected_count >= selected_count ||
        raw.continuation.selected_begin < 0 ||
        raw.continuation.selected_count <= 0 ||
        static_cast<std::size_t>(raw.continuation.selected_begin) >
            working.size() ||
        static_cast<std::size_t>(raw.continuation.selected_count) >
            working.size() -
                static_cast<std::size_t>(raw.continuation.selected_begin)) {
      hypothesis.refused = true;
      hypothesis.refusal = RnaControllerRefusal::InvalidContinuation;
      hypothesis.metrics.observe_refusal(hypothesis.refusal);
      break;
    }
    selected_begin = raw.continuation.selected_begin;
    selected_count = raw.continuation.selected_count;
    chain_score = raw.continuation.chain_score;
    initial_chain_score = raw.continuation.initial_chain_score;
  }
  return hypothesis;
}

void apply_transcript_strand(RnaSpliceHypothesisResult& hypothesis,
                             int transcript_strand,
                             const SpliceControllerOptions& options) {
  const char strand = transcript_strand == 1 ? '+'
                    : transcript_strand == 2 ? '-' : '\0';
  for (std::size_t index = 0; index < hypothesis.segments.size(); ++index) {
    hypothesis.segments[index].transcript_strand = strand;
    SpliceRegionEvidence& evidence = hypothesis.region_evidence[index];
    if (!evidence.is_spliced) continue;
    const int pair = options.match + options.mismatch;
    if (transcript_strand == 1 || transcript_strand == 2)
      evidence.dp_maximum += pair + (pair >> 1);
    else if (transcript_strand == 3)
      evidence.dp_maximum -= pair;
    if (evidence.dp_maximum < 0) evidence.dp_maximum = 0;
    if (index == 0)
      hypothesis.primary_dp_maximum = evidence.dp_maximum;
  }
}

}  // namespace

RnaSpliceRealizationResult realize_exact_anchor_path(
    const RnaSpliceRealizationRequest& request,
    const ExactAnchorPath& bundle) noexcept {
  RnaSpliceRealizationResult result;
  result.anchor_path_hash = bundle.anchor_path_hash;
  result.metrics.top_level_controller_calls = 1;
  try {
    if (!request.query_forward || !request.query_reverse ||
        !request.reference || request.query_reverse->size() !=
            request.query_forward->size() || request.reference->empty()) {
      refuse(result, RnaControllerRefusal::InvalidRequest);
      return result;
    }
    const std::size_t maximum = static_cast<std::size_t>(
        std::numeric_limits<int>::max());
    if (request.query_forward->size() > maximum ||
        request.reference->size() > maximum) {
      refuse(result, RnaControllerRefusal::UnrepresentableDomain);
      return result;
    }
    if (static_cast<int>(request.query_forward->size()) !=
        bundle.query_length) {
      refuse(result, RnaControllerRefusal::InvalidRequest);
      return result;
    }
    if (request.options.k != request.index_k) {
      refuse(result, RnaControllerRefusal::UnsupportedProfile);
      return result;
    }
    std::string boundary_error;
    const SpliceAnchorView view = make_splice_anchor_view(
        bundle, request.index_k, &boundary_error);
    if (!boundary_error.empty()) {
      refuse(result, RnaControllerRefusal::InvalidBundle);
      return result;
    }
    if (view.anchors.empty()) {
      refuse(result, RnaControllerRefusal::EmptyAnchorView);
      return result;
    }
    result.controller_view_hash = view.content_hash;
    // Prepared segments are keyed by selected range over this anchor geometry only.
    if (request.scratch) request.scratch->prepared.clear();

    const bool forced = request.forced_transcript_orientation != 0;
    if (forced && request.forced_transcript_orientation != 1 &&
        request.forced_transcript_orientation != 2 &&
        request.forced_transcript_orientation != 3) {
      refuse(result, RnaControllerRefusal::InvalidRequest);
      return result;
    }
    if (forced) {
      const auto orientation = static_cast<TranscriptOrientation>(
          request.forced_transcript_orientation);
      result.hypotheses.push_back(
          run_hypothesis(request, bundle, view, orientation));
    } else {
      result.hypotheses.push_back(run_hypothesis(
          request, bundle, view, TranscriptOrientation::Forward));
      result.hypotheses.push_back(run_hypothesis(
          request, bundle, view, TranscriptOrientation::Reverse));
    }
    for (const RnaSpliceHypothesisResult& hypothesis : result.hypotheses)
      result.metrics.merge(hypothesis.metrics);
    for (const RnaSpliceHypothesisResult& hypothesis : result.hypotheses) {
      if (hypothesis.refused || hypothesis.segments.empty()) {
        result.refused = true;
        result.refusal = hypothesis.refused
            ? hypothesis.refusal : RnaControllerRefusal::IncompleteHypothesis;
        if (!hypothesis.refused)
          result.metrics.observe_refusal(result.refusal);
        return result;
      }
    }

    std::size_t winner = 0;
    int transcript_strand = static_cast<int>(
        result.hypotheses.front().orientation);
    // -un scores no motif, so there is no transcript strand: 0 stamps no ts:A and applies
    // neither the motif bonus nor the tie penalty (3 is the exact-tie code).
    if (result.hypotheses.front().orientation == TranscriptOrientation::None)
      transcript_strand = 0;
    if (!forced) {
      const int forward_score = result.hypotheses[0].primary_dp_score;
      const int reverse_score = result.hypotheses[1].primary_dp_score;
      if (forward_score > reverse_score) {
        winner = 0;
        transcript_strand = 1;
      } else if (forward_score < reverse_score) {
        winner = 1;
        transcript_strand = 2;
      } else {
        transcript_strand = 3;
        winner = static_cast<std::size_t>(
            (bundle.query_length + forward_score) & 1);
      }
    }
    result.transcript_strand = transcript_strand;
    result.winning_hypothesis = result.hypotheses[winner].orientation;
    const SpliceControllerOptions& options = request.options;
    apply_transcript_strand(
        result.hypotheses[winner], transcript_strand, options);
    // Quality verdicts are kept as evidence without dropping a valid placement; only
    // invalid query geometry refuses.
    for (const SpliceRegionEvidence& evidence :
         result.hypotheses[winner].region_evidence) {
      const SpliceRegionReject verdict =
          classify_splice_region(evidence, bundle.query_length, options);
      if (verdict == SpliceRegionReject::None)
        continue;
      if (result.region_reject == SpliceRegionReject::None) {
        result.region_reject = verdict;
        result.rejected_region = evidence;
      }
      if (verdict == SpliceRegionReject::QueryRange) {
        refuse(result, RnaControllerRefusal::RegionFilter);
        result.region_reject = verdict;
        result.rejected_region = evidence;
        return result;
      }
    }
    // Only the winner's segments become records, so only they pay for CIGAR, cs:Z and
    // MD:Z strings.
    RnaSpliceHypothesisResult& elected = result.hypotheses[winner];
    for (std::size_t index = 0; index < elected.segments.size(); ++index) {
      if (materialize_segment_cigar(bundle, elected.segment_cigars[index],
                                    elected.segments[index]) &&
          attach_difference_strings(request, bundle, elected.segments[index]))
        continue;
      refuse(result, RnaControllerRefusal::CigarReconciliation);
      return result;
    }
    result.segments = elected.segments;
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

/**** known_junctions ****/

namespace fa::cpu::lr::rna {
namespace {
void strip_carriage_return(std::string& line) {
  if (!line.empty() && line.back() == '\r')
    line.pop_back();
}

std::vector<std::string_view> fields(std::string_view line,
                                     char delimiter = '\t') {
  std::vector<std::string_view> out;
  while (true) {
    const auto pos = line.find(delimiter);
    out.push_back(line.substr(0, pos));
    if (pos == std::string_view::npos)
      break;
    line.remove_prefix(pos + 1);
  }
  return out;
}
std::uint32_t number(std::string_view text, const std::string& path,
                     std::uint64_t line) {
  std::uint64_t value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (text.empty() || result.ec != std::errc{} ||
      result.ptr != text.data() + text.size() ||
      value > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument(path + ":" + std::to_string(line) +
                                ": invalid BED coordinate");
  return static_cast<std::uint32_t>(value);
}
std::vector<std::string> read_lines(const std::string& path) {
  std::vector<std::string> out;
  if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
    gzFile file = gzopen(path.c_str(), "rb");
    if (!file)
      throw std::invalid_argument("cannot open RNA junction BED: " + path);
    std::string line;
    char buffer[65536];
    while (gzgets(file, buffer, sizeof buffer)) {
      line += buffer;
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
        strip_carriage_return(line);
        out.push_back(line);
        line.clear();
      }
    }
    int gz_error = Z_OK;
    (void)gzerror(file, &gz_error);
    const int error = gzclose(file);
    if (!line.empty()) {
      strip_carriage_return(line);
      out.push_back(line);
    }
    if ((gz_error != Z_OK && gz_error != Z_STREAM_END) || error != Z_OK)
      throw std::invalid_argument("cannot decompress RNA junction BED: " +
                                  path);
  } else {
    std::ifstream file(path);
    if (!file)
      throw std::invalid_argument("cannot open RNA junction BED: " + path);
    std::string line;
    while (std::getline(file, line)) {
      strip_carriage_return(line);
      out.push_back(line);
    }
    if (file.bad())
      throw std::invalid_argument("cannot read RNA junction BED: " + path);
  }
  return out;
}
} // namespace

KnownJunctionStore KnownJunctionStore::load_bed(
    const std::string& path, const std::vector<std::string>& names,
    const std::vector<std::vector<std::uint8_t>>& references) {
  if (names.size() != references.size())
    throw std::invalid_argument(
        "RNA junction BED: invalid reference dictionary");
  std::unordered_map<std::string, std::size_t> ids;
  for (std::size_t i = 0; i < names.size(); ++i)
    if (!ids.emplace(names[i], i).second)
      throw std::invalid_argument(
          "RNA junction BED: duplicate reference contig: " + names[i]);
  KnownJunctionStore store;
  store.contigs_.resize(names.size());
  std::uint64_t line_number = 0;
  for (const std::string& line : read_lines(path)) {
    ++line_number;
    if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 ||
        line.rfind("browser", 0) == 0)
      continue;
    ++store.counters_.rows;
    const auto f = fields(line);
    if (f.size() != 6 && f.size() != 12)
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": expected BED6 (strand in column 6) or BED12");
    const auto found = ids.find(std::string(f[0]));
    if (found == ids.end())
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": unknown contig: " + std::string(f[0]));
    if (f[5] != "+" && f[5] != "-")
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": strand must be + or -");
    const std::uint8_t strand = f[5] == "+" ? 1 : 2;
    const auto chrom_start = number(f[1], path, line_number);
    const auto chrom_end = number(f[2], path, line_number);
    if (chrom_end <= chrom_start ||
        chrom_end > references[found->second].size())
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": BED interval outside reference");
    auto add = [&](std::uint64_t start, std::uint64_t end) {
      if (end <= start || end > references[found->second].size() ||
          end > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                    ": intron outside reference");
      store.contigs_[found->second].push_back(
          {static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(end),
           strand});
      ++store.counters_.introns;
    };
    if (f.size() == 6)
      add(chrom_start, chrom_end);
    else {
      const auto count = number(f[9], path, line_number);
      const auto sizes = fields(f[10], ',');
      const auto starts = fields(f[11], ',');
      auto usable = [](const std::vector<std::string_view>& v) {
        return !v.empty() && v.back().empty() ? v.size() - 1 : v.size();
      };
      if (count == 0 || usable(sizes) != count || usable(starts) != count)
        throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                    ": malformed BED12 blocks");
      if (count == 1)
        ++store.counters_.single_block_rows;
      std::uint64_t previous_end = 0;
      for (std::uint32_t i = 0; i < count; ++i) {
        const auto block_start = number(starts[i], path, line_number);
        const auto block_size = number(sizes[i], path, line_number);
        const std::uint64_t exon_start =
            std::uint64_t(chrom_start) + block_start;
        const std::uint64_t exon_end = exon_start + block_size;
        if (block_size == 0 || exon_end > chrom_end ||
            (i && exon_start < previous_end))
          throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                      ": invalid BED12 block geometry");
        if (i)
          add(previous_end, exon_start);
        previous_end = exon_end;
      }
    }
  }
  for (auto& contig : store.contigs_) {
    std::sort(contig.begin(), contig.end(), [](const auto& a, const auto& b) {
      return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    std::vector<KnownIntron> merged;
    for (const auto& intron : contig) {
      if (!merged.empty() && merged.back().start == intron.start &&
          merged.back().end == intron.end) {
        merged.back().strand_mask |= intron.strand_mask;
        ++store.counters_.duplicates_merged;
      } else
        merged.push_back(intron);
    }
    contig = std::move(merged);
  }
  return store;
}
KnownJunctionContigView
KnownJunctionStore::contig(std::size_t id) const noexcept {
  if (id >= contigs_.size())
    return {};
  return {contigs_[id].data(), contigs_[id].data() + contigs_[id].size()};
}
std::uint64_t materialize_known_junction_mask(KnownJunctionContigView view,
                                              int begin, int end,
                                              std::vector<std::uint8_t>& out) {
  if (begin < 0 || end <= begin) {
    out.clear();
    return 0;
  }
  out.assign(static_cast<std::size_t>(end - begin), 0);
  if (view.empty())
    return 0;
  auto it = std::lower_bound(view.begin, view.end, begin,
                             [](const KnownIntron& i, int p) {
                               return i.start < static_cast<std::uint32_t>(p);
                             });
  for (; it != view.end && it->start < static_cast<std::uint32_t>(end); ++it) {
    if (it->end > static_cast<std::uint32_t>(end))
      continue;
    auto& donor = out[it->start - begin];
    auto& acceptor = out[it->end - 1 - begin];
    if (it->strand_mask & 1) {
      donor |= 1;
      acceptor |= 2;
    }
    if (it->strand_mask & 2) {
      donor |= 8;
      acceptor |= 4;
    }
  }
  return static_cast<std::uint64_t>(std::count_if(
      out.begin(), out.end(), [](std::uint8_t value) { return value != 0; }));
}
} // namespace fa::cpu::lr::rna
