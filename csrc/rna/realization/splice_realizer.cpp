#include "splice_realizer.h"
#include "splice_anchor_view.h"

#include "../../core/cigar.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <new>

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
  // itself (rival_lifecycle.cpp n_mis).
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
// segment, by core/cigar.cpp's replay as in DNA (`~` intron token included). The replay's
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
        raw.matches, raw.dp_maximum, raw.is_spliced});

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
    if (request.index_k != 15 || request.options.k != request.index_k) {
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
