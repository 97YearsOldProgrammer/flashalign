// Realization of a fixed exact anchor path under one or both transcript hypotheses.
#pragma once

#include "../../core/cigar.h"   // output::CigarReplayRequest (cs/MD)
#include "../../core/types.h"
#include "controller_metrics.h"
#include "../anchoring/exact_anchor_path.h"
#include "splice_controller.h"

#include <cstdint>
#include <string>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct RnaSpliceRealizationRequest {
  const std::vector<std::uint8_t>* query_forward = nullptr;
  const std::vector<std::uint8_t>* query_reverse = nullptr;
  const std::vector<std::uint8_t>* reference = nullptr;
  std::string reference_name;
  KnownJunctionContigView known_junctions;
  int index_k = 15;
  SpliceControllerOptions options = ordinary_long_read_splice_options();
  // cs:Z / MD:Z serialization (minimap2 --cs / --MD) of the elected segments, by the
  // same core/cigar.cpp replay as DNA, `~` intron token included. A default request
  // costs nothing.
  ::fa::cpu::output::CigarReplayRequest cigar_replay_request;
  // 0 runs both transcript hypotheses and elects one by DP score (minimap2 -ub); 1, 2
  // and 3 run one. The orientation is read-relative: 1 is the sense strand (-uf), 2
  // antisense (-ur), 3 scores no splice motif (-un) and the record gets no ts:A tag. The
  // controller maps it to the reference-frame motif per chain, so both mapping strands
  // are still needed upstream.
  int forced_transcript_orientation = 0;
  // Optional per-worker scratch shared by both hypotheses and every segment; null keeps
  // buffers function-local.
  SpliceRealizationScratch* scratch = nullptr;
};

// Rebinds every reference-owned view (sequence, name, known junctions) at once, so a
// request reused for another contig cannot keep the previous contig's annotation view.
bool bind_splice_realization_reference(
    RnaSpliceRealizationRequest& request,
    const std::vector<std::vector<std::uint8_t>>& references,
    const std::vector<std::string>& reference_names,
    const KnownJunctionStore* known_junctions, int reference_id) noexcept;

struct SpliceRegionEvidence {
  int query_begin = 0;
  int query_end = 0;
  int region_anchor_count = 0;
  int matches = 0;
  int dp_maximum = 0;
  bool is_spliced = false;
};

struct RnaSpliceHypothesisResult {
  TranscriptOrientation orientation = TranscriptOrientation::Forward;
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;
  std::vector<AlignResult> segments;
  std::vector<SpliceRegionEvidence> region_evidence;
  int primary_dp_score = 0;
  int primary_dp_maximum = 0;
  // One numeric CIGAR per segment. Public CIGAR strings are built for the winning
  // hypothesis only; the pick reads DP scores, never a CIGAR.
  std::vector<std::vector<std::uint32_t>> segment_cigars;
  RnaControllerMetrics metrics;
};

enum class SpliceRegionReject : std::uint8_t {
  None = 0,
  AnchorCount,
  MatchBases,
  DpMaximum,
  QueryRange,
};

inline SpliceRegionReject
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

struct RnaSpliceRealizationResult {
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;
  SpliceRegionReject region_reject = SpliceRegionReject::None;
  // Evidence for the first region verdict, including non-refusing quality.
  SpliceRegionEvidence rejected_region;
  std::uint64_t anchor_path_hash = 0;
  std::uint64_t controller_view_hash = 0;
  int transcript_strand = 0;  // 1 '+', 2 '-', 3 exact tie, 0 unavailable
  TranscriptOrientation winning_hypothesis = TranscriptOrientation::Forward;
  std::vector<AlignResult> segments;
  std::vector<RnaSpliceHypothesisResult> hypotheses;
  RnaControllerMetrics metrics;
};

RnaSpliceRealizationResult realize_exact_anchor_path(
    const RnaSpliceRealizationRequest& request,
    const ExactAnchorPath& bundle) noexcept;

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
