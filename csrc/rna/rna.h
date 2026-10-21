// The RNA lane's context, config, scratch, result, anchor-path and metric types, the query-span geometry, and the declarations of its stages.
#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include RNA context"
#endif

#include "../chain/chain.h"
#include "../core/core.h" // CoordinateDomainRefusal, output::CigarReplayRequest (cs/MD)
#include "../dp/dp.h"
#include "../index/format.h" // KmerPostingView
#include "../index/index.h" // engine::ReferenceContext
#include "../index/seed.h" // QuerySeed
#include "../seed/seed.h" // LongReadSeedContext, select_spread_from_evidence, VoteHit, SeedDrop

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**** controller_metrics ****/
// Saturating work counters for the RNA splice controller.

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

/**** config ****/
// Runtime configuration of the RNA pipeline, built once per configure by
// LongReadEngine::compose_rna_runtime and read by RnaBackend::map_read.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Transcript-strand policy, read-relative like minimap2's -u: Unknown runs both
// hypotheses and elects one (-ub); Forward treats every read as the sense strand (-uf),
// Reverse as antisense (-ur). It never fixes the mapping strand: seeding, voting and
// chaining stay two-stranded and only the second realization pass is skipped.
// None is minimap2's -un: no splice motif is scored, one hypothesis runs and no ts:A tag
// is emitted.
enum class StrandMode : uint8_t {
  Unknown = 0,
  Forward = 1,
  Reverse = 2,
  None = 3,
};

StrandMode parse_strand_mode(const char* text, StrandMode fallback);

struct RnaConfig {
  // Seeding and intron geometry, mirrored from the shared long-read config.
  int k = 15;
  int min_support = 2;
  int min_intron = 20;
  int max_intron = 200000;
  int max_locus_chains = 6;
  int max_chain_predecessors = 64;

  // Unknown runs minimap2's two transcript hypotheses; a forward/reverse hint runs one.
  StrandMode strand_mode = StrandMode::Unknown;
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** result ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct Result : AlignResult {};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** refusal ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class AnchoringRefusal : std::uint8_t {
  None = 0,
  RegionReferenceSpanOverBudget,
  RegionVolumeOverBudget,
  RegionCountDomain,
  RegionCoordinateDomain,
  EmptyOrInvalidInput,
  InvalidSkeletonNode,
  NonMonotoneSkeleton,
  SkeletonSkipOverBudget,
  HarvestCountDomain,
  InvalidHarvestInput,
  InvalidIndexGeometry,
  RawCandidateSafetyLimit,
  DuplicateQueryStart,
  MirrorProjectionDomain,
  Rank1SkeletonCountDomain,
  InvalidRank1Skeleton,
  Rank1SkeletonCoordinateDomain,
  PoolRefused,
  PoolCountDomain,
  StableIdDomain,
  EmptyPool,
  NoAcceptedChain,
  RawSizeDrift,
  DuplicateStableId,
  TracebackRange,
  TracebackIdDomain,
  TracebackId,
  SelectedCountDomain,
  BundleSchema,
  BundleIdentity,
  BundleEmpty,
  BundleCountDomain,
  BundleStableIdDomain,
  BundleAnchorRange,
  BundleRawOrder,
  BundleSelectedCountDomain,
  BundleCount,
  BundleSelectedRange,
  BundleSelectedOrder,
  BundleSelectedSpan,
  BundleHash,
};

constexpr const char* anchoring_refusal_name(AnchoringRefusal refusal) noexcept {
  switch (refusal) {
    case AnchoringRefusal::None: return "";
    case AnchoringRefusal::RegionReferenceSpanOverBudget:
      return "REGION_REFERENCE_SPAN_OVER_BUDGET";
    case AnchoringRefusal::RegionVolumeOverBudget:
      return "REGION_VOLUME_OVER_BUDGET";
    case AnchoringRefusal::RegionCountDomain: return "REGION_COUNT_DOMAIN";
    case AnchoringRefusal::RegionCoordinateDomain:
      return "REGION_COORDINATE_DOMAIN";
    case AnchoringRefusal::EmptyOrInvalidInput:
      return "EMPTY_OR_INVALID_INPUT";
    case AnchoringRefusal::InvalidSkeletonNode:
      return "INVALID_SKELETON_NODE";
    case AnchoringRefusal::NonMonotoneSkeleton:
      return "NON_MONOTONE_SKELETON";
    case AnchoringRefusal::SkeletonSkipOverBudget:
      return "SKELETON_SKIP_OVER_BUDGET";
    case AnchoringRefusal::HarvestCountDomain: return "HARVEST_COUNT_DOMAIN";
    case AnchoringRefusal::InvalidHarvestInput: return "INVALID_HARVEST_INPUT";
    case AnchoringRefusal::InvalidIndexGeometry:
      return "INVALID_INDEX_GEOMETRY";
    case AnchoringRefusal::RawCandidateSafetyLimit:
      return "RAW_CANDIDATE_SAFETY_LIMIT";
    case AnchoringRefusal::DuplicateQueryStart:
      return "DUPLICATE_QUERY_START";
    case AnchoringRefusal::MirrorProjectionDomain:
      return "MIRROR_PROJECTION_DOMAIN";
    case AnchoringRefusal::Rank1SkeletonCountDomain:
      return "RANK1_SKELETON_COUNT_DOMAIN";
    case AnchoringRefusal::InvalidRank1Skeleton:
      return "INVALID_RANK1_SKELETON";
    case AnchoringRefusal::Rank1SkeletonCoordinateDomain:
      return "RANK1_SKELETON_COORDINATE_DOMAIN";
    case AnchoringRefusal::PoolRefused: return "M0_POOL_REFUSED";
    case AnchoringRefusal::PoolCountDomain: return "M0_POOL_COUNT_DOMAIN";
    case AnchoringRefusal::StableIdDomain: return "M0_STABLE_ID_DOMAIN";
    case AnchoringRefusal::EmptyPool: return "M0_EMPTY_POOL";
    case AnchoringRefusal::NoAcceptedChain: return "M0_NO_ACCEPTED_CHAIN";
    case AnchoringRefusal::RawSizeDrift: return "M0_RAW_SIZE_DRIFT";
    case AnchoringRefusal::DuplicateStableId:
      return "M0_DUPLICATE_STABLE_ID";
    case AnchoringRefusal::TracebackRange: return "M0_TRACEBACK_RANGE";
    case AnchoringRefusal::TracebackIdDomain:
      return "M0_TRACEBACK_ID_DOMAIN";
    case AnchoringRefusal::TracebackId: return "M0_TRACEBACK_ID";
    case AnchoringRefusal::SelectedCountDomain:
      return "M0_SELECTED_COUNT_DOMAIN";
    case AnchoringRefusal::BundleSchema: return "M0_BUNDLE_SCHEMA";
    case AnchoringRefusal::BundleIdentity: return "M0_BUNDLE_IDENTITY";
    case AnchoringRefusal::BundleEmpty: return "M0_BUNDLE_EMPTY";
    case AnchoringRefusal::BundleCountDomain:
      return "M0_BUNDLE_COUNT_DOMAIN";
    case AnchoringRefusal::BundleStableIdDomain:
      return "M0_BUNDLE_STABLE_ID_DOMAIN";
    case AnchoringRefusal::BundleAnchorRange:
      return "M0_BUNDLE_ANCHOR_RANGE";
    case AnchoringRefusal::BundleRawOrder: return "M0_BUNDLE_RAW_ORDER";
    case AnchoringRefusal::BundleSelectedCountDomain:
      return "M0_BUNDLE_SELECTED_COUNT_DOMAIN";
    case AnchoringRefusal::BundleCount: return "M0_BUNDLE_COUNT";
    case AnchoringRefusal::BundleSelectedRange:
      return "M0_BUNDLE_SELECTED_RANGE";
    case AnchoringRefusal::BundleSelectedOrder:
      return "M0_BUNDLE_SELECTED_ORDER";
    case AnchoringRefusal::BundleSelectedSpan:
      return "M0_BUNDLE_SELECTED_SPAN";
    case AnchoringRefusal::BundleHash: return "M0_BUNDLE_HASH";
  }
  return "";
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** query_span ****/
// Forward-frame query spans and minimap2's mask_level competition test (mm_set_parent
// with uncov_len = 0), shared by the fine chain's sibling scan and the chain MAPQ.

namespace fa::cpu::lr::rna {

// minimap2's mask_level: a rival overlapping the committed hypothesis over
// more than half the shorter query span competes with it.
inline constexpr double kRnaChainMapqMaskLevel = 0.5;

// A query span in its chain's own mapping frame: the read for a forward chain, its
// reverse complement for a reverse one. Spans from opposite strands are comparable only
// after rna_forward_query_span.
struct RnaQuerySpan {
  int begin = -1;
  int end = -1;
};

// minimap2's mask_level test on forward-frame spans: an invalid span on either side
// competes (censored evidence is not disjoint evidence), a non-positive overlap does not,
// otherwise the overlap must exceed half the shorter span.
bool rna_query_spans_compete(RnaQuerySpan winner,
                             RnaQuerySpan rival) noexcept;

// Projects a mapping-oriented span into the read's forward frame; a reverse span mirrors
// to [read_len - end, read_len - begin). An invalid span, or a reverse span leaving the
// read, projects to the invalid span.
RnaQuerySpan rna_forward_query_span(int oriented_begin, int oriented_end,
                                    bool reverse, int read_len) noexcept;

// Rival query bases the winner's span does not cover; 0 when either span is invalid.
int rna_uncovered_query_bases(RnaQuerySpan winner, RnaQuerySpan rival) noexcept;

} // namespace fa::cpu::lr::rna

/**** placement/types ****/
// Value types produced by RNA coarse locus placement.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

struct CoarseDiagonalPeak {
  int reference_id = -1;
  bool reverse = false;
  std::uint32_t oriented_query_begin = 0;
  std::uint32_t oriented_query_end = 0;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  std::int64_t median_diagonal = 0;
  std::uint32_t query_covered_bases = 0;
  std::uint32_t distinct_seed_support = 0;
};

struct CoarseLocus {
  std::vector<std::uint32_t> peak_indices;
  // Bounded intron-aware colinear score over the envelope's peaks; the RNA MAPQ reads it
  // as coarse confidence.
  std::int64_t score = 0;
  // Catalogue order key: how much of the read the envelope explains with seed matches.
  std::int64_t rank_score = 0;
  int reference_id = -1;
  bool reverse = false;
  std::uint32_t oriented_query_begin = 0;
  std::uint32_t oriented_query_end = 0;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  std::uint32_t intron_count = 0;
  // Union of the peaks' oriented query intervals, in read bases; a selector tie-break.
  std::int64_t union_query_coverage_bases = 0;
};

}  // namespace placement
}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** splice_anchor ****/
// The selected RNA anchor path as realization consumes it.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum SpliceAnchorFlag : std::uint8_t {
  kSpliceAnchorLongJoin = 1u << 0,
  kSpliceAnchorIgnore = 1u << 1,
  kSpliceAnchorTandem = 1u << 2,
};

struct SpliceAnchor {
  std::uint32_t stable_id = 0;
  std::int32_t query_begin = 0;
  std::int32_t reference_begin = 0;
  std::int32_t span = 0;
  std::uint8_t flags = 0;
};

struct SpliceAnchorView {
  std::vector<SpliceAnchor> anchors;
  int selected_begin = 0;
  int selected_count = 0;
  bool mapping_reverse = false;
  std::vector<std::uint32_t> view_to_raw_index;
  std::vector<std::uint32_t> excluded_interior_raw_indices;
  std::uint64_t content_hash = 0;
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** rival_lifecycle_config ****/
// Rival lifecycle options, kept apart so the RNA context can carry them without
// including the realizer.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Defaults follow minimap2: pri_ratio and min_diff are its secondary retention band.
// minimap2's min_diff is 2k; 30 is its value at k = 15 and is not rescaled with the
// index's k.

struct RnaRivalLifecycleConfig {
  // -p / --rival-min-diff: a competitor is admitted when its chain score is within
  // either band of the primary's (inclusive).
  double pri_ratio = 0.8;
  int min_diff = 30;
  // -N: bound on rival realizations per read (minimap2's best_n). The scan
  // itself covers the whole catalogue.
  int realize_max = 5;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** known_junctions ****/

namespace fa::cpu::lr::rna {

struct KnownIntron {
  std::uint32_t start = 0;
  std::uint32_t end = 0;
  std::uint8_t strand_mask = 0; // 1 plus, 2 minus
};

struct KnownJunctionLoadCounters {
  std::uint64_t rows = 0;
  std::uint64_t introns = 0;
  std::uint64_t duplicates_merged = 0;
  std::uint64_t single_block_rows = 0;
};

struct KnownJunctionContigView {
  const KnownIntron* begin = nullptr;
  const KnownIntron* end = nullptr;
  bool empty() const noexcept { return begin == end; }
};

class KnownJunctionStore {
public:
  static KnownJunctionStore
  load_bed(const std::string& path, const std::vector<std::string>& names,
           const std::vector<std::vector<std::uint8_t>>& references);
  KnownJunctionContigView contig(std::size_t id) const noexcept;
  const KnownJunctionLoadCounters& counters() const noexcept {
    return counters_;
  }

private:
  std::vector<std::vector<KnownIntron>> contigs_;
  KnownJunctionLoadCounters counters_;
};

std::uint64_t
materialize_known_junction_mask(KnownJunctionContigView junctions,
                                int reference_begin, int reference_end,
                                std::vector<std::uint8_t>& output);

} // namespace fa::cpu::lr::rna

/**** skeleton_regions ****/
// Rank-1 search regions. Coordinates are zero-based, half-open, chromosome-local and in
// the selected strand's query frame.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class SearchRegionKind : uint8_t { kCore = 0 };

struct SkeletonNode {
  uint32_t node_id = 0;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int64_t reference_begin = 0;
  int64_t reference_end = 0;
  int64_t median_diagonal = 0;
  uint32_t support = 0;
};

struct SearchRegion {
  uint32_t region_id = 0;
  SearchRegionKind kind = SearchRegionKind::kCore;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int64_t reference_begin = 0;
  int64_t reference_end = 0;
  int64_t diagonal_low = 0;
  int64_t diagonal_high = 0;  // inclusive
  uint32_t per_region_budget = 0;
};

struct SkeletonRegionParams {
  // Sanity limit on one gene's genomic footprint; posting enumeration is restricted to
  // this interval.
  int64_t maximum_locus_reference_span = 4000000;
  uint32_t locus_budget = 262144;
};

struct SearchRegionBuildResult {
  std::vector<SearchRegion> regions;
  uint64_t total_volume = 0;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** skeleton_harvest ****/
// Fine-anchor candidate pools for the RNA locus harvest: value types, the posting memo
// and the harvest entry points.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class CandidatePoolKind : uint8_t {
  kSkeletonUnionUncapped = 1,
  kSkeletonUnionBounded = 2,
};

// 3 is unused; the values are part of the pool's provenance record, so they keep their
// numbers.
enum class FineAnchorDropReason : uint8_t {
  kNone = 0,
  kOccurrence = 1,
  kGeometry = 2,
  kFeatureCap = 4,
  kRegionCap = 5,
  kGlobalCap = 6,
};

struct FineAnchor {
  uint32_t anchor_id = 0;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int32_t reference_begin = 0;
  int32_t reference_end = 0;
  int32_t span = 0;
  uint32_t source_seed_id = 0;
  uint32_t source_occurrence = 0;
  uint32_t evidence_weight = 0;
  std::vector<uint32_t> region_ids;
  std::vector<uint32_t> compatible_core_nodes;
  uint32_t region_kind_mask = 0;
  uint32_t raw_multiplicity = 0;
  // Realization flags carried by the pool; region tags never set them.
  bool ignore_flag = false;
  bool tandem_flag = false;
  bool long_join_flag = false;
  bool raw_seen = false;
  bool survived_occurrence = false;
  bool survived_geometry = false;
  bool survived_feature_cap = false;
  bool survived_region_cap = false;
  bool survived_global_cap = false;
  bool final_pool = false;
  FineAnchorDropReason drop_reason = FineAnchorDropReason::kNone;
};

struct CandidatePoolLedger {
  uint64_t raw_seen = 0;
  uint64_t occurrence_filtered = 0;
  uint64_t geometry_filtered = 0;
  uint64_t feature_capped = 0;
  uint64_t region_capped = 0;
  uint64_t global_capped = 0;
  uint64_t final_pool = 0;
};

struct CandidatePool {
  CandidatePoolKind kind = CandidatePoolKind::kSkeletonUnionUncapped;
  std::vector<FineAnchor> anchors; // includes filtered/capped records
  CandidatePoolLedger ledger;
  uint64_t content_hash = 0;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

struct SkeletonHarvestParams {
  // Skip enumerating the postings of a key the occurrence cap will reject anyway. The cap
  // reads the key's genome-wide count, which lookup_interval already reports, so those
  // postings would all fail survived_occurrence. They stay counted in the ledger, the
  // safety limit and the anchor identifiers, so the pool, ledger, hash and refusals match
  // full enumeration exactly; false enumerates every posting.
  bool pre_enumeration_occurrence_skip = true;
  // Read against the key's genome-wide count, so a key common anywhere is rejected even
  // where it is unique.
  uint32_t occurrence_cap = 200;
  uint32_t per_feature_budget = 16;
  uint32_t global_budget = 32768;
  // Refuse a seed whose canonical k-mer is low complexity before its postings are read:
  // >= 11/15 of one base, or <= 4 distinct 3-mers among the 13 overlapping ones (scaled
  // from k = 15). Under the log-priced intron transition a lone homopolymer anchor
  // anywhere within max_intron pays for itself, and such k-mers occur everywhere, so the
  // occurrence cap alone does not stop them. A masked seed is bookkept like an absent
  // key; the mirrored frame carries the same keys and masks the same seeds.
  bool mask_low_complexity = true;
  uint64_t raw_candidate_safety_limit = 250000;
  int64_t core_diagonal_tolerance = 64;
  uint32_t minimum_core_evidence_weight = 1;

  // Strand routing: every inspected posting goes to the lane its stored orientation bit
  // is compatible with (packed_ref_orientation_compatible against frame_reverse_lane or
  // its negation); exactly one holds, since no palindromic key is a seed. One raw
  // enumeration yields the nominated pool, in this frame, and the opposite pool,
  // mirrored into the other lane's coordinates. `frame_reverse_lane` is the lane
  // `strand_seeds` belongs to: false forward, true reverse complement.
  bool frame_reverse_lane = false;
  // The read's seeds re-admitted past the cap (seed/seed.h mark_readmitted_seeds),
  // one byte per forward read position, or null. Such a seed passes the occurrence
  // cap in both lanes.
  const std::uint8_t* readmitted = nullptr;
};

// Strand-routing work counters. `postings_bulk_skipped` restates
// HarvestWorkCounters::postings_skipped_occurrence: those postings are never read, so
// their orientation is unknown and they enter neither pool.
struct RoutingCounters {
  uint64_t postings_inspected = 0;
  uint64_t postings_bulk_skipped = 0;
  uint64_t routed_nominated = 0;
  uint64_t routed_opposite = 0;
  uint64_t admitted_nominated = 0;
  uint64_t admitted_opposite = 0;
};

// Harvest work counters, kept out of CandidatePool so they never affect the pool's
// ledger or content hash. Postings the harvest declines to enumerate still
// count as available:
// postings_available == postings_enumerated + postings_skipped_occurrence.
struct HarvestWorkCounters {
  uint64_t seeds_visited = 0;
  uint64_t seeds_without_region = 0;  // no applicable search region
  uint64_t seeds_key_absent = 0;      // key not in the index at all
  uint64_t seeds_over_occurrence = 0; // key present, over the occurrence cap
  // Seeds refused by the low-complexity mask; disjoint from the counters below.
  uint64_t seeds_masked_low_complexity = 0;
  uint64_t interval_queries = 0;       // lookup_interval() calls
  uint64_t interval_queries_empty = 0; // calls returning no posting in-window
  uint64_t postings_available = 0;     // postings inside the queried intervals
  uint64_t postings_enumerated = 0;    // postings whose position was read
  // `over_occurrence`: available postings in intervals the occurrence cap rejects;
  // `skipped_occurrence`: those the harvest declined to enumerate (equal when the
  // pre-enumeration skip is on).
  uint64_t postings_over_occurrence = 0;
  uint64_t postings_skipped_occurrence = 0;
  uint64_t canonical_anchors = 0; // distinct (query, reference) records
  uint64_t regions_queried = 0;   // sum over seeds of applicable regions
  uint64_t merged_intervals = 0;  // sum over seeds of merged intervals
};

// Per-read memo of full posting views, shared by both strand streams. One read runs
// several harvests over its two seed vectors, each probing every key; lookup() is pure
// over an immutable index, so caching it changes only CPU. The reverse stream is the
// forward one walked backwards with the same canonical keys, so reverse seed j is
// forward seed n-1-j and storage is in the forward frame. `views` borrow the index's
// posting array; a memo is valid for one index.
//
// Two guards, because a wrong view is silent corruption: bind_reverse refuses a reverse
// stream whose size differs from the forward one's (a dropped seed breaks n-1-j), and
// every hit re-checks the stored key, falling back to a real probe and counting
// `mirror_key_mismatches`, which must stay zero. Freshness is a generation stamp that
// new_read() bumps in O(1); the buffers are never re-zeroed per read.
struct SeedPostingMemo {
  static constexpr std::size_t kNoSlot = static_cast<std::size_t>(-1);

  // Storage frame: the forward stream.
  const void* forward_identity = nullptr;
  std::size_t forward_size = 0;
  // Mirrored frame, bound only when it is provably a mirror of the forward one.
  const void* reverse_identity = nullptr;
  std::size_t reverse_size = 0;

  // 0 means "no read in progress": nothing hits and nothing is stored, so a
  // memo that never saw new_read() degrades to plain lookups.
  std::uint32_t generation = 0;
  std::vector<KmerPostingView> views;
  std::vector<std::uint64_t> keys;
  std::vector<std::uint32_t> stamps;

  // Guard 2's tripwire. Always zero unless the mirroring premise is false.
  std::uint64_t mirror_key_mismatches = 0;

  // Invalidate every slot in O(1). Called once per read by the caller.
  void new_read();

  // Install the forward stream and size the buffers. Idempotent. Buffers only
  // ever grow: shrinking would cost the O(n) write the stamps exist to avoid.
  void bind_forward(const std::vector<QuerySeed>& seeds);

  // Install the reverse stream, but only if it is provably forward's mirror.
  // Idempotent and optional; an unbound reverse stream is simply not memoized.
  void bind_reverse(const std::vector<QuerySeed>& seeds);

  // Which forward-frame slot backs `stream`'s seed `i`, or kNoSlot when this
  // stream is not a bound frame.
  std::size_t slot_for(const std::vector<QuerySeed>& stream,
                       std::size_t i) const;

  // The seed's full posting view, probing `index` only on a miss.
  KmerPostingView resolve(const SeedIndex& index,
                          const std::vector<QuerySeed>& stream, std::size_t i);
};

// A block of anchors full enumeration would have created and dropped on occurrence. The
// harvest visits seeds in ascending query order with distinct query starts, so a skipped
// seed's anchors occupy one contiguous identifier run. The block keeps the seed's query
// interval so the frame projection can mirror it (q -> length - q_end).
struct SkippedBlock {
  int32_t query_begin = 0;
  int32_t query_end = 0;
  uint64_t count = 0;
};

struct SkeletonCandidatePools {
  CandidatePool uncapped;
  CandidatePool bounded;
  // `nominated` is in the harvested frame's coordinates and carries the skeleton's core
  // tags. `opposite` is in the other lane's coordinates (collected from the raw postings,
  // then mirrored) and carries none. Each is capped on its own. `opposite` is refused
  // with MirrorProjectionDomain when the window has more than one search region.
  CandidatePool nominated;
  CandidatePool opposite;
  bool routed = false;
  RoutingCounters routing;
  HarvestWorkCounters work;
  // Ascending in query; kept so the frame projection can replay the identifier walk.
  std::vector<SkippedBlock> skipped_blocks;
  // The same blocks in the opposite lane's coordinates, produced by the same
  // mirror that produced `opposite`.
  std::vector<SkippedBlock> opposite_skipped_blocks;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** exact_anchor_path ****/
// Immutable handoff from the RNA exact-anchor-path chain to realization. Coordinates are
// in the mapping-oriented query frame: the read for forward mappings, its reverse
// complement for reverse ones. Only final output intervals are mirrored back.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct SelectedPoolAnchor {
  uint32_t stable_id = 0;
  int32_t query_begin = 0;
  int32_t reference_begin = 0;
  int32_t span = 0;
  uint32_t source_seed_id = 0;
  uint32_t source_occurrence = 0;
  bool ignore_flag = false;
  bool tandem_flag = false;
  bool long_join_flag = false;
  // Provenance: all three fields enter anchor_path_hash and are checked by validation.
  uint32_t region_kind_mask = 0;
  std::vector<uint32_t> region_ids;
  std::vector<uint32_t> core_node_ids;
};

struct ExactAnchorPath {
  static constexpr uint32_t kSchemaVersion = 1;

  uint32_t schema_version = kSchemaVersion;
  uint64_t read_key = 0;
  int32_t query_length = 0;
  int32_t reference_id = -1;
  bool mapping_reverse = false;
  std::vector<SelectedPoolAnchor> selected_pool_anchors;
  std::vector<uint32_t> selected_raw_indices;
  int32_t anchor_path_score = 0;
  uint32_t anchor_path_anchor_count = 0;
  uint64_t candidate_pool_hash = 0;
  uint64_t anchor_path_hash = 0;
};

struct ExactAnchorPathResult {
  ExactAnchorPath bundle;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

// Endpoint facts about the selected path, computed with validation of the selected
// indices and without building chaining::Anchor copies.
struct ExactAnchorPathSummary {
  uint64_t selected_count = 0;
  int32_t first_query_begin = 0;
  int32_t last_query_end = 0;
  uint64_t outer_query_span = 0;
  int32_t first_reference_begin = 0;
  int32_t last_reference_end = 0;
  uint64_t outer_reference_span = 0;
};

bool validate_exact_anchor_path(const ExactAnchorPath& bundle,
                                AnchoringRefusal* refusal = nullptr);

bool summarize_exact_anchor_path(const ExactAnchorPath& bundle,
                                 ExactAnchorPathSummary& summary,
                                 AnchoringRefusal* refusal = nullptr) noexcept;

// Query bases covered by the union of the selected anchors' query intervals, read by the
// PAF matches field and the MAPQ coverage brake. False, with covered_bases 0, when an
// interval leaves the query or the union is empty.
bool anchor_union_query_coverage(const ExactAnchorPath& bundle,
                                 uint64_t& covered_bases) noexcept;

// Stable read identity: FNV-1a over the upper-cased read sequence.
uint64_t hash_rna_read_key(const std::string& read) noexcept;

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** chain_mapq ****/
// The specification of the RNA MAPQ is the comment under hit.cpp's chain_mapq banner.

namespace fa::cpu::lr::rna {

// minimap2's q_coef.
inline constexpr double kRnaChainMapqCoef = 40.0;
// The spliced chain's min_sc, standing in for minimap2's min_chain_score in the subsc
// floor.
inline constexpr int kRnaChainMapqMinChainScore = 20;
// Default tau of the query-coverage brake: pen is discounted linearly below an
// anchor-union coverage of tau * read_len. Presets set RnaChainMapqEvidence::qcov_tau.
inline constexpr double kRnaChainMapqQcovTau = 0.7;
// Scale on minimap2's 6.02 in the near-tie cap; the cap's shape, operands and order
// stay minimap2's.
inline constexpr double kRnaChainMapqBwaCapScale = 0.6;

// One catalogue rival of the committed hypothesis. A chained rival carries its fine-chain
// measurement; an unchained one only its vote and screening span (censored evidence is
// not zero evidence).
struct RnaChainMapqRival {
  // Catalogue index; not used by the formula.
  int candidate = -1;
  // Fine-chain measurement (chained rivals only).
  int chain_score = 0;
  int chain_anchors = 0;
  // The rival's chain query span in its own mapping frame (see `reverse`); -1 when
  // censored.
  int q_begin = -1;
  int q_end = -1;
  // The coarse catalogue's rank_score (the RNA vote evidence).
  std::int64_t rank_score = 0;
  bool chained = false;
  // The strand the rival's chain was taken on, the frame of q_begin/q_end: the
  // post-repair locus's `reverse` (a repaired rival's locus is already flipped, so do not
  // XOR again). Not read by rna_chain_mapq.
  bool reverse = false;
  // This rival's realized DP score is evidence.dp2; lets n_sub keep minimap2's per-rival
  // OR.
  bool owns_dp2 = false;
  // Realized with a positive DP maximum, so rho already bounds this rival and it stays
  // out of the undamped x_floor. Several rivals may be realized while one owns dp2; the
  // formula reads `realized_dp || owns_dp2`.
  bool realized_dp = false;
  // This rival's fine chain is the committed chain's own placement seen from a
  // neighbouring catalogue window (same contig and strand, overlapping or continuing
  // reference spans; rna_same_place). The opposite-lane postings of the winner's envelope
  // can form a small window of their own whose harvest re-chains the committed anchors;
  // priced as a rival that is an exact realized tie and a false Q0. A shadow stays in the
  // lifecycle but never supplies dp2 and is inadmissible to the formula (x, x_floor, f2,
  // n_sub, tie). Set in map.cpp, hit.cpp and leftover.cpp against the committed and the
  // elected chain.
  bool shadow = false;
  // Under RnaChainMapqEvidence::parent_rule: a region of this rival was a child of the
  // priced region (rna.h).
  bool parented = false;
};

struct RnaChainMapqEvidence {
  int f1 = 0;  // committed fine-chain score            (mm2 score0 / score)
  int cnt = 0; // committed fine-chain anchors          (mm2 cnt)
  // The committed window's second-best chain from the same recurrence, a free lower
  // bound on subsc; <= 0 when absent.
  int sib_f2 = 0;
  int winner_q_begin = -1;
  int winner_q_end = -1;
  // The committed locus's coarse rank_score (rs1), the vote denominator.
  std::int64_t rank1_score = 0;
  // rs2: the largest rank_score over catalogue loci other than the committed one (on an
  // adopted read, the demoted incumbent's). Taken from the catalogue rather than
  // `rivals` so the floor v counts a locus that was never chained. 0 when the catalogue
  // held no other locus.
  std::int64_t rank2_score = 0;
  int read_len = 0;
  // The committed hypothesis's strand, with RnaChainMapqRival::reverse's semantics; the
  // frame of winner_q_begin/winner_q_end. Not read by the formula.
  bool winner_reverse = false;
  // Anchor-union query bases covered by the committed fine chain, the brake's
  // numerator (the map-only PAF matches field).
  std::uint64_t selected_query_covered_bases = 0;
  // The brake's denominator: 0 means read_len (the primary); a piece passes its region's
  // query span.
  int qcov_bases = 0;
  std::vector<RnaChainMapqRival> rivals;
  // Realized DP maxima with CIGAR output; 0 for map-only.
  double dp1 = 0.0;
  double dp2 = 0.0;
  // matches/block_len of the realized primary; 1.0 for map-only output, which makes the
  // identity factor inert.
  double identity = 1.0;
  int match_sc = 1; // minimap2 opt->a (cigar_dp_match; 1 in both presets)
  // minimap2 opt->a * 2 + opt->b: 4 (splice) / 6 (splice:hq). The threshold
  // under which the realized rival still counts toward n_sub.
  int sub_diff = 0;
  // With CIGAR output, which chains every catalogue rival, an unchained (vote-only)
  // rival contributes x strength 0 instead of its vote ratio. False for map-only
  // output, which chains no rival. x is deliberately not clamped at 1: x > 1 from a
  // demoted incumbent brakes an adopted read.
  bool unchained_zero = false;
  // The preset calibration values travel with the evidence so the formula stays a pure
  // function of this struct; map.cpp and hit.cpp fill them from ResolvedOptions.
  //
  // Coverage brake threshold, in (0, 1]: pen's third limb is
  // clip((qcov/read_len)/qcov_tau, 0, 1).
  double qcov_tau = kRnaChainMapqQcovTau;
  // Coefficient of the damped disjoint-vote floor v, in [0, 1]; 0 disables it with no
  // float operation.
  double disjoint_vote_damp = 0.0;
  // minimap2's parent rule in place of the span test: a chained rival is admissible
  // only when `parented`, and n_sub is parent_n_sub, the priced region's own count.
  bool parent_rule = false;
  int parent_n_sub = 0;
};

// The formula's verdict on one rival, for the caller.
struct RnaChainMapqRivalVerdict {
  // minimap2's mask_level test on the query spans.
  bool admissible = false;
  // Chained exact tie (f_i == f1 && cnt_i == cnt): undecidable, x = 1.
  bool tie = false;
  // The rival's contribution to x (0 when it contributes nothing).
  double strength = 0.0;
  // The rival counted toward n_sub. The standalone verdict applies only the
  // anchor rule; rna_chain_mapq adds minimap2's DP rule for the dp2 owner.
  bool counted = false;
};

// Every intermediate the formula forms. Filled only when the caller passes a
// non-null pointer.
struct RnaChainMapqBreakdown {
  int rivals_total = 0;
  int rivals_admissible = 0;
  int rivals_tie = 0;
  int f2 = 0; // best rival or sibling chain score
  int n_sub = 0;
  double vote_ratio = 0.0; // largest rs_i / rs1 over the admissible rivals
  double x_used = 0.0;     // base x, before branch C's damping
  double pen = 0.0;        // after every discount, as consumed
  double qcov_brake = 1.0; // the brake term alone, before the min
  double identity = 1.0;   // as consumed (1.0 for map-only)
  bool dp_branch = false;  // the emitted q read ln(dp1) instead of ln(f1)
  // `branch_c`: the emitted q used (1 - x_c^2). rho = dp2/dp1, x_floor the undamped
  // floor over the unrealized rivals (sibling included), x_c = max(rho * x_used,
  // x_floor). Off branch C: rho = 0, x_floor = the sibling share, x_c = x_used.
  bool branch_c = false;
  double rho = 0.0;
  double x_floor = 0.0;
  double x_c = 0.0;
  // minimap2's BWA-MEM-like near-tie cap, or -1 when not evaluated (no realized dp1/dp2
  // pair); `bwa_cap_bound` says whether it lowered q.
  int bwa_cap = -1;
  bool bwa_cap_bound = false;
};

// minimap2's mask_level test (uncov_len = 0) on the raw oriented spans: a rival competes
// when it overlaps the committed hypothesis over more than half the shorter span. A
// missing span on either side is censored geometry and stays admissible; a shadow never
// is. The forward-frame helpers in rna.h (query_span) serve the lifecycle, the chimeric floor
// and rival pricing.
bool rna_chain_mapq_admissible(const RnaChainMapqEvidence& evidence,
                               const RnaChainMapqRival& rival) noexcept;

// The RNA MAPQ for the committed hypothesis, in [0,60]. Returns 0 when there
// is no committed fine chain (f1 <= 0). `verdicts`, when non-null, receives
// one entry per evidence.rivals element, in order.
int rna_chain_mapq(const RnaChainMapqEvidence& evidence,
                   RnaChainMapqBreakdown* breakdown = nullptr,
                   std::vector<RnaChainMapqRivalVerdict>* verdicts = nullptr);

} // namespace fa::cpu::lr::rna

/**** splice_scratch ****/
// Caller-owned per-worker realization scratch, like minimap2's mm_tbuf_t: grown on demand,
// never shrunk, and reused across packets, segments, hypotheses and reads. A null pointer
// keeps every buffer function-local.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// One ignore/long-join bit set the cluster filters added to one anchor.
struct SpliceAnchorFlagMutation {
  std::uint32_t index = 0;
  std::uint8_t added = 0;
};

// One segment's orientation-blind preparation: the terminal-anchor filter, the two
// cluster filters, the pivots and the terminal envelope. None of it reads the transcript
// orientation or the flags the cluster filters add, so the second hypothesis replays the
// first's.
struct SplicePreparedSegment {
  int selected_begin = 0;
  int selected_count = 0;
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;
  // The segment-result fields the preparation wrote before any refusal, so a replay
  // writes the same ones.
  bool wrote_filtered = false;
  bool wrote_envelope = false;
  int filtered_begin = 0;
  int filtered_count = 0;
  int first_reference_pivot = 0;
  int first_query_pivot = 0;
  int last_reference_pivot = 0;
  int last_query_pivot = 0;
  int envelope_reference_begin = 0;
  int envelope_query_begin = 0;
  int envelope_reference_end = 0;
  int envelope_query_end = 0;
  int left_neighbor = -1;
  int right_neighbor = -1;
  std::vector<SpliceAnchorFlagMutation> anchor_flags;
};

struct SpliceRealizationScratch {
  // ksw2 arena and the reused ez.cigar.
  dp::SpliceKernelWorkspace kernel;
  // The kernel's last result, reused so its CIGAR keeps its capacity.
  dp::SpliceKernelResult aligned;
  // One packet's known-junction mask.
  std::vector<std::uint8_t> junction;
  // The left packet's reversed query and reference flanks.
  std::vector<std::uint8_t> reversed_query;
  std::vector<std::uint8_t> reversed_reference;
  // The hypothesis' mutable copy of the frozen anchor view.
  std::vector<SpliceAnchor> working;
  // Segments prepared during this realization, keyed by selected range. Anchor geometry
  // is fixed for a realization, so entries stay valid for later segments and the second
  // hypothesis; realize_exact_anchor_path clears them for each new anchor path.
  std::vector<SplicePreparedSegment> prepared;
  // Flag bytes of the selected range before the cluster filters ran.
  std::vector<std::uint8_t> flags_before;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** splice_controller ****/
// Splice realization of one selected anchor segment. The controller works on a mutable
// copy of the anchor view and holds no global state; the DP kernel is dp/dp.h.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The transcript hypothesis a realization pass scores. None is minimap2's -un: no splice
// motif, so there is a single orientation-free pass.
enum class TranscriptOrientation : std::uint8_t {
  Forward = 1,
  Reverse = 2,
  None = 3,
};

enum class SplicePacketRole : std::uint8_t {
  Left = 1,
  Internal = 2,
  Right = 3,
};

struct SpliceControllerOptions {
  int k = 15;
  int match = 1;
  int mismatch = 2;
  int gap_open = 2;
  int gap_extend = 1;
  int long_gap_open = 32;
  int long_gap_extend = 0;
  int transition = 0;
  int ambiguous = 1;
  int junction_bonus = 9;
  int junction_penalty = 5;
  int zdrop = 200;
  int inversion_zdrop = 100;
  int end_bonus = -1;
  int maximum_gap = 2000;
  int maximum_reference_gap = 200000;
  int minimum_anchor_count = 3;
  int minimum_match_bases = 40;
  int minimum_dp_maximum = 80;
  int minimum_packet_length = 200;
  // Shortest reference gap called an intron, used only by the packet-shape accounting;
  // the kernel derives its own threshold from the long-gap scores.
  int minimum_intron = 20;
  int anchor_extension_length = 20;
  int anchor_extension_shift = 6;
  std::int64_t maximum_dp_cells = 0;
};

SpliceControllerOptions ordinary_long_read_splice_options() noexcept;
bool splice_controller_options_supported(
    const SpliceControllerOptions& options) noexcept;

struct SpliceControllerRequest {
  const std::vector<std::uint8_t>* query_forward = nullptr;
  const std::vector<std::uint8_t>* query_reverse = nullptr;
  const std::vector<std::uint8_t>* reference = nullptr;
  std::vector<SpliceAnchor>* working_anchors = nullptr;
  int selected_begin = 0;
  int selected_count = 0;
  int chain_score = 0;
  int initial_chain_score = 0;
  bool mapping_reverse = false;
  TranscriptOrientation transcript_orientation =
      TranscriptOrientation::Forward;
  SpliceControllerOptions options = ordinary_long_read_splice_options();
  KnownJunctionContigView known_junctions;
  // Optional per-worker scratch. Null keeps every buffer function-local.
  SpliceRealizationScratch* scratch = nullptr;
};

struct SpliceContinuation {
  int selected_begin = 0;
  int selected_count = 0;
  int chain_score = 0;
  int initial_chain_score = 0;

  bool present() const noexcept { return selected_count > 0; }
};

struct SpliceSegmentResult {
  bool refused = false;
  RnaControllerRefusal refusal =
      RnaControllerRefusal::None;

  int filtered_begin = 0;
  int filtered_count = 0;
  int first_reference_pivot = 0;
  int first_query_pivot = 0;
  int last_reference_pivot = 0;
  int last_query_pivot = 0;
  int envelope_reference_begin = 0;
  int envelope_query_begin = 0;
  int envelope_reference_end = 0;
  int envelope_query_end = 0;
  int left_neighbor_index = -1;
  int right_neighbor_index = -1;

  int query_begin = 0;
  int query_end = 0;
  int reference_begin = 0;
  int reference_end = 0;
  int region_anchor_count = 0;
  int chain_score = 0;
  int matches = 0;
  int block_length = 0;
  bool is_spliced = false;
  int dp_score = 0;
  int dp_maximum = 0;
  int dp_maximum_before_strand = 0;
  int ambiguous_bases = 0;
  std::vector<std::uint32_t> cigar;

  SpliceContinuation continuation;
  std::uint32_t split_flags = 0;
  RnaControllerMetrics metrics;
};

SpliceSegmentResult realize_splice_segment(
    const SpliceControllerRequest& request) noexcept;

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** context ****/
// Dependency bundle handed to RnaBackend::map_read (the RNA counterpart of dna/dna.h):
// reference views plus an options snapshot composed once per configure.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Execution snapshot composed once from the engine's resolved options.
struct ResolvedOptions {
  // The RNA config, with the entry-point clamps applied.
  rna::RnaConfig cfg;
  // The resolved splice-controller scoring row (defaults: minimap2's `splice`), used
  // unchanged by every hypothesis and continuation.
  SpliceControllerOptions splice_controller =
      ordinary_long_read_splice_options();

  // Rank-1 harvest and exact-anchor-path chain options: the fine chain's query gap, -g,
  // and its band, -r's first value.
  int cigar_local_global_occ = 200;
  int fine_chain_query_gap = 5000;
  int fine_chain_band = 200000;
  // The vote's seeds per strand on a read of length L: max(the seed context's,
  // ceil(L / vote_seed_spacing)); 0 keeps the seed context's.
  int vote_seed_spacing = 0;
  // The within-read repeat filter's fraction (--q-occ-frac); 0 is off.
  float query_occ_frac = 0.0f;
  // minimap2's occ_dist (-e): an over-cap streak re-admits one seed per this many
  // bases; 0 is off.
  int occ_dist = 0;

  // Rival realization lifecycle options; the lifecycle is the requested-CIGAR path's
  // only arbitration.
  RnaRivalLifecycleConfig rival_lifecycle;
  // The MAPQ formula's two preset calibration constants (rna.h), copied onto
  // every RnaChainMapqEvidence by map.cpp and hit.cpp. The defaults reproduce the formula when
  // nothing was composed.
  double mapq_qcov_tau = kRnaChainMapqQcovTau;
  double mapq_disjoint_vote_damp = 0.0;
  // Gates the placement-only return versus splice realization in map.cpp.
  bool enable_full_read_cigar = true;
  // cs:Z / MD:Z request (minimap2 --cs / --MD), copied onto every realization request.
  ::fa::cpu::output::CigarReplayRequest cigar_replay_request;
  std::shared_ptr<const KnownJunctionStore> known_junctions;
};

// Non-owning views plus the options snapshot, composed once per config; the engine owns
// the storage.
struct Context {
  engine::ReferenceContext ref;
  ResolvedOptions opts;
  // Computed once when the RNA runtime is composed. Per-read code refuses a reference
  // whose counts or contig lengths do not fit the signed coordinate domain.
  mapping::CoordinateDomainRefusal coordinate_domain_refusal =
      mapping::CoordinateDomainRefusal::None;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** fused_capture ****/
// Single-pass two-strand RNA vote capture. The reverse seed stream is the forward stream
// with projected query coordinates and the same canonical keys, so each forward seed's
// posting view is resolved once and reused by both strands. Postings append to a flat
// per-strand log while a reused counter tracks each bin's distinct-seed support, so bins
// below min_support are dropped before any per-bin container exists. Peaks are emitted
// in ascending bin key, forward strand first.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// One captured posting with its diagonal-bin key and the bin's stable ordinal in the
// read's DiagonalBinCounter, so the drain reads bin support by index.
struct CapturedVote {
  std::int64_t key = 0;
  VoteHit hit{};
  std::uint32_t bin = 0;
};

// Reused open-addressing counter of distinct-seed support and postings per bin key. Keys
// are vote_pack_key values (>= 0), so -1 marks an empty slot; reset() once per read.
// Slots map key -> bin ordinal and the counts live in `bins_`. Ordinals are assigned in
// first-touch order and survive grow(), so add()'s return value can be stored.
class DiagonalBinCounter {
public:
  void reset(std::size_t expected_bins);
  // `first_posting_of_seed_in_bin` counts a new distinct seed; a seed's postings in one
  // bin are contiguous, so the caller passes key != previous key. Returns the key's
  // stable bin ordinal.
  std::uint32_t add(std::int64_t key, bool first_posting_of_seed_in_bin);
  std::uint32_t support(std::int64_t key) const;
  // support(key) by bin ordinal, with no hash probe.
  std::uint32_t support_at(std::uint32_t bin) const {
    return bins_[bin].support;
  }
  std::size_t bin_count() const { return bins_.size(); }
  std::uint64_t max_bin_postings() const { return max_bin_postings_; }

private:
  // Generation-stamped slots make reset() O(1): a slot is live only while its
  // generation matches the counter's.
  struct Slot {
    std::int64_t key = -1;
    std::uint32_t generation = 0;
    std::uint32_t bin = 0;
  };
  struct Bin {
    std::uint32_t support = 0;
    std::uint32_t postings = 0;
  };
  bool live(const Slot& slot) const { return slot.generation == generation_; }
  void grow();
  std::vector<Slot> slots_;
  std::vector<Bin> bins_;
  std::size_t mask_ = 0;
  std::uint32_t generation_ = 0;
  std::uint64_t max_bin_postings_ = 0;
};

struct FusedCaptureScratch {
  std::vector<CapturedVote> fwd_log;
  std::vector<CapturedVote> rc_log;
  DiagonalBinCounter fwd_bins;
  DiagonalBinCounter rc_bins;
  // One posting-list view per forward seed, resolved once and reused by the
  // per-strand occurrence evidence, the tile-spread selection, and the vote.
  std::vector<KmerPostingView> views;
  std::vector<::fa::cpu::QuerySeedOccurrenceEvidence> evidence;
  std::vector<QuerySeed> selected_fwd;
  std::vector<QuerySeed> selected_rc;
  // Aggregation scratch (survivor filter + per-bin peak math).
  std::vector<CapturedVote> survivors;
  std::vector<int> read_pos_scratch;
  std::vector<int> ref_start_scratch;
  std::vector<int> seed_id_scratch;
  // The aggregation's peak p has its distinct seeds' read positions, in its lane's
  // frame, at peak_read_pos[peak_read_pos_offset[p] .. peak_read_pos_offset[p + 1])
  // (select_coarse_loci's query_cover).
  std::vector<std::uint32_t> peak_read_pos_offset;
  std::vector<int> peak_read_pos;
};

// Post-selection seed counts per strand.
struct FusedCaptureSummary {
  std::size_t fwd_voted_seeds = 0;
  std::size_t rc_voted_seeds = 0;
  // Strand-compatibility counts per lane: `inspected` postings reached the test,
  // `compatible` passed it and were logged, `rejected` failed it.
  std::uint64_t fwd_postings_inspected = 0;
  std::uint64_t fwd_postings_compatible = 0;
  std::uint64_t fwd_postings_rejected_strand = 0;
  std::uint64_t rc_postings_inspected = 0;
  std::uint64_t rc_postings_compatible = 0;
  std::uint64_t rc_postings_rejected_strand = 0;
};

// Captures votes over the forward canonical seed stream and its reverse projection with
// one posting lookup per forward seed: each strand's tile-spread selection of
// max_seeds (0: every seed) is replayed from the cached evidence, then both
// strands vote under the per-seed gates (index hit, posting count, occurrence policy).
// `readmitted`, when non-null, is refilled by mark_readmitted_seeds (seed/seed.h) over
// the forward stream at the occurrence cap, one seed per `occ_dist` bases of an
// over-cap streak; a seed it marks passes the occurrence gate.
FusedCaptureSummary capture_votes_fused(
    const LongReadSeedContext& ctx, const std::vector<QuerySeed>& fwd_seeds,
    const std::vector<QuerySeed>& rc_seeds, int read_len, int max_seeds,
    FusedCaptureScratch& scratch,
    int occ_dist, std::vector<std::uint8_t>* readmitted);

// Aggregates the captured logs into diagonal peaks: forward-strand bins in ascending key
// order, then reverse-strand bins.
std::vector<CoarseDiagonalPeak>
aggregate_fused_capture_to_peaks(FusedCaptureScratch& scratch, int k,
                                 int min_support);

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** coarse_chain ****/
// Coarse locus-envelope selection over RNA diagonal peaks.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

struct CoarseLocusOptions {
  std::uint32_t diagonal_band = 24;
  std::uint32_t max_query_gap = 200;
  std::uint32_t max_query_overlap = 64;
  std::uint32_t min_intron = 20;
  std::uint32_t max_intron = 200000;
  // Reference gap across which two envelopes may still be coalesced: several
  // consecutive exons below the vote threshold can push the next peak beyond one legal
  // intron. Defaults to 4 * max_intron; 0 disables coalescing.
  std::uint32_t coalesce_gap = 800000;
  // Ceiling on a coalesced envelope's reference span, so the envelope plus the
  // +/- max_intron harvest pad fits skeleton_regions' 4 Mb locus bound. Defaults to
  // 16 * max_intron.
  std::uint32_t max_coalesced_span = 3200000;
  // At least one of two coalesced envelopes must carry this much distinct-seed support,
  // so two threshold-level peaks cannot fabricate a transcript. Set by the caller from
  // the vote threshold.
  std::uint32_t coalesce_min_support = 8;
  // Minimum oriented query extent of either piece of a coalesce (default 2 seed spans):
  // a narrower scrap is not an exon, and its few novel bases could lift a rival's
  // rank_score past the true locus.
  std::uint32_t coalesce_min_piece_span = 30;
  std::uint32_t max_predecessors = 64;
  std::uint32_t max_locus_chains = 6;
  // How many ranked envelopes beyond the cut the selector may additionally report in
  // `overflow_out`; 0 leaves the selection unchanged.
  std::uint32_t max_overflow_loci = 0;
  std::int32_t splice_open = 12;
  std::int32_t splice_log_scale = 2;
  std::int32_t query_gap_scale = 1;
  std::int32_t diagonal_drift_scale = 1;
  std::int32_t overlap_scale = 1;
};

CoarseLocusOptions make_coarse_locus_options(std::uint32_t seed_span,
                                             std::uint32_t min_intron,
                                             std::uint32_t max_intron);

// Caller-owned reusable buffers for one worker's coarse locus selection: cleared per
// call, never shrunk, holding no state any decision reads.
struct CoarseLocusScratch {
  // Per-call.
  std::vector<std::uint32_t> order;
  std::vector<CoarseLocus> envelopes;
  std::vector<std::uint64_t> support;
  std::vector<std::int64_t> coverage;
  // One flag per surviving envelope slot: did it absorb another envelope in the
  // coalescing pass? Only those slots pay the evidence recompute.
  std::vector<std::uint8_t> coalesce_absorbed;
  std::vector<std::uint32_t> ranked;
  // Per-envelope temporaries; distinct buffers, so no helper can alias another.
  std::vector<std::uint32_t> query_order;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> spans;
  std::vector<std::int64_t> chain_score;
  std::vector<std::uint32_t> by_reference;
  // Per read base: covered by a kept locus's seeds (query_cover admissions).
  std::vector<std::uint8_t> query_cover;
};

// Each peak's seeds, for select_coarse_loci's admissions past the cut.
struct CoarseLocusQueryCover {
  // Peak p's distinct seeds' read positions, in its lane's frame, are
  // read_pos[offset[p] .. offset[p + 1]).
  const std::vector<std::uint32_t>* offset = nullptr;
  const std::vector<int>* read_pos = nullptr;
  int read_len = 0;
  int k = 0;
  int min_new_bases = 0;
};

// Selects candidate locus envelopes for the downstream harvest. A reference-ordered sweep
// merges same-contig, same-strand peaks separated by at most one legal intron (covering
// one unpeaked exon); a second pass coalesces envelopes that are colinear pieces of one
// transcript within `coalesce_gap` (covering a run of unpeaked exons). The catalogue is
// ranked by seed evidence, then query coverage, support and peak count, then
// coarse_locus_tie_hash under `tie_seed`; contig slot and start settle only a hash
// collision. At most one envelope is kept per reference window, strand ignored.
//
// `coalesced`, when non-null, receives how many envelopes the coalescing absorbed.
// `overflow_out`, when non-null and max_overflow_loci > 0, receives the next ranked
// distinct envelopes past the cut, unfinished (no colinear score, no intron estimate,
// peaks in reference order); the kept prefix is unaffected.
// `query_cover`, when non-null and `overflow_out` is null, also keeps, past the cut and
// in rank order, each distinct envelope whose seeds' [pos, pos + k), in the forward
// frame, hold at least min_new_bases bases that no kept envelope's seeds cover.
std::vector<CoarseLocus> select_coarse_loci(
    const std::vector<CoarseDiagonalPeak>& peaks,
    const CoarseLocusOptions& options, CoarseLocusScratch& scratch,
    std::uint32_t tie_seed, std::size_t* coalesced = nullptr,
    std::vector<CoarseLocus>* overflow_out = nullptr,
    const CoarseLocusQueryCover* query_cover = nullptr);

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** query_partition ****/
// The re-vote's posting budget, scratch and stats (the re-vote itself, which the second
// coarse pass runs over the whole read, is in vote.cpp's second_pass section), the
// unexplained pieces of a read, the leftover harvest over the query no printed record
// covers, and the segment geometry the shadow rule and the pieces read.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Posting budget per read of the re-vote and the leftover harvest: past it they decline
// the read.
inline constexpr std::int64_t kRnaExplainMaxPostings = 65536;

// The geometry every colinearity decision reads.
struct RnaSegmentGeometry {
  int reference_id = -1;
  bool reverse = false;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  RnaQuerySpan forward_span; // forward-frame query hull
};

// Same contig and strand, and reference spans that overlap or continue each other: two
// descriptions of one segment, never competing (not a rival or a demoted vote).
bool rna_same_place(const RnaSegmentGeometry& a, const RnaSegmentGeometry& b,
                    int max_intron) noexcept;

// Caller-owned reusable buffers: cleared per call, never shrunk.
struct RnaExplainScratch {
  // The re-vote's own capture area (the two strand logs and bin counters the
  // fused capture's drain reads; its view/selection members stay empty).
  placement::FusedCaptureScratch revote_capture;
};

// The re-vote's own accounting: how many extracted seeds were re-voted, how
// many postings they hold, and whether that exceeded kRnaExplainMaxPostings
// (in which case nothing at all was voted and both counters are 0).
struct RnaRevoteStats {
  int seeds_revoted = 0;
  std::int64_t postings = 0;
  bool over_budget = false;
  // Records logged per strand, read off the logs. `postings_inspected` counts postings
  // that reached the strand test, the same for both strands since each posting is
  // tested for both.
  std::int64_t emitted_forward = 0;
  std::int64_t emitted_reverse = 0;
  std::int64_t postings_inspected = 0;
  std::int64_t rejected_forward = 0;
  std::int64_t rejected_reverse = 0;
};

// The leftover's intervals: every forward-frame interval of at least `min_bp` that
// neither `primary_segments` nor `explained` covers, the clips outside their joint hull
// with `primary_hull` included.
void rna_leftover_intervals(RnaQuerySpan primary_hull,
                            const std::vector<RnaQuerySpan>& primary_segments,
                            const std::vector<RnaQuerySpan>& explained,
                            int read_len, int min_bp,
                            std::vector<RnaQuerySpan>& out);

// minimap2's chain floor (min_cnt, min_chain_score), on the leftover's chains only.
inline constexpr int kRnaLeftoverChainMinCount = 3;
inline constexpr int kRnaLeftoverChainMinScore = 40;

// One contig and lane's anchors from rna_leftover_pools, in that lane's query frame.
struct RnaLeftoverPool {
  int reference_id = -1;
  bool reverse = false;
  CandidatePool pool;
};

// The leftover harvest: rna_revote_intervals's walk over `intervals` (every extracted
// seed whose k-mer lies inside one, under the vote's occurrence gate, each posting in the
// lane its orientation bit allows), each posting an anchor in its lane's frame instead of
// a vote. The anchors are pooled per contig and lane in reference, then query order,
// duplicates removed. Over kRnaExplainMaxPostings postings nothing is pooled and false is
// returned. A seed marked in `readmitted` (by forward read position, or null) passes the
// occurrence gate.
bool rna_leftover_pools(const LongReadSeedContext& seed_ctx,
                        const std::vector<QuerySeed>& forward_seeds,
                        const std::vector<KmerPostingView>& views,
                        const std::vector<RnaQuerySpan>& intervals, int read_len,
                        int k, int n_chr, const std::uint8_t* readmitted,
                        std::vector<RnaLeftoverPool>& pools);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** worker_scratch ****/
// Caller-owned per-worker RNA mapping scratch, one per worker thread (minimap2's
// mm_tbuf_t). Buffers are reused across reads.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct WorkerScratch {
  placement::FusedCaptureScratch capture;
  placement::CoarseLocusScratch coarse;
  // The re-vote's capture.
  RnaExplainScratch explain;
  // Per-read posting memo shared by both strands (the reverse projection keeps canonical
  // keys and indexes it through n-1-j). Scoped to a read by new_read()'s generation bump:
  // the stream identity check cannot tell two reads apart when the allocator recycles a
  // block.
  SeedPostingMemo postings;
  // Realization's DP arena and packet buffers, shared by rank 1, every rival and both
  // transcript hypotheses, which run one after another on this thread.
  SpliceRealizationScratch realization;
  // seed/seed.h drop_query_repeats: the read's sorted seed keys, and its forward
  // stream without the repeated ones. Each read that uses them refills them.
  std::vector<std::uint64_t> query_occ_keys;
  std::vector<QuerySeed> query_occ_seeds;
  // The read's seeds re-admitted past the cap (ReadFacts::readmitted). Each read that
  // uses it refills it.
  std::vector<std::uint8_t> readmitted;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** backend ****/
// The RNA spliced backend. The shared engine dispatches each read as
// Backend::map_read(...); this header must not include engine/engine.h.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Maps one read: vote capture, exon peaks, coarse locus catalogue, exact anchor path,
// then map-only projection or splice realization with the rival lifecycle. rctx and
// seed_ctx are composed by the engine (compose_rna_runtime); worker_scratch is the
// caller's per-worker scratch. `read_name_hash` is tie_name_hash of the read name (0
// without one, as minimap2 without a qname); with the read length it seeds the tie-break
// between two equally scored loci.
struct RnaBackend {
  using WorkerScratch = ::fa::cpu::lr::rna::WorkerScratch;

  static AlignResult map_read(const Context &rctx,
                           const LongReadSeedContext &seed_ctx,
                           WorkerScratch &worker_scratch,
                           const std::string &read,
                           const std::vector<uint8_t> &fwd_enc,
                           std::vector<uint8_t> &rc_enc,
                           const std::vector<QuerySeed> *shared_fwd_syncmer_seeds,
                           std::uint32_t read_name_hash);
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** aggregate ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// Read bases covered by the union of [pos, pos + k) over ascending seed positions.
std::uint32_t union_seed_coverage(
    const std::vector<int>& sorted_read_pos, int k);

}  // namespace placement
}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** frame_diagonal_support ****/
// Which query frame a harvested candidate pool is diagonally coherent in, read off the
// pool before any chain runs. In the right frame an exon's anchors share a diagonal r - q.
// In the mirrored frame q' = length - span - q, and every anchor of a pool has the same
// span, so the mirror's diagonal is the anti-diagonal r + q shifted by a constant, which
// run and band clustering cannot see. One pool therefore answers for both frames. The
// fine stage elects each window's frame this way and chains both frames only when the
// verdict is ambiguous.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Diagonal coherence of one pool in both query frames at one tolerance, over the
// final-pool anchors only. A run is a maximal group of sorted keys whose consecutive
// differences are within the tolerance: `mass` counts the anchors in runs of two or more,
// `best_run` is the largest run and `runs_ge3` counts runs of three or more. `band` is
// the most keys inside any window of width `tolerance`; unlike run mass it stays
// discriminating on a dense pool, where the wrong frame's keys chain into long runs.
struct FrameDiagonalSupport {
  int32_t tolerance = 0;
  uint32_t anchors = 0;
  uint32_t mass_nominated = 0;
  uint32_t best_run_nominated = 0;
  uint32_t runs_ge3_nominated = 0;
  uint32_t mass_mirror = 0;
  uint32_t best_run_mirror = 0;
  uint32_t runs_ge3_mirror = 0;
  uint32_t band_nominated = 0;
  uint32_t band_mirror = 0;
};

// Clusters the pool's final anchors in both key spaces. `tolerance` is in reference
// bases; a negative value is read as 0. A refused or empty pool reports zeros.
FrameDiagonalSupport compute_frame_diagonal_support(const CandidatePool& pool,
                                                    int32_t tolerance);

// kAmbiguous means the bands are too close to settle the frame cheaply; the window is
// then chained in both frames and settled by the chain scores.
enum class FrameElection { kNominated, kMirror, kAmbiguous };

// With margin = |band_n - band_m|, the window is decided when
//
//     margin >= min_margin  and  margin >= min_ratio * (band_n + band_m)
//
// and a decided window elects the mirror iff band_m > band_n, the tie rule of the
// chain-score comparison. The absolute bound keeps a tiny pool from deciding on noise;
// the relative bound keeps a dense pool from deciding on a small margin.
FrameElection elect_query_frame(const FrameDiagonalSupport& support);

// The same reading off the two strand-routed pools: `nominated` is the harvested lane's
// pool and `opposite` the other lane's, already mirrored to q' = length - span - q. The
// opposite pool's own diagonal key r - q' is therefore the mirror key, so both pools are
// read with r - q. Unlike the one-pool reading, the two masses come from different anchor
// sets, each lane capped on its own. A refused pool contributes zeros.
FrameDiagonalSupport
routed_frame_diagonal_support(const CandidatePool& nominated,
                              const CandidatePool& opposite, int32_t tolerance);

// Tolerance at which the fine stage computes the support, the same for rank-1 and rival
// windows.
constexpr int32_t kFrameElectionTolerance = 8;

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** selected_locus_anchoring ****/
// Candidate pool harvest around the selected (rank-1) coarse locus.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct Rank1HarvestOptions {
  int global_occurrence_cap = 200;
  int chain_query_gap = 5000;
  // The chain's diagonal band.
  int chain_band = 200000;
  int chain_min_count = 2;
  int chain_min_score = 20;
  float chain_gap_penalty = 0.0f;
  // Also the reference pad the harvest adds on each side of the coarse envelope.
  int max_intron = 200000;
  // SkeletonHarvestParams::readmitted.
  const std::uint8_t* readmitted = nullptr;
};

struct Rank1HarvestResult {
  CandidatePool bounded_pool;
  // `nominated_pool` is the locus's own strand lane and `opposite_pool` the other lane,
  // already in its query coordinates; each is capped on its own. A refused harvest never
  // sets `routed`.
  CandidatePool nominated_pool;
  CandidatePool opposite_pool;
  bool routed = false;
  RoutingCounters routing;
  chaining::SplicedChainParams anchor_path_params;
  // Harvest work counters.
  HarvestWorkCounters work;
  // What the pool was built against, kept so mirror_candidate_pool_query_frame can
  // project the opposite query frame out of this pool instead of harvesting again. Set
  // whenever the region build succeeded.
  std::vector<SearchRegion> regions;
  SkeletonHarvestParams harvest_params;
  std::vector<SkippedBlock> skipped_blocks;
  // The same blocks in the opposite lane's coordinates.
  std::vector<SkippedBlock> opposite_skipped_blocks;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

// The fine chain's parameters for a harvest under `options`; k is the index's seed
// length.
chaining::SplicedChainParams rank1_anchor_path_params(
    int k, const Rank1HarvestOptions& options);

// `posting_memo` is the caller's per-read memo, shared by both strand streams. It only
// saves posting probes across the read's harvests; the result is the same without it.
Rank1HarvestResult harvest_rank1(
    const SeedIndex& index, const std::vector<QuerySeed>& strand_seeds,
    const placement::CoarseLocus& rank1_locus,
    const std::vector<placement::CoarseDiagonalPeak>& retained_peaks,
    int read_length, int reference_length, const Rank1HarvestOptions& options,
    SeedPostingMemo* posting_memo = nullptr);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** select_exact_anchor_path ****/
// Runs the spliced colinear chain over an RNA candidate pool and converts the best chain
// into an ExactAnchorPath.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The chains the recurrence built but the selection did not take; the best competing one
// feeds the chain MAPQ. Every field stays -1 until the recurrence runs.
struct SiblingChainStats {
  // Best score among the non-selected chains that compete with the selected one under
  // minimap2's mask_level (rna_query_spans_compete); a query-disjoint piece of the same
  // read is not an alternative.
  int32_t sib_score = -1;
  int32_t sib_count = -1;     // that chain's anchor count
  int32_t chains_total = -1;  // how many chains the recurrence returned
};

// A chain of the window other than the selected one that does not compete with it under
// mask_level: its own path over the same pool, and its own competing sibling.
struct DisjointAnchorPath {
  ExactAnchorPath bundle;
  SiblingChainStats siblings;
};

// `disjoint`, when non-null, is cleared and, on an accepted selection, receives every
// other chain that does not compete with the selected one, in the recurrence's order; a
// chain whose path fails validation is left out.
ExactAnchorPathResult select_exact_anchor_path(
    uint64_t read_key, int32_t query_length, int32_t reference_id,
    bool mapping_reverse, const CandidatePool& pool,
    const chaining::SplicedChainParams& params,
    SiblingChainStats* siblings = nullptr,
    std::vector<DisjointAnchorPath>* disjoint = nullptr);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** plain_elect ****/
// The plain PAF lane's catalogue election. An intronless copy of a gene keeps the k-mers
// spanning each junction, which the spliced gene cannot have, so on anchor count alone
// the copy can outscore the gene. The score gives each junction back what it costs.

namespace fa::cpu::lr::rna {

// Scores are integers in thousandths of an anchor.
inline constexpr std::int64_t kRnaPlainElectScale = 1000;

// The credit per junction, in thousandths of an anchor. A junction removes about k - 1
// anchor positions, (k - 1) * 2 / (k - s + 1) = 14/3 anchors at the splice presets'
// default seeding (k = 15, s = 10). The credit is 4.7, the measured value, not 14/3,
// and is not rescaled for an index of another k or s.
inline constexpr std::int64_t kRnaPlainElectJunctionCredit = 4700;

} // namespace fa::cpu::lr::rna

/**** splice_realizer ****/
// Realization of a fixed exact anchor path under one or both transcript hypotheses.

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
  // same core/core.cpp replay as DNA, `~` intron token included. A default request
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
  // The segment's own chain score: the score it was handed less what its
  // continuation took (minimap2's mm_split_reg).
  int chain_score = 0;
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

SpliceRegionReject
classify_splice_region(const SpliceRegionEvidence& evidence, int query_length,
                       const SpliceControllerOptions& options) noexcept;

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

/**** output ****/
// Commits the selected realized RNA family to the output record.

namespace fa { namespace cpu { namespace lr { namespace rna {

// A family's chain tags, each -1 when none: cm:i and s1:i, the anchors and score of its
// own chain, and s2:i, the best competing chain score its MAPQ was weighed against.
struct RnaChainTags {
  int chain_anchors = -1;
  int chain_score = -1;
  int secondary_chain_score = -1;
};

// A family emitted beside the committed one: a supplementary piece of the parent rule
// (hit.cpp), with its own chain MAPQ.
// `result` is borrowed from the caller's realization store and must outlive the commit.
struct RnaChimericFamily {
  const RnaSpliceRealizationResult* result = nullptr;
  int mapq = 0;
  // This family's own chain tags.
  RnaChainTags chain_tags;
};

// Attaches one runner-up family to `out`'s alternative hypotheses, moving its segments
// out. Call only after commit_realized_alignment, which overwrites `out`. `chain_tags`,
// the family's own chain (s2:i left -1), go on every record of it.
void attach_runner_up_family(RnaSpliceRealizationResult& runner_up,
                             const RnaChainTags& chain_tags,
                             AlignResult& out);

}}}}  // namespace fa::cpu::lr::rna

/**** rival_lifecycle ****/
// The RNA rival realization lifecycle, in minimap2's shape:
//   1. Classify every retained chain against the primary by query geometry (the
//      mask_level test of mm_set_parent): an overlapping chain competes for the same
//      bases, a disjoint one is a co-primary, a different part of the read.
//   2. Admit competitors by chain score (pri_ratio / min_diff), co-primaries by their
//      own chain-quality floor.
//   3. Realize the admitted competitors, bounded, and the chain-stage co-primaries.
//   4. Elect the primary on the realized DP maximum (mm_hit_sort on dp_max).
//   5. Feed the demoted maxima to the MAPQ as dp2 / n_sub evidence.
// This header is the pure policy; the backend owns storage, harvests, realizations and
// the adoption.
//
// The election fails closed: a competitor takes the primary slot only when its own and
// the incumbent's realized DP maxima are both present and positive and its own is
// strictly greater. Exact ties (queue keys, both elections, the dp2 owner, and in the
// backend the recalibration leader and runner-up order) go to the lower catalogue index,
// the incumbent's (0) first. The catalogue is ordered by vote, then by minimap2's
// read-seeded hash, so loci tied on votes too are decided by the hash.
//
// A co-primary never takes part in the election and never supplies dp2 or n_sub. A
// realized one may print as a supplementary piece with its own chain MAPQ; that writes
// nothing back onto the hypotheses.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Borrowed by pointer only; declared here to avoid a dependency cycle with the coarse
// catalogue.
namespace placement {
struct CoarseLocus;
}

// What a hypothesis IS relative to the rank-1 incumbent, decided once before
// any realization on forward-frame query geometry.
enum class RnaRivalClass : std::uint8_t {
  // The rank-1 chain, already realized by the caller. Exactly one hypothesis carries
  // this, and it is always element 0 of the hypothesis vector.
  Incumbent = 0,
  // Competes for the same query bases (mask_level overlap, or censored geometry): it may
  // win the primary slot, and if it loses it is what dp2 / n_sub measure.
  Competitor = 1,
  // Explains a different part of the read: never elected, never dp2 or n_sub.
  CoPrimary = 2,
};

// One hypothesis in the lifecycle. Its bundle, summary and locus are borrowed from the
// backend's per-read storage, which outlives the vector, so reordering the vector moves
// only pointers. `result` is the one owned member, the realization run for a
// non-incumbent; it stays empty for element 0, whose realization the caller holds.
struct RnaRealizedHypothesis {
  int catalogue_index = -1;
  int reference_id = -1;
  const placement::CoarseLocus* locus = nullptr;
  const ExactAnchorPath* bundle = nullptr;
  const ExactAnchorPathSummary* summary = nullptr;
  SiblingChainStats siblings;
  // The catalogue's order key, so the committed-frame MAPQ rivals can be rebuilt from
  // this vector.
  std::int64_t rank_score = 0;
  // The chain's own values, in its own mapping frame; `reverse` is the post-repair
  // strand, read straight.
  bool chained = false;
  int chain_score = 0;
  int chain_anchors = 0;
  int q_begin = -1;
  int q_end = -1;
  bool reverse = false;
  bool repaired = false;
  // Anchor-union query coverage of `bundle`, the commit's qcov input.
  std::uint64_t query_covered_bases = 0;
  bool coverage_ok = false;
  // Lifecycle verdicts.
  RnaRivalClass klass = RnaRivalClass::Competitor;
  bool admitted = false;
  bool realized = false;
  // This chain is the committed placement re-found from a neighbouring catalogue window
  // (RnaChainMapqRival::shadow). It stays in the lifecycle but never supplies dp2 and is
  // inadmissible to the MAPQ.
  bool shadow = false;
  // The election price: segment 0's realized DP maximum, minimap2's dp_max currency.
  // The election and the dp2 selection read only this field.
  std::optional<int> dp_maximum;
  // Segment 0 alone, and the realized segments' forward-frame query hull
  // [aq_begin, aq_end), -1/-1 when unmeasured (rna.h).
  std::optional<int> dp_segment0;
  int aq_begin = -1;
  int aq_end = -1;
  // The rank-recalibration price, set only when the recalibration fired for this read
  // and this hypothesis's accounting could be formed. dp_maximum keeps the raw price.
  std::optional<int> dp_recal;
  RnaSpliceRealizationResult result;
};

// Step 1: minimap2's mask_level verdict on the forward-frame spans
// (rna_query_spans_compete, evaluated by the caller). Censored geometry competes: an
// unmeasurable span is not evidence of disjointness.
RnaRivalClass rna_classify_rival(bool spans_compete) noexcept;

// Step 2a: minimap2's secondary retention band, as an OR: within pri_ratio of the
// primary chain score, or within min_diff of it. Both inclusive; an unchained rival is
// never admitted.
bool rna_competitor_admitted(int chain_score, int primary_chain_score,
                             const RnaRivalLifecycleConfig& config) noexcept;

// Step 2b: a co-primary's floor, on chain score, anchor count and the query bases it adds
// over the incumbent (rna_uncovered_query_bases). minimap2's numbers, fixed.
inline constexpr int kRnaCoprimaryMinChainScore = 40;
inline constexpr int kRnaCoprimaryMinAnchors = 3;
inline constexpr int kRnaCoprimaryMinQueryBases = 50;

bool rna_coprimary_admitted(int chain_score, int chain_anchors,
                            int uncovered_query_bases) noexcept;

// The query a second-pass placement or rival must explain. Fixed, not tunable.
inline constexpr int kRnaChimeraMinQueryBases = 100;

// Step 3: the realization queue's total order: chain score, then anchor count, then
// catalogue index. Total, so the bounded budget never depends on the sort.
bool rna_rival_queue_before(const RnaRealizedHypothesis& a,
                            const RnaRealizedHypothesis& b) noexcept;

// Step 3b: minimap2's rank recalibration (mm_update_dp_max / mm_recal_max_dp), between
// realization and the election. On raw DP a spliced alignment pays the long-gap open at
// every junction while an intronless retrocopy pays nothing, so the retrocopy can
// outscore the gene. When two aligned hits both cover most of the read, minimap2 ranks
// every hit by an intron-free, divergence-amplified score instead:
//
//     match_sc * (mlen - b2*n_mis - gap_cost)
//
// Introns contribute nothing, a mismatch costs b2 = 0.5/divergence match units, and a
// gap costs b2 + log2(1+len) per operation with no per-base cost. A hypothesis whose
// accounting cannot be formed keeps its raw price.

// minimap2's rank_min_len, rank_frac and divergence floor.
inline constexpr int kRnaRankRecalMinReadLen = 500;
inline constexpr double kRnaRankRecalFrac = 0.9;
inline constexpr double kRnaRankRecalMinDivergence = 0.02;

// One realized hypothesis's primary-segment accounting and numeric CIGAR, borrowed from
// RnaSpliceHypothesisResult::segment_cigars. The counters follow minimap2's
// mm_update_extra: `matches` is mlen (unambiguous equal residues), `ambiguities` is
// n_ambi (excluded from the block length), and `mismatches` is blen - mlen - n_gap. The
// CIGAR uses the realizer's encoding, op = value & 0xf over "MIDNSHP=XB".
struct RnaRankRecalAccounting {
  int matches = 0;
  int mismatches = 0;
  int ambiguities = 0;
  const std::uint32_t* cigar = nullptr;
  std::size_t cigar_len = 0;
};

// minimap2's mm_event_identity, mlen / (blen + n_ambi - n_gap + n_gapo): identity over
// events, each gap counting once however long it is. blen is matches + mismatches +
// n_gap; introns are neither gap nor block. Returns -1 on an empty or malformed CIGAR,
// negative counters or a non-positive denominator, and the caller then skips the
// recalibration.
double rna_event_identity(const RnaRankRecalAccounting& accounting) noexcept;

// minimap2's divergence-to-mismatch weight:
//
//     div = 1 - identity;  if (div < 0.02) div = 0.02;
//     b2 = .5 / div;  if (b2 * a < b) b2 = (double)a / b;
//
// The second line is transcribed as minimap2 writes it, reciprocal included, so the same
// hypothesis is elected. `mismatch_sc` is minimap2's positive `b`; a non-positive a or b
// skips that line.
double rna_rank_recal_b2(double top_event_identity, int match_sc,
                         int mismatch_sc) noexcept;

// minimap2's mm_recal_max_dp. The caller floors the result at 0, as mm_update_dp_max
// does.
int rna_recal_max_dp(const RnaRankRecalAccounting& accounting, double b2,
                     int match_sc) noexcept;

// minimap2's gate: the read is long enough to rank, the top hit is aligned across most of
// it (its realized forward-frame hull [aq_begin, aq_end)), and the runner-up is within
// `frac` of the top. Both DP prices must be positive.
bool rna_rank_recal_gate(int read_len, int top_aq_begin, int top_aq_end,
                         int top_dp, int second_dp) noexcept;

// Raw-price veto band for the combined election, in match units. Not a minimap2 number:
// electing on the recalibrated price alone flips reads the raw price placed correctly,
// so a flip may overturn at most this much raw margin.
inline constexpr int kRnaRankRecalRawVetoMargin = 20;

// Step 4: the elected primary's index; element 0 must be the incumbent. Only realized
// incumbent and competitor entries with a positive DP maximum take part. Returns 0 unless
// a competitor strictly beat a present incumbent maximum.
std::size_t rna_elect_primary_hypothesis(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept;

// Step 4 on a recalibrated read: argmax over raw segment-0 plus recalibrated price with
// the same tie discipline. A winner that differs from the raw election's takes the slot
// only when the raw margin it overturns is within kRnaRankRecalRawVetoMargin
// (inclusive). Missing evidence falls back to the raw verdict.
std::size_t rna_elect_primary_rank_recal(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept;

// Step 5: the hypothesis supplying the MAPQ's dp2: the largest realized DP maximum among
// the non-elected competitors (the demoted incumbent included), ties to the lower
// catalogue index. Co-primaries and shadows are excluded.
std::optional<std::size_t>
rna_select_dp2_owner(const std::vector<RnaRealizedHypothesis>& hypotheses,
                     std::size_t elected) noexcept;

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** rna_parent ****/
// minimap2's region bookkeeping (hit.c, align.c) over a small hit struct:
// mm_set_parent, mm_update_dp_max and mm_hit_sort's order. Pure; the caller builds the
// hits and reads the verdicts back by `owner`.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// minimap2's MM_PARENT_UNSET.
inline constexpr int kRnaParentUnset = -1;

// One region (mm_reg1_t). Query coordinates are on the forward read.
struct RnaParentHit {
  int qs = 0;
  int qe = 0;
  int score = 0; // the region's own chain score
  int cnt = 0;
  int rid = -1;
  std::int64_t rs = 0;
  std::int64_t re = 0;
  // minimap2's r->p != 0: the region is aligned and dp_max / dp_max2 are set.
  bool has_dp = false;
  int dp_max = 0;
  int dp_max2 = 0;
  // mm_update_dp_max's input for an aligned region; accounting_ok false when it could
  // not be formed.
  RnaRankRecalAccounting accounting;
  bool accounting_ok = false;
  // minimap2's bookkeeping, written by rna_parent_set.
  int parent = kRnaParentUnset;
  int subsc = 0;
  int n_sub = 0;
  // The caller's keys: which hypothesis the region belongs to, its catalogue index (the
  // declared tie-break) and its segment within the hypothesis.
  int owner = -1;
  int catalogue = 0;
  int segment = 0;
};

// mm_set_parent with hard_mask_level 0 and no mask_len bound, without the alt-contig
// terms. Sets parent, subsc, n_sub and, on aligned pairs, dp_max2; r[0]'s n_sub is never
// reset (minimap2 hit.c:182).
void rna_parent_set(std::vector<RnaParentHit>& hits, float mask_level, int sub_diff);

// mm_update_dp_max with rank_min_len: reprices every aligned region's dp_max by
// mm_recal_max_dp when the read is at least min_len long, there are two aligned regions,
// the top spans at least frac of the read and the second is within frac of the top.
// Atomic: a region whose accounting is missing, or a top region without an identity,
// leaves every price raw. Returns whether the prices were rewritten.
bool rna_parent_update_dp_max(int qlen, int min_len, std::vector<RnaParentHit>& hits,
                              float frac, int a, int b);

// mm_hit_sort's order: dp_max for aligned regions, score otherwise, descending; ties to
// the lower catalogue index, then the lower segment. Stable.
void rna_parent_hit_sort(std::vector<RnaParentHit>& hits);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** rank2_arbitration ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** rival_pricing ****/
// Prices a realized RNA hypothesis for the rival lifecycle. A realized hypothesis is a
// family (a primary plus its continuation segments), while minimap2's dp_max is per hit.
// Every segment's dp_maximum is already in region_evidence, in the read's forward frame,
// so pricing runs no DP.

namespace fa::cpu::lr::rna {

// The winning transcript-orientation hypothesis, or null when the realization was refused
// or produced no segment.
const RnaSpliceHypothesisResult*
rna_winning_hypothesis(const RnaSpliceRealizationResult& result) noexcept;

// Segment 0's dp_maximum, the promotion statistic of minimap2's post-DP dp_max sort; null
// without a winning hypothesis.
std::optional<int>
rna_segment0_dp_maximum(const RnaSpliceRealizationResult& result) noexcept;

// The family's realized query footprint in the read's forward frame: the
// half-open hull [min begin, max end) over its segments. Null when there is
// no winning hypothesis or no segment carries a valid span.
std::optional<RnaQuerySpan>
rna_realized_query_span(const RnaSpliceRealizationResult& result) noexcept;

} // namespace fa::cpu::lr::rna

/**** read_state ****/
// map_read's per-read state in five structs, by lifetime: the read's facts, its state, the
// attempt's, the rival lifecycle's and the committed chain's. map_read binds each field's
// name to it, and the stage functions of map.cpp, vote.cpp, hit.cpp and leftover.cpp take
// the structs by reference.


namespace fa::cpu::lr::rna {

// The read's facts: map_read's arguments and what the prologue and Stages 1-3 fix for
// the whole read. Read-only once set, except the lazily encoded rc_enc.
struct ReadFacts {
  ReadFacts(const Context& rctx, const LongReadSeedContext& seed_ctx,
            WorkerScratch& worker_scratch, const std::string& read,
            const std::vector<uint8_t>& fwd_enc, std::vector<uint8_t>& rc_enc,
            const std::vector<QuerySeed>* shared_fwd_syncmer_seeds)
      : rctx(rctx), seed_ctx(seed_ctx), worker_scratch(worker_scratch),
        read(read), fwd_enc(fwd_enc), rc_enc(rc_enc),
        shared_fwd_syncmer_seeds(shared_fwd_syncmer_seeds), cfg(rctx.opts.cfg),
        postings(worker_scratch.postings), capture(worker_scratch.capture) {}
  ReadFacts(const ReadFacts&) = delete;
  ReadFacts& operator=(const ReadFacts&) = delete;

  const Context& rctx;
  const LongReadSeedContext& seed_ctx;
  WorkerScratch& worker_scratch;
  const std::string& read;
  const std::vector<uint8_t>& fwd_enc;
  std::vector<uint8_t>& rc_enc;
  const std::vector<QuerySeed>* shared_fwd_syncmer_seeds;
  int read_len = 0;
  int reference_count = 0;
  std::uint32_t tie_seed = 0;
  const RnaConfig& cfg;
  std::uint64_t read_key = 0;
  int k = 0;
  int min_support = 0;
  int min_intron = 0;
  int max_intron = 0;
  int max_loci = 0;
  std::vector<QuerySeed> projected_shared_reverse;
  SeedPostingMemo& postings;
  placement::FusedCaptureScratch& capture;
  // Stage 1's seeds re-admitted past the cap, by forward read position, or null; every
  // later occurrence gate of the read reads them.
  const std::uint8_t* readmitted = nullptr;
  placement::CoarseLocusOptions copt;
  const std::vector<QuerySeed>* fwd_seeds = nullptr;
  Rank1HarvestOptions harvest_options;
};

// The catalogue and the read's record, mutated across the stages: the second pass
// replaces peaks and loci, so readers take them from here at call time.
struct ReadState {
  explicit ReadState(Result& out) : out(out) {}
  ReadState(const ReadState&) = delete;
  ReadState& operator=(const ReadState&) = delete;

  Result& out;
  std::vector<placement::CoarseDiagonalPeak> peaks;
  std::vector<placement::CoarseLocus> loci;
  int coarse_pass = 1;
  std::int64_t catalogue_rank1_score = 0;
  std::int64_t catalogue_rank2_score = 0;
};

// One Stage-4 attempt's state, reset at the top of each attempt; the accepted
// attempt's is what Stage 5 reads. selected_locus may point at repaired_locus and
// committed_geometry_ptr at committed_geometry.
struct AttemptState {
  AttemptState() = default;
  AttemptState(const AttemptState&) = delete;
  AttemptState& operator=(const AttemptState&) = delete;

  // Orientation repair may replace the nomination with its mirror locus, which needs
  // storage the catalogue does not own; both output modes project the same one.
  placement::CoarseLocus repaired_locus;
  const placement::CoarseLocus* selected_locus = nullptr;
  SiblingChainStats rank1_siblings;
  SiblingChainStats mirror_siblings;
  int chr = -1;
  bool fine_accepted = false;
  bool selected_reverse = false;
  bool repaired_orientation = false;
  bool path_coverage_ok = false;
  std::uint64_t selected_query_covered_bases = 0;
  ExactAnchorPathResult anchor_path_built;
  ExactAnchorPathSummary selected_path_summary;
  // The committed chain's geometry for the shadow rule; null until this attempt's chain
  // and summary are fixed.
  RnaSegmentGeometry committed_geometry;
  const RnaSegmentGeometry* committed_geometry_ptr = nullptr;

  // The per-attempt reset at the top of each Stage-4 attempt: selected_locus back at
  // the catalogue's rank 1 and the frame and trace fields cleared. chr,
  // anchor_path_built and the committed geometry are set later in the attempt.
  void reset_for_attempt(const placement::CoarseLocus& rank1_locus) {
    selected_locus = &rank1_locus;
    selected_reverse = false;
    repaired_orientation = false;
    path_coverage_ok = false;
    selected_query_covered_bases = 0;
    rank1_siblings = rna::SiblingChainStats{};
    mirror_siblings = rna::SiblingChainStats{};
  }
};
// A catalogue rival's chain evidence for the MAPQ formula. A rival that could not be
// chained still competes through its vote ratio (censored evidence is not zero
// evidence) and never counts toward n_sub.
struct ProductionRival {
  rna::RnaChainMapqRival rival;
  const char* refusal = "none";
  int retry = 0;
  // The rival's chain, kept so the rival can be realized without
  // re-harvesting its window; `owns_path` says whether `chained_path` holds one.
  rna::ExactAnchorPath chained_path;
  bool owns_path = false;
  // The locus the chain was taken on (the catalogue entry, or its mirror when
  // `repaired`), with that window's summary and sibling evidence; `reference_id` is the
  // harvested contig.
  rna::placement::CoarseLocus locus;
  rna::ExactAnchorPathSummary summary;
  rna::SiblingChainStats siblings;
  int reference_id = -1;
  bool repaired = false;
};

// The supplementary pieces: an after-DP primary outside the elected family, its count
// of children and dp_max2 without its shadows, the hypotheses owning its children at
// either stage, and the first child at dp_max2.
struct PieceRegion {
  rna::RnaParentHit region;
  std::vector<std::size_t> children;
  std::optional<std::size_t> dp2_owner;
  // Its hypothesis and every hypothesis at its place, kept off the runner-ups.
  std::vector<std::size_t> place;
};

// The rival lifecycle's state (Stage 5): the rival chains and their hypotheses, the
// election's verdicts and the parent rule's pieces, which the MAPQ, the leftover, the
// pieces' records and the runner-ups read after it. hypotheses borrow
// production_rivals' and leftover_rivals' chains.
struct Lifecycle {
  Lifecycle() = default;
  Lifecycle(const Lifecycle&) = delete;
  Lifecycle& operator=(const Lifecycle&) = delete;

  // The rival set, chained once for the lifecycle and reused by the MAPQ.
  std::vector<ProductionRival> production_rivals;
  // The leftover's chains, which their hypotheses borrow.
  std::vector<ProductionRival> leftover_rivals;
  std::vector<rna::RnaRealizedHypothesis> hypotheses;
  // The demoted incumbent's realization, parked here on an adoption. Element 0's
  // `result` stays empty by the borrowed-storage contract (rna.h). Read by
  // the pieces and the runner-up attach.
  rna::RnaSpliceRealizationResult demoted_incumbent;
  std::vector<char> locus_on_wire;
  std::size_t lifecycle_elected = 0;
  std::optional<std::size_t> lifecycle_dp2_owner;
  // Whether the rank recalibration repriced this read; selects the election rule.
  bool rank_recal_fired = false;
  // The parent rule's chain-stage hits (owner = hypothesis index), in the queue's order.
  std::vector<rna::RnaParentHit> parent_chain_hits;
  std::vector<PieceRegion> pieces;
  int parent_sub_diff = 0;
  // The after-election parent rule's input, kept for the leftover's rerun: its regions
  // in order (the elected family's first), the hypotheses it leaves out, and the chain
  // stage's regions.
  std::vector<rna::RnaParentHit> pieces_dp_input;
  std::vector<char> pieces_left_out;
  std::vector<rna::RnaParentHit> pieces_chain;
};

// Stage 5's committed frame: the record's locus, path, strand, contig and catalogue rank,
// and its own MAPQ inputs. The incumbent's until the lifecycle adopts a rival.
struct Committed {
  Committed() = default;
  Committed(const Committed&) = delete;
  Committed& operator=(const Committed&) = delete;

  const placement::CoarseLocus* committed_locus = nullptr;
  const ExactAnchorPath* committed_path = nullptr;
  const ExactAnchorPathSummary* committed_path_summary = nullptr;
  bool committed_reverse = false;
  int committed_chr = -1;
  std::size_t committed_rank = 1;
  bool commit_coverage_ok = false;
  std::uint64_t commit_qcov = 0;
  const SiblingChainStats* committed_siblings = nullptr;
};

} // namespace fa::cpu::lr::rna

/**** second_pass ****/

namespace fa::cpu::lr::rna {

// The second coarse pass (vote.cpp): every extracted seed of the whole read re-voted and
// the catalogue rebuilt, replacing rs.peaks and rs.loci and setting rs.coarse_pass to 2.
// False, with peaks and loci untouched, when the re-vote is over its posting budget or
// the rebuilt catalogue is empty.
bool build_second_pass_catalogue(const ReadFacts& facts, ReadState& rs);

} // namespace fa::cpu::lr::rna

/**** plain_election ****/

namespace fa::cpu::lr::rna {

// The plain PAF lane's election (hit.cpp) between the incumbent chain `anchor_path` and
// the chained production rivals. On an adoption it repoints the committed_* outputs at
// the elected rival and rewrites chain_evidence for the elected chain; otherwise it
// leaves both as the caller set them (the incumbent).
void rna_plain_election(const ReadFacts& facts, const ReadState& rs,
                        const AttemptState& as, const ExactAnchorPath& anchor_path,
                        const std::vector<ProductionRival>& production_rivals,
                        RnaChainMapqEvidence& chain_evidence,
                        const placement::CoarseLocus*& committed_locus,
                        const ExactAnchorPathSummary*& committed_summary,
                        std::uint64_t& committed_qcov);

} // namespace fa::cpu::lr::rna

/**** parent_rule ****/
// The rival lifecycle's helpers (hit.cpp), shared by the lifecycle and the leftover.

namespace fa::cpu::lr::rna {

// Binds a realization request to contig reference_id; false when it will not bind.
bool bind_realization_reference(const ReadFacts& facts,
                                RnaSpliceRealizationRequest& target,
                                int reference_id);

// Prices one realization for the lifecycle (segment 0's dp_maximum, the family price and
// the realized hull).
void price_hypothesis(RnaRealizedHypothesis& h,
                      const RnaSpliceRealizationResult& result);

// mm_set_parent over hits in their stage's order, less the hits whose owner `left_out`
// marks.
std::vector<RnaParentHit> parent_regions(const Lifecycle& lc,
                                         const std::vector<RnaParentHit>& ordered,
                                         const std::vector<char>* left_out);

// A chained hypothesis's place for the shadow rule; false without a chain.
bool shadow_geometry(const ReadFacts& facts, const RnaRealizedHypothesis& h,
                     RnaSegmentGeometry& g);

// Whether lc.hypotheses a and b are chained at the same place.
bool same_place(const ReadFacts& facts, const Lifecycle& lc, std::size_t a,
                std::size_t b);

// Appends lc.hypotheses[i]'s realized regions, one per segment of `result`.
void append_dp_hits(const Lifecycle& lc, std::size_t i,
                    const RnaSpliceRealizationResult& result,
                    std::vector<RnaParentHit>& hits);

// Rebuilds lc.pieces from the after-election regions `dp` and the chain stage's `chain`.
void decide_pieces(const ReadFacts& facts, Lifecycle& lc,
                   const std::vector<RnaParentHit>& dp,
                   const std::vector<RnaParentHit>& chain);

} // namespace fa::cpu::lr::rna

/**** lifecycle ****/
// Stage 5's rival lifecycle stages (hit.cpp), called by map.cpp's run_lifecycle in this
// order.

namespace fa::cpu::lr::rna {

// (2) and (3): the hypotheses, the incumbent first, classified and admitted.
void build_hypotheses(const ReadFacts& facts, const ReadState& rs,
                      const AttemptState& as, Lifecycle& lc,
                      const ExactAnchorPath& anchor_path,
                      const RnaRivalLifecycleConfig& lifecycle_cfg,
                      const RnaSpliceRealizationResult& realized);

// (4): the realization queue, competitors within the budget then the chain-stage
// co-primaries.
std::vector<std::size_t> queue_realizations(const ReadFacts& facts, ReadState& rs,
                                            Lifecycle& lc,
                                            const RnaRivalLifecycleConfig& lifecycle_cfg);

// The rank recalibration (minimap2's mm_update_dp_max) of the realized hypotheses.
void rank_recalibrate(const ReadFacts& facts, Lifecycle& lc,
                      const RnaSpliceRealizationResult& realized);

// (5) up to the adoption: the election, the shadow rule, the dp2 owner, the pieces.
void elect_and_decide_pieces(const ReadFacts& facts, Lifecycle& lc,
                             const RnaSpliceRealizationResult& realized);

} // namespace fa::cpu::lr::rna

/**** committed_evidence ****/

namespace fa::cpu::lr::rna {

// The committed record's chain MAPQ evidence in the committed frame: its chain, vote,
// coverage, identity and realized DP maxima, and every other hypothesis as a rival.
void committed_chain_evidence(const ReadFacts& facts, const ReadState& rs,
                              const Lifecycle& lc, const Committed& cf,
                              const std::optional<int>& committed_dp_max,
                              const std::optional<int>& rival_dp_max,
                              const RnaSpliceRealizationResult& realized,
                              RnaChainMapqEvidence& chain_evidence);

} // namespace fa::cpu::lr::rna

/**** leftover ****/

namespace fa::cpu::lr::rna {

// The leftover, after the primary's MAPQ: the forward intervals no printed region
// covers, each contig and lane's pool chained and realized as a co-primary, and the
// after-election parent rule rerun with their regions behind the elected family's.
void run_leftover(const ReadFacts& facts, ReadState& rs, Lifecycle& lc,
                  const RnaSpliceRealizationResult& realized,
                  const RnaSpliceRealizationRequest& request);

} // namespace fa::cpu::lr::rna

/**** pieces ****/

namespace fa::cpu::lr::rna {

// The pieces' records: each piece's region segment priced by the chain MAPQ on its own
// region against its children and stored as one flat supplementary family; the
// places the pieces print go on the wire.
void price_pieces(const ReadFacts& facts, Lifecycle& lc,
                  const RnaSpliceRealizationResult& realized,
                  const RnaChainMapqEvidence& chain_evidence,
                  std::vector<RnaSpliceRealizationResult>& chimeric_results,
                  std::vector<RnaChimericFamily>& chimeric_families);

} // namespace fa::cpu::lr::rna

/**** runner_ups ****/

namespace fa::cpu::lr::rna {

// The runner-ups: one realized family per surviving competitor of the committed record,
// the best priced first and one per place, moved out of the hypotheses.
void attach_runner_ups(const ReadFacts& facts, ReadState& rs, Lifecycle& lc);

} // namespace fa::cpu::lr::rna
