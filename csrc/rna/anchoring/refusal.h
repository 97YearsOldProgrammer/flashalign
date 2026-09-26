#pragma once

#include <cstdint>

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
