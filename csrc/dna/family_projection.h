#pragma once

#include "context.h"
#include "placement_family_adapter.h"
#include "result.h"

#include <cstdint>

namespace fa::cpu::lr {

struct DnaPlacementChainingResult;
struct DnaPlacementCandidateChain;

// Projects the part of an accepted chain inside `block` (forward-query tiles)
// to a record without CIGAR. Returns false on invalid geometry.
bool dna_project_chained_block(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block, AlignResult& output);

enum class DnaProjectionFallbackReason : std::uint8_t {
  None,
  ChainingUnavailableOrRefused,
  NoProjectedRecords,
  MissingOrInvalidSelectedChain,
  InvalidBlockGeometry,
  OverlappingProjectedQueryBlocks,
  TransactionValidationFailure,
  InvalidPlacementFamily,
};

const char* dna_projection_fallback_name(
    DnaProjectionFallbackReason reason) noexcept;

struct DnaFamilyProjectionResult {
  dna::Result output;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  bool committed = false;
  DnaProjectionFallbackReason fallback =
      DnaProjectionFallbackReason::None;

  bool fell_back() const noexcept {
    return fallback != DnaProjectionFallbackReason::None;
  }
};

// Returns the committed projection, or the incumbent with score 0 and the
// refusal reason recorded.
dna::Result settle_map_only_projection(DnaFamilyProjectionResult result);

DnaFamilyProjectionResult project_map_only_placement_family(
    const DnaContext& context,
    const DnaPlacementChainingResult& placement,
    const dna::Result& incumbent);

}  // namespace fa::cpu::lr
