#pragma once

#include "context.h"
#include "placement_family_adapter.h"
#include "result.h"
#include "../chaining/anchor.h"

#include <cstdint>
#include <vector>

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
  NoOwnerChain,
  NoProjectedRecords,
  MissingOrInvalidSelectedChain,
  InvalidBlockGeometry,
  OverlappingProjectedQueryBlocks,
  TransactionValidationFailure,
  InvalidPlacementFamily,
};

const char* dna_projection_fallback_name(
    DnaProjectionFallbackReason reason) noexcept;

// A projected block record with what the join below reads: its owner, contig
// and strand, the first and last anchors it was projected from (oriented
// order), and the fuzzy walk over those anchors before the record's trims.
struct DnaJoinPiece {
  AlignResult record;
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  int chromosome = -1;
  bool reverse = false;
  chaining::Anchor first;
  chaining::Anchor last;
  int walk_matches = 0;
  int walk_length = 0;
  // The family.block_parts entry, -1 for none.
  int part = -1;
};

struct DnaFamilyProjectionResult {
  dna::Result output;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  bool committed = false;
  DnaProjectionFallbackReason fallback =
      DnaProjectionFallbackReason::None;
  // Every projected block record, in partition order.
  std::vector<DnaJoinPiece> pieces;

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

// Joins consecutive block records of the committed map-only family, in
// forward-query order, into one record wherever they pass the bridge test (no
// duplicate seam anchor, then plan_dna_family_join) with the centers of the
// records' corner anchors as the cuts. A record that owns no block never
// joins. The joined record spans both, its columns 10 and 11 are the walk over
// both anchor runs, and it is owned by the piece of the higher whole chain
// score, the leftmost on a tie. No DP and no reference base. `realized` and `primary_candidate`
// change only when two records join.
void dna_join_map_only_family(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    std::vector<DnaJoinPiece> pieces, dna::Result& realized,
    ::fa::cpu::voting::CandidateId& primary_candidate);

}  // namespace fa::cpu::lr
