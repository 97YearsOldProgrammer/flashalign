#pragma once

#include "placement_family_adapter.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

struct DnaContext;
struct DnaPlacementCandidateChain;
struct DnaPlacementChainingResult;

enum class DnaAlternativeRefusal : std::uint8_t {
  None,
  NoStableFamily,
  NotSingleBlock,
  NoCandidate,
  NoLocus,
  NoChain,
  ShadowOnly,
  DisjointOnly,
  IncredibleOnly,
  ExactRestoreFailed,
};

struct DnaAlternativeSelection {
  ::fa::cpu::voting::CandidateId incumbent =
      ::fa::cpu::voting::kNullCandidate;
  ::fa::cpu::voting::CandidateId candidate =
      ::fa::cpu::voting::kNullCandidate;
  DnaAlternativeRefusal refusal = DnaAlternativeRefusal::NoStableFamily;
};

// Picks at most one alternative to the owner of a single-block family, from
// the screening chains. A candidate qualifies when it is not selected, is not
// a shadow of a selected locus (same contig and strand, vote start within
// max(2 kb, read length / 10)), overlaps the block by at least half the
// shorter query span, and scores at least -p (0.8) times the owner or within
// 2k of it. Best by (score, anchors, overlap, lowest id).
DnaAlternativeSelection select_dna_alternative_hypothesis(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable);

// Every candidate select_dna_alternative_hypothesis would qualify, best first
// in its order, at most `limit`: the first is its choice.
std::vector<::fa::cpu::voting::CandidateId> rank_dna_alternative_hypotheses(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable, std::size_t limit);

}  // namespace fa::cpu::lr
