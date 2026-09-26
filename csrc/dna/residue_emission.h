#pragma once

#include "geometry_work.h"
#include "context.h"
#include "placement_chaining.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

// Terminal-clip geometry: a chain on the primary's contig whose nearest
// reference edge lies within 10 kb of the primary's.
inline constexpr std::int64_t kDnaResidueEmissionLocalBp = 10000;
inline constexpr int kDnaResidueEmissionMinAnchors = 10;
inline constexpr int kDnaResidueEmissionMinSpanBp = 100;

std::int64_t dna_residue_reference_gap(std::int64_t left_begin,
                                       std::int64_t left_end,
                                       std::int64_t right_begin,
                                       std::int64_t right_end) noexcept;

struct DnaResidueBlockOutcome {
  bool accepted = false;
  int supporting_tiles = 0;
  explicit operator bool() const noexcept { return accepted; }
};

// Converts a post-commit residue chain to the forward-query tile block it
// covers. A single-tile block must also reach `min_anchors` anchors and
// kDnaResidueEmissionMinSpanBp query bases; a caller with a scaled admission
// bar passes the same floor here.
DnaResidueBlockOutcome dna_residue_chain_block(
    const DnaPlacementFamily& family, const DnaPlacementCandidate& candidate,
    const DnaPlacementCandidateChain& chain,
    ::fa::cpu::voting::QueryBlock& block,
    int min_anchors = kDnaResidueEmissionMinAnchors);

// Whether post-commit residue emission may run. It never runs after an
// alternative hypothesis was promoted.
bool dna_residue_emission_boundary(bool mapped, bool placement_chaining_ran,
                                   bool alternative_promoted,
                                   bool cigar_lane,
                                   bool primary_has_cigar) noexcept;

// Emits one residue chain as a supplementary of `primary`, realized with
// CIGAR output and projected otherwise. `window` is caller-owned scratch whose
// selected blocks are overwritten.
bool dna_emit_residue_record(
    const DnaContext& context, DnaPlacementFamily& window,
    const DnaPlacementChainingResult& placement,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query, bool cigar_lane,
    dna::Result& primary, int& ksw2_attempts, std::int64_t& estimated_cells,
    DnaGeometryWork& geometry);

}  // namespace fa::cpu::lr
