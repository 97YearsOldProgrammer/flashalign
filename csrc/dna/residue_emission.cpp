#include "residue_emission.h"

#include "family_projection.h"
#include "family_realization.h"
#include "placement_family_adapter.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <limits>
#include <utility>

namespace fa::cpu::lr {
namespace {

// A terminal-clip record owns no block: it carries kNullCandidate and
// inherits the primary's MAPQ.
void append_emission(AlignResult record, dna::Result& primary) {
  std::vector<AlignResult> children = std::move(record.supplementary);
  record.supplementary.clear();
  record.secondary.clear();
  primary.supplementary.push_back(std::move(record));
  primary.supplementary_candidates.push_back(::fa::cpu::voting::kNullCandidate);
  for (AlignResult& child : children) {
    child.supplementary.clear();
    child.secondary.clear();
    primary.supplementary.push_back(std::move(child));
    primary.supplementary_candidates.push_back(
        ::fa::cpu::voting::kNullCandidate);
  }
}

}  // namespace

bool dna_emit_residue_record(
    const DnaContext& context, DnaPlacementFamily& window,
    const DnaPlacementChainingResult& placement,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query, bool cigar_lane,
    dna::Result& primary, int& ksw2_attempts, std::int64_t& estimated_cells,
    DnaGeometryWork& geometry) {
  if (cigar_lane) {
    window.partition.selected.blocks.assign(1, block);
    DnaFamilyRealizationRequest request;
    request.family = &window;
    request.placement = &placement;
    request.forward_query = &forward_query;
    request.reverse_query = &reverse_query;
    DnaFamilyRealizationOutcome outcome =
        realize_full_cigar_family(context, request);
    ksw2_attempts += outcome.ksw2_attempts;
    estimated_cells += outcome.estimated_cells;
    geometry.add(outcome.geometry);
    if (!outcome.accepted() || outcome.output.cigar.empty()) return false;
    append_emission(std::move(outcome.output), primary);
    return true;
  }
  AlignResult emission;
  if (!dna_project_chained_block(context, window, chain, block, emission))
    return false;
  append_emission(std::move(emission), primary);
  return true;
}

std::int64_t dna_residue_reference_gap(std::int64_t left_begin,
                                       std::int64_t left_end,
                                       std::int64_t right_begin,
                                       std::int64_t right_end) noexcept {
  return std::max<std::int64_t>(0, std::max(left_begin, right_begin) -
                                      std::min(left_end, right_end));
}

DnaResidueBlockOutcome dna_residue_chain_block(
    const DnaPlacementFamily& family, const DnaPlacementCandidate& candidate,
    const DnaPlacementCandidateChain& chain,
    ::fa::cpu::voting::QueryBlock& block, int min_anchors) {
  DnaResidueBlockOutcome outcome;
  ::fa::cpu::voting::QueryTileMask tiles;
  int begin = std::numeric_limits<int>::max();
  int end = std::numeric_limits<int>::min();
  for (const chaining::Anchor& anchor : chain.primary) {
    const int tile = dna_forward_query_tile(
        anchor.q, candidate.peak.is_rc, family.read_length,
        family.seed_length, family.tile_count);
    if (tile < 0 || tile >= family.tile_count)
      return outcome;
    tiles.set(tile);
    begin = std::min(begin, tile);
    end = std::max(end, tile + 1);
  }
  if (begin >= end)
    return outcome;
  const int supporting = static_cast<int>(tiles.count());
  outcome.supporting_tiles = supporting;
  if (supporting < kDnaPostCommitRecordMinBlockTiles)
    return outcome;
  if (supporting < 2) {
    const int anchors = static_cast<int>(chain.primary.size());
    const int span =
        std::max(0, chain.forward_query_end - chain.forward_query_begin);
    if (anchors < min_anchors || span < kDnaResidueEmissionMinSpanBp)
      return outcome;
  }
  block = {};
  block.candidate = candidate.id;
  block.query_tile_begin = begin;
  block.query_tile_end = end;
  block.supporting_tiles = supporting;
  outcome.accepted = true;
  return outcome;
}

bool dna_residue_emission_boundary(bool mapped, bool placement_chaining_ran,
                                   bool alternative_promoted,
                                   bool cigar_lane,
                                   bool primary_has_cigar) noexcept {
  return mapped && placement_chaining_ran && !alternative_promoted &&
         cigar_lane == primary_has_cigar;
}

}  // namespace fa::cpu::lr
