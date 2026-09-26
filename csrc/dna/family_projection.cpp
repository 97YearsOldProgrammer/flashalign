#include "family_projection.h"

#include "placement_chaining.h"
#include "record_family.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

bool project_chained_block_impl(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block,
    AlignResult& output) {
  const DnaPlacementCandidate* candidate = family.find(block.candidate);
  if (candidate == nullptr ||
      chain.status != DnaPlacementChainStatus::Accepted ||
      chain.primary.empty() || context.ref.names == nullptr ||
      candidate->peak.chr < 0 ||
      candidate->peak.chr >=
          static_cast<int>(context.ref.names->size()) ||
      candidate->peak.chr >= context.ref.contig_count() ||
      block.query_tile_begin < 0 ||
      block.query_tile_end <= block.query_tile_begin ||
      block.query_tile_end > ::fa::cpu::voting::kQueryTileCount)
    return false;
  if (block.supporting_tiles < kDnaPostCommitRecordMinBlockTiles)
    return false;

  const int forward_begin = ::fa::cpu::voting::query_tile_begin(
      block.query_tile_begin, family.read_length, family.seed_length);
  const int forward_end =
      block.query_tile_end == ::fa::cpu::voting::kQueryTileCount
          ? family.read_length
          : ::fa::cpu::voting::query_tile_begin(
                block.query_tile_end, family.read_length,
                family.seed_length);
  const int oriented_begin =
      candidate->peak.is_rc ? family.read_length - forward_end
                            : forward_begin;
  const int oriented_end =
      candidate->peak.is_rc ? family.read_length - forward_begin
                            : forward_end;
  if (forward_begin < 0 || forward_end <= forward_begin ||
      forward_end > family.read_length || oriented_begin < 0 ||
      oriented_end <= oriented_begin)
    return false;

  std::vector<const chaining::Anchor*> anchors;
  anchors.reserve(chain.primary.size());
  for (const chaining::Anchor& anchor : chain.primary) {
    const int tile = dna_forward_query_tile(
        anchor.q, candidate->peak.is_rc, family.read_length,
        family.seed_length);
    if (tile >= block.query_tile_begin && tile < block.query_tile_end)
      anchors.push_back(&anchor);
  }
  if (anchors.empty())
    return false;

  const chaining::Anchor& first = *anchors.front();
  const chaining::Anchor& last = *anchors.back();
  const int query_begin = std::max(oriented_begin, static_cast<int>(first.q));
  const int query_end = std::min(oriented_end, static_cast<int>(last.q_end()));
  if (query_end <= query_begin)
    return false;
  const int left_trim = query_begin - first.q;
  const int right_trim = last.q_end() - query_end;
  const int reference_begin = first.r + left_trim;
  const int reference_end = last.r_end() - right_trim;
  // Contig lengths come from the index offsets: map-only projection never
  // reads a reference base.
  const int reference_length =
      static_cast<int>(context.ref.contig_length(candidate->peak.chr));
  if (reference_begin < 0 || reference_end <= reference_begin ||
      reference_end > reference_length)
    return false;

  int fuzzy_matches = first.span;
  int fuzzy_block_length = first.span;
  for (std::size_t index = 1; index < anchors.size(); ++index) {
    const chaining::Anchor& previous = *anchors[index - 1];
    const chaining::Anchor& current = *anchors[index];
    const int target_length = current.r - previous.r;
    const int query_length = current.q - previous.q;
    if (target_length <= 0 || query_length <= 0)
      return false;
    fuzzy_block_length += std::max(target_length, query_length);
    fuzzy_matches +=
        target_length > current.span && query_length > current.span
            ? current.span
            : std::min(target_length, query_length);
  }
  fuzzy_block_length =
      std::max({1, fuzzy_block_length - left_trim - right_trim,
                query_end - query_begin,
                reference_end - reference_begin});
  fuzzy_matches = std::clamp(
      fuzzy_matches - left_trim - right_trim, 1, fuzzy_block_length);

  output = {};
  output.read_len = family.read_length;
  output.chromosome =
      (*context.ref.names)[static_cast<std::size_t>(candidate->peak.chr)];
  output.pos = reference_begin;
  output.target_end = reference_end;
  output.query_start =
      candidate->peak.is_rc ? family.read_length - query_end : query_begin;
  output.query_end =
      candidate->peak.is_rc ? family.read_length - query_begin : query_end;
  output.is_reverse = candidate->peak.is_rc;
  output.score = std::max(1, chain.chain_score);
  output.matches = fuzzy_matches;
  output.block_len = fuzzy_block_length;
  output.mapq = 0;
  output.origin = AlignmentOrigin::DnaQueryPartition;
  output.target_regions.emplace_back(reference_begin, reference_end);
  if (!(output.mapped() && output.cigar.empty() &&
        !output.alignment_accounting_valid && output.edit_distance < 0 &&
        output.query_start >= forward_begin &&
        output.query_end <= forward_end))
    return false;
  return true;
}

void refuse(DnaFamilyProjectionResult& result,
            DnaProjectionFallbackReason reason) {
  result.fallback = reason;
}

}  // namespace

bool dna_project_chained_block(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block, AlignResult& output) {
  return project_chained_block_impl(context, family, chain, block, output);
}

const char* dna_projection_fallback_name(
    DnaProjectionFallbackReason reason) noexcept {
  switch (reason) {
    case DnaProjectionFallbackReason::None:
      return "none";
    case DnaProjectionFallbackReason::ChainingUnavailableOrRefused:
      return "chaining_unavailable_or_refused";
    case DnaProjectionFallbackReason::NoProjectedRecords:
      return "no_projected_records";
    case DnaProjectionFallbackReason::MissingOrInvalidSelectedChain:
      return "missing_or_invalid_selected_chain";
    case DnaProjectionFallbackReason::InvalidBlockGeometry:
      return "invalid_block_geometry";
    case DnaProjectionFallbackReason::OverlappingProjectedQueryBlocks:
      return "overlapping_projected_query_blocks";
    case DnaProjectionFallbackReason::TransactionValidationFailure:
      return "transaction_validation_failure";
    case DnaProjectionFallbackReason::InvalidPlacementFamily:
      return "invalid_placement_family";
    }
  return "unknown";
}

DnaFamilyProjectionResult project_map_only_placement_family(
    const DnaContext& context,
    const DnaPlacementChainingResult& placement,
    const dna::Result& incumbent) {
  DnaFamilyProjectionResult result;
  result.output = incumbent;
  const DnaPlacementFamily& family = placement.family;
  if (!placement.accepted) {
    refuse(
        result,
        DnaProjectionFallbackReason::ChainingUnavailableOrRefused);
    return result;
  }
  if (!family.valid) {
    refuse(result, DnaProjectionFallbackReason::InvalidPlacementFamily);
    return result;
  }
  if (!incumbent.mapped() || !incumbent.cigar.empty() ||
      incumbent.alignment_accounting_valid) {
    refuse(result,
           DnaProjectionFallbackReason::TransactionValidationFailure);
    return result;
  }

  std::vector<DnaSegmentRecord> records;
  for (const auto& block : family.partition.selected.blocks) {
    if (block.candidate == ::fa::cpu::voting::kNullCandidate) continue;
    const DnaPlacementCandidateChain* chain =
        placement.find(block.candidate);
    AlignResult projected;
    if (chain == nullptr || !chain->exact ||
        chain->status != DnaPlacementChainStatus::Accepted ||
        chain->primary.empty()) {
      refuse(
          result,
          DnaProjectionFallbackReason::MissingOrInvalidSelectedChain);
      return result;
    }
    if (!project_chained_block_impl(
            context, family, *chain, block, projected)) {
      refuse(result, DnaProjectionFallbackReason::InvalidBlockGeometry);
      return result;
    }
    records.push_back({std::move(projected), block.candidate});
  }
  if (records.empty()) {
    refuse(result, DnaProjectionFallbackReason::NoProjectedRecords);
    return result;
  }
  DnaRecordFamily assembled = assemble_dna_record_family(std::move(records));
  if (!assembled.valid) {
    refuse(result,
           DnaProjectionFallbackReason::OverlappingProjectedQueryBlocks);
    return result;
  }

  result.output = incumbent;
  static_cast<AlignResult&>(result.output) = std::move(assembled.primary);
  result.primary_candidate = assembled.primary_candidate;
  result.output.mapq = incumbent.mapq;
  result.output.supplementary.clear();
  result.output.supplementary_candidates.clear();
  for (DnaSegmentRecord& record : assembled.supplementary) {
    result.output.supplementary.push_back(std::move(record.alignment));
    result.output.supplementary_candidates.push_back(record.candidate);
  }
  result.committed = true;
  return result;
}

dna::Result settle_map_only_projection(DnaFamilyProjectionResult result) {
  dna::Result output = std::move(result.output);
  if (!result.committed)
    output.score = 0;
  return output;
}

}  // namespace fa::cpu::lr
