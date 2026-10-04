#include "family_projection.h"

#include "chain_mapq.h"  // kDnaChainMapqStudyBridgeOwner
#include "family_realization.h"  // dna_selected_sibling_path, the seam test
#include "placement_chaining.h"
#include "record_family.h"
#include "../dp/params.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

// One step of the fuzzy walk that gives a projected record its PAF columns 10
// and 11; false when the step does not advance on both axes.
bool walk_step(const chaining::Anchor& previous,
               const chaining::Anchor& current, int& matches, int& length) {
  const int target_length = current.r - previous.r;
  const int query_length = current.q - previous.q;
  if (target_length <= 0 || query_length <= 0)
    return false;
  length += std::max(target_length, query_length);
  matches += target_length > current.span && query_length > current.span
                 ? current.span
                 : std::min(target_length, query_length);
  return true;
}

// Columns 10 and 11 from the walk less the record's trims, floored by both
// spans.
void walk_columns(int walk_matches, int walk_length, int trims,
                  int query_span, int reference_span, int& matches,
                  int& block_length) {
  block_length =
      std::max({1, walk_length - trims, query_span, reference_span});
  matches = std::clamp(walk_matches - trims, 1, block_length);
}

// A slice of a sibling's path scores its share of the sibling's chain score by
// anchor count, as minimap2 splits a chain's score (mm_split_reg).
int slice_score(const DnaPlacementCandidateChain& chain, int sibling,
                std::size_t slice_anchors) {
  const std::size_t which = static_cast<std::size_t>(sibling);
  const float share = static_cast<float>(slice_anchors) /
                      static_cast<float>(chain.sibling_paths[which].size());
  return static_cast<int>(chain.sibling_scores[which] * share + .499);
}

bool project_chained_block_impl(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block,
    AlignResult& output, DnaJoinPiece* piece = nullptr) {
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
      block.query_tile_end > family.tile_count)
    return false;
  if (block.supporting_tiles < kDnaPostCommitRecordMinBlockTiles)
    return false;

  const int forward_begin = ::fa::cpu::voting::query_tile_begin(
      block.query_tile_begin, family.read_length, family.seed_length,
      family.tile_count);
  const int forward_end =
      block.query_tile_end == family.tile_count
          ? family.read_length
          : ::fa::cpu::voting::query_tile_begin(
                block.query_tile_end, family.read_length,
                family.seed_length, family.tile_count);
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
        family.seed_length, family.tile_count);
    if (tile >= block.query_tile_begin && tile < block.query_tile_end)
      anchors.push_back(&anchor);
  }
  // A block that holds none of the primary's anchors is projected from a
  // sibling path, as realization does (dna_selected_sibling_path); the record
  // then carries the slice's share of that sibling's chain score.
  std::vector<chaining::Anchor> sibling_slice;
  int sibling = -1;
  if (anchors.empty()) {
    sibling_slice = dna_selected_sibling_path(chain, oriented_begin,
                                              oriented_end, &sibling);
    for (const chaining::Anchor& anchor : sibling_slice)
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

  int walk_matches = first.span;
  int walk_length = first.span;
  for (std::size_t index = 1; index < anchors.size(); ++index) {
    if (!walk_step(*anchors[index - 1], *anchors[index], walk_matches,
                   walk_length))
      return false;
  }
  int fuzzy_matches = 0;
  int fuzzy_block_length = 0;
  walk_columns(walk_matches, walk_length, left_trim + right_trim,
               query_end - query_begin, reference_end - reference_begin,
               fuzzy_matches, fuzzy_block_length);

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
  // A block from the owner's own chain keeps the whole chain score.
  output.score = std::max(
      1, sibling < 0 ? chain.chain_score
                     : slice_score(chain, sibling, sibling_slice.size()));
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
  if (piece != nullptr) {
    piece->record = output;
    piece->candidate = block.candidate;
    piece->chromosome = candidate->peak.chr;
    piece->reverse = candidate->peak.is_rc;
    piece->first = first;
    piece->last = last;
    piece->walk_matches = walk_matches;
    piece->walk_length = walk_length;
  }
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
    DnaJoinPiece piece;
    if (!project_chained_block_impl(
            context, family, *chain, block, projected, &piece)) {
      refuse(result, DnaProjectionFallbackReason::InvalidBlockGeometry);
      return result;
    }
    result.pieces.push_back(std::move(piece));
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

namespace {

// A record's first and last anchors as the block path the seam test reads.
ordered_anchor::OrderedAnchorPath corner_path(const DnaJoinPiece& piece,
                                              int k) {
  ordered_anchor::OrderedAnchorPath path;
  path.minimizer_k = k;
  path.selected = {piece.first, piece.last};
  path.realization_end = path.selected.size();
  return path;
}

// Whether the -c lane would bridge `left` and `right`, adjacent in
// forward-query order (bridge_geometry in family_realization.cpp).
bool seam_joins(const DnaContext& context, const DnaJoinPiece& left,
                const DnaJoinPiece& right) {
  if (left.chromosome != right.chromosome || left.reverse != right.reverse)
    return false;
  const int k = context.ref.index->k();
  const ordered_anchor::OrderedAnchorPath left_path = corner_path(left, k);
  const ordered_anchor::OrderedAnchorPath right_path = corner_path(right, k);
  if (dna_family_seam_has_duplicate_anchor(left_path, right_path,
                                           left.reverse))
    return false;
  const ordered_anchor::AdjustedEndpoint left_begin =
      ordered_anchor::adjusted_endpoint(left_path, left.first);
  const ordered_anchor::AdjustedEndpoint left_end =
      ordered_anchor::adjusted_endpoint(left_path, left.last);
  const ordered_anchor::AdjustedEndpoint right_begin =
      ordered_anchor::adjusted_endpoint(right_path, right.first);
  const ordered_anchor::AdjustedEndpoint right_end =
      ordered_anchor::adjusted_endpoint(right_path, right.last);
  const ::fa::cpu::DpMapOpt dp;
  const DnaFamilyJoinPlan plan = plan_dna_family_join(DnaFamilyJoinGeometry{
      left.chromosome, right.chromosome, left.reverse, right.reverse,
      left.record.query_start, left.record.query_end,
      right.record.query_start, right.record.query_end, left_begin.query,
      left_end.query, right_begin.query, right_end.query, left_begin.target,
      left_end.target, right_begin.target, right_end.target,
      context.opts.cigar_dp_max_gap, context.opts.cigar_dp_max_gap,
      dp.max_sw_mat});
  return plan.kind != DnaFamilyJoinKind::NotEligible;
}

// Joins `next`, the record after `into` in forward-query order, onto it.
// False, `into` untouched, when the walk does not advance across the seam.
bool merge_pieces(DnaJoinPiece& into, const DnaJoinPiece& next) {
  const DnaJoinPiece& low = into.reverse ? next : into;
  const DnaJoinPiece& high = into.reverse ? into : next;
  int walk_matches = low.walk_matches + high.walk_matches - high.first.span;
  int walk_length = low.walk_length + high.walk_length - high.first.span;
  if (!walk_step(low.last, high.first, walk_matches, walk_length))
    return false;
  const int read_length = into.record.read_len;
  const int query_begin = into.reverse ? read_length - low.record.query_end
                                       : low.record.query_start;
  const int query_end = into.reverse ? read_length - high.record.query_start
                                     : high.record.query_end;
  const int trims =
      (query_begin - low.first.q) + (high.last.q_end() - query_end);
  const int reference_begin = low.record.pos;
  const int reference_end = high.record.target_end;
  const chaining::Anchor first = low.first;
  const chaining::Anchor last = high.last;
  int matches = 0;
  int block_length = 0;
  walk_columns(walk_matches, walk_length, trims, query_end - query_begin,
               reference_end - reference_begin, matches, block_length);
  into.record.query_end = next.record.query_end;
  into.record.pos = reference_begin;
  into.record.target_end = reference_end;
  into.record.target_regions.assign(1, {reference_begin, reference_end});
  into.record.matches = matches;
  into.record.block_len = block_length;
  into.first = first;
  into.last = last;
  into.walk_matches = walk_matches;
  into.walk_length = walk_length;
  return true;
}

}  // namespace

void dna_join_map_only_family(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    std::vector<DnaJoinPiece> pieces, dna::Result& realized,
    ::fa::cpu::voting::CandidateId& primary_candidate) {
  namespace voting = ::fa::cpu::voting;
  if (pieces.size() < 2 || !realized.mapped())
    return;
  // Forward-query order. The terminal-clip records sort in after the blocks,
  // so a pair joins only when no clip record lies between them.
  std::vector<std::pair<int, std::size_t>> order;
  for (std::size_t index = 0; index < pieces.size(); ++index)
    order.emplace_back(pieces[index].record.query_start, index);
  for (std::size_t index = 0; index < realized.supplementary.size(); ++index)
    if (realized.supplementary_candidates[index] == voting::kNullCandidate)
      order.emplace_back(realized.supplementary[index].query_start,
                         pieces.size() + index);
  std::sort(order.begin(), order.end());

  const bool bridge_owner =
      (dna_chain_mapq_study_bits(context.opts.chain_mapq_hifi_margin) &
       kDnaChainMapqStudyBridgeOwner) != 0;
  const auto chain_score = [&placement](voting::CandidateId id) {
    const DnaPlacementCandidateChain* chain = placement.find(id);
    return chain == nullptr ? 0 : chain->chain_score;
  };
  // The first piece of each joined record, which accumulates the rest.
  std::vector<std::size_t> runs;
  for (std::size_t at = 0; at < order.size(); ++at) {
    const std::size_t index = order[at].second;
    if (index >= pieces.size())
      continue;
    const std::size_t before = at == 0 ? pieces.size() : order[at - 1].second;
    const DnaJoinPiece& current = pieces[index];
    if (before < pieces.size() &&
        seam_joins(context, pieces[before], current)) {
      DnaJoinPiece& into = pieces[runs.back()];
      if (merge_pieces(into, current)) {
        // The -c lane's owner of a bridged unit: the leftmost block's, or
        // under the HiFi rule the higher chain score. AS is the owner's.
        if (bridge_owner &&
            chain_score(current.candidate) > chain_score(into.candidate)) {
          into.candidate = current.candidate;
          into.record.score = current.record.score;
        }
        continue;
      }
    }
    runs.push_back(index);
  }
  if (runs.size() == pieces.size())
    return;

  std::vector<DnaSegmentRecord> records;
  records.reserve(runs.size());
  for (const std::size_t index : runs)
    records.push_back(
        {std::move(pieces[index].record), pieces[index].candidate});
  std::vector<AlignResult> clips;
  for (std::size_t index = 0; index < realized.supplementary.size(); ++index)
    if (realized.supplementary_candidates[index] == voting::kNullCandidate)
      clips.push_back(std::move(realized.supplementary[index]));
  DnaRecordFamily assembled = assemble_dna_record_family(std::move(records));
  std::vector<AlignResult> secondary = std::move(realized.secondary);
  const int mapq = realized.mapq;
  static_cast<AlignResult&>(realized) = std::move(assembled.primary);
  realized.mapq = mapq;
  realized.secondary = std::move(secondary);
  realized.supplementary.clear();
  realized.supplementary_candidates.clear();
  for (DnaSegmentRecord& record : assembled.supplementary) {
    realized.supplementary.push_back(std::move(record.alignment));
    realized.supplementary_candidates.push_back(record.candidate);
  }
  for (AlignResult& clip : clips) {
    realized.supplementary.push_back(std::move(clip));
    realized.supplementary_candidates.push_back(voting::kNullCandidate);
  }
  primary_candidate = assembled.primary_candidate;
}

}  // namespace fa::cpu::lr
