#include "select_exact_anchor_path.h"

#include "../../chaining/colinear_chain.h"
#include "../../core/checked_range.h"
#include "query_span.h" // RnaQuerySpan, rna_query_spans_compete

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

// One chain's query span [min q, max q_end) in the recurrence's mapping frame. Anchors
// are query-ordered but spans vary, so the maximum end is not always the last anchor's.
// An invalid traceback yields the invalid span, which rna_query_spans_compete treats as
// competing.
RnaQuerySpan chain_query_span(const chaining::ChainResult& chained,
                              const chaining::Chain& chain) {
  RnaQuerySpan span;
  for (const int32_t index : chain.idx) {
    if (index < 0 || static_cast<std::size_t>(index) >= chained.anchors.size())
      return RnaQuerySpan();
    const chaining::Anchor& anchor =
        chained.anchors[static_cast<std::size_t>(index)];
    if (span.begin < 0 || anchor.q < span.begin)
      span.begin = anchor.q;
    if (anchor.q_end() > span.end)
      span.end = anchor.q_end();
  }
  return span;
}

} // namespace

ExactAnchorPathResult select_exact_anchor_path(
    uint64_t read_key, int32_t query_length, int32_t reference_id,
    bool mapping_reverse, const CandidatePool& pool,
    const chaining::SplicedChainParams& params, SiblingChainStats* siblings,
    std::vector<DisjointAnchorPath>* disjoint) {
  if (disjoint)
    disjoint->clear();
  ExactAnchorPathResult result;
  ExactAnchorPath& bundle = result.bundle;
  bundle.read_key = read_key;
  bundle.query_length = query_length;
  bundle.reference_id = reference_id;
  bundle.mapping_reverse = mapping_reverse;
  bundle.candidate_pool_hash = pool.content_hash;
  if (pool.refused) {
    result.refused = true;
    result.refusal = pool.refusal == AnchoringRefusal::None
                         ? AnchoringRefusal::PoolRefused
                         : pool.refusal;
    return result;
  }
  int pool_count = 0;
  if (!mapping::checked_size_to_int(pool.anchors.size(), pool_count)) {
    result.refused = true;
    result.refusal = AnchoringRefusal::PoolCountDomain;
    return result;
  }
  for (const FineAnchor& anchor : pool.anchors) {
    if (!anchor.final_pool) continue;
    if (anchor.anchor_id >
        static_cast<uint32_t>(mapping::kSignedCoordinateMaximum)) {
      result.refused = true;
      result.refusal = AnchoringRefusal::StableIdDomain;
      return result;
    }
    SelectedPoolAnchor raw;
    raw.stable_id = anchor.anchor_id;
    raw.query_begin = anchor.query_begin;
    raw.reference_begin = anchor.reference_begin;
    raw.span = anchor.span;
    raw.source_seed_id = anchor.source_seed_id;
    raw.source_occurrence = anchor.source_occurrence;
    raw.ignore_flag = anchor.ignore_flag;
    raw.tandem_flag = anchor.tandem_flag;
    raw.long_join_flag = anchor.long_join_flag;
    raw.region_kind_mask = anchor.region_kind_mask;
    raw.region_ids = anchor.region_ids;
    raw.core_node_ids = anchor.compatible_core_nodes;
    bundle.selected_pool_anchors.push_back(std::move(raw));
  }
  std::sort(bundle.selected_pool_anchors.begin(), bundle.selected_pool_anchors.end(),
            selected_pool_anchor_precedes);
  if (bundle.selected_pool_anchors.empty()) {
    result.refused = true;
    result.refusal = AnchoringRefusal::EmptyPool;
    return result;
  }

  std::vector<chaining::Anchor> input;
  input.reserve(bundle.selected_pool_anchors.size());
  for (const SelectedPoolAnchor& anchor : bundle.selected_pool_anchors) {
    uint32_t flags = 0;
    if (anchor.ignore_flag) flags |= chaining::ANCHOR_IGNORE;
    if (anchor.tandem_flag) flags |= chaining::ANCHOR_TANDEM;
    if (anchor.long_join_flag) flags |= chaining::ANCHOR_LONG_JOIN;
    input.push_back({anchor.reference_begin, anchor.query_begin, anchor.span,
                     static_cast<int32_t>(anchor.stable_id), flags});
  }

  const chaining::ChainResult chained =
      chaining::chain_spliced_colinear(std::move(input), params);
  const int best = chaining::best_chain_index(chained);
  // Filled before the `best < 0` refusal so a refused window still reports its chain
  // count.
  if (siblings) {
    siblings->chains_total = static_cast<int32_t>(chained.chains.size());
    // minimap2's mask_level rule (mm_set_parent): a chain competes with the selected
    // one only when they overlap over more than half the shorter query span. A colinear
    // piece the recurrence could not bridge explains different bases and stays out of
    // sib_score. Both chains share the mapping frame, so the spans compare directly.
    const RnaQuerySpan winner_span =
        best >= 0 ? chain_query_span(chained,
                                     chained.chains[static_cast<size_t>(best)])
                  : RnaQuerySpan();
    for (size_t i = 0; i < chained.chains.size(); ++i) {
      if (best >= 0 && i == static_cast<size_t>(best)) continue;
      const chaining::Chain& sibling = chained.chains[i];
      // A chain that cannot beat the running maximum cannot change the answer.
      if (sibling.score <= siblings->sib_score)
        continue;
      if (best >= 0 && !rna_query_spans_compete(
                           winner_span, chain_query_span(chained, sibling)))
        continue;
      siblings->sib_score = sibling.score;
      siblings->sib_count = static_cast<int32_t>(sibling.idx.size());
    }
  }
  if (best < 0) {
    result.refused = true;
    result.refusal = AnchoringRefusal::NoAcceptedChain;
    return result;
  }
  if (chained.anchors.size() != bundle.selected_pool_anchors.size()) {
    result.refused = true;
    result.refusal = AnchoringRefusal::RawSizeDrift;
    return result;
  }
  std::unordered_map<uint32_t, uint32_t> raw_by_id;
  for (uint32_t i = 0; i < bundle.selected_pool_anchors.size(); ++i) {
    if (!raw_by_id.emplace(bundle.selected_pool_anchors[i].stable_id, i).second) {
      result.refused = true;
      result.refusal = AnchoringRefusal::DuplicateStableId;
      return result;
    }
  }
  const chaining::Chain& selected =
      chained.chains[static_cast<size_t>(best)];
  bundle.anchor_path_score = selected.score;
  for (int32_t chained_index : selected.idx) {
    if (chained_index < 0 ||
        static_cast<size_t>(chained_index) >= chained.anchors.size()) {
      result.refused = true;
      result.refusal = AnchoringRefusal::TracebackRange;
      return result;
    }
    const int32_t chained_id =
        chained.anchors[static_cast<size_t>(chained_index)].id;
    if (chained_id < 0) {
      result.refused = true;
      result.refusal = AnchoringRefusal::TracebackIdDomain;
      return result;
    }
    const uint32_t stable_id = static_cast<uint32_t>(chained_id);
    const auto found = raw_by_id.find(stable_id);
    if (found == raw_by_id.end()) {
      result.refused = true;
      result.refusal = AnchoringRefusal::TracebackId;
      return result;
    }
    bundle.selected_raw_indices.push_back(found->second);
  }
  int selected_count = 0;
  if (!mapping::checked_size_to_int(bundle.selected_raw_indices.size(),
                                    selected_count)) {
    result.refused = true;
    result.refusal = AnchoringRefusal::SelectedCountDomain;
    return result;
  }
  bundle.anchor_path_anchor_count = static_cast<uint32_t>(selected_count);
  bundle.anchor_path_hash = hash_exact_anchor_path(bundle);
  AnchoringRefusal validation_refusal = AnchoringRefusal::None;
  if (!validate_exact_anchor_path(bundle, &validation_refusal)) {
    result.refused = true;
    result.refusal = validation_refusal;
  }
  if (disjoint && !result.refused) {
    const RnaQuerySpan winner_span = chain_query_span(chained, selected);
    std::vector<RnaQuerySpan> spans(chained.chains.size());
    for (size_t i = 0; i < chained.chains.size(); ++i)
      spans[i] = chain_query_span(chained, chained.chains[i]);
    for (size_t i = 0; i < chained.chains.size(); ++i) {
      if (i == static_cast<size_t>(best) ||
          rna_query_spans_compete(winner_span, spans[i]))
        continue;
      const chaining::Chain& chain = chained.chains[i];
      DisjointAnchorPath piece;
      ExactAnchorPath& path = piece.bundle;
      path.read_key = bundle.read_key;
      path.query_length = bundle.query_length;
      path.reference_id = bundle.reference_id;
      path.mapping_reverse = bundle.mapping_reverse;
      path.candidate_pool_hash = bundle.candidate_pool_hash;
      path.selected_pool_anchors = bundle.selected_pool_anchors;
      path.anchor_path_score = chain.score;
      bool ok = true;
      for (int32_t chained_index : chain.idx) {
        if (chained_index < 0 ||
            static_cast<size_t>(chained_index) >= chained.anchors.size() ||
            chained.anchors[static_cast<size_t>(chained_index)].id < 0) {
          ok = false;
          break;
        }
        const auto found = raw_by_id.find(static_cast<uint32_t>(
            chained.anchors[static_cast<size_t>(chained_index)].id));
        if (found == raw_by_id.end()) {
          ok = false;
          break;
        }
        path.selected_raw_indices.push_back(found->second);
      }
      int count = 0;
      if (!ok ||
          !mapping::checked_size_to_int(path.selected_raw_indices.size(), count))
        continue;
      path.anchor_path_anchor_count = static_cast<uint32_t>(count);
      path.anchor_path_hash = hash_exact_anchor_path(path);
      if (!validate_exact_anchor_path(path, nullptr))
        continue;
      // Its own competing sibling, by the rule the selected chain's uses.
      piece.siblings.chains_total = static_cast<int32_t>(chained.chains.size());
      for (size_t j = 0; j < chained.chains.size(); ++j) {
        if (j == i || chained.chains[j].score <= piece.siblings.sib_score ||
            !rna_query_spans_compete(spans[i], spans[j]))
          continue;
        piece.siblings.sib_score = chained.chains[j].score;
        piece.siblings.sib_count =
            static_cast<int32_t>(chained.chains[j].idx.size());
      }
      disjoint->push_back(std::move(piece));
    }
  }
  return result;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
