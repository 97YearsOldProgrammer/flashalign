#include "alternative_hypothesis.h"

#include "context.h"
#include "placement_chaining.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

constexpr int kShadowFloorBp = 2000;
constexpr int kShadowFraction = 10;

// minimap2's mm_select_sub: a rival is credible at pri_ratio times the
// owner's score or within credible_gap of it. Exact for int scores at
// p = 0.8, where it is 5 * score >= 4 * owner.
bool credible_rival(int score, int owner, double pri_ratio, int credible_gap) {
  return static_cast<double>(score) >=
             pri_ratio * static_cast<double>(owner) ||
         owner - score <= credible_gap;
}

int query_overlap(int left_begin, int left_end, int right_begin,
                  int right_end) {
  return std::max(0, std::min(left_end, right_end) -
                         std::max(left_begin, right_begin));
}

bool same_locus_shadow(const DnaPlacementCandidate& candidate,
                       const DnaPlacementCandidate& selected,
                       int read_length) {
  if (candidate.peak.chr != selected.peak.chr ||
      candidate.peak.is_rc != selected.peak.is_rc)
    return false;
  const std::int64_t window = std::max<std::int64_t>(
      kShadowFloorBp, read_length / kShadowFraction);
  return std::llabs(candidate.peak.raw_ref_start -
                    selected.peak.raw_ref_start) <= window;
}

// A qualifying candidate's sort key: (-score, -anchors, -overlap, id).
using RankKey = std::tuple<int, int, int, ::fa::cpu::voting::CandidateId>;

// select_dna_alternative_hypothesis, which also appends every qualifying
// candidate's key to `qualified` when given.
DnaAlternativeSelection select_alternative(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable,
    std::vector<RankKey>* qualified) {
  DnaAlternativeSelection out;
  if (!stable.accepted) return out;

  int non_null_blocks = 0;
  int block_begin = 0;
  int block_end = 0;
  for (const auto& block : family.partition.selected.blocks) {
    if (block.candidate == ::fa::cpu::voting::kNullCandidate) continue;
    ++non_null_blocks;
    out.incumbent = block.candidate;
    block_begin = ::fa::cpu::voting::query_tile_begin(
        block.query_tile_begin, family.read_length, family.seed_length,
        family.tile_count);
    block_end = ::fa::cpu::voting::query_tile_end(
        block.query_tile_end, family.read_length, family.seed_length,
        family.tile_count);
  }
  if (non_null_blocks != 1) {
    out.refusal = DnaAlternativeRefusal::NotSingleBlock;
    return out;
  }
  const DnaPlacementCandidateChain* incumbent = stable.find(out.incumbent);
  if (incumbent == nullptr) return out;

  std::vector<::fa::cpu::voting::CandidateId> selected;
  for (const auto id : family.partition.selected.assignment) {
    if (id != ::fa::cpu::voting::kNullCandidate &&
        std::find(selected.begin(), selected.end(), id) == selected.end())
      selected.push_back(id);
  }
  const int contig_count = context.ref.names == nullptr
      ? 0 : static_cast<int>(context.ref.names->size());
  const int credible_gap = 2 * std::max(1, context.opts.k);
  DnaAlternativeRefusal worst = DnaAlternativeRefusal::NoCandidate;
  const auto note = [&](DnaAlternativeRefusal refusal) {
    if (static_cast<int>(refusal) > static_cast<int>(worst)) worst = refusal;
  };

  const DnaPlacementCandidate* best = nullptr;
  const DnaPlacementCandidateChain* best_record = nullptr;
  int best_overlap = 0;
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    if (std::find(selected.begin(), selected.end(), candidate.id) !=
        selected.end())
      continue;
    if (candidate.peak.chr < 0 || candidate.peak.chr >= contig_count) {
      note(DnaAlternativeRefusal::NoLocus);
      continue;
    }
    const DnaPlacementCandidateChain* record = stable.find(candidate.id);
    if (record == nullptr || record->screening_chain_score <= 0 ||
        record->screening_forward_query_begin < 0 ||
        record->screening_forward_query_end <=
            record->screening_forward_query_begin) {
      note(DnaAlternativeRefusal::NoChain);
      continue;
    }
    bool shadow = false;
    for (const auto id : selected) {
      const DnaPlacementCandidate* chosen = family.find(id);
      if (chosen != nullptr &&
          same_locus_shadow(candidate, *chosen, family.read_length)) {
        shadow = true;
        break;
      }
    }
    if (shadow) {
      note(DnaAlternativeRefusal::ShadowOnly);
      continue;
    }
    const int span = record->screening_forward_query_end -
                     record->screening_forward_query_begin;
    const int overlap = query_overlap(
        record->screening_forward_query_begin,
        record->screening_forward_query_end, block_begin, block_end);
    const int shorter = std::min(span, std::max(0, block_end - block_begin));
    if (static_cast<std::int64_t>(overlap) * 2 < shorter) {
      note(DnaAlternativeRefusal::DisjointOnly);
      continue;
    }
    const int score = record->screening_chain_score;
    const int incumbent_score = incumbent->screening_chain_score;
    if (!credible_rival(score, incumbent_score, context.opts.pri_ratio,
                        credible_gap)) {
      note(DnaAlternativeRefusal::IncredibleOnly);
      continue;
    }
    if (qualified != nullptr)
      qualified->emplace_back(-score, -record->screening_chain_anchors,
                              -overlap, candidate.id);
    if (best == nullptr ||
        std::tuple(-score, -record->screening_chain_anchors, -overlap,
                   candidate.id) <
            std::tuple(-best_record->screening_chain_score,
                       -best_record->screening_chain_anchors, -best_overlap,
                       best->id)) {
      best = &candidate;
      best_record = record;
      best_overlap = overlap;
    }
  }
  if (best == nullptr) {
    out.refusal = worst;
    return out;
  }
  out.candidate = best->id;
  out.refusal = DnaAlternativeRefusal::None;
  return out;
}

}  // namespace

DnaAlternativeSelection select_dna_alternative_hypothesis(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable) {
  return select_alternative(context, family, stable, nullptr);
}

std::vector<::fa::cpu::voting::CandidateId> rank_dna_alternative_hypotheses(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable, std::size_t limit) {
  std::vector<RankKey> qualified;
  select_alternative(context, family, stable, &qualified);
  std::sort(qualified.begin(), qualified.end());
  if (qualified.size() > limit) qualified.resize(limit);
  std::vector<::fa::cpu::voting::CandidateId> ranked;
  ranked.reserve(qualified.size());
  for (const RankKey& key : qualified) ranked.push_back(std::get<3>(key));
  return ranked;
}

}  // namespace fa::cpu::lr
