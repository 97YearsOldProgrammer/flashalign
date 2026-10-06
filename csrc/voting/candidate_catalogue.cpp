#include "candidate_catalogue.h"
#include "query_partition.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <tuple>

namespace fa::cpu::voting {

CandidateCatalogue build_candidate_catalogue(
    std::vector<CandidateInput> inputs, int max_per_lane,
    double admission_ratio, CandidateMaskSource* mask_source, int tile_count,
    bool unpartitioned) {
  const int lane_bound_ceiling =
      unpartitioned ? kAllChainsLaneBound : kMaxCatalogueLaneBound;
  if (max_per_lane < 1 || max_per_lane > lane_bound_ceiling) {
    throw std::invalid_argument(
        "candidate catalogue per-lane bound must be within [1," +
        std::to_string(lane_bound_ceiling) + "]");
  }
  if (admission_ratio < 0.0 || admission_ratio > 1.0)
    throw std::invalid_argument(
        "candidate catalogue admission ratio must be within [0,1]");
  if (mask_source) {
    for (const CandidateInput& input : inputs) {
      if (input.mask_slot < 0)
        throw std::invalid_argument(
            "candidate catalogue mask source requires a mask slot per input");
    }
  }

  // Catalogue order is the caller's ranking; inputs stay in place so admission decisions
  // and mask slots keep their identity.
  std::vector<std::size_t> order(inputs.size());
  for (std::size_t index = 0; index < order.size(); ++index)
    order[index] = index;
  std::stable_sort(order.begin(), order.end(),
                   [&inputs](std::size_t left, std::size_t right) {
    return std::tuple(inputs[left].catalogue_rank,
                      -inputs[left].vote_evidence, inputs[left].lane,
                      inputs[left].equivalence_key) <
           std::tuple(inputs[right].catalogue_rank,
                      -inputs[right].vote_evidence, inputs[right].lane,
                      inputs[right].equivalence_key);
  });

  const auto mask_of = [&](std::size_t index) -> const QueryTileMask& {
    return mask_source ? mask_source->mask(inputs[index].mask_slot)
                       : inputs[index].support;
  };

  // Ratio admission judges a candidate against the best vote among the candidates that
  // compete for its query interval (as minimap2's mask_level: their tile masks overlap by at
  // least half of the smaller one), not against the read's best. A chimeric segment
  // competes with nothing and is admitted on its own vote; a paralog of the primary is held
  // to the ratio. A candidate with no tiles competes with everything. Only strictly
  // stronger candidates can drop one, so the decisions are order-independent.
  std::vector<char> admitted(inputs.size(), 1);
  if (admission_ratio > 0.0 && !inputs.empty()) {
    int vote_max = 0;
    for (const CandidateInput& input : inputs)
      vote_max = std::max(vote_max, input.vote_evidence);
    const auto ratio_floor = [admission_ratio](int vote) {
      return static_cast<int>(
          std::ceil(admission_ratio * static_cast<double>(vote)));
    };
    const auto coarse_valid = [&inputs](std::size_t index) {
      return inputs[index].coarse_tile_lo <= inputs[index].coarse_tile_hi;
    };
    const auto coarse_range = [&inputs, tile_count](std::size_t index) {
      QueryTileMask range;
      const int lo = std::max(0, inputs[index].coarse_tile_lo);
      const int hi = std::min(tile_count - 1, inputs[index].coarse_tile_hi);
      for (int tile = lo; tile <= hi; ++tile) range.set(tile);
      return range;
    };
    // Judge in descending vote order, materializing the unconditional
    // admissions' masks first, so every mask that can decide a weaker
    // candidate is available to the containment shortcut when it fires.
    std::vector<std::size_t> by_vote(order);
    std::stable_sort(by_vote.begin(), by_vote.end(),
                     [&inputs](std::size_t left, std::size_t right) {
      return inputs[left].vote_evidence > inputs[right].vote_evidence;
    });
    std::vector<char> materialized(inputs.size(), 0);
    const auto materialize = [&](std::size_t index) -> const QueryTileMask& {
      materialized[index] = 1;
      return mask_of(index);
    };
    const int strong_floor = ratio_floor(vote_max);
    for (const std::size_t index : by_vote) {
      const CandidateInput& input = inputs[index];
      if (input.vote_evidence >= strong_floor) {
        // No competitor can raise the floor above the global best: admitted
        // with no competition test. The mask is needed downstream.
        materialize(index);
        continue;
      }
      bool dropped = false;
      // Containment shortcut: a coarse range inside a materialized stronger mask proves the
      // competition without this candidate's mask, which then overlaps it fully or is empty.
      if (coarse_valid(index)) {
        const QueryTileMask range = coarse_range(index);
        for (const std::size_t other_index : by_vote) {
          if (other_index == index) continue;
          if (input.vote_evidence >=
              ratio_floor(inputs[other_index].vote_evidence))
            break;
          if (!materialized[other_index]) continue;
          const QueryTileMask& other = mask_of(other_index);
          if (other.count() == 0) continue;
          if ((range & other) == range) {
            dropped = true;
            break;
          }
        }
      }
      if (!dropped) {
        const QueryTileMask& own = materialize(index);
        const int own_tiles = own.count();
        if (own_tiles == 0) {
          dropped = true;  // competes with everything, and vote_max is above
        } else {
          for (const std::size_t other_index : by_vote) {
            if (other_index == index) continue;
            if (input.vote_evidence >=
                ratio_floor(inputs[other_index].vote_evidence))
              break;
            // Disjoint coarse ranges bound disjoint masks: not competing.
            if (coarse_valid(index) && coarse_valid(other_index) &&
                (inputs[index].coarse_tile_hi <
                     inputs[other_index].coarse_tile_lo ||
                 inputs[other_index].coarse_tile_hi <
                     inputs[index].coarse_tile_lo))
              continue;
            const QueryTileMask& other = materialize(other_index);
            const int other_tiles = other.count();
            if (other_tiles == 0) continue;
            const int shared = (own & other).count();
            if (2 * shared >= std::min(own_tiles, other_tiles)) {
              dropped = true;
              break;
            }
          }
        }
      }
      admitted[index] = dropped ? 0 : 1;
    }
  }

  CandidateCatalogue result;
  std::set<std::uint64_t> seen;
  std::map<int, int> lane_counts;
  for (const std::size_t index : order) {
    const CandidateInput& input = inputs[index];
    if (!seen.insert(input.equivalence_key).second) continue;
    if (lane_counts[input.lane] >= max_per_lane) continue;
    if (!admitted[index]) continue;
    ++lane_counts[input.lane];
    QueryCandidate candidate;
    candidate.id = static_cast<CandidateId>(result.candidates.size());
    candidate.equivalence_key = input.equivalence_key;
    candidate.lane = input.lane;
    candidate.catalogue_rank = input.catalogue_rank;
    candidate.vote_evidence = input.vote_evidence;
    candidate.chain_evidence = input.chain_evidence;
    candidate.support = mask_of(index);
    result.candidates.push_back(candidate);
  }
  return result;
}

}  // namespace fa::cpu::voting
