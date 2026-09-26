#include "partition.h"

#include "colinear_chain.h"
#include "geometry.h"

#include <algorithm>

namespace fa {
namespace cpu {
namespace chaining {

ChainPartition partition_chains(
    const ChainResult& result, ::fa::cpu::split::QueryMaskPolicy policy,
    int max_segments) {
  using ::fa::cpu::split::classify_query_relation;
  using ::fa::cpu::split::QueryInterval;
  using ::fa::cpu::split::QueryRelation;

  ChainPartition partition;
  const int best = best_chain_index(result);
  if (best < 0) return partition;
  partition.primary = best;

  auto span = [&](int chain_index) {
    int32_t low = 0;
    int32_t high = 0;
    chain_query_span(result,
                     result.chains[static_cast<size_t>(chain_index)],
                     low, high);
    return QueryInterval{low, high};
  };
  const QueryInterval primary_span = span(best);

  std::vector<int> order;
  order.reserve(result.chains.size());
  for (size_t i = 0; i < result.chains.size(); ++i) {
    if (static_cast<int>(i) != best) order.push_back(static_cast<int>(i));
  }
  std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
    return result.chains[static_cast<size_t>(left)].score >
           result.chains[static_cast<size_t>(right)].score;
  });

  std::vector<QueryInterval> selected;
  selected.push_back(primary_span);
  for (int chain_index : order) {
    const QueryInterval candidate_span = span(chain_index);
    if (candidate_span.length() <= 0) continue;

    if (classify_query_relation(candidate_span, primary_span, policy) ==
        QueryRelation::Rival) {
      ++partition.n_sub;
      if (result.chains[static_cast<size_t>(chain_index)].score >
          partition.f2) {
        partition.f2 =
            result.chains[static_cast<size_t>(chain_index)].score;
        partition.f2_index = chain_index;
      }
    }

    if (static_cast<int>(selected.size()) >= max_segments) continue;
    bool is_segment = true;
    for (const QueryInterval& already_selected : selected) {
      if (classify_query_relation(candidate_span, already_selected, policy) ==
          QueryRelation::Rival) {
        is_segment = false;
        break;
      }
    }
    if (is_segment) {
      partition.supplementary.push_back(chain_index);
      selected.push_back(candidate_span);
    }
  }
  return partition;
}

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
