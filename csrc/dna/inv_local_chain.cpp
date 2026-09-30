#include "inv_local_chain.h"

#include "retained_seed_density.h"
#include "../index/format.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <utility>
#include <vector>

namespace fa::cpu::lr {

int dna_inv_local_chain(const void* opaque, int query_begin, int query_end,
                        int target_begin, int target_end) {
  const DnaInvLocalChainSource& source =
      *static_cast<const DnaInvLocalChainSource*>(opaque);
  const int k = source.seed_length;
  // The window in the opposite lane's frame, the reverse complement of the
  // block's.
  const int lane_begin =
      source.read_length - (source.query_offset + query_end);
  const int lane_end =
      source.read_length - (source.query_offset + query_begin);
  const int reference_begin = source.target_offset + target_begin;
  const int reference_end = source.target_offset + target_end;
  if (lane_end - lane_begin < k || reference_end - reference_begin < k)
    return 0;
  const bool lane_reverse = !source.block_reverse;
  const std::vector<RetainedSeedRef>& lane =
      lane_reverse ? source.seeds->fine_reverse()
                   : source.seeds->fine_forward();
  // Both lanes are in ascending read_pos.
  auto seed = std::lower_bound(
      lane.begin(), lane.end(), lane_begin,
      [](const RetainedSeedRef& ref, int position) {
        return ref.seed.read_pos < position;
      });
  std::vector<std::pair<int, int>> anchors;  // (q', r)
  for (; seed != lane.end() && seed->seed.read_pos <= lane_end - k; ++seed) {
    // Postings starting in [reference_begin, reference_end - k].
    const KmerPostingIntervalView interval = source.seeds->slice(
        *source.index, seed->entry, source.chromosome,
        static_cast<std::uint32_t>(reference_begin),
        static_cast<std::uint32_t>(reference_end - k + 1));
    if (!interval.found() || interval.count == 0 ||
        (source.occurrence_cap != 0 &&
         interval.global_count > source.occurrence_cap))
      continue;
    for (std::uint32_t posting = 0; posting < interval.count; ++posting) {
      const PackedRefPos packed = interval.positions.packed_at(posting);
      if (!packed_ref_orientation_compatible(seed->seed.z, packed,
                                             lane_reverse))
        continue;
      anchors.emplace_back(seed->seed.read_pos,
                           static_cast<int>(packed_ref_local(packed)));
    }
  }
  if (anchors.empty()) return 0;
  std::sort(anchors.begin(), anchors.end());
  // Longest chain, exact over the predecessors within kDnaInvLocalGap in q'.
  std::vector<int> length(anchors.size(), 1);
  int best = 1;
  std::size_t low = 0;
  for (std::size_t j = 0; j < anchors.size(); ++j) {
    const int qj = anchors[j].first;
    const int rj = anchors[j].second;
    while (anchors[low].first < qj - kDnaInvLocalGap) ++low;
    const int dj = rj - qj;
    for (std::size_t i = low; i < j; ++i) {
      const int qi = anchors[i].first;
      const int ri = anchors[i].second;
      if (qi >= qj || ri >= rj || rj - ri > kDnaInvLocalGap ||
          std::abs((ri - qi) - dj) > kDnaInvLocalDiagonal)
        continue;
      length[j] = std::max(length[j], length[i] + 1);
    }
    best = std::max(best, length[j]);
  }
  return best;
}

}  // namespace fa::cpu::lr
