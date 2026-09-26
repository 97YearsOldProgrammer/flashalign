// Aggregates one diagonal bin's captured postings into a CoarseDiagonalPeak, templated
// over how a hit is fetched.
#pragma once

#include "../../seeding/types.h" // VoteHit
#include "aggregate.h" // union_seed_coverage
#include "types.h"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// Returns false when the bin is empty, carries a negative reference id, or
// fails the distinct-seed min_support gate. Scratch vectors are caller-owned
// and reused; their prior contents are discarded.
template <class HitAt>
inline bool
build_bin_peak(std::int64_t key, std::size_t hit_count, HitAt&& hit_at, int k,
               bool reverse, int min_support, std::vector<int>& read_pos,
               std::vector<int>& ref_starts, std::vector<int>& seed_ids,
               CoarseDiagonalPeak& out) {
  if (hit_count == 0)
    return false;
  const int chr = static_cast<int>(key >> 32);
  if (chr < 0)
    return false;

  read_pos.clear();
  ref_starts.clear();
  seed_ids.clear();
  seed_ids.reserve(hit_count);
  std::int64_t ref_lo = INT64_MAX;
  std::int64_t ref_hi = INT64_MIN;
  int q_lo = INT32_MAX;
  int q_hi = INT32_MIN;
  for (std::size_t i = 0; i < hit_count; ++i) {
    const VoteHit& h = hit_at(i);
    read_pos.push_back(h.read_pos);
    ref_starts.push_back(h.ref_start);
    const std::int64_t genomic = static_cast<std::int64_t>(h.ref_start) +
                                 static_cast<std::int64_t>(h.read_pos);
    ref_lo = std::min(ref_lo, genomic);
    ref_hi = std::max(ref_hi, genomic + k);
    q_lo = std::min(q_lo, h.read_pos);
    q_hi = std::max(q_hi, h.read_pos + k);
    seed_ids.push_back(h.seed_idx);
  }
  std::sort(seed_ids.begin(), seed_ids.end());
  int distinct = 0;
  int last_seed = INT32_MIN;
  for (int s : seed_ids) {
    if (s != last_seed) {
      ++distinct;
      last_seed = s;
    }
  }
  if (distinct < min_support)
    return false;
  if (ref_lo < 0)
    ref_lo = 0;

  std::sort(read_pos.begin(), read_pos.end());
  std::nth_element(ref_starts.begin(),
                   ref_starts.begin() +
                       static_cast<long>(ref_starts.size() / 2),
                   ref_starts.end());
  const std::int64_t intercept = ref_starts[ref_starts.size() / 2];

  out.reference_id = chr;
  out.reverse = reverse;
  out.oriented_query_begin = static_cast<std::uint32_t>(std::max(0, q_lo));
  out.oriented_query_end = static_cast<std::uint32_t>(std::max(q_lo, q_hi));
  out.reference_begin = static_cast<std::uint64_t>(ref_lo);
  out.reference_end = static_cast<std::uint64_t>(std::max(ref_lo, ref_hi));
  out.median_diagonal = intercept;
  out.query_covered_bases = union_seed_coverage(read_pos, k);
  out.distinct_seed_support = static_cast<std::uint32_t>(distinct);
  return true;
}

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
