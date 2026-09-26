// Occurrence distribution of a built index, and minimap2's mid-occurrence quantile over
// it. Every key's exact occurrence is stored: an inline singleton is 1, otherwise the
// payload's high word is the posting count.
#pragma once

#include "faix.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

namespace fa {
namespace cpu {

// Occurrences at or below this bound get a dense bucket; the rare tail above it is kept
// exactly in a std::map.
inline constexpr std::uint32_t kFaixOccHistogramBound = 1u << 16;

struct FaixOccDistribution {
  // counts[v - 1] == number of distinct keys whose occurrence is exactly v,
  // for 1 <= v <= kFaixOccHistogramBound.
  std::vector<std::uint64_t> counts;
  // Exact key counts for occurrence values above kFaixOccHistogramBound.
  std::map<std::uint32_t, std::uint64_t> tail;
  std::uint64_t n_distinct = 0;
  // Every key's occurrences summed: all the reference seeds.
  std::uint64_t total_postings = 0;
};

// One streaming pass over the compact record table (the shards laid end to end).
inline FaixOccDistribution faix_occ_distribution(const FaixIndex& index) {
  FaixOccDistribution dist;
  dist.counts.assign(kFaixOccHistogramBound, 0);
  const std::uint32_t* tags = index.tags_data();
  const std::uint64_t* payloads = index.payloads_data();
  if (index.empty() || tags == nullptr || payloads == nullptr)
    return dist;
  const std::uint64_t records = index.distinct_keys();
  for (std::uint64_t r = 0; r < records; ++r) {
    const std::uint32_t occ =
        (tags[r] & kFaixTagMultiBit) == 0
            ? 1u
            : static_cast<std::uint32_t>(payloads[r] >> 32);
    if (occ == 0)
      continue; // structurally impossible; never counted
    ++dist.n_distinct;
    dist.total_postings += occ;
    if (occ <= kFaixOccHistogramBound) {
      ++dist.counts[occ - 1];
    } else {
      ++dist.tail[occ];
    }
  }
  return dist;
}

// minimap2's mm_idx_cal_max_occ(mi, f), evaluated from the cumulative counts. With the n
// distinct occurrences sorted ascending, minimap2 returns occ[(size_t)((1 - f) * n)] + 1:
// the smallest v whose cumulative count exceeds that 0-based rank. f <= 0 disables the
// quantile (INT32_MAX), as in minimap2.
inline std::int64_t faix_cal_max_occ(const FaixOccDistribution& dist,
                                     double f) {
  if (!(f > 0.0))
    return static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max());
  if (dist.n_distinct == 0)
    return 0;
  const double position = (1.0 - f) * static_cast<double>(dist.n_distinct);
  std::uint64_t rank =
      position <= 0.0 ? 0 : static_cast<std::uint64_t>(position);
  if (rank >= dist.n_distinct)
    rank = dist.n_distinct - 1;
  std::uint64_t cumulative = 0;
  for (std::uint32_t v = 1; v <= kFaixOccHistogramBound; ++v) {
    cumulative += dist.counts[v - 1];
    if (cumulative > rank)
      return static_cast<std::int64_t>(v) + 1;
  }
  for (const auto& entry : dist.tail) {
    cumulative += entry.second;
    if (cumulative > rank)
      return static_cast<std::int64_t>(entry.first) + 1;
  }
  return 0; // unreachable: the ranks are bounded by n_distinct
}

} // namespace cpu
} // namespace fa
