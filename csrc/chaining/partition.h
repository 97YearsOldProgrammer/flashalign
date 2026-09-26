// Splits chains into a primary and non-rival supplementaries, and counts the primary's
// rivals (n_sub, f2) for MAPQ.
#pragma once

#include "result.h"
#include "../split/query_geometry.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace chaining {

inline constexpr int kDefaultMaximumSegments = 8;

struct ChainPartition {
  int primary = -1;
  std::vector<int> supplementary;
  int32_t f2 = 0;
  // Index in ChainResult::chains of the best-scoring rival (the chain f2 comes from);
  // -1 when no chain rivals the primary.
  int f2_index = -1;
  int n_sub = 0;
};

ChainPartition partition_chains(
    const ChainResult& result,
    ::fa::cpu::split::QueryMaskPolicy policy =
        ::fa::cpu::split::kHalfFloor,
    int max_segments = kDefaultMaximumSegments);

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
