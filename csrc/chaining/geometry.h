#pragma once

#include "result.h"

#include <cstddef>
#include <cstdint>

namespace fa {
namespace cpu {
namespace chaining {

inline void chain_query_span(const ChainResult& result, const Chain& chain,
                             int32_t& query_low, int32_t& query_high) {
  query_low = INT32_MAX;
  query_high = INT32_MIN;
  for (int32_t index : chain.idx) {
    const Anchor& anchor = result.anchors[static_cast<size_t>(index)];
    if (anchor.q < query_low) query_low = anchor.q;
    if (anchor.q + anchor.span > query_high)
      query_high = anchor.q + anchor.span;
  }
}

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
