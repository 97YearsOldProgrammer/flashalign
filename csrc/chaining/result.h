#pragma once

#include "anchor.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace chaining {

struct Chain {
  std::vector<int32_t> idx;
  int32_t score = 0;
};

struct ChainResult {
  std::vector<Anchor> anchors;
  std::vector<Chain> chains;
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
