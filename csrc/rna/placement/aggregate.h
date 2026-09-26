#pragma once

#include "types.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// Read bases covered by the union of [pos, pos + k) over ascending seed positions.
std::uint32_t union_seed_coverage(
    const std::vector<int>& sorted_read_pos, int k);

}  // namespace placement
}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
