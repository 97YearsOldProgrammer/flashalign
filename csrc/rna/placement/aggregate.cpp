#include "aggregate.h"

#include <algorithm>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

std::uint32_t union_seed_coverage(const std::vector<int>& sorted_read_pos,
                                  int k) {
  std::uint32_t covered = 0;
  int cur_lo = 0;
  int cur_hi = 0;
  bool open = false;
  for (int p : sorted_read_pos) {
    const int lo = p;
    const int hi = p + k;
    if (!open) {
      cur_lo = lo;
      cur_hi = hi;
      open = true;
    } else if (lo <= cur_hi) {
      cur_hi = std::max(cur_hi, hi);
    } else {
      covered += static_cast<std::uint32_t>(cur_hi - cur_lo);
      cur_lo = lo;
      cur_hi = hi;
    }
  }
  if (open)
    covered += static_cast<std::uint32_t>(cur_hi - cur_lo);
  return covered;
}

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
