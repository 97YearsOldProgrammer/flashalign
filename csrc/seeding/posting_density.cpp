#include "posting_density.h"

namespace fa::cpu {

bool posting_density_boundaries_valid(
    const std::uint64_t* offsets,
    std::uint64_t chromosome_count) noexcept {
  if (offsets == nullptr || chromosome_count == 0 || offsets[0] != 0)
    return false;
  for (std::uint64_t chromosome = 0; chromosome < chromosome_count;
       ++chromosome) {
    if (offsets[chromosome + 1] < offsets[chromosome])
      return false;
  }
  return true;
}

}  // namespace fa::cpu
