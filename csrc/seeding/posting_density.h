#pragma once

#include <cstdint>

namespace fa::cpu {

// Validate the chromosome directory used to interpret global posting
// coordinates. Empty chromosomes are permitted; the global coordinate space
// must start at zero and remain nondecreasing.
bool posting_density_boundaries_valid(
    const std::uint64_t* chromosome_offsets,
    std::uint64_t chromosome_count) noexcept;

}  // namespace fa::cpu
