#pragma once

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

// Reverse complement of a 0-4 encoded sequence (A C G T N); N stays N.
std::vector<uint8_t>
reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc);

} // namespace lr
} // namespace cpu
} // namespace fa
