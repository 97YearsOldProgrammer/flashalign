#include "sequence.h"

namespace fa {
namespace cpu {
namespace lr {

std::vector<uint8_t>
reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc) {
  std::vector<uint8_t> rc_enc(fwd_enc.size());
  for (size_t i = 0; i < fwd_enc.size(); i++) {
    const uint8_t b = fwd_enc[fwd_enc.size() - 1 - i];
    rc_enc[i] = (b > 3) ? 4 : static_cast<uint8_t>(3 - b);
  }
  return rc_enc;
}

} // namespace lr
} // namespace cpu
} // namespace fa
