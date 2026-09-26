// Anchor type shared by the DNA and RNA chainers.
#pragma once

#include <cstdint>

namespace fa {
namespace cpu {
namespace chaining {

// As minimap2's MM_SEED_IGNORE, MM_SEED_TANDEM and MM_SEED_LONG_JOIN.
enum AnchorFlag : uint32_t {
  ANCHOR_IGNORE = 1u << 0,
  ANCHOR_TANDEM = 1u << 1,
  ANCHOR_LONG_JOIN = 1u << 2,
};

// `r` and `q` are reference/query coordinates. The recurrence uses pairwise
// differences, so uniform-span START and END coordinates are equivalent.
struct Anchor {
  int32_t r = 0;
  int32_t q = 0;
  int32_t span = 0;
  int32_t id = -1;
  uint32_t flags = 0;

  int32_t q_end() const { return q + span; }
  int32_t r_end() const { return r + span; }
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
