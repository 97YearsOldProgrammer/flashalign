// The selected RNA anchor path as realization consumes it.
#pragma once

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum SpliceAnchorFlag : std::uint8_t {
  kSpliceAnchorLongJoin = 1u << 0,
  kSpliceAnchorIgnore = 1u << 1,
  kSpliceAnchorTandem = 1u << 2,
};

struct SpliceAnchor {
  std::uint32_t stable_id = 0;
  std::int32_t query_begin = 0;
  std::int32_t reference_begin = 0;
  std::int32_t span = 0;
  std::uint8_t flags = 0;
};

struct SpliceAnchorView {
  std::vector<SpliceAnchor> anchors;
  int selected_begin = 0;
  int selected_count = 0;
  bool mapping_reverse = false;
  std::vector<std::uint32_t> view_to_raw_index;
  std::vector<std::uint32_t> excluded_interior_raw_indices;
  std::uint64_t content_hash = 0;
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
