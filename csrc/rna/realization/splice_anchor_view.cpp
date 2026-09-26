#include "splice_anchor_view.h"

#include "../../core/checked_range.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

constexpr uint64_t kFnvOffset = 1469598103934665603ull;
constexpr uint64_t kFnvPrime = 1099511628211ull;

uint64_t mix_u64(uint64_t hash, uint64_t value) {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (byte * 8)) & 0xffu;
    hash *= kFnvPrime;
  }
  return hash;
}

uint64_t mix_vector(uint64_t hash, const std::vector<uint32_t>& values) {
  hash = mix_u64(hash, values.size());
  for (uint32_t value : values) hash = mix_u64(hash, value);
  return hash;
}

uint64_t hash_view(const SpliceAnchorView& view) {
  uint64_t hash = kFnvOffset;
  hash = mix_u64(hash, static_cast<uint32_t>(view.selected_begin));
  hash = mix_u64(hash, static_cast<uint32_t>(view.selected_count));
  for (size_t i = 0; i < view.anchors.size(); ++i) {
    const SpliceAnchor& anchor = view.anchors[i];
    const int64_t query_end = static_cast<int64_t>(anchor.query_begin) +
                              anchor.span - 1;
    const int64_t reference_end =
        static_cast<int64_t>(anchor.reference_begin) + anchor.span - 1;
    const uint64_t conceptual_x =
        (view.mapping_reverse ? (1ull << 63) : 0ull) |
        static_cast<uint32_t>(reference_end);
    uint64_t conceptual_y = static_cast<uint64_t>(anchor.span) << 32 |
                            static_cast<uint32_t>(query_end);
    if ((anchor.flags & kSpliceAnchorLongJoin) != 0)
      conceptual_y |= 1ull << 40;
    if ((anchor.flags & kSpliceAnchorIgnore) != 0)
      conceptual_y |= 1ull << 41;
    if ((anchor.flags & kSpliceAnchorTandem) != 0)
      conceptual_y |= 1ull << 42;
    hash = mix_u64(hash, conceptual_x);
    hash = mix_u64(hash, conceptual_y);
    hash = mix_u64(hash, anchor.stable_id);
    hash = mix_u64(hash, view.view_to_raw_index[i]);
  }
  return mix_vector(hash, view.excluded_interior_raw_indices);
}

}  // namespace

SpliceAnchorView make_splice_anchor_view(
    const ExactAnchorPath& bundle, int index_k, std::string* error) {
  SpliceAnchorView view;
  AnchoringRefusal validation_refusal = AnchoringRefusal::None;
  if (!validate_exact_anchor_path(bundle, &validation_refusal)) {
    if (error) *error = anchoring_refusal_name(validation_refusal);
    return view;
  }
  if (index_k <= 0 || index_k > 255) {
    if (error) *error = "SPLICE_VIEW_INDEX_K";
    return view;
  }
  int raw_count = 0;
  int selected_count = 0;
  if (!mapping::checked_size_to_int(bundle.selected_pool_anchors.size(), raw_count) ||
      !mapping::checked_size_to_int(
          bundle.selected_raw_indices.size(), selected_count)) {
    if (error) *error = "SPLICE_VIEW_COUNT_DOMAIN";
    return view;
  }
  view.anchors.reserve(static_cast<std::size_t>(raw_count));
  view.view_to_raw_index.reserve(static_cast<std::size_t>(raw_count));

  const uint32_t first = bundle.selected_raw_indices.front();
  const uint32_t last = bundle.selected_raw_indices.back();
  std::vector<uint8_t> selected(bundle.selected_pool_anchors.size(), 0);
  for (uint32_t index : bundle.selected_raw_indices) selected[index] = 1;

  view.mapping_reverse = bundle.mapping_reverse;
  auto append = [&](uint32_t raw_index) {
    const SelectedPoolAnchor& raw = bundle.selected_pool_anchors[raw_index];
    if (raw.span != index_k) return false;
    const int64_t query_end = static_cast<int64_t>(raw.query_begin) +
                              raw.span - 1;
    const int64_t reference_end =
        static_cast<int64_t>(raw.reference_begin) + raw.span - 1;
    if (query_end < 0 || query_end > UINT32_MAX || reference_end < 0 ||
        reference_end > UINT32_MAX)
      return false;

    SpliceAnchor anchor;
    anchor.stable_id = raw.stable_id;
    anchor.query_begin = raw.query_begin;
    anchor.reference_begin = raw.reference_begin;
    anchor.span = raw.span;
    if (raw.long_join_flag) anchor.flags |= kSpliceAnchorLongJoin;
    if (raw.ignore_flag) anchor.flags |= kSpliceAnchorIgnore;
    if (raw.tandem_flag) anchor.flags |= kSpliceAnchorTandem;
    view.anchors.push_back(anchor);
    view.view_to_raw_index.push_back(raw_index);
    return true;
  };

  for (uint32_t index = 0; index < first; ++index) {
    if (!append(index)) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  if (!mapping::checked_size_to_int(view.anchors.size(), view.selected_begin)) {
    if (error) *error = "SPLICE_VIEW_COUNT_DOMAIN";
    return {};
  }
  for (uint32_t index : bundle.selected_raw_indices) {
    if (!append(index)) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  view.selected_count = selected_count;
  for (uint32_t index = first; index <= last; ++index) {
    if (!selected[index]) view.excluded_interior_raw_indices.push_back(index);
  }
  for (uint64_t index = static_cast<uint64_t>(last) + 1;
       index < bundle.selected_pool_anchors.size(); ++index) {
    if (!append(static_cast<uint32_t>(index))) {
      if (error) *error = "SPLICE_VIEW_ANCHOR_ENCODING";
      return {};
    }
  }
  view.content_hash = hash_view(view);
  if (error) error->clear();
  return view;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
