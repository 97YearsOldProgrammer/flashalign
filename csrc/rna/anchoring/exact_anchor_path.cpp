#include "exact_anchor_path.h"

#include "../../core/checked_range.h"

#include <algorithm>
#include <limits>
#include <utility>
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

}  // namespace

bool selected_pool_anchor_precedes(const SelectedPoolAnchor& left,
                         const SelectedPoolAnchor& right) noexcept {
  if (left.reference_begin != right.reference_begin)
    return left.reference_begin < right.reference_begin;
  if (left.query_begin != right.query_begin)
    return left.query_begin < right.query_begin;
  if (left.span != right.span) return left.span < right.span;
  return left.stable_id < right.stable_id;
}

uint64_t hash_rna_read_key(const std::string& read) noexcept {
  uint64_t hash = kFnvOffset;
  for (unsigned char value : read) {
    const unsigned char upper =
        value >= 'a' && value <= 'z'
            ? static_cast<unsigned char>(value - ('a' - 'A'))
            : value;
    hash ^= static_cast<uint64_t>(upper);
    hash *= kFnvPrime;
  }
  return hash;
}

uint64_t hash_exact_anchor_path(const ExactAnchorPath& bundle) {
  uint64_t hash = kFnvOffset;
  hash = mix_u64(hash, bundle.schema_version);
  hash = mix_u64(hash, bundle.read_key);
  hash = mix_u64(hash, static_cast<uint32_t>(bundle.query_length));
  hash = mix_u64(hash, static_cast<uint32_t>(bundle.reference_id));
  hash = mix_u64(hash, bundle.mapping_reverse ? 1u : 0u);
  hash = mix_u64(hash, bundle.candidate_pool_hash);
  hash = mix_u64(hash, static_cast<uint32_t>(bundle.anchor_path_score));
  hash = mix_u64(hash, bundle.anchor_path_anchor_count);
  hash = mix_u64(hash, bundle.selected_pool_anchors.size());
  for (const SelectedPoolAnchor& anchor : bundle.selected_pool_anchors) {
    hash = mix_u64(hash, anchor.stable_id);
    hash = mix_u64(hash, static_cast<uint32_t>(anchor.query_begin));
    hash = mix_u64(hash, static_cast<uint32_t>(anchor.reference_begin));
    hash = mix_u64(hash, static_cast<uint32_t>(anchor.span));
    hash = mix_u64(hash, anchor.source_seed_id);
    hash = mix_u64(hash, anchor.source_occurrence);
    hash = mix_u64(hash, anchor.ignore_flag ? 1u : 0u);
    hash = mix_u64(hash, anchor.tandem_flag ? 1u : 0u);
    hash = mix_u64(hash, anchor.long_join_flag ? 1u : 0u);
    hash = mix_u64(hash, anchor.region_kind_mask);
    hash = mix_vector(hash, anchor.region_ids);
    hash = mix_vector(hash, anchor.core_node_ids);
  }
  return mix_vector(hash, bundle.selected_raw_indices);
}

bool summarize_exact_anchor_path(const ExactAnchorPath& bundle,
                                ExactAnchorPathSummary& summary,
                                AnchoringRefusal* refusal) noexcept {
  summary = {};
  auto fail = [&](AnchoringRefusal reason) {
    if (refusal) *refusal = reason;
    return false;
  };
  int selected_count = 0;
  if (!mapping::checked_size_to_int(bundle.selected_raw_indices.size(),
                                    selected_count))
    return fail(AnchoringRefusal::BundleSelectedCountDomain);
  if (selected_count <= 0) return fail(AnchoringRefusal::BundleEmpty);
  if (bundle.anchor_path_anchor_count != bundle.selected_raw_indices.size())
    return fail(AnchoringRefusal::BundleCount);

  uint32_t previous = 0;
  bool first = true;
  int first_query_begin = 0;
  int last_query_end = 0;
  int first_reference_begin = 0;
  int last_reference_end = 0;
  for (uint32_t raw_index : bundle.selected_raw_indices) {
    if (raw_index >= bundle.selected_pool_anchors.size())
      return fail(AnchoringRefusal::BundleSelectedRange);
    if (!first && raw_index <= previous)
      return fail(AnchoringRefusal::BundleSelectedOrder);
    const SelectedPoolAnchor& anchor = bundle.selected_pool_anchors[raw_index];
    int query_end = 0;
    int reference_end = 0;
    if (anchor.span <= 0 ||
        !mapping::checked_signed_endpoint(
            anchor.query_begin, anchor.span, bundle.query_length, query_end) ||
        !mapping::checked_signed_endpoint(
            anchor.reference_begin, anchor.span,
            mapping::kSignedCoordinateMaximum, reference_end))
      return fail(AnchoringRefusal::BundleAnchorRange);
    if (first) {
      first_query_begin = anchor.query_begin;
      first_reference_begin = anchor.reference_begin;
    }
    last_query_end = query_end;
    last_reference_end = reference_end;
    previous = raw_index;
    first = false;
  }
  if (last_query_end < first_query_begin ||
      last_reference_end < first_reference_begin)
    return fail(AnchoringRefusal::BundleSelectedSpan);

  summary.selected_count = static_cast<uint64_t>(selected_count);
  summary.first_query_begin = first_query_begin;
  summary.last_query_end = last_query_end;
  summary.outer_query_span = static_cast<uint64_t>(
      static_cast<int64_t>(last_query_end) - first_query_begin);
  summary.first_reference_begin = first_reference_begin;
  summary.last_reference_end = last_reference_end;
  summary.outer_reference_span = static_cast<uint64_t>(
      static_cast<int64_t>(last_reference_end) - first_reference_begin);
  if (refusal) *refusal = AnchoringRefusal::None;
  return true;
}

bool anchor_union_query_coverage(const ExactAnchorPath& bundle,
                                 uint64_t& covered_bases) noexcept {
  covered_bases = 0;
  if (bundle.query_length <= 0 || bundle.selected_raw_indices.empty())
    return false;
  std::vector<std::pair<int32_t, int32_t>> intervals;
  try {
    intervals.reserve(bundle.selected_raw_indices.size());
    for (uint32_t raw_index : bundle.selected_raw_indices) {
      if (raw_index >= bundle.selected_pool_anchors.size()) return false;
      const SelectedPoolAnchor& anchor = bundle.selected_pool_anchors[raw_index];
      const int64_t end = static_cast<int64_t>(anchor.query_begin) +
                          static_cast<int64_t>(anchor.span);
      if (anchor.query_begin < 0 || anchor.span <= 0 ||
          end > bundle.query_length)
        return false;
      intervals.emplace_back(anchor.query_begin, static_cast<int32_t>(end));
    }
    std::sort(intervals.begin(), intervals.end());
  } catch (...) {
    return false;
  }
  uint64_t covered = 0;
  int32_t union_begin = intervals.front().first;
  int32_t union_end = intervals.front().second;
  for (size_t i = 1; i < intervals.size(); ++i) {
    if (intervals[i].first > union_end) {
      covered += static_cast<uint64_t>(union_end - union_begin);
      union_begin = intervals[i].first;
      union_end = intervals[i].second;
    } else {
      union_end = std::max(union_end, intervals[i].second);
    }
  }
  covered += static_cast<uint64_t>(union_end - union_begin);
  if (covered == 0 || covered > static_cast<uint64_t>(bundle.query_length))
    return false;
  covered_bases = covered;
  return true;
}

bool validate_exact_anchor_path(const ExactAnchorPath& bundle,
                                AnchoringRefusal* refusal) {
  auto fail = [&](AnchoringRefusal reason) {
    if (refusal) *refusal = reason;
    return false;
  };
  if (bundle.schema_version != ExactAnchorPath::kSchemaVersion)
    return fail(AnchoringRefusal::BundleSchema);
  if (bundle.query_length <= 0 || bundle.reference_id < 0)
    return fail(AnchoringRefusal::BundleIdentity);
  if (bundle.selected_pool_anchors.empty() || bundle.selected_raw_indices.empty())
    return fail(AnchoringRefusal::BundleEmpty);
  int raw_count = 0;
  int selected_count = 0;
  if (!mapping::checked_size_to_int(bundle.selected_pool_anchors.size(), raw_count) ||
      !mapping::checked_size_to_int(bundle.selected_raw_indices.size(),
                                    selected_count))
    return fail(AnchoringRefusal::BundleCountDomain);
  for (size_t i = 0; i < bundle.selected_pool_anchors.size(); ++i) {
    const SelectedPoolAnchor& anchor = bundle.selected_pool_anchors[i];
    int query_end = 0;
    int reference_end = 0;
    if (anchor.stable_id >
            static_cast<uint32_t>(mapping::kSignedCoordinateMaximum))
      return fail(AnchoringRefusal::BundleStableIdDomain);
    if (anchor.span <= 0 ||
        !mapping::checked_signed_endpoint(
            anchor.query_begin, anchor.span, bundle.query_length, query_end) ||
        !mapping::checked_signed_endpoint(
            anchor.reference_begin, anchor.span,
            mapping::kSignedCoordinateMaximum, reference_end))
      return fail(AnchoringRefusal::BundleAnchorRange);
    if (i > 0 &&
        selected_pool_anchor_precedes(anchor, bundle.selected_pool_anchors[i - 1]))
      return fail(AnchoringRefusal::BundleRawOrder);
  }
  ExactAnchorPathSummary summary;
  AnchoringRefusal summary_refusal = AnchoringRefusal::None;
  if (!summarize_exact_anchor_path(bundle, summary, &summary_refusal)) {
    if (refusal) *refusal = summary_refusal;
    return false;
  }
  if (hash_exact_anchor_path(bundle) != bundle.anchor_path_hash)
    return fail(AnchoringRefusal::BundleHash);
  if (refusal) *refusal = AnchoringRefusal::None;
  return true;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
