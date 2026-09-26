// Immutable handoff from the RNA exact-anchor-path chain to realization. Coordinates are
// in the mapping-oriented query frame: the read for forward mappings, its reverse
// complement for reverse ones. Only final output intervals are mirrored back.
#pragma once

#include "refusal.h"

#include <cstdint>
#include <string>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct SelectedPoolAnchor {
  uint32_t stable_id = 0;
  int32_t query_begin = 0;
  int32_t reference_begin = 0;
  int32_t span = 0;
  uint32_t source_seed_id = 0;
  uint32_t source_occurrence = 0;
  bool ignore_flag = false;
  bool tandem_flag = false;
  bool long_join_flag = false;
  // Provenance: all three fields enter anchor_path_hash and are checked by validation.
  uint32_t region_kind_mask = 0;
  std::vector<uint32_t> region_ids;
  std::vector<uint32_t> core_node_ids;
};

struct ExactAnchorPath {
  static constexpr uint32_t kSchemaVersion = 1;

  uint32_t schema_version = kSchemaVersion;
  uint64_t read_key = 0;
  int32_t query_length = 0;
  int32_t reference_id = -1;
  bool mapping_reverse = false;
  std::vector<SelectedPoolAnchor> selected_pool_anchors;
  std::vector<uint32_t> selected_raw_indices;
  int32_t anchor_path_score = 0;
  uint32_t anchor_path_anchor_count = 0;
  uint64_t candidate_pool_hash = 0;
  uint64_t anchor_path_hash = 0;
};

struct ExactAnchorPathResult {
  ExactAnchorPath bundle;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

// Endpoint facts about the selected path, computed with validation of the selected
// indices and without building chaining::Anchor copies.
struct ExactAnchorPathSummary {
  uint64_t selected_count = 0;
  int32_t first_query_begin = 0;
  int32_t last_query_end = 0;
  uint64_t outer_query_span = 0;
  int32_t first_reference_begin = 0;
  int32_t last_reference_end = 0;
  uint64_t outer_reference_span = 0;
};

bool validate_exact_anchor_path(const ExactAnchorPath& bundle,
                                AnchoringRefusal* refusal = nullptr);

bool summarize_exact_anchor_path(const ExactAnchorPath& bundle,
                                 ExactAnchorPathSummary& summary,
                                 AnchoringRefusal* refusal = nullptr) noexcept;

// Query bases covered by the union of the selected anchors' query intervals, read by the
// PAF matches field and the MAPQ coverage brake. False, with covered_bases 0, when an
// interval leaves the query or the union is empty.
bool anchor_union_query_coverage(const ExactAnchorPath& bundle,
                                 uint64_t& covered_bases) noexcept;

uint64_t hash_exact_anchor_path(const ExactAnchorPath& bundle);

// Stable read identity: FNV-1a over the upper-cased read sequence.
uint64_t hash_rna_read_key(const std::string& read) noexcept;

// Canonical anchor order, shared by bundle validation and the anchor-path pipeline.
bool selected_pool_anchor_precedes(const SelectedPoolAnchor& left,
                                   const SelectedPoolAnchor& right) noexcept;

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
