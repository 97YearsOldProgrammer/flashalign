// The plain PAF lane's catalogue election. An intronless copy of a gene keeps the k-mers
// spanning each junction, which the spliced gene cannot have, so on anchor count alone
// the copy can outscore the gene. The score gives each junction back what it costs.
#pragma once

#include "anchoring/exact_anchor_path.h"

#include <cstddef>
#include <cstdint>

namespace fa::cpu::lr::rna {

// Scores are integers in thousandths of an anchor.
inline constexpr std::int64_t kRnaPlainElectScale = 1000;

// The credit per junction, in thousandths of an anchor. A junction removes about k - 1
// anchor positions, (k - 1) * 2 / (k - s + 1) = 14/3 anchors at the splice presets'
// default seeding (k = 15, s = 10). The credit is 4.7, the measured value, not 14/3,
// and is not rescaled for an index of another k or s.
inline constexpr std::int64_t kRnaPlainElectJunctionCredit = 4700;

// The chain's junctions: consecutive chained anchors whose reference gap exceeds their
// query gap by at least min_intron. `path` must be a bundle summarize_exact_anchor_path
// accepted, so its selected indices are in range and in chain order.
inline int rna_plain_elect_junctions(const ExactAnchorPath& path,
                                     int min_intron) noexcept {
  int junctions = 0;
  for (std::size_t i = 1; i < path.selected_raw_indices.size(); ++i) {
    const SelectedPoolAnchor& left =
        path.selected_pool_anchors[path.selected_raw_indices[i - 1]];
    const SelectedPoolAnchor& right =
        path.selected_pool_anchors[path.selected_raw_indices[i]];
    const std::int64_t reference_gap =
        static_cast<std::int64_t>(right.reference_begin) - left.reference_begin;
    const std::int64_t query_gap =
        static_cast<std::int64_t>(right.query_begin) - left.query_begin;
    if (reference_gap - query_gap >= min_intron)
      ++junctions;
  }
  return junctions;
}

// Anchor count plus the junction credit, in thousandths of an anchor.
inline std::int64_t rna_plain_elect_score(const ExactAnchorPath& path,
                                          int min_intron) noexcept {
  return static_cast<std::int64_t>(path.anchor_path_anchor_count) *
             kRnaPlainElectScale +
         rna_plain_elect_junctions(path, min_intron) *
             kRnaPlainElectJunctionCredit;
}

} // namespace fa::cpu::lr::rna
