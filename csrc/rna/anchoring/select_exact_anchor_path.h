// Runs the spliced colinear chain over an RNA candidate pool and converts the best chain
// into an ExactAnchorPath.
#pragma once

#include "../../chaining/colinear_params.h"
#include "exact_anchor_path.h"
#include "skeleton_harvest.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The chains the recurrence built but the selection did not take; the best competing one
// feeds the chain MAPQ. Every field stays -1 until the recurrence runs.
struct SiblingChainStats {
  // Best score among the non-selected chains that compete with the selected one under
  // minimap2's mask_level (rna_query_spans_compete); a query-disjoint piece of the same
  // read is not an alternative.
  int32_t sib_score = -1;
  int32_t sib_count = -1;     // that chain's anchor count
  int32_t chains_total = -1;  // how many chains the recurrence returned
};

// A chain of the window other than the selected one that does not compete with it under
// mask_level: its own path over the same pool, and its own competing sibling.
struct DisjointAnchorPath {
  ExactAnchorPath bundle;
  SiblingChainStats siblings;
};

// `disjoint`, when non-null, is cleared and, on an accepted selection, receives every
// other chain that does not compete with the selected one, in the recurrence's order; a
// chain whose path fails validation is left out.
ExactAnchorPathResult select_exact_anchor_path(
    uint64_t read_key, int32_t query_length, int32_t reference_id,
    bool mapping_reverse, const CandidatePool& pool,
    const chaining::SplicedChainParams& params,
    SiblingChainStats* siblings = nullptr,
    std::vector<DisjointAnchorPath>* disjoint = nullptr);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
