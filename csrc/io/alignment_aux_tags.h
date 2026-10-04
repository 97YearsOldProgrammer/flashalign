// Auxiliary tag values shared by the PAF and SAM writers.
#pragma once

#include "../core/types.h"

#include <algorithm>
#include <optional>
#include <stdexcept>

namespace fa::cpu::output {

struct AlignmentAuxiliaryTags {
  std::optional<int> alignment_score;
  std::optional<int> edit_distance;
  std::optional<int> ambiguities;
  // Absent when the engine did not compute the value; no tag is written.
  std::optional<int> max_segment_score;     // ms:i
  std::optional<int> max_score_margin;      // md:i
  std::optional<int> chain_anchors;         // cm:i
  std::optional<int> chain_score;           // s1:i
  std::optional<int> secondary_chain_score; // s2:i
};

inline AlignmentAuxiliaryTags alignment_auxiliary_tags(
    const AlignResult& result) {
  AlignmentAuxiliaryTags tags;
  if (!result.mapped()) return tags;

  tags.alignment_score = result.score;
  if (result.alignment_accounting_valid) {
    tags.edit_distance = result.edit_distance;
    tags.ambiguities = result.ambiguities;
  }
  // md is minibwa's margin: this record's ms minus the best ms among its
  // mapped alternatives (0 when there are none). It may be negative.
  if (result.dp_max_segment >= 0) {
    tags.max_segment_score = result.dp_max_segment;
    int best_alternative = 0;
    for (const AlignResult& alternative : result.secondary) {
      if (!alternative.mapped() || alternative.dp_max_segment < 0)
        continue;
      best_alternative = std::max(best_alternative, alternative.dp_max_segment);
    }
    tags.max_score_margin = result.dp_max_segment - best_alternative;
  }
  if (result.chain_anchors >= 0)
    tags.chain_anchors = result.chain_anchors;
  if (result.chain_score >= 0)
    tags.chain_score = result.chain_score;
  if (result.secondary_chain_score >= 0)
    tags.secondary_chain_score = result.secondary_chain_score;
  return tags;
}

inline int factual_sa_edit_distance(const AlignResult& result) {
  const AlignmentAuxiliaryTags tags = alignment_auxiliary_tags(result);
  if (!tags.edit_distance) {
    throw std::logic_error(
        "mapped split-family member lacks exact alignment accounting");
  }
  return *tags.edit_distance;
}

}  // namespace fa::cpu::output
