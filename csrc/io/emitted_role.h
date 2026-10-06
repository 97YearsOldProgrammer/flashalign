#pragma once

#include "../core/types.h"

#include <stdexcept>

namespace fa::cpu::output {

// A record's role in the output, from its place in the result tree.
enum class EmittedRole {
  Primary,
  Supplementary,
  Secondary,
  SecondarySupplementary,
};

inline bool emitted_role_is_secondary(EmittedRole role) noexcept {
  return role == EmittedRole::Secondary ||
         role == EmittedRole::SecondarySupplementary;
}

inline bool emitted_role_is_supplementary(EmittedRole role) noexcept {
  return role == EmittedRole::Supplementary ||
         role == EmittedRole::SecondarySupplementary;
}

inline int emitted_role_sam_flag(EmittedRole role) noexcept {
  return (emitted_role_is_secondary(role) ? 0x100 : 0) |
         (emitted_role_is_supplementary(role) ? 0x800 : 0);
}

inline char emitted_role_paf_type(EmittedRole role) noexcept {
  return emitted_role_is_secondary(role) ? 'S' : 'P';
}

// The tp:A type as minimap2's write_tags sets it:
//   type = r->id == r->parent ? (r->inv? 'I' : 'P') : (r->inv? 'i' : 'S');
// Only DNA local-inversion records (the mm_align1_inv port) get 'I' or 'i'.
// An all-chains record is 'S' in any role: minimap2's -P sets no parent.
inline char emitted_role_paf_type(EmittedRole role,
                                  const AlignResult& result) noexcept {
  if (result.origin == AlignmentOrigin::DnaLocalInversion)
    return emitted_role_is_secondary(role) ? 'i' : 'I';
  if (result.origin == AlignmentOrigin::DnaAllChains)
    return 'S';
  return emitted_role_paf_type(role);
}

inline void validate_one_level_hypotheses(const AlignResult& result) {
  for (const AlignResult& secondary : result.secondary) {
    if (!secondary.secondary.empty()) {
      throw std::invalid_argument(
          "alignment secondary hypotheses may not be nested");
    }
    for (const AlignResult& supplementary : secondary.supplementary) {
      if (!supplementary.secondary.empty() ||
          !supplementary.supplementary.empty()) {
        throw std::invalid_argument(
            "alignment secondary supplementary segments must be leaves");
      }
    }
  }
  for (const AlignResult& supplementary : result.supplementary) {
    if (!supplementary.secondary.empty() ||
        !supplementary.supplementary.empty()) {
      throw std::invalid_argument(
          "alignment supplementary segments must be leaves");
    }
  }
}

}  // namespace fa::cpu::output
