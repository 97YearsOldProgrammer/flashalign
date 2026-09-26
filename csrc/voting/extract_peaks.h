#pragma once

#include "../seeding/types.h" // VotePeak, chain_peak_better

#include <vector>

namespace fa {
namespace cpu {
namespace lr {

// The top peak of a drained vote by chain_peak_better; nullptr iff peaks is empty. The
// first wins on ties.
inline const VotePeak *best_peak(const std::vector<VotePeak> &peaks) {
  const VotePeak *best = nullptr;
  for (const auto &cp : peaks) {
    if (!best || chain_peak_better(cp, *best))
      best = &cp;
  }
  return best;
}

} // namespace lr
} // namespace cpu
} // namespace fa
