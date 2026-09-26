#pragma once

#include "splice_realizer.h"

#include <cstdint>
#include <optional>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The promotion statistic corresponding to minimap2's post-DP dp_max sort.
// Unavailable or refused winning hypotheses have no comparable score.
std::optional<int>
realized_primary_dp_maximum(const RnaSpliceRealizationResult& result) noexcept;

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
