#include "rank2_arbitration.h"

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

std::optional<int>
realized_primary_dp_maximum(const RnaSpliceRealizationResult& result) noexcept {
  if (result.refused || result.segments.empty())
    return std::nullopt;
  for (const auto& hypothesis : result.hypotheses) {
    if (hypothesis.orientation == result.winning_hypothesis &&
        !hypothesis.refused && !hypothesis.segments.empty())
      return hypothesis.primary_dp_maximum;
  }
  return std::nullopt;
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
