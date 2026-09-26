#include "splice_kernel.h"

#include "ksw2.h"
#include "ksw_arena.h"

#include <cstdlib>
#include <memory>
#include <new>
#include <utility>

namespace fa {
namespace cpu {
namespace dp {

namespace {

// Every scalar back to its declared default, keeping the CIGAR buffer's
// capacity. Must stay in step with SpliceKernelResult's member initializers.
void reset_result(SpliceKernelResult& result) noexcept {
  result.succeeded = false;
  result.failure = SpliceKernelFailure::None;
  result.score = 0;
  result.maximum = 0;
  result.maximum_query = -1;
  result.maximum_reference = -1;
  result.query_end_score = 0;
  result.query_end_reference = -1;
  result.reference_end_score = 0;
  result.reference_end_query = -1;
  result.reached_end = false;
  result.zdropped = false;
  result.cigar.clear();
}

}  // namespace

std::array<std::int8_t, 25> make_splice_score_matrix(
    const SpliceKernelScoring& scoring, int transition) noexcept {
  std::array<std::int8_t, 25> matrix{};
  const std::int8_t match = static_cast<std::int8_t>(
      scoring.match < 0 ? -scoring.match : scoring.match);
  const std::int8_t mismatch = static_cast<std::int8_t>(
      scoring.mismatch > 0 ? -scoring.mismatch : scoring.mismatch);
  const std::int8_t ambiguous = static_cast<std::int8_t>(
      scoring.ambiguous > 0 ? -scoring.ambiguous : scoring.ambiguous);
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column)
      matrix[static_cast<std::size_t>(row * 5 + column)] =
          row == column ? match : mismatch;
    matrix[static_cast<std::size_t>(row * 5 + 4)] = ambiguous;
  }
  for (int column = 0; column < 5; ++column)
    matrix[static_cast<std::size_t>(20 + column)] = ambiguous;
  if (transition != 0 && transition != scoring.mismatch) {
    const std::int8_t value = static_cast<std::int8_t>(
        transition > 0 ? -transition : transition);
    matrix[2] = value;
    matrix[8] = value;
    matrix[10] = value;
    matrix[16] = value;
  }
  return matrix;
}

struct SpliceKernelWorkspace::State {
  ksw_extz_t extension{};
  // One DP arena per workspace. ez.cigar grows by krealloc from NULL, so it lives outside
  // the arena and is reused across resets.
  ksw_arena_t* arena = ksw_arena_create();

  ~State() {
    std::free(extension.cigar);
    ksw_arena_destroy(arena);
  }
};

SpliceKernelWorkspace::SpliceKernelWorkspace()
    : state_(std::make_unique<State>()) {}

SpliceKernelWorkspace::~SpliceKernelWorkspace() = default;
SpliceKernelWorkspace::SpliceKernelWorkspace(SpliceKernelWorkspace&&) noexcept =
    default;
SpliceKernelWorkspace& SpliceKernelWorkspace::operator=(
    SpliceKernelWorkspace&&) noexcept = default;

void SpliceKernelWorkspace::align(const SpliceKernelRequest& request,
                                  SpliceKernelResult& result) noexcept {
  reset_result(result);
  if (!state_ || !request.query || !request.reference ||
      !request.score_matrix || request.query_length <= 0 ||
      request.reference_length <= 0) {
    result.failure = SpliceKernelFailure::InvalidRequest;
    return;
  }
  try {
    ksw_arena_reset(state_->arena);
    ksw_exts2_sse(
        state_->arena, request.query_length, request.query,
        request.reference_length, request.reference, 5, request.score_matrix,
        static_cast<std::int8_t>(request.scoring.gap_open),
        static_cast<std::int8_t>(request.scoring.gap_extend),
        static_cast<std::int8_t>(request.scoring.long_gap_open),
        // `noncan`: unused, since kSpliceComplex is always set below.
        0,
        request.zdrop, request.end_bonus,
        static_cast<std::int8_t>(request.scoring.junction_bonus),
        static_cast<std::int8_t>(request.scoring.junction_penalty),
        request.flags | kSpliceComplex, request.junction, &state_->extension);
    const ksw_extz_t& extension = state_->extension;
    result.score = extension.score;
    result.maximum = static_cast<int>(extension.max);
    result.maximum_query = extension.max_q;
    result.maximum_reference = extension.max_t;
    result.query_end_score = extension.mqe;
    result.query_end_reference = extension.mqe_t;
    result.reference_end_score = extension.mte;
    result.reference_end_query = extension.mte_q;
    result.reached_end = extension.reach_end != 0;
    result.zdropped = extension.zdropped != 0;
    if (extension.n_cigar > 0 && extension.cigar != nullptr)
      result.cigar.assign(extension.cigar,
                          extension.cigar + extension.n_cigar);
    result.succeeded = extension.n_cigar >= 0 &&
        (extension.n_cigar == 0 || extension.cigar != nullptr);
    if (!result.succeeded) result.failure = SpliceKernelFailure::KernelFailure;
  } catch (const std::bad_alloc&) {
    reset_result(result);
    result.failure = SpliceKernelFailure::AllocationFailure;
  } catch (...) {
    reset_result(result);
    result.failure = SpliceKernelFailure::KernelFailure;
  }
}

LocalAlignmentScore local_alignment_score(
    const std::uint8_t* query, int query_length,
    const std::uint8_t* reference, int reference_length,
    const std::array<std::int8_t, 25>& score_matrix,
    int gap_open, int gap_extend) noexcept {
  LocalAlignmentScore result;
  if (!query || !reference || query_length <= 0 || reference_length <= 0) {
    result.failure = SpliceKernelFailure::InvalidRequest;
    return result;
  }
  using Profile = std::unique_ptr<void, decltype(&std::free)>;
  Profile profile(ksw_ll_qinit(
      nullptr, 2, query_length, query, 5, score_matrix.data()), &std::free);
  if (!profile) {
    result.failure = SpliceKernelFailure::AllocationFailure;
    return result;
  }
  result.score = ksw_ll_i16(profile.get(), reference_length, reference,
                            gap_open, gap_extend, &result.query_end,
                            &result.reference_end);
  result.succeeded = true;
  return result;
}

}  // namespace dp
}  // namespace cpu
}  // namespace fa
