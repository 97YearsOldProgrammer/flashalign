// C++ wrapper around ksw2's splice kernel (ksw_exts2_sse) and its local-score kernel.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace fa {
namespace cpu {
namespace dp {

enum SpliceKernelFlag : int {
  kSpliceRightAlign = 0x02,
  kSpliceApproximateMaximum = 0x08,
  kSpliceExtensionOnly = 0x40,
  kSpliceReverseCigar = 0x80,
  kSpliceForward = 0x100,
  kSpliceReverse = 0x200,
  kSpliceFlank = 0x400,
  kSpliceComplex = 0x800,
};

struct SpliceKernelScoring {
  int match = 1;
  int mismatch = 2;
  int gap_open = 2;
  int gap_extend = 1;
  int long_gap_open = 32;
  int ambiguous = 1;
  int junction_bonus = 9;
  // Inert unless KSW_EZ_SPLICE_SCORE, which SpliceKernelFlag does not expose, is set.
  int junction_penalty = 5;
};

struct SpliceKernelRequest {
  const std::uint8_t* query = nullptr;
  int query_length = 0;
  const std::uint8_t* reference = nullptr;
  int reference_length = 0;
  // Optional annotated-junction mask over the reference slice; null means no junctions.
  const std::uint8_t* junction = nullptr;
  // Required 5x5 substitution matrix, owned by the caller.
  const std::int8_t* score_matrix = nullptr;
  int zdrop = 200;
  int end_bonus = -1;
  int flags = 0;
  SpliceKernelScoring scoring;
};

enum class SpliceKernelFailure : std::uint8_t {
  None = 0,
  InvalidRequest,
  AllocationFailure,
  KernelFailure,
};

struct SpliceKernelResult {
  bool succeeded = false;
  SpliceKernelFailure failure = SpliceKernelFailure::None;
  int score = 0;
  int maximum = 0;
  int maximum_query = -1;
  int maximum_reference = -1;
  int query_end_score = 0;
  int query_end_reference = -1;
  int reference_end_score = 0;
  int reference_end_query = -1;
  bool reached_end = false;
  bool zdropped = false;
  std::vector<std::uint32_t> cigar;
};

std::array<std::int8_t, 25> make_splice_score_matrix(
    const SpliceKernelScoring& scoring, int transition = 0) noexcept;

class SpliceKernelWorkspace {
 public:
  SpliceKernelWorkspace();
  ~SpliceKernelWorkspace();
  SpliceKernelWorkspace(const SpliceKernelWorkspace&) = delete;
  SpliceKernelWorkspace& operator=(const SpliceKernelWorkspace&) = delete;
  SpliceKernelWorkspace(SpliceKernelWorkspace&&) noexcept;
  SpliceKernelWorkspace& operator=(SpliceKernelWorkspace&&) noexcept;

  // Fills `result` in place. Keep one workspace per worker thread: its ksw2 arena and
  // `result`'s CIGAR buffer are reused across calls.
  void align(const SpliceKernelRequest& request,
             SpliceKernelResult& result) noexcept;

 private:
  struct State;
  std::unique_ptr<State> state_;
};

struct LocalAlignmentScore {
  bool succeeded = false;
  SpliceKernelFailure failure = SpliceKernelFailure::None;
  int score = 0;
  int query_end = -1;
  int reference_end = -1;
};

// Precondition, shared with SpliceKernelWorkspace::align: the sequences are nt4 (0..4);
// they are not re-validated here.
LocalAlignmentScore local_alignment_score(
    const std::uint8_t* query, int query_length,
    const std::uint8_t* reference, int reference_length,
    const std::array<std::int8_t, 25>& score_matrix,
    int gap_open, int gap_extend) noexcept;

}  // namespace dp
}  // namespace cpu
}  // namespace fa
