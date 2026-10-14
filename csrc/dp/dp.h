// DP options and result types, and the declarations of the ksw2 controller and the splice kernel wrapper.
#pragma once

#include "ksw2.h" // ksw2 kernels, ksw_extz_t and flags

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <vector>

/**** params ****/

namespace fa { namespace cpu {

// DP options for the controller (dp.h), as minimap2's: scoring (-A/-B/-O/-E/-O2/-E2/
// --score-N) and ksw2 control (-r bands, Z-drop, end bonus, max_sw_mat). Penalties are stored
// positive. -g max_gap is not here: as in minimap2 it sizes the DP window around the kernel.
// The defaults are minimap2's map-ont row; presets override them per call.
struct DpMapOpt {
    // scoring (minimap2 -A/-B/--score-N and the two gap systems -O/-E,-O2/-E2)
    int a = 2;            // -A  match            (>0)
    int b = 4;            // -B  mismatch penalty (stored positive)
    int sc_ambi = 1;      // --score-N            (stored positive)
    int q = 4,  e = 2;    // -O/-E   gap system 1
    int q2 = 24, e2 = 1;  // -O2/-E2 gap system 2

    // control
    int zdrop = 400;      // --z-drop : off-diagonal drop-off that stops/splits DP
    int zdrop_inv = 200;  // minimap2 opt->zdrop_inv: the local-inversion Z-drop.
                          // Gates mm_test_zdrop's reverse-complement probe and
                          // becomes the retry Z-drop on a code-2 verdict.
    int end_bonus = -1;   // flank bonus for reaching the query end (EXTZ_ONLY only)
    int junc_bonus = 9;   // --junc-bonus; inert while junc == nullptr
    int junc_pen   = 5;   // --junc-pen; inert while junc == nullptr
    int bw = 500;         // -r flank band      (narrow; pre ×1.5)
    int bw_long = 20000;  // -r internal band   (wide;   pre ×1.5)
    int64_t max_sw_mat = 100000000; // per-call tlen*qlen cap; over it -> skip (zdropped)

    // Effective bands, as mm_align1 (bw * 1.5 + 1); the internal band is never below the
    // flank band.
    int bw_eff()      const { return static_cast<int>(bw * 1.5 + 1.0); }
    int bw_long_eff() const {
        const int v = static_cast<int>(bw_long * 1.5 + 1.0);
        const int f = bw_eff();
        return v < f ? f : v;
    }
    // Kernel selection, mirroring mm_align_pair(): equal gap systems -> single
    // affine (ksw_extz2_sse); unequal -> dual affine (ksw_extd2_sse).
    bool single_affine() const { return q == q2 && e == e2; }
};

}}  // namespace fa::cpu

/**** result ****/

namespace fa { namespace cpu {

struct GotohCigarResult {
    std::string cigar;
    int ref_offset;      // ref bases skipped at start (semi-global)
    int ref_consumed;    // ref bases consumed (M + D count)
    int matches;         // exact character matches
    int score = 0;       // terminal DP score for accepted tracebacks
    // On a Z-drop the traceback ends at the best-scoring cell and the controller splits
    // there; max_q/max_t are ksw2's zero-based coordinates of that cell.
    bool zdropped = false;
    int max_q = -1;      // raw zero-based query max-cell coordinate
    int max_t = -1;      // raw zero-based target max-cell coordinate
};

}}  // namespace fa::cpu

/**** ksw2_align ****/
// ksw2 helpers for the DP controller (dp.h): the score matrix.

namespace fa { namespace cpu {

// minimap2's m=5 score matrix (A,C,G,T,N), as ksw_gen_simple_mat: +|match| on the diagonal,
// -|mismatch| off it, -|ambi| on the N row and column.
void ksw2_simple_mat(int8_t* mat, int match, int mismatch, int ambi);

}}  // namespace fa::cpu

/**** control ****/
// minimap2-style DP controller over the ksw2 kernels, as align.c: kernel selection, the
// per-call matrix cap, flank extension and gap fill with minimap2's bands, end bonus and
// Z-drop. It reports a Z-drop and its max cell; splitting the chain is the caller's job.

namespace fa {
namespace cpu {
namespace dp {

enum class DpInversionProbeStatus : std::uint8_t {
  Disabled,
  DropBelowInversionThreshold,
  InvalidBounds,
  WeakReverseScore,
  Accepted,
  // The control's gate kept the probe from running.
  Gated,
};

// A count over a drop window (query_begin, query_end, target_begin,
// target_end; fill-local, half-open), from the caller's source.
using DpInversionProbeGate = int (*)(const void* source, int query_begin,
                                     int query_end, int target_begin,
                                     int target_end);

struct DpInversionProbeControl {
  bool enabled = false;
  int max_gap = 0;
  int min_chain_score = 0;
  int min_dp_max = 0;
  // When set, the probe runs only where the gate counts at least
  // gate_min_count over the drop window.
  DpInversionProbeGate gate = nullptr;
  const void* gate_source = nullptr;
  int gate_min_count = 0;
};

struct DpLocalScoreResult {
  bool attempted = false;
  int score = 0;
  int query_end = -1;
  int target_end = -1;
};

// Result of minimap2's mm_test_zdrop. Coordinates are half-open target/query offsets of
// the largest drop. Code 2 (inversion) is possible only with a DpInversionProbeControl.
struct DpZdropResult {
  bool tested = false;
  int code = 0;
  int max_drop = 0;
  int target_begin = -1;
  int target_end = -1;
  int query_begin = -1;
  int query_end = -1;
  DpInversionProbeStatus inversion_status = DpInversionProbeStatus::Disabled;
  DpLocalScoreResult reverse_probe;
};

// The ksw_extz_t fields a caller needs to soft-clip at the max cell, snap to the query end
// on reach_end, or split on a Z-drop. The CIGAR is BAM-packed (len << 4 | op, with
// 0=M 1=I 2=D 3=N); dp_to_gotoh() renders it as text.
struct DpAlnResult {
  bool has_cigar = false;
  bool zdropped = false;  // Z-drop fired -> caller may split the chain
  bool reach_end = false; // extension reached the query end (flank snapped)
  int max = 0, max_q = -1, max_t = -1; // best cell  (soft-clip / split point)
  int mqe = 0, mqe_t = -1;             // best cell reaching the query end
  int mte = 0, mte_q = -1;             // best cell reaching the target end
  int corner_score = 0;                // ez.score (both ends consumed; global)
  std::vector<uint32_t> cigar;         // BAM-packed
  DpZdropResult zdrop_test;
  bool exact_retry = false;
  int exact_retry_zdrop = -1;
};

enum class DpKernelId : uint8_t {
  Extz2,
  Extd2,
};

struct DpAttemptObserver {
  using Callback = void (*)(void* user_data, DpKernelId kernel, int flags,
                            int band, int zdrop, bool approximate,
                            const DpAlnResult& raw);

  Callback callback = nullptr;
  void* user_data = nullptr;
};

// Local Smith-Waterman score and end points (ksw_ll_i16) for the inversion checks.
DpLocalScoreResult dp_local_score(const uint8_t* qseq, int qlen,
                                  const uint8_t* tseq, int tlen,
                                  const int8_t* mat, int gap_open,
                                  int gap_extend);

// Flank extension (5' if to_left, else 3'), as mm_align1: narrow band, EXTZ_ONLY, end bonus
// and Z-drop. ksw2 backtracks from the max cell unless mqe + end_bonus beats it
// (reach_end). The left flank is aligned reversed with RIGHT | REV_CIGAR, so its CIGAR
// comes back in forward orientation.
DpAlnResult dp_extend(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                      int tlen, const uint8_t* tseq, bool to_left,
                      const DpAttemptObserver* observer = nullptr);

// Opposite-strand middle of a local inversion, as minimap2's mm_align1_inv: an open-ended
// extension over caller-bounded slices at band int(opt.bw * 1.5), with no end bonus.
DpAlnResult
dp_align_inversion_middle(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                          int tlen, const uint8_t* tseq,
                          const DpAttemptObserver* observer = nullptr);

// Gap fill between anchors, as mm_align1: global alignment over the wide band (or
// max(qlen, tlen) for a long join) with Z-drop. The first pass is approximate
// (KSW_EZ_APPROX_MAX); when dp_test_zdrop confirms a drop, one exact pass reruns. The
// result carries zdropped and the max cell so the caller can split.
DpAlnResult
dp_fill_gap(const DpMapOpt& opt, int qlen, const uint8_t* qseq, int tlen,
            const uint8_t* tseq, bool long_join = false,
            const DpAttemptObserver* observer = nullptr,
            const DpInversionProbeControl* inversion_control = nullptr);

// Text CIGAR and ref/match accounting for a result, over the same (qseq, tseq). The score
// follows minimap2: the corner score for a gap fill (`global`), mqe for a flank on
// reach_end, else the max cell. ksw2 emits M, so matches are counted by comparing bases.
// With render_text false, g.cigar stays empty and only the accounting is computed.
GotohCigarResult dp_to_gotoh(const DpAlnResult& r, const uint8_t* qseq,
                             int qlen, const uint8_t* tseq, int tlen,
                             bool global, bool render_text = true);

} // namespace dp
} // namespace cpu
} // namespace fa

/**** splice_kernel ****/
// C++ wrapper around ksw2's splice kernel (ksw_exts2_sse) and its local-score kernel.

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
