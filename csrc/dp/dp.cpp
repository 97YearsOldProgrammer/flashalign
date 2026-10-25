// The minimap2-style controller over the ksw2 kernels: kernel selection, the band, Z-drop, the per-thread arena, the splice kernel wrapper.
#include "dp.h"

#include "ksw2.h"
#include "ksw_arena.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

/**** ksw2_align ****/

namespace fa { namespace cpu {

void ksw2_simple_mat(int8_t* mat, int match, int mismatch, int ambi) {
    const int8_t a = static_cast<int8_t>(match < 0 ? -match : match);
    const int8_t b = static_cast<int8_t>(mismatch > 0 ? -mismatch : mismatch);
    const int8_t sc_ambi = static_cast<int8_t>(ambi > 0 ? -ambi : ambi);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) mat[i * 5 + j] = (i == j) ? a : b;
        mat[i * 5 + 4] = sc_ambi;
    }
    for (int j = 0; j < 5; ++j) mat[4 * 5 + j] = sc_ambi;
}

}}  // namespace fa::cpu

/**** ksw_arena_tls ****/
// Per-thread ksw arena for the DNA DP controller, created on first use and reset once per
// kernel call, like minimap2's per-thread kalloc pool. ez->cigar grows by krealloc from NULL,
// so it stays on the system heap and the controller's free(ez.cigar) is correct.
// The RNA splice kernel owns its own arena.

namespace fa {
namespace cpu {
namespace dp {
namespace detail {

struct ArenaHolder {
  ksw_arena_t* arena = ksw_arena_create();
  ~ArenaHolder() { ksw_arena_destroy(arena); }
  ArenaHolder() = default;
  ArenaHolder(const ArenaHolder&) = delete;
  ArenaHolder& operator=(const ArenaHolder&) = delete;
};

static ArenaHolder& thread_arena_holder() {
  static thread_local ArenaHolder holder;
  return holder;
}

static ksw_arena_t* thread_arena() {
  return thread_arena_holder().arena; // may be NULL: the system allocator is used
}

} // namespace detail

std::size_t thread_arena_bytes() {
  return ksw_arena_capacity(detail::thread_arena());
}

void release_thread_arena() {
  detail::ArenaHolder& holder = detail::thread_arena_holder();
  ksw_arena_destroy(holder.arena);
  holder.arena = ksw_arena_create();
}

} // namespace dp
} // namespace cpu
} // namespace fa

/**** control ****/

namespace fa {
namespace cpu {
namespace dp {

namespace detail {

static void fill_result(DpAlnResult& out, const ksw_extz_t& ez,
                        bool copy_cigar = true) {
  out.zdropped = ez.zdropped != 0;
  out.reach_end = ez.reach_end != 0;
  out.max = static_cast<int>(ez.max);
  out.max_q = ez.max_q;
  out.max_t = ez.max_t;
  out.mqe = ez.mqe;
  out.mqe_t = ez.mqe_t;
  out.mte = ez.mte;
  out.mte_q = ez.mte_q;
  out.corner_score = ez.score;
  if (ez.n_cigar > 0 && ez.cigar != nullptr) {
    if (copy_cigar)
      out.cigar.assign(ez.cigar, ez.cigar + ez.n_cigar);
    out.has_cigar = true;
  }
}

static void notify_attempt(const DpAttemptObserver* observer, DpKernelId kernel,
                           int flags, int band, int zdrop, const ksw_extz_t& ez) {
  if (observer == nullptr || observer->callback == nullptr)
    return;
  DpAlnResult raw;
  // Observers read only metadata, so the CIGAR is not copied (has_cigar is still set).
  constexpr bool kAttemptCigarCopy = false;
  fill_result(raw, ez, kAttemptCigarCopy);
  const bool approximate =
      (flags & (KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP)) != 0;
  observer->callback(observer->user_data, kernel, flags, band, zdrop,
                     approximate, raw);
}

} // namespace detail

// Kernel dispatch, as minimap2's mm_align_pair: single or dual affine, with the per-call
// max_sw_mat guard. `ez` must be value-initialized by the caller (ksw_reset_extz leaves
// ez.cigar and ez.m_cigar untouched), and so must `zstop`, which the dual-affine kernel
// fills (ksw_extd2_sse_zstop). Splice alignment goes through dp.h.
static void dp_align_pair(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                          int tlen, const uint8_t* tseq, const int8_t* mat,
                          int w, int end_bonus, int zdrop, int ksw_flag,
                          ksw_extz_t* ez,
                          const DpAttemptObserver* observer = nullptr,
                          ksw_zstop_t* zstop = nullptr) {
  if (opt.max_sw_mat > 0 &&
      static_cast<int64_t>(tlen) * qlen > opt.max_sw_mat) {
    ksw_reset_extz(ez); // too big: skip and report Z-dropped, as minimap2
    ez->zdropped = 1;
    return;
  }
  // Per-thread arena for the kernel's DP rows, traceback and offsets, reset per call.
  // ez->cigar is allocated outside it (see dp.cpp).
  ksw_arena_t* const km = detail::thread_arena();
  ksw_arena_reset(km);
  DpKernelId kernel = DpKernelId::Extd2;
  if (opt.single_affine()) {
    kernel = DpKernelId::Extz2;
    ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat,
                  static_cast<int8_t>(opt.q), static_cast<int8_t>(opt.e), w,
                  zdrop, end_bonus, ksw_flag, ez);
  } else if (zstop != nullptr) {
    ksw_extd2_sse_zstop(km, qlen, qseq, tlen, tseq, 5, mat,
                        static_cast<int8_t>(opt.q), static_cast<int8_t>(opt.e),
                        static_cast<int8_t>(opt.q2),
                        static_cast<int8_t>(opt.e2), w, zdrop, end_bonus,
                        ksw_flag, ez, zstop);
  } else {
    ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat,
                  static_cast<int8_t>(opt.q), static_cast<int8_t>(opt.e),
                  static_cast<int8_t>(opt.q2), static_cast<int8_t>(opt.e2), w,
                  zdrop, end_bonus, ksw_flag, ez);
  }
  detail::notify_attempt(observer, kernel, ksw_flag, w, zdrop, *ez);
}

DpLocalScoreResult dp_local_score(const uint8_t* qseq, int qlen,
                                  const uint8_t* tseq, int tlen,
                                  const int8_t* mat, int gap_open,
                                  int gap_extend) {
  DpLocalScoreResult result;
  if (qseq == nullptr || tseq == nullptr || mat == nullptr || qlen <= 0 ||
      tlen <= 0) {
    return result;
  }
  void* profile = ksw_ll_qinit(nullptr, 2, qlen, qseq, 5, mat);
  if (profile == nullptr)
    return result;
  result.attempted = true;
  result.score = ksw_ll_i16(profile, tlen, tseq, gap_open, gap_extend,
                            &result.query_end, &result.target_end);
  std::free(profile);
  return result;
}

// Port of minimap2's mm_test_zdrop: re-scores the CIGAR over (qseq, tseq), records the
// region with the largest diagonal-corrected score drop and returns code 1 when it exceeds
// opt.zdrop. With an enabled inversion control it also runs the reverse-complement probe
// over that region, unless the control's gate is closed there, and returns code 2 when
// both reverse-score floors pass.
static DpZdropResult
dp_test_zdrop(const DpMapOpt& opt, const uint8_t* qseq, const uint8_t* tseq,
              const std::vector<uint32_t>& cigar, const int8_t* mat,
              const DpInversionProbeControl* inversion_control = nullptr) {
  DpZdropResult result;
  result.tested = true;
  int32_t score = 0, max = std::numeric_limits<int32_t>::min();
  int32_t max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
  const int e = opt.e;
  auto upd = [&](int32_t sc, int ii, int jj) {
    if (sc < max) {
      const int li = ii - max_i, lj = jj - max_j;
      const int diff = li > lj ? li - lj : lj - li;
      const int z = max - sc - diff * e;
      if (z > max_zdrop) {
        max_zdrop = z;
        result.target_begin = max_i;
        result.target_end = ii;
        result.query_begin = max_j;
        result.query_end = jj;
      }
    } else {
      max = sc;
      max_i = ii;
      max_j = jj;
    }
  };
  for (uint32_t k = 0; k < cigar.size(); ++k) {
    const uint32_t op = cigar[k] & 0xf, len = cigar[k] >> 4;
    if (op == KSW_CIGAR_MATCH) {
      for (uint32_t l = 0; l < len; ++l) {
        score += mat[tseq[i + l] * 5 + qseq[j + l]];
        upd(score, i + l, j + l);
      }
      i += len;
      j += len;
    } else if (op == KSW_CIGAR_INS || op == KSW_CIGAR_DEL ||
               op == KSW_CIGAR_N_SKIP) {
      score -= opt.q + opt.e * static_cast<int>(len);
      if (op == KSW_CIGAR_INS)
        j += len;
      else
        i += len;
      upd(score, i, j);
    }
  }
  result.max_drop = max_zdrop;
  if (inversion_control != nullptr && inversion_control->enabled) {
    if (max_zdrop <= opt.zdrop_inv) {
      result.inversion_status =
          DpInversionProbeStatus::DropBelowInversionThreshold;
    } else {
      const int query_length = result.query_end - result.query_begin;
      const int target_length = result.target_end - result.target_begin;
      const bool valid_bounds =
          inversion_control->max_gap > 0 && query_length > 0 &&
          target_length > 0 && query_length < inversion_control->max_gap &&
          target_length < inversion_control->max_gap &&
          result.query_begin >= 0 && result.target_begin >= 0 &&
          result.query_end <= j && result.target_end <= i;
      if (!valid_bounds) {
        result.inversion_status = DpInversionProbeStatus::InvalidBounds;
      } else if (inversion_control->gate != nullptr &&
                 inversion_control->gate(inversion_control->gate_source,
                                         result.query_begin, result.query_end,
                                         result.target_begin,
                                         result.target_end) <
                     inversion_control->gate_min_count) {
        result.inversion_status = DpInversionProbeStatus::Gated;
      } else {
        std::vector<uint8_t> reverse_query(
            static_cast<std::size_t>(query_length));
        for (int offset = 0; offset < query_length; ++offset) {
          const int base = qseq[result.query_end - offset - 1];
          reverse_query[static_cast<std::size_t>(offset)] =
              static_cast<uint8_t>(base >= 4 ? 4 : 3 - base);
        }
        result.reverse_probe = dp_local_score(
            reverse_query.data(), query_length, tseq + result.target_begin,
            target_length, mat, opt.q, opt.e);
        const std::int64_t chain_floor =
            static_cast<std::int64_t>(inversion_control->min_chain_score) *
            static_cast<std::int64_t>(opt.a);
        if (result.reverse_probe.attempted &&
            result.reverse_probe.query_end >= 0 &&
            result.reverse_probe.target_end >= 0 &&
            result.reverse_probe.score >= chain_floor &&
            result.reverse_probe.score >= inversion_control->min_dp_max) {
          result.inversion_status = DpInversionProbeStatus::Accepted;
          result.code = 2;
          return result;
        }
        result.inversion_status = DpInversionProbeStatus::WeakReverseScore;
      }
    }
  }
  result.code = max_zdrop > opt.zdrop ? 1 : 0;
  return result;
}

DpAlnResult dp_extend(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                      int tlen, const uint8_t* tseq, bool to_left,
                      const DpAttemptObserver* observer) {
  DpAlnResult out;
  if (qlen <= 0 || tlen <= 0)
    return out;

  int8_t mat[25];
  ::fa::cpu::ksw2_simple_mat(mat, opt.a, opt.b, opt.sc_ambi);

  int flag = KSW_EZ_EXTZ_ONLY;
  const uint8_t* q = qseq;
  const uint8_t* t = tseq;
  std::vector<uint8_t> qr, tr;
  if (to_left) {
    qr.assign(qseq, qseq + qlen);
    std::reverse(qr.begin(), qr.end());
    tr.assign(tseq, tseq + tlen);
    std::reverse(tr.begin(), tr.end());
    q = qr.data();
    t = tr.data();
    flag |= KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR;
  }

  ksw_extz_t ez{};
  dp_align_pair(opt, qlen, q, tlen, t, mat, opt.bw_eff(), opt.end_bonus,
                opt.zdrop, flag, &ez, observer);
  detail::fill_result(out, ez);
  if (ez.cigar)
    free(ez.cigar);
  return out;
}

DpAlnResult
dp_align_inversion_middle(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                          int tlen, const uint8_t* tseq,
                          const DpAttemptObserver* observer) {
  DpAlnResult out;
  if (qlen <= 0 || tlen <= 0)
    return out;
  int8_t mat[25];
  ::fa::cpu::ksw2_simple_mat(mat, opt.a, opt.b, opt.sc_ambi);
  ksw_extz_t ez{};
  const int band = static_cast<int>(opt.bw * 1.5);
  dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, band, -1, opt.zdrop,
                KSW_EZ_EXTZ_ONLY, &ez, observer);
  detail::fill_result(out, ez);
  if (ez.cigar)
    free(ez.cigar);
  return out;
}

static void free_zstop(ksw_zstop_t& zstop) {
  for (ksw_extz_t& stop : zstop.ez)
    free(stop.cigar);
}

DpAlnResult
dp_fill_gap(const DpMapOpt& opt, int qlen, const uint8_t* qseq, int tlen,
            const uint8_t* tseq, bool long_join,
            const DpAttemptObserver* observer,
            const DpInversionProbeControl* inversion_control,
            bool exact_first) {
  DpAlnResult out;
  if (qlen <= 0 || tlen <= 0)
    return out;

  int8_t mat[25];
  ::fa::cpu::ksw2_simple_mat(mat, opt.a, opt.b, opt.sc_ambi);
  const int w = long_join ? std::max(qlen, tlen) : opt.bw_long_eff();

  ksw_extz_t ez{};
  ksw_zstop_t zstop{};
  // Pass 1 is minimap2's approximate pass or, where the caller asks for it (the DNA
  // realization does when its window holds kDnaExactFirstUnconfirmed unconfirmed retries),
  // the exact pass run to the corner, recording where it would stop at zdrop and zdrop_inv.
  // Only the dual-affine kernel records, and only a fill with Z-drop retries. An exact
  // pass 1 that Z-drops without stopping (the band ran out or the matrix cap held) is
  // replaced by the approximate one.
  exact_first = exact_first && opt.zdrop >= 0 && !opt.single_affine();
  if (exact_first) {
    zstop.zdrop[0] = opt.zdrop;
    zstop.zdrop[1] = opt.zdrop_inv;
    dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, w, -1, -1, 0, &ez,
                  observer, &zstop);
    if (ez.zdropped) {
      free(ez.cigar);
      ez = ksw_extz_t{};
      exact_first = false;
    }
  }
  if (!exact_first)
    dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, w, -1, opt.zdrop,
                  KSW_EZ_APPROX_MAX, &ez, observer);
  // As minimap2, re-score pass 1's traceback even when the pass itself Z-dropped, and
  // confirm a drop with the exact pass.
  if (opt.zdrop >= 0) {
    // O(ops) bound before the O(bases) re-score. dp_test_zdrop charges gaps at system 1,
    // so its final score is T' = ez.score + sum over gaps of [min(s1, s2) - s1]. Every
    // prefix and suffix scores at most a*min(qlen,tlen), so no drop exceeds
    // 2*a*min(qlen,tlen) - T'; at or under the threshold the re-score returns code 0.
    // A Z-dropped or CIGAR-less pass has no corner score and is always re-scored.
    bool provably_no_drop = false;
    if (!ez.zdropped && ez.n_cigar > 0 && ez.cigar != nullptr) {
      const int threshold =
          inversion_control != nullptr && inversion_control->enabled
              ? std::min(opt.zdrop, opt.zdrop_inv)
              : opt.zdrop;
      std::int64_t rescore_total = ez.score;
      for (int k = 0; k < ez.n_cigar; ++k) {
        const uint32_t op = ez.cigar[k] & 0xf;
        if (op == KSW_CIGAR_INS || op == KSW_CIGAR_DEL) {
          const std::int64_t len = ez.cigar[k] >> 4;
          const std::int64_t s1 = opt.q + opt.e * len;
          const std::int64_t s2 = opt.q2 + opt.e2 * len;
          if (s2 < s1)
            rescore_total += s2 - s1;
        }
      }
      const std::int64_t half =
          static_cast<std::int64_t>(opt.a) * (qlen < tlen ? qlen : tlen);
      provably_no_drop = 2 * half - rescore_total <= threshold;
    }
    if (provably_no_drop) {
      out.zdrop_test.tested = true; // proven code 0; max_drop not computed
      detail::fill_result(out, ez);
      if (ez.cigar)
        free(ez.cigar);
      free_zstop(zstop);
      return out;
    }
    std::vector<uint32_t> tmp;
    if (ez.n_cigar > 0 && ez.cigar != nullptr)
      tmp.assign(ez.cigar, ez.cigar + ez.n_cigar);
    out.zdrop_test =
        dp_test_zdrop(opt, qseq, tseq, tmp, mat, inversion_control);
    if (out.zdrop_test.code != 0) {
      const int retry_zdrop =
          out.zdrop_test.code == 2 ? opt.zdrop_inv : opt.zdrop;
      if (exact_first) {
        // The exact pass at retry_zdrop is pass 1's recorded stop there, or pass 1 itself
        // where it records none.
        const int k = out.zdrop_test.code == 2 ? 1 : 0;
        if (zstop.hit[k]) {
          free(ez.cigar);
          ez = zstop.ez[k];
          zstop.ez[k].cigar = nullptr;
        }
      } else {
        free(ez.cigar);
        ez = ksw_extz_t{};
        dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, w, -1, retry_zdrop, 0,
                      &ez, observer);
      }
      out.exact_retry = true;
      out.exact_retry_zdrop = retry_zdrop;
    }
  }
  detail::fill_result(out, ez);
  if (ez.cigar)
    free(ez.cigar);
  free_zstop(zstop);
  return out;
}

GotohCigarResult dp_to_gotoh(const DpAlnResult& r, const uint8_t* qseq,
                             int qlen, const uint8_t* tseq, int tlen,
                             bool global, bool render_text) {
  GotohCigarResult g;
  g.cigar.clear();
  g.ref_offset = 0;
  g.ref_consumed = 0;
  g.matches = 0;
  g.score = 0;
  if (!r.has_cigar)
    return g;

  static const char OP[4] = {'M', 'I', 'D', 'N'};
  std::string text;
  int qi = 0, ti = 0, ref_consumed = 0, matches = 0;
  for (uint32_t c : r.cigar) {
    int op = static_cast<int>(c & 0xf);
    const int len = static_cast<int>(c >> 4);
    if (op < 0 || op > 3)
      op = 0;
    if (render_text) {
      text += std::to_string(len);
      text += OP[op];
    }
    if (op == 0) { // M: query + ref
      for (int x = 0; x < len; ++x)
        if (qi + x < qlen && ti + x < tlen && qseq[qi + x] == tseq[ti + x])
          ++matches;
      qi += len;
      ti += len;
      ref_consumed += len;
    } else if (op == 1) { // I: query only
      qi += len;
    } else { // D / N: ref only
      ti += len;
      ref_consumed += len;
    }
  }
  g.cigar = std::move(text);
  g.ref_consumed = ref_consumed;
  g.matches = matches;
  // A Z-dropped fill never reached the corner: report the max-cell score.
  g.score = r.zdropped
                ? r.max
                : (global ? r.corner_score : (r.reach_end ? r.mqe : r.max));
  g.zdropped = r.zdropped;
  g.max_q = r.max_q;
  g.max_t = r.max_t;
  return g;
}

} // namespace dp
} // namespace cpu
} // namespace fa

/**** splice_kernel ****/

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
