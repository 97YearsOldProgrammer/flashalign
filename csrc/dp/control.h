// minimap2-style DP controller over the ksw2 kernels, as align.c: kernel selection, the
// per-call matrix cap, flank extension and gap fill with minimap2's bands, end bonus and
// Z-drop. It reports a Z-drop and its max cell; splitting the chain is the caller's job.
#pragma once

#include "ksw2.h"          // ksw2 kernels, ksw_extz_t and flags
#include "ksw_arena_tls.h" // per-thread DP arena for every kernel call here
#include "params.h"        // DpMapOpt
#include "result.h"        // GotohCigarResult
#include "ksw2_align.h" // ksw2_simple_mat() and the call counters

#include <algorithm>
#include <cstdint>
#include <cstdlib> // free
#include <limits>
#include <string>
#include <vector>

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

namespace detail {

inline void fill_result(DpAlnResult& out, const ksw_extz_t& ez,
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

inline void notify_attempt(const DpAttemptObserver* observer, DpKernelId kernel,
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
// ez.cigar and ez.m_cigar untouched). Splice alignment goes through splice_kernel.h.
inline void dp_align_pair(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                          int tlen, const uint8_t* tseq, const int8_t* mat,
                          int w, int end_bonus, int zdrop, int ksw_flag,
                          ksw_extz_t* ez,
                          const DpAttemptObserver* observer = nullptr) {
  if (opt.max_sw_mat > 0 &&
      static_cast<int64_t>(tlen) * qlen > opt.max_sw_mat) {
    ksw_reset_extz(ez); // too big: skip and report Z-dropped, as minimap2
    ez->zdropped = 1;
    return;
  }
  // Per-thread arena for the kernel's DP rows, traceback and offsets, reset per call.
  // ez->cigar is allocated outside it (see ksw_arena_tls.h).
  ksw_arena_t* const km = detail::thread_arena();
  ksw_arena_reset(km);
  DpKernelId kernel = DpKernelId::Extd2;
  if (opt.single_affine()) {
    kernel = DpKernelId::Extz2;
    ::fa::cpu::ksw2_extz2_call_counter().fetch_add(1u,
                                                   std::memory_order_relaxed);
    ksw_extz2_sse(km, qlen, qseq, tlen, tseq, 5, mat,
                  static_cast<int8_t>(opt.q), static_cast<int8_t>(opt.e), w,
                  zdrop, end_bonus, ksw_flag, ez);
  } else {
    ::fa::cpu::ksw2_extd2_call_counter().fetch_add(1u,
                                                   std::memory_order_relaxed);
    ksw_extd2_sse(km, qlen, qseq, tlen, tseq, 5, mat,
                  static_cast<int8_t>(opt.q), static_cast<int8_t>(opt.e),
                  static_cast<int8_t>(opt.q2), static_cast<int8_t>(opt.e2), w,
                  zdrop, end_bonus, ksw_flag, ez);
  }
  detail::notify_attempt(observer, kernel, ksw_flag, w, zdrop, *ez);
}

// Local Smith-Waterman score and end points (ksw_ll_i16) for the inversion checks.
inline DpLocalScoreResult dp_local_score(const uint8_t* qseq, int qlen,
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
inline DpZdropResult
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

// Flank extension (5' if to_left, else 3'), as mm_align1: narrow band, EXTZ_ONLY, end bonus
// and Z-drop. ksw2 backtracks from the max cell unless mqe + end_bonus beats it
// (reach_end). The left flank is aligned reversed with RIGHT | REV_CIGAR, so its CIGAR
// comes back in forward orientation.
inline DpAlnResult dp_extend(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                             int tlen, const uint8_t* tseq, bool to_left,
                             const DpAttemptObserver* observer = nullptr) {
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

// Opposite-strand middle of a local inversion, as minimap2's mm_align1_inv: an open-ended
// extension over caller-bounded slices at band int(opt.bw * 1.5), with no end bonus.
inline DpAlnResult
dp_align_inversion_middle(const DpMapOpt& opt, int qlen, const uint8_t* qseq,
                          int tlen, const uint8_t* tseq,
                          const DpAttemptObserver* observer = nullptr) {
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

// Gap fill between anchors, as mm_align1: global alignment over the wide band (or
// max(qlen, tlen) for a long join) with Z-drop. The first pass is approximate
// (KSW_EZ_APPROX_MAX); when dp_test_zdrop confirms a drop, one exact pass reruns. The
// result carries zdropped and the max cell so the caller can split.
inline DpAlnResult
dp_fill_gap(const DpMapOpt& opt, int qlen, const uint8_t* qseq, int tlen,
            const uint8_t* tseq, bool long_join = false,
            const DpAttemptObserver* observer = nullptr,
            const DpInversionProbeControl* inversion_control = nullptr) {
  DpAlnResult out;
  if (qlen <= 0 || tlen <= 0)
    return out;

  int8_t mat[25];
  ::fa::cpu::ksw2_simple_mat(mat, opt.a, opt.b, opt.sc_ambi);
  const int w = long_join ? std::max(qlen, tlen) : opt.bw_long_eff();

  ksw_extz_t ez{};
  dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, w, -1, opt.zdrop,
                KSW_EZ_APPROX_MAX, &ez, observer);
  // As minimap2, re-score the approximate traceback even when the pass itself Z-dropped,
  // and confirm a drop with one exact pass.
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
      return out;
    }
    std::vector<uint32_t> tmp;
    if (ez.n_cigar > 0 && ez.cigar != nullptr)
      tmp.assign(ez.cigar, ez.cigar + ez.n_cigar);
    out.zdrop_test =
        dp_test_zdrop(opt, qseq, tseq, tmp, mat, inversion_control);
    if (out.zdrop_test.code != 0) {
      free(ez.cigar);
      ez = ksw_extz_t{};
      const int retry_zdrop =
          out.zdrop_test.code == 2 ? opt.zdrop_inv : opt.zdrop;
      dp_align_pair(opt, qlen, qseq, tlen, tseq, mat, w, -1, retry_zdrop, 0,
                    &ez, observer);
      out.exact_retry = true;
      out.exact_retry_zdrop = retry_zdrop;
    }
  }
  detail::fill_result(out, ez);
  if (ez.cigar)
    free(ez.cigar);
  return out;
}

// Text CIGAR and ref/match accounting for a result, over the same (qseq, tseq). The score
// follows minimap2: the corner score for a gap fill (`global`), mqe for a flank on
// reach_end, else the max cell. ksw2 emits M, so matches are counted by comparing bases.
// With render_text false, g.cigar stays empty and only the accounting is computed.
inline GotohCigarResult dp_to_gotoh(const DpAlnResult& r, const uint8_t* qseq,
                                    int qlen, const uint8_t* tseq, int tlen,
                                    bool global, bool render_text = true) {
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
