// The whole-read winner's per-read line d(q) = a + b q, with d = local - q. A read whose
// anchors drift by b per query base covers L + b L bases of reference, not L, so the
// harvest window [raw_ref_start - pad, raw_ref_start + L + pad] can leave its ends
// outside. The line is fitted on the winner's own seeds, and dna/backend.cpp widens the
// winner's window to its union with the line's span [a, a + L + stretch(L, b)]
// (VotePeak::harvest_below / harvest_above). The vote itself is unchanged.
//
// The fit is integer only (int64, __int128), so no target's float contraction can move a
// window edge:
//   S0    one (q, d) pair per seed of the winner lane's exact-refine set: its
//         strand-compatible posting on the winner's contig within one bin (width W) of
//         the winner's bin, nearest raw_ref_start. Least squares gives line 0.
//   S1    the same along line 0: per seed, the posting within max(W, 256) of the line
//         and nearest it. Least squares, one trim (residuals beyond max(3 MAD, 64)
//         dropped), least squares again: line 1, of slope b1.
//   gate  line 1 rests on at least 8 pairs spanning at least L / 4 of the query; b1 is
//         clamped to [-6 %, +10 %] and carried as round(b1 2^20). A read that fails the
//         gate keeps the plain window.
//   a     floor(d at q = 0) on the line of the clamped slope through line 1's centroid.
#pragma once

#include "../index/format.h"    // KmerPostingView, packed_ref_orientation_compatible
#include "../seeding/types.h"   // VotePeak, kVoteSlopeFractionBits

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// One anchor of the fit: q the seed's oriented read_pos, d = local - q.
struct VoteSlopePair {
  int q = 0;
  std::int64_t d = 0;
};

// Floor division by a positive divisor, correct for a negative numerator.
template <class T>
inline T vote_slope_floor_div(T numerator, T divisor) {
  T quotient = numerator / divisor;
  if (numerator % divisor != 0 && numerator < 0) --quotient;
  return quotient;
}

// A least-squares line through its pairs' exact centroid, slope in fixed point. Sums are
// taken about the integer centre (q0, d0) to keep products small; sx and sy are the
// centre's remainders, in [0, n).
struct VoteSlopeLine {
  bool valid = false;
  std::int64_t n = 0;
  std::int64_t q0 = 0;
  std::int64_t d0 = 0;
  std::int64_t sx = 0;
  std::int64_t sy = 0;
  // round(b 2^20), saturated to int32 so at() stays inside __int128.
  std::int64_t b_q20 = 0;

  // floor(d at q) on the line: d0 + (sy + b (n (q - q0) - sx)) / n.
  std::int64_t at(std::int64_t q) const {
    const __int128 one = __int128{1} << kVoteSlopeFractionBits;
    const __int128 numerator =
        static_cast<__int128>(sy) * one +
        static_cast<__int128>(b_q20) *
            (static_cast<__int128>(n) * (q - q0) - sx);
    return d0 + static_cast<std::int64_t>(vote_slope_floor_div<__int128>(
                    numerator, static_cast<__int128>(n) * one));
  }
};

inline VoteSlopeLine vote_slope_fit(const std::vector<VoteSlopePair>& pairs) {
  VoteSlopeLine line;
  const std::int64_t n = static_cast<std::int64_t>(pairs.size());
  if (n < 2) return line;
  std::int64_t sum_q = 0;
  std::int64_t sum_d = 0;
  for (const VoteSlopePair& pair : pairs) {
    sum_q += pair.q;
    sum_d += pair.d;
  }
  line.n = n;
  line.q0 = vote_slope_floor_div<std::int64_t>(sum_q, n);
  line.d0 = vote_slope_floor_div<std::int64_t>(sum_d, n);
  __int128 sxx = 0;
  __int128 sxy = 0;
  for (const VoteSlopePair& pair : pairs) {
    const std::int64_t x = pair.q - line.q0;
    const std::int64_t y = pair.d - line.d0;
    line.sx += x;
    line.sy += y;
    sxx += static_cast<__int128>(x) * x;
    sxy += static_cast<__int128>(x) * y;
  }
  // n^2 times the variance of q and the covariance: both exact.
  const __int128 den =
      static_cast<__int128>(n) * sxx - static_cast<__int128>(line.sx) * line.sx;
  if (den <= 0) return line;  // every q equal: no slope
  const __int128 num =
      static_cast<__int128>(n) * sxy - static_cast<__int128>(line.sx) * line.sy;
  // round half up: floor((2 num 2^20 + den) / (2 den)).
  __int128 slope = vote_slope_floor_div<__int128>(
      num * (__int128{2} << kVoteSlopeFractionBits) + den, 2 * den);
  const __int128 limit = std::numeric_limits<std::int32_t>::max();
  slope = std::max(-limit, std::min(limit, slope));
  line.b_q20 = static_cast<std::int64_t>(slope);
  line.valid = true;
  return line;
}

// The d range [lo, hi] a seed's posting must lie in, and the centre it is ranked by.
struct VoteSlopeWindow {
  std::int64_t centre = 0;
  std::int64_t lo = 0;
  std::int64_t hi = 0;
};

// One pair per seed of `seeds` (DnaLongSeedView or ChainWindowRetainedSeed): the
// strand-compatible posting on contig `chr` whose d lies in window_of(q), nearest the
// window's centre, the first on ties.
template <class SeedViews, class WindowOf>
inline void vote_slope_collect(const SeedViews& seeds,
                               const std::uint64_t* chr_bounds, int chr,
                               bool is_rc, WindowOf window_of,
                               std::vector<VoteSlopePair>& out) {
  out.clear();
  const std::uint64_t chr_lo = chr_bounds[static_cast<std::size_t>(chr)];
  const std::uint64_t chr_hi = chr_bounds[static_cast<std::size_t>(chr + 1)];
  for (const auto& item : seeds) {
    const QuerySeed& seed = item.seed;
    const KmerPostingView& view = item.view;
    if (!view.found() || view.count == 0) continue;
    const VoteSlopeWindow window = window_of(seed.read_pos);
    bool found = false;
    std::int64_t best_d = 0;
    std::int64_t best_cost = 0;
    for (std::uint32_t i = 0; i < view.count; ++i) {
      const std::uint64_t g = static_cast<std::uint64_t>(view.positions[i]);
      if (g < chr_lo) continue;
      // Postings are stored ascending, so d only grows from here.
      if (g >= chr_hi) break;
      const std::int64_t d =
          static_cast<std::int64_t>(g - chr_lo) - seed.read_pos;
      if (d > window.hi) break;
      if (d < window.lo) continue;
      // The vote's strand test: a posting counts only in its lane.
      if (!packed_ref_orientation_compatible(
              seed.z, view.positions.packed_at(i), is_rc))
        continue;
      const std::int64_t cost = std::llabs(d - window.centre);
      if (!found || cost < best_cost) {
        found = true;
        best_cost = cost;
        best_d = d;
      }
    }
    if (found) out.push_back({seed.read_pos, best_d});
  }
}

// The clamp on b1, the gate's pair minimum, and the floors of S1's reach and the trim.
inline constexpr std::int64_t kVoteSlopeMinQ20 =
    -(std::int64_t{6} << kVoteSlopeFractionBits) / 100;  // -6 %
inline constexpr std::int64_t kVoteSlopeMaxQ20 =
    (std::int64_t{10} << kVoteSlopeFractionBits) / 100;  // +10 %
inline constexpr int kVoteSlopeMinPairs = 8;
inline constexpr std::int64_t kVoteSlopeExtendFloor = 256;
inline constexpr std::int64_t kVoteSlopeTrimFloor = 64;

struct VoteSlopeFit {
  int n_s1 = 0;     // pairs line 1 rests on (S1 after its trim)
  int span_s1 = 0;  // their query span, max q - min q
  std::int64_t b1_q20 = 0;  // clamped to [kVoteSlopeMinQ20, kVoteSlopeMaxQ20]
  std::int64_t a = 0;       // line 1's d at q = 0 under the clamped b1
  bool gate = false;
};

// Fits the line of `winner` on `seeds`, its lane's exact-refine seeds. `pairs` and
// `scratch` are caller-owned buffers.
template <class SeedViews>
inline VoteSlopeFit vote_slope_fit_winner(const SeedViews& seeds,
                                          const std::uint64_t* chr_bounds,
                                          const VotePeak& winner,
                                          int read_len,
                                          std::vector<VoteSlopePair>& pairs,
                                          std::vector<std::int64_t>& scratch) {
  VoteSlopeFit fit;
  const std::int64_t w = std::max(1, winner.anchor.ref_start_bin_width);
  const std::int64_t bin = winner.anchor.ref_start_bin;
  // S0: +-1 bin of the winner's, d in [(bin - 1) W, (bin + 2) W).
  const VoteSlopeWindow fold{winner.raw_ref_start, (bin - 1) * w,
                             (bin + 2) * w - 1};
  vote_slope_collect(seeds, chr_bounds, winner.chr, winner.is_rc,
                     [&fold](int) { return fold; }, pairs);
  const VoteSlopeLine line0 = vote_slope_fit(pairs);
  if (!line0.valid) return fit;
  // S1: the same seeds along line 0.
  const std::int64_t reach = std::max(w, kVoteSlopeExtendFloor);
  vote_slope_collect(seeds, chr_bounds, winner.chr, winner.is_rc,
                     [&line0, reach](int q) {
                       const std::int64_t centre = line0.at(q);
                       return VoteSlopeWindow{centre, centre - reach,
                                              centre + reach};
                     },
                     pairs);
  const VoteSlopeLine extended = vote_slope_fit(pairs);
  if (!extended.valid) return fit;
  // One trim at max(3 MAD, 64) about the extended fit. A median is the upper middle
  // element, as in the vote.
  scratch.clear();
  for (const VoteSlopePair& pair : pairs)
    scratch.push_back(pair.d - extended.at(pair.q));
  const std::size_t middle = scratch.size() / 2;
  std::nth_element(scratch.begin(),
                   scratch.begin() + static_cast<std::ptrdiff_t>(middle),
                   scratch.end());
  const std::int64_t median = scratch[middle];
  for (std::int64_t& value : scratch)
    value = std::llabs(value - median);
  std::nth_element(scratch.begin(),
                   scratch.begin() + static_cast<std::ptrdiff_t>(middle),
                   scratch.end());
  const std::int64_t trim = std::max(3 * scratch[middle], kVoteSlopeTrimFloor);
  std::size_t kept = 0;
  for (std::size_t i = 0; i < pairs.size(); ++i)
    if (std::llabs(pairs[i].d - extended.at(pairs[i].q)) <= trim)
      pairs[kept++] = pairs[i];
  pairs.resize(kept);
  const VoteSlopeLine line1 = vote_slope_fit(pairs);
  fit.n_s1 = static_cast<int>(pairs.size());
  if (!line1.valid) return fit;
  int q_lo = pairs.front().q;
  int q_hi = pairs.front().q;
  for (const VoteSlopePair& pair : pairs) {
    q_lo = std::min(q_lo, pair.q);
    q_hi = std::max(q_hi, pair.q);
  }
  fit.span_s1 = q_hi - q_lo;
  fit.b1_q20 =
      std::max(kVoteSlopeMinQ20, std::min(kVoteSlopeMaxQ20, line1.b_q20));
  VoteSlopeLine clamped = line1;
  clamped.b_q20 = fit.b1_q20;
  fit.a = clamped.at(0);
  fit.gate = fit.n_s1 >= kVoteSlopeMinPairs &&
             4 * static_cast<std::int64_t>(fit.span_s1) >= read_len;
  return fit;
}

}}}  // namespace fa::cpu::lr
