#include "rna_parent.h"

#include <algorithm>
#include <utility>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

void rna_parent_set(std::vector<RnaParentHit>& r, float mask_level,
                    int sub_diff) {
  const int n = static_cast<int>(r.size());
  if (n <= 0)
    return;
  std::vector<std::pair<int, int>> cov;
  cov.reserve(static_cast<std::size_t>(n));
  std::vector<int> w(static_cast<std::size_t>(n));
  w[0] = 0;
  r[0].parent = 0;
  int k = 1;
  for (int i = 1; i < n; ++i) {
    RnaParentHit& ri = r[i];
    const int si = ri.qs;
    const int ei = ri.qe;
    cov.clear();
    for (int j = 0; j < k; ++j) {
      const RnaParentHit& rp = r[w[j]];
      int sj = rp.qs;
      int ej = rp.qe;
      if (ej <= si || sj >= ei)
        continue;
      if (sj < si)
        sj = si;
      if (ej > ei)
        ej = ei;
      cov.emplace_back(sj, ej);
    }
    int j = k;
    if (!cov.empty()) {
      // The query length no existing primary covers.
      int uncov_len = 0;
      std::sort(cov.begin(), cov.end());
      int x = si;
      for (const auto& [cs, ce] : cov) {
        if (cs > x)
          uncov_len += cs - x;
        x = ce > x ? ce : x;
      }
      if (ei > x)
        uncov_len += ei - x;
      for (j = 0; j < k; ++j) {
        RnaParentHit& rp = r[w[j]];
        const int sj = rp.qs;
        const int ej = rp.qe;
        if (ej <= si || sj >= ei)
          continue;
        const int min = ej - sj < ei - si ? ej - sj : ei - si;
        const int max = ej - sj > ei - si ? ej - sj : ei - si;
        const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                               : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
        // minimap2 divides by min unguarded; a zero min gives NaN there and fails the
        // test, as the guard does here.
        if (min > 0 &&
            static_cast<float>(ol) / static_cast<float>(min) -
                    static_cast<float>(uncov_len) / static_cast<float>(max) >
                mask_level) {
          bool cnt_sub = false;
          ri.parent = rp.parent;
          rp.subsc = rp.subsc > ri.score ? rp.subsc : ri.score;
          if (ri.cnt >= rp.cnt)
            cnt_sub = true;
          // The last condition excludes identical hits after DP.
          if (rp.has_dp && ri.has_dp &&
              (rp.rid != ri.rid || rp.rs != ri.rs || rp.re != ri.re ||
               ol != min)) {
            rp.dp_max2 = rp.dp_max2 > ri.dp_max ? rp.dp_max2 : ri.dp_max;
            if (rp.dp_max - ri.dp_max <= sub_diff)
              cnt_sub = true;
          }
          if (cnt_sub)
            ++rp.n_sub;
          break;
        }
      }
    }
    if (j == k) {
      w[static_cast<std::size_t>(k++)] = i;
      ri.parent = i;
      ri.n_sub = 0;
    }
  }
}

bool rna_parent_update_dp_max(int qlen, int min_len,
                              std::vector<RnaParentHit>& r, float frac, int a,
                              int b) {
  if (qlen < min_len || r.size() < 2)
    return false;
  int max = -1;
  int max2 = -1;
  int max_i = -1;
  for (std::size_t i = 0; i < r.size(); ++i) {
    if (!r[i].has_dp)
      continue;
    if (r[i].dp_max > max) {
      max2 = max;
      max = r[i].dp_max;
      max_i = static_cast<int>(i);
    } else if (r[i].dp_max > max2) {
      max2 = r[i].dp_max;
    }
  }
  if (max_i < 0 || max < 0 || max2 < 0)
    return false;
  const RnaParentHit& top = r[static_cast<std::size_t>(max_i)];
  if (top.qe - top.qs < static_cast<double>(qlen) * frac)
    return false;
  if (max2 < static_cast<double>(max) * frac)
    return false;
  for (const RnaParentHit& h : r)
    if (h.has_dp && !h.accounting_ok)
      return false;
  const double identity = rna_event_identity(top.accounting);
  if (identity < 0.0)
    return false;
  const double b2 = rna_rank_recal_b2(identity, a, b);
  for (RnaParentHit& h : r) {
    if (!h.has_dp)
      continue;
    const int dp = rna_recal_max_dp(h.accounting, b2, a);
    h.dp_max = dp < 0 ? 0 : dp;
  }
  return true;
}

void rna_parent_hit_sort(std::vector<RnaParentHit>& r) {
  std::stable_sort(r.begin(), r.end(),
                   [](const RnaParentHit& x, const RnaParentHit& y) {
                     const int kx = x.has_dp ? x.dp_max : x.score;
                     const int ky = y.has_dp ? y.dp_max : y.score;
                     if (kx != ky)
                       return kx > ky;
                     if (x.catalogue != y.catalogue)
                       return x.catalogue < y.catalogue;
                     return x.segment < y.segment;
                   });
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
