// Groups RNA diagonal peaks into candidate locus envelopes. An RNA read spans several
// exons, each on its own diagonal, so peaks are grouped rather than collapsed to one
// placement. This stage reads no reference bases. Envelopes are ranked by their seed
// evidence; the colinear score is computed for the kept loci only and is the coarse
// confidence the RNA MAPQ reads.

#include "coarse_chain.h"
#include "coarse_transition.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// Derives the internal score units from the seed span; only the intron bounds come from
// the user.
CoarseLocusOptions make_coarse_locus_options(uint32_t seed_span,
                                             uint32_t min_intron,
                                             uint32_t max_intron) {
  CoarseLocusOptions o;
  const uint32_t s = seed_span == 0 ? 15u : seed_span;
  o.diagonal_band = std::max<uint32_t>(16, s + s / 2); // ~1.5 seeds of drift
  o.max_query_gap = std::max<uint32_t>(64, s * 8);
  o.max_query_overlap = std::max<uint32_t>(16, s * 2);
  o.coalesce_min_piece_span = std::max<uint32_t>(16, 2 * s);
  o.min_intron = min_intron == 0 ? 20u : min_intron;
  o.max_intron = max_intron == 0 ? 200000u : max_intron;
  // Coalescing reach and ceiling in units of the maximum intron (see CoarseLocusOptions);
  // saturating, so an extreme -G cannot wrap them.
  const uint64_t mi = o.max_intron;
  o.coalesce_gap =
      static_cast<uint32_t>(std::min<uint64_t>(4u * mi, 0xffffffffull));
  o.max_coalesced_span =
      static_cast<uint32_t>(std::min<uint64_t>(16u * mi, 0xffffffffull));
  o.splice_open = static_cast<int32_t>(std::max<uint32_t>(4, s / 2));
  o.splice_log_scale = 2;
  o.query_gap_scale = 1;
  o.diagonal_drift_scale = 1;
  o.overlap_scale = 1;
  return o;
}

namespace {

// Union of the peaks' oriented query intervals, in bases. Only a tie-break: as a primary
// key it would reward a paralogous envelope for accumulating scattered peaks.
int64_t
union_query_coverage(const std::vector<CoarseDiagonalPeak>& peaks,
                     const std::vector<uint32_t>& indices,
                     std::vector<std::pair<uint32_t, uint32_t>>& spans) {
  spans.clear();
  spans.reserve(indices.size());
  for (uint32_t i : indices) {
    const CoarseDiagonalPeak& p = peaks[i];
    if (p.oriented_query_end > p.oriented_query_begin)
      spans.push_back({p.oriented_query_begin, p.oriented_query_end});
  }
  if (spans.empty())
    return 0;
  std::sort(spans.begin(), spans.end());
  int64_t total = 0;
  uint32_t lo = spans.front().first;
  uint32_t hi = spans.front().second;
  for (size_t i = 1; i < spans.size(); ++i) {
    if (spans[i].first > hi) {
      total += static_cast<int64_t>(hi) - lo;
      lo = spans[i].first;
      hi = spans[i].second;
    } else if (spans[i].second > hi) {
      hi = spans[i].second;
    }
  }
  return total + (static_cast<int64_t>(hi) - lo);
}

// How much of the read the envelope explains with seed matches: peaks in query order,
// each credited the smaller of its own seed coverage and the query it newly explains, so
// no query base counts twice. Unlike best_colinear_score's span increments, a sparse
// peak straddling a wide query range earns only the bases it matched. `indices` may be
// in any order; `scratch` is caller-owned.
int64_t seed_evidence_score(const std::vector<CoarseDiagonalPeak>& peaks,
                            const std::vector<uint32_t>& indices,
                            std::vector<uint32_t>& scratch) {
  scratch.assign(indices.begin(), indices.end());
  std::stable_sort(scratch.begin(), scratch.end(), [&](uint32_t a, uint32_t b) {
    const CoarseDiagonalPeak& pa = peaks[a];
    const CoarseDiagonalPeak& pb = peaks[b];
    if (pa.oriented_query_begin != pb.oriented_query_begin)
      return pa.oriented_query_begin < pb.oriented_query_begin;
    return pa.oriented_query_end < pb.oriented_query_end;
  });
  int64_t total = 0;
  int64_t reach = -1;
  for (uint32_t i : scratch) {
    const CoarseDiagonalPeak& p = peaks[i];
    const int64_t begin =
        std::max<int64_t>(reach, static_cast<int64_t>(p.oriented_query_begin));
    const int64_t novel = static_cast<int64_t>(p.oriented_query_end) - begin;
    if (novel <= 0)
      continue;
    total += std::min<int64_t>(p.query_covered_bases, novel);
    reach =
        std::max<int64_t>(reach, static_cast<int64_t>(p.oriented_query_end));
  }
  return total;
}

// Best colinear chain score over one complete envelope's peaks. It ranks envelopes
// (colinear evidence separates the true locus from a paralog of scattered peaks) but
// never bounds them. `indices` must be sorted by reference_begin.
int64_t best_colinear_score(const std::vector<CoarseDiagonalPeak>& peaks,
                            const std::vector<uint32_t>& indices,
                            const CoarseLocusOptions& o,
                            std::vector<int64_t>& score) {
  const size_t m = indices.size();
  if (m == 0)
    return 0;
  score.assign(m, 0);
  int64_t best = 0;
  for (size_t i = 0; i < m; ++i) {
    const CoarseDiagonalPeak& pi = peaks[indices[i]];
    score[i] = static_cast<int64_t>(pi.query_covered_bases);
    uint32_t scanned = 0;
    for (size_t j = i; j-- > 0 && scanned < o.max_predecessors;) {
      ++scanned;
      const CoarseDiagonalPeak& pj = peaks[indices[j]];
      const detail::Transition t = detail::score_transition(pj, pi, o);
      if (t.kind == detail::TransitionKind::kInvalid)
        continue;
      const int64_t candidate =
          score[j] + detail::incremental_query_reward(pj, pi) - t.penalty;
      if (candidate > score[i])
        score[i] = candidate;
    }
    if (score[i] > best)
      best = score[i];
  }
  return best;
}

// Reference-ordered count of transitions that look like an intron; context only, the
// realized junction count is the splice controller's.
uint32_t estimate_intron_count(const std::vector<CoarseDiagonalPeak>& peaks,
                               const std::vector<uint32_t>& indices,
                               const CoarseLocusOptions& o,
                               std::vector<uint32_t>& by_reference) {
  by_reference.assign(indices.begin(), indices.end());
  std::stable_sort(by_reference.begin(), by_reference.end(),
                   [&](uint32_t a, uint32_t b) {
                     return peaks[a].reference_begin < peaks[b].reference_begin;
                   });
  uint32_t introns = 0;
  for (size_t i = 1; i < by_reference.size(); ++i) {
    const CoarseDiagonalPeak& left = peaks[by_reference[i - 1]];
    const CoarseDiagonalPeak& right = peaks[by_reference[i]];
    if (right.reference_begin <= left.reference_end)
      continue;
    const int64_t dr = static_cast<int64_t>(right.reference_begin) -
                       static_cast<int64_t>(left.reference_end);
    const int64_t dq = static_cast<int64_t>(right.oriented_query_begin) -
                       static_cast<int64_t>(left.oriented_query_end);
    const int64_t intron = dr - dq;
    if (intron >= static_cast<int64_t>(o.min_intron) &&
        intron <= static_cast<int64_t>(o.max_intron)) {
      ++introns;
    }
  }
  return introns;
}

} // namespace

std::vector<CoarseLocus>
select_coarse_loci(const std::vector<CoarseDiagonalPeak>& peaks,
                   const CoarseLocusOptions& o, CoarseLocusScratch& scratch,
                   uint32_t tie_seed, size_t* coalesced,
                   std::vector<CoarseLocus>* overflow_out) {
  std::vector<CoarseLocus> out;
  if (overflow_out)
    overflow_out->clear();
  // Zeroed up front so the early returns report no coalescing.
  if (coalesced)
    *coalesced = 0;
  const size_t n = peaks.size();
  if (n == 0)
    return out;

  std::vector<uint32_t>& order = scratch.order;
  order.assign(n, 0);
  for (size_t i = 0; i < n; ++i)
    order[i] = static_cast<uint32_t>(i);
  std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
    const CoarseDiagonalPeak& pa = peaks[a];
    const CoarseDiagonalPeak& pb = peaks[b];
    if (pa.reference_id != pb.reference_id)
      return pa.reference_id < pb.reference_id;
    if (pa.reverse != pb.reverse)
      return pa.reverse < pb.reverse;
    if (pa.reference_begin != pb.reference_begin)
      return pa.reference_begin < pb.reference_begin;
    if (pa.reference_end != pb.reference_end)
      return pa.reference_end < pb.reference_end;
    return a < b; // original index: total order
  });

  // One linear sweep: a new envelope starts whenever the strand/contig changes
  // or the reference gap exceeds one legal intron.
  std::vector<CoarseLocus>& envelopes = scratch.envelopes;
  std::vector<uint64_t>& support = scratch.support;
  std::vector<int64_t>& coverage = scratch.coverage;
  std::vector<uint32_t>& query_order = scratch.query_order;
  envelopes.clear();
  support.clear();
  coverage.clear();
  size_t i = 0;
  while (i < n) {
    const CoarseDiagonalPeak& first = peaks[order[i]];
    if (first.reference_id < 0) {
      ++i;
      continue;
    }
    CoarseLocus locus;
    locus.reference_id = first.reference_id;
    locus.reverse = first.reverse;
    locus.reference_begin = first.reference_begin;
    locus.reference_end = first.reference_end;
    locus.oriented_query_begin = first.oriented_query_begin;
    locus.oriented_query_end = first.oriented_query_end;
    locus.peak_indices.push_back(order[i]);
    uint64_t locus_support = first.distinct_seed_support;
    size_t j = i + 1;
    for (; j < n; ++j) {
      const CoarseDiagonalPeak& p = peaks[order[j]];
      if (p.reference_id != locus.reference_id || p.reverse != locus.reverse) {
        break;
      }
      if (p.reference_begin >
          locus.reference_end + static_cast<uint64_t>(o.max_intron)) {
        break;
      }
      locus.reference_end = std::max(locus.reference_end, p.reference_end);
      locus.oriented_query_begin =
          std::min(locus.oriented_query_begin, p.oriented_query_begin);
      locus.oriented_query_end =
          std::max(locus.oriented_query_end, p.oriented_query_end);
      locus.peak_indices.push_back(order[j]);
      locus_support += p.distinct_seed_support;
    }
    // Rank by seed evidence; coverage and support are tie-breaks. `score` and
    // `intron_count` are not ranking terms and are computed for the kept loci only.
    locus.rank_score =
        seed_evidence_score(peaks, locus.peak_indices, query_order);
    locus.union_query_coverage_bases =
        union_query_coverage(peaks, locus.peak_indices, scratch.spans);
    coverage.push_back(locus.union_query_coverage_bases);
    envelopes.push_back(std::move(locus));
    support.push_back(locus_support);
    i = j;
  }
  if (envelopes.empty())
    return out;

  // Second pass: coalesce envelopes that are pieces of one transcript. The sweep splits
  // an envelope when consecutive small exons each fall below min_support and the next
  // peak lies beyond max_intron; the trailing piece would be spent as a rival and its
  // junctions fall outside the harvest window. A pair merges only if it is colinear in
  // the read (the later piece starts after the earlier one ends in oriented query
  // coordinates); a paralog explains the same read interval and never merges. Each
  // envelope is offered to every surviving envelope of its group, so a spurious upstream
  // envelope neither blocks nor captures real pieces; the qualifying target reaching
  // furthest right wins.
  size_t coalesced_count = 0;
  if (o.coalesce_gap != 0 && envelopes.size() > 1) {
    std::vector<uint8_t>& absorbed = scratch.coalesce_absorbed;
    absorbed.assign(envelopes.size(), 0); // every slot starts clean
    size_t write = 0;                     // next surviving slot
    size_t group_begin = 0; // first surviving slot of the current group
    // A target's evidence terms are re-derived once, when its group closes: both depend
    // on the whole peak set.
    auto close_group = [&]() {
      for (size_t t = group_begin; t < write; ++t) {
        if (!absorbed[t])
          continue;
        CoarseLocus& m = envelopes[t];
        m.rank_score = seed_evidence_score(peaks, m.peak_indices, query_order);
        m.union_query_coverage_bases =
            union_query_coverage(peaks, m.peak_indices, scratch.spans);
        coverage[t] = m.union_query_coverage_bases;
      }
    };
    for (size_t e = 0; e < envelopes.size(); ++e) {
      const CoarseLocus& c = envelopes[e];
      if (write > group_begin &&
          (envelopes[group_begin].reference_id != c.reference_id ||
           envelopes[group_begin].reverse != c.reverse)) {
        close_group();
        group_begin = write;
      }
      // The best earlier target: of this group's surviving envelopes that could
      // absorb `c`, the one whose reference_end is largest.
      size_t best = write; // == write means "none qualified"
      for (size_t t = group_begin; t < write; ++t) {
        const CoarseLocus& m = envelopes[t];
        // A group is strictly ordered and disjoint after the sweep; the test also keeps
        // the unsigned subtraction safe.
        if (c.reference_begin <= m.reference_end ||
            c.reference_begin - m.reference_end >
                static_cast<uint64_t>(o.coalesce_gap)) {
          continue;
        }
        // Both pieces must be exon-shaped (see coalesce_min_piece_span).
        if (m.oriented_query_end - m.oriented_query_begin <
                o.coalesce_min_piece_span ||
            c.oriented_query_end - c.oriented_query_begin <
                o.coalesce_min_piece_span) {
          continue;
        }
        if (static_cast<int64_t>(m.oriented_query_end) >
            static_cast<int64_t>(c.oriented_query_begin) +
                static_cast<int64_t>(o.max_query_overlap)) {
          continue; // not colinear in the read
        }
        const uint64_t span =
            std::max(m.reference_end, c.reference_end) - m.reference_begin;
        if (span > static_cast<uint64_t>(o.max_coalesced_span))
          continue;
        // Evidence gate: two scattered repeat peaks of a short junk read can land
        // colinear by chance, so one credible piece is required.
        if (std::max(support[t], support[e]) <
            static_cast<uint64_t>(o.coalesce_min_support)) {
          continue;
        }
        if (best == write || envelopes[best].reference_end < m.reference_end)
          best = t;
      }
      if (best != write) {
        CoarseLocus& m = envelopes[best];
        m.reference_end = std::max(m.reference_end, c.reference_end);
        m.oriented_query_begin =
            std::min(m.oriented_query_begin, c.oriented_query_begin);
        m.oriented_query_end =
            std::max(m.oriented_query_end, c.oriented_query_end);
        // Every peak of `c` begins past the target's reference end, so appending keeps
        // peak_indices reference-ordered, as best_colinear_score requires.
        m.peak_indices.insert(m.peak_indices.end(), c.peak_indices.begin(),
                              c.peak_indices.end());
        support[best] += support[e];
        absorbed[best] = 1;
        ++coalesced_count;
        continue;
      }
      // `e` survives as a target for what follows; compact it over the absorbed slots
      // together with its parallel entries.
      if (write != e) {
        envelopes[write] = std::move(envelopes[e]);
        support[write] = support[e];
        coverage[write] = coverage[e];
      }
      ++write;
    }
    close_group();
    envelopes.resize(write);
    support.resize(write);
    coverage.resize(write);
  }
  if (coalesced)
    *coalesced = coalesced_count;

  std::vector<uint32_t>& ranked = scratch.ranked;
  ranked.assign(envelopes.size(), 0);
  for (size_t e = 0; e < envelopes.size(); ++e)
    ranked[e] = static_cast<uint32_t>(e);
  std::stable_sort(ranked.begin(), ranked.end(), [&](uint32_t a, uint32_t b) {
    if (envelopes[a].rank_score != envelopes[b].rank_score)
      return envelopes[a].rank_score > envelopes[b].rank_score;
    if (coverage[a] != coverage[b])
      return coverage[a] > coverage[b];
    if (support[a] != support[b])
      return support[a] > support[b];
    if (envelopes[a].peak_indices.size() != envelopes[b].peak_indices.size())
      return envelopes[a].peak_indices.size() >
             envelopes[b].peak_indices.size();
    // Equal on every evidence key: minimap2's read-seeded locus hash decides rather
    // than contig order; contig slot and start settle only a 64-bit collision.
    const uint64_t a_hash = coarse_locus_tie_hash(tie_seed, envelopes[a]);
    const uint64_t b_hash = coarse_locus_tie_hash(tie_seed, envelopes[b]);
    if (a_hash != b_hash)
      return a_hash < b_hash;
    if (envelopes[a].reference_id != envelopes[b].reference_id)
      return envelopes[a].reference_id < envelopes[b].reference_id;
    return envelopes[a].reference_begin < envelopes[b].reference_begin;
  });

  const size_t keep = std::min<size_t>(o.max_locus_chains, ranked.size());
  // With overflow requested, up to `admit` loci are taken and everything past `keep`
  // goes to `overflow_out` unfinished.
  const size_t admit =
      overflow_out ? std::min<size_t>(keep + o.max_overflow_loci, ranked.size())
                   : keep;
  out.reserve(keep);
  // Walk in rank order, keeping at most one envelope per distinct reference window: a
  // candidate overlapping a kept locus on the same contig, on either strand, names the
  // same place and is dropped, and the next ranked entry takes the slot. Ignoring strand
  // removes the winner's own opposite-strand mirror; same-strand envelopes are disjoint
  // by construction.
  for (size_t r = 0; r < ranked.size() &&
                     out.size() + (overflow_out ? overflow_out->size() : 0) <
                         admit;
       ++r) {
    const CoarseLocus& candidate = envelopes[ranked[r]];
    bool duplicate_window = false;
    for (const CoarseLocus& kept : out) {
      if (kept.reference_id != candidate.reference_id)
        continue;
      // Half-open intersection, strand ignored: the mirror is the same window.
      if (candidate.reference_begin < kept.reference_end &&
          kept.reference_begin < candidate.reference_end) {
        duplicate_window = true;
        break;
      }
    }
    // Overflow loci must be distinct places too.
    if (!duplicate_window && overflow_out) {
      for (const CoarseLocus& kept : *overflow_out) {
        if (kept.reference_id != candidate.reference_id)
          continue;
        if (candidate.reference_begin < kept.reference_end &&
            kept.reference_begin < candidate.reference_end) {
          duplicate_window = true;
          break;
        }
      }
    }
    if (duplicate_window)
      continue;
    CoarseLocus locus = std::move(envelopes[ranked[r]]);
    if (out.size() >= keep) {
      overflow_out->push_back(std::move(locus));
      continue;
    }
    // Non-ranking statistics, for the kept loci only (best_colinear_score is an
    // O(m * max_predecessors) DP). Both need the sweep's reference-ordered peak_indices,
    // so they run before the query-order sort. `score` is the coarse confidence the RNA
    // MAPQ reads.
    locus.score =
        best_colinear_score(peaks, locus.peak_indices, o, scratch.chain_score);
    locus.intron_count = estimate_intron_count(peaks, locus.peak_indices, o,
                                               scratch.by_reference);
    // Query order is what the skeleton and the downstream anchor priority read.
    std::stable_sort(locus.peak_indices.begin(), locus.peak_indices.end(),
                     [&](uint32_t a, uint32_t b) {
                       const CoarseDiagonalPeak& pa = peaks[a];
                       const CoarseDiagonalPeak& pb = peaks[b];
                       if (pa.oriented_query_begin != pb.oriented_query_begin)
                         return pa.oriented_query_begin <
                                pb.oriented_query_begin;
                       if (pa.reference_begin != pb.reference_begin)
                         return pa.reference_begin < pb.reference_begin;
                       return a < b;
                     });
    out.push_back(std::move(locus));
  }
  return out;
}

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
