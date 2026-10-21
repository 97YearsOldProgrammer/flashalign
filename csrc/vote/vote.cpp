// The DNA whole-read vote over one window: postings accumulated into diagonal bins, ranked, refined by re-walking the best bins, emitted as peaks; the candidate catalogue those peaks form.
#include "vote.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

/**** accumulate ****/
// The DNA vote's accumulate stage: one distinct-query vote per seed into each reference
// start (diagonal) bin, into scratch.buckets, with the used-seed bookkeeping.

namespace fa {
namespace cpu {
namespace lr {

// Accumulates distinct-query diagonal votes into scratch.buckets.
static void
vote_accumulate(const VoteWindowState &st, ChainWindowPeakScratch &scratch,
                const std::vector<QuerySeed> &seeds,
                const std::vector<DnaLongSeedView> *seed_view_override,
                bool retain_used_views, ChainSeedLookupCache *lookup_cache) {
  const LongReadSeedContext &ctx = st.ctx;
  auto &seed_views = scratch.seed_views;
  auto &seed_used = scratch.seed_used;
  auto &retained_seeds = scratch.retained_seeds;
  auto &buckets = scratch.buckets;
  uint32_t seed_epoch = 0; // per-window, ++ once per voting seed (stamp)
  // A rescued seed (DnaLongSeedView::rescued) passes the occurrence cap.
  auto visit_seed_view_for_vote = [&](const QuerySeed &seed,
                                      const KmerPostingView &v, bool rescued) {
    if (!retain_used_views) {
      seed_views.push_back(v);
    }
    uint8_t used = 0;
    if (!v.found() || v.count == 0) {
      if (!retain_used_views)
        seed_used.push_back(used);
      return;
    }
    if (!rescued && !vote_seed_occ_allowed(st, v)) {
      // Rejected by the occurrence cap.
      if (!retain_used_views)
        seed_used.push_back(used);
      return;
    }
    used = 1;
    if (!retain_used_views) {
      seed_used.push_back(used);
    } else {
      retained_seeds.push_back({seed, v, rescued});
    }
    int chr_idx = chromosome_index_for_global_pos(st.chr_bounds, st.n_chr,
                                                  v.positions[0]);
    if (chr_idx < 0 || chr_idx >= st.n_chr)
      return;
    uint64_t chr_lo = st.chr_bounds[static_cast<size_t>(chr_idx)];
    uint64_t chr_hi = st.chr_bounds[static_cast<size_t>(chr_idx + 1)];
    const int seed_weight = 1; // unweighted vote: one unit per distinct query seed
    // Stamp the bucket with this seed's epoch so only the seed's first posting in a diagonal
    // adds support. Postings ascend, so a seed's hits in one bin are contiguous and the
    // bucket pointer is re-probed only when the bin changes.
    const uint32_t epoch = ++seed_epoch;
    int64_t cached_key = std::numeric_limits<int64_t>::min();
    ChainWindowBucketVote *cached = nullptr;
    for (uint32_t i = 0; i < v.count; ++i) {
      const uint64_t g = static_cast<uint64_t>(v.positions[i]);
      advance_contig_cursor(st.chr_bounds, st.n_chr, g, chr_idx, chr_lo,
                            chr_hi);
      if (g < chr_lo || g >= chr_hi)
        continue;
      // A posting votes only in the lane it is compatible with: the oriented query k-mer
      // must be the reference k-mer there. Tested after the cursor advance, which every
      // posting must make.
      if (!packed_ref_orientation_compatible(seed.z, v.positions.packed_at(i),
                                             st.is_rc))
        continue;
      const int local = static_cast<int>(g - chr_lo);
      // The read's own k-mer casts no vote, as minimap2's -D drops the anchors on the
      // exact diagonal; a self hit off that diagonal still votes.
      if (chr_idx == ctx.self_contig && !st.is_rc && local == seed.read_pos)
        continue;
      const int ref_start = local - seed.read_pos;
      const int64_t bucket_key =
          vote_pack_key(chr_idx, vote_floor_div(ref_start, st.W));
      if (bucket_key != cached_key) {
        cached = &buckets[bucket_key];
        cached_key = bucket_key;
      }
      if (cached->voter_epoch != epoch) {
        cached->support++;
        cached->center++;
        cached->score += seed_weight;
        cached->voter_epoch = epoch;
      }
    }
  };
  if (seed_view_override) {
    // Views were resolved upstream, so this loop only accumulates.
    for (const auto &sv : *seed_view_override) {
      visit_seed_view_for_vote(sv.seed, sv.view, sv.rescued);
    }
  } else {
    for (const auto &seed : seeds) {
      const KmerPostingView v = lookup_cache
                                    ? lookup_cache->lookup(*ctx.index, seed.key)
                                    : ctx.index->lookup(seed.key);
      visit_seed_view_for_vote(seed, v, false);
    }
  }
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** rank ****/

namespace fa { namespace cpu { namespace lr {

// The vote's rank stage: folds each bucket of scratch.buckets with its vote_radius
// neighbours into scratch.ranked, drops diagonals below min_support, and heapifies by
// chain_window_ranked_better under the read's tie seed.
static void vote_rank(const VoteWindowState& st,
                      ChainWindowPeakScratch& scratch) {
    const int vote_radius = st.vote_radius;
    const int min_support = st.min_support;
    auto& buckets = scratch.buckets;
    auto& ranked = scratch.ranked;
    ranked.reserve(buckets.size());
    for (const auto& kv : buckets) {
        const int key_chr =
            static_cast<int>(static_cast<uint64_t>(kv.first) >> 32);
        const int key_bin = static_cast<int>(
            static_cast<int32_t>(
                static_cast<uint32_t>(kv.first & 0xFFFFFFFF)));
        int score = kv.second.score;
        int support = kv.second.support;
        for (int delta = -vote_radius; delta <= vote_radius; ++delta) {
            if (delta == 0) continue;
            const auto it = buckets.find(vote_pack_key(key_chr, key_bin + delta));
            if (it != buckets.end()) {
                score += it->second.score;
                support += it->second.support;
            }
        }
        ranked.push_back({key_chr, key_bin, score, support, kv.second.center});
    }
    // Gate before ranking, so a high-scoring bucket below min_support cannot hide a
    // well-supported diagonal.
    ranked.erase(
        std::remove_if(ranked.begin(), ranked.end(),
            [min_support](const ChainWindowRankedBucket& rb) {
                return rb.support < min_support;
            }),
        ranked.end());
    std::make_heap(ranked.begin(), ranked.end(), VoteHeapBelow{&st});
}

}}}  // namespace fa::cpu::lr

/**** emit ****/

namespace fa { namespace cpu { namespace lr {

// Appends cp unless an emitted peak on its contig lies within W * (vote_radius + 1)
// of its raw median diagonal (raw_ref_start, not the clamped ref_pos), so peaks that
// clamp to one start at a contig end stay distinct. Returns true iff appended.
static bool vote_append_peak_if_distinct(const VoteWindowState& st,
                                         std::vector<VotePeak>& out,
                                         VotePeak& cp) {
    const std::int64_t reach =
        static_cast<std::int64_t>(st.W) * (st.vote_radius + 1);
    bool distinct = true;
    for (const auto& prev : out) {
        if (prev.chr != cp.chr) continue;
        if (std::llabs(prev.raw_ref_start - cp.raw_ref_start) <= reach) {
            distinct = false;
            break;
        }
    }
    if (distinct) out.push_back(cp);
    return distinct;
}

}}}  // namespace fa::cpu::lr

/**** refine ****/
// The vote's refine and emit stage (the heap drain): pops diagonals best first, re-walks
// each bucket's used seeds (the exact-refine pass) to recompute the median reference
// start and median occurrence, and emits one VotePeak per bucket.

namespace fa { namespace cpu { namespace lr {

// Batched drain (ctx.vote_batched_refine): the same output as the sequential drain in
// vote_emit_peaks, but the exact-refine walk is paid once per chunk of drained buckets
// instead of once per bucket. Buckets are processed in heap order, so support counts,
// medians, the min_support skip and distinctness replay the sequential path, and a chunk
// never pops more buckets than emissions still needed.
static void vote_emit_peaks_batched(const VoteWindowState& st,
                        ChainWindowPeakScratch& scratch,
                        std::vector<VotePeak>& out,
                        const std::vector<QuerySeed>& seeds,
                        bool retain_used_views,
                        size_t active_seed_count) {
    auto& ranked = scratch.ranked;
    auto& seed_views = scratch.seed_views;
    auto& seed_used = scratch.seed_used;
    auto& retained_seeds = scratch.retained_seeds;
    const std::vector<DnaLongSeedView>* exact_refine_seed_views =
        !scratch.pending_exact_refine_views.empty()
            ? &scratch.pending_exact_refine_views
            : (!scratch.seed_bundle.exact_refine_views.empty()
                   ? &scratch.seed_bundle.exact_refine_views
                   : nullptr);
    const bool evidence_interval_valid =
        exact_refine_seed_views == nullptr || retain_used_views;
    const int coarse_radius = std::max(1, st.vote_radius);
    std::vector<ChainWindowRankedBucket> chunk;
    std::vector<std::vector<int>> chunk_ref_starts;
    std::vector<std::vector<uint32_t>> chunk_occurrences;
    std::vector<int> chunk_seed_count;
    std::vector<int> chunk_evidence_lo;
    std::vector<int> chunk_evidence_hi;
    std::vector<uint8_t> seed_hit;
    std::vector<size_t> seed_hit_entries;
    // The chunk's entries in (chr, bin) order, so a posting visits only the entries of its
    // chromosome within coarse_radius bins of it. Each entry still sees its postings in
    // order, so the emitted peaks match a full scan.
    std::vector<size_t> by_chr_bin;
    // The chunk's contigs in ascending order: flattened bounds and their entries'
    // range in by_chr_bin.
    struct ChunkContig {
        uint64_t lo;
        uint64_t hi;
        int chr;
        size_t first;
        size_t last;
    };
    std::vector<ChunkContig> chunk_contigs;
    while (!ranked.empty() && static_cast<int>(out.size()) < st.limit) {
        chunk.clear();
        const int needed = st.limit - static_cast<int>(out.size());
        while (!ranked.empty() && static_cast<int>(chunk.size()) < needed) {
            // The heap order vote_rank built.
            std::pop_heap(ranked.begin(), ranked.end(), VoteHeapBelow{&st});
            const ChainWindowRankedBucket rb = ranked.back();
            ranked.pop_back();
            if (rb.support < st.min_support) continue;
            chunk.push_back(rb);
        }
        if (chunk.empty()) break;
        const size_t chunk_size = chunk.size();
        chunk_ref_starts.assign(chunk_size, {});
        chunk_occurrences.assign(chunk_size, {});
        chunk_seed_count.assign(chunk_size, 0);
        chunk_evidence_lo.assign(chunk_size, std::numeric_limits<int>::max());
        chunk_evidence_hi.assign(chunk_size, std::numeric_limits<int>::min());
        seed_hit.assign(chunk_size, 0);
        for (size_t c = 0; c < chunk_size; ++c) {
            const size_t support =
                static_cast<size_t>(std::max(0, chunk[c].support));
            chunk_ref_starts[c].reserve(support);
            chunk_occurrences[c].reserve(std::min(active_seed_count, support));
        }
        by_chr_bin.resize(chunk_size);
        for (size_t c = 0; c < chunk_size; ++c) by_chr_bin[c] = c;
        std::sort(by_chr_bin.begin(), by_chr_bin.end(),
                  [&chunk](size_t a, size_t b) {
                      if (chunk[a].chr != chunk[b].chr)
                          return chunk[a].chr < chunk[b].chr;
                      if (chunk[a].bin != chunk[b].bin)
                          return chunk[a].bin < chunk[b].bin;
                      return a < b;
                  });
        chunk_contigs.clear();
        for (size_t o = 0; o < chunk_size; ++o) {
            const int chr = chunk[by_chr_bin[o]].chr;
            if (!chunk_contigs.empty() && chunk_contigs.back().chr == chr) {
                chunk_contigs.back().last = o + 1;
                continue;
            }
            chunk_contigs.push_back({st.chr_bounds[static_cast<size_t>(chr)],
                                     st.chr_bounds[static_cast<size_t>(chr + 1)],
                                     chr, o, o + 1});
        }
        // A posting list ascends: k stays on the first chunk contig ending past the
        // postings read so far, and a posting off the chunk's contigs reaches no entry.
        auto visit_seed = [&](const QuerySeed& seed,
                              const KmerPostingView& v,
                              int /*seed_id*/) {
            if (!v.found() || v.count == 0) return;
            size_t k = 0;
            for (uint32_t i = 0; i < v.count; ++i) {
                const uint64_t g = static_cast<uint64_t>(v.positions[i]);
                if (g >= chunk_contigs[k].hi) {
                    k = static_cast<size_t>(
                        std::upper_bound(chunk_contigs.begin() +
                                             static_cast<ptrdiff_t>(k + 1),
                                         chunk_contigs.end(), g,
                                         [](uint64_t x, const ChunkContig& cc) {
                                             return x < cc.hi;
                                         }) -
                        chunk_contigs.begin());
                    if (k == chunk_contigs.size()) break;
                }
                const ChunkContig& cc = chunk_contigs[k];
                if (g < cc.lo) continue;
                // The accumulate stage's strand test, so every statistic here counts
                // compatible postings only.
                if (!packed_ref_orientation_compatible(
                        seed.z, v.positions.packed_at(i), st.is_rc))
                    continue;
                const int local = static_cast<int>(g - cc.lo);
                const int ref_start = local - seed.read_pos;
                const int bin = vote_floor_div(ref_start, st.W);
                // First entry of this chromosome at or past bin - coarse_radius.
                const int bin_lo = bin - coarse_radius;
                size_t lo = cc.first, hi = cc.last;
                while (lo < hi) {
                    const size_t mid = lo + (hi - lo) / 2;
                    if (chunk[by_chr_bin[mid]].bin < bin_lo)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                for (size_t o = lo; o < cc.last; ++o) {
                    const size_t c = by_chr_bin[o];
                    const int bin_distance = std::abs(bin - chunk[c].bin);
                    if (bin - chunk[c].bin < -coarse_radius) break;
                    if (bin_distance <= coarse_radius) {
                        chunk_evidence_lo[c] =
                            std::min(chunk_evidence_lo[c], seed.read_pos);
                        chunk_evidence_hi[c] =
                            std::max(chunk_evidence_hi[c], seed.read_pos);
                    }
                    if (bin_distance <= st.vote_radius) {
                        chunk_ref_starts[c].push_back(ref_start);
                        if (!seed_hit[c]) {
                            seed_hit[c] = 1;
                            seed_hit_entries.push_back(c);
                        }
                    }
                }
            }
            // Unweighted vote: one unit per distinct query seed.
            for (const size_t c : seed_hit_entries) {
                seed_hit[c] = 0;
                ++chunk_seed_count[c];
                chunk_occurrences[c].push_back(v.occurrence);
            }
            seed_hit_entries.clear();
        };
        if (exact_refine_seed_views) {
          for (size_t seed_i = 0; seed_i < exact_refine_seed_views->size();
               ++seed_i) {
            const auto& selected = (*exact_refine_seed_views)[seed_i];
            visit_seed(selected.seed, selected.view, static_cast<int>(seed_i));
          }
        } else if (retain_used_views) {
          for (size_t seed_i = 0; seed_i < retained_seeds.size(); ++seed_i) {
            const auto& retained = retained_seeds[seed_i];
            visit_seed(retained.seed, retained.view, static_cast<int>(seed_i));
          }
        } else {
          for (size_t seed_i = 0; seed_i < seeds.size(); ++seed_i) {
            if (seed_i >= seed_used.size() || !seed_used[seed_i])
              continue;
            visit_seed(seeds[seed_i], seed_views[seed_i],
                       static_cast<int>(seed_i));
          }
        }
        for (size_t c = 0; c < chunk_size; ++c) {
            auto& ref_starts = chunk_ref_starts[c];
            if (ref_starts.empty()) continue;
            const ChainWindowRankedBucket& rb = chunk[c];
            auto& occurrence_values = chunk_occurrences[c];
            uint32_t median_occurrence = 0;
            if (!occurrence_values.empty()) {
                const size_t median_idx = occurrence_values.size() / 2;
                std::nth_element(
                    occurrence_values.begin(),
                    occurrence_values.begin() +
                        static_cast<ptrdiff_t>(median_idx),
                    occurrence_values.end());
                median_occurrence = occurrence_values[median_idx];
            }
            std::nth_element(
                ref_starts.begin(),
                ref_starts.begin() +
                    static_cast<ptrdiff_t>(ref_starts.size() / 2),
                ref_starts.end());
            int64_t ref_pos = ref_starts[ref_starts.size() / 2];
            const int64_t raw_ref_start = ref_pos;
            const int64_t chr_len =
                static_cast<int64_t>(st.chr_bounds[rb.chr + 1]) -
                static_cast<int64_t>(st.chr_bounds[rb.chr]);
            if (ref_pos < 0) ref_pos = 0;
            if (ref_pos > chr_len - static_cast<int64_t>(st.span)) {
                ref_pos = std::max<int64_t>(
                    0, chr_len - static_cast<int64_t>(st.span));
            }
            VotePeak cp;
            cp.read_lo = 0;
            cp.read_hi = st.span;
            cp.ref_pos = static_cast<int>(ref_pos);
            cp.chr = rb.chr;
            cp.is_rc = st.is_rc;
            cp.support = chunk_seed_count[c];
            cp.center_support = std::min(rb.center, chunk_seed_count[c]);
            cp.vote_score = chunk_seed_count[c];
            cp.anchor.ref_start_bin = rb.bin;
            cp.anchor.ref_start_bin_width = st.W;
            cp.anchor.median_occurrence = median_occurrence;
            cp.raw_ref_start = raw_ref_start;
            if (evidence_interval_valid &&
                chunk_evidence_lo[c] <= chunk_evidence_hi[c]) {
                cp.evidence_read_lo = chunk_evidence_lo[c];
                cp.evidence_read_hi = chunk_evidence_hi[c];
            }
            vote_append_peak_if_distinct(st, out, cp);
        }
    }
    std::sort(out.begin(), out.end(), chain_peak_better);
}

// Drains the ranked heap into `out`. active_seed_count sizes the scratch reservations.
static void vote_emit_peaks(const VoteWindowState& st,
                        ChainWindowPeakScratch& scratch,
                        std::vector<VotePeak>& out,
                        const std::vector<QuerySeed>& seeds,
                        bool retain_used_views,
                        size_t active_seed_count) {
    if (st.ctx.vote_batched_refine) {
        vote_emit_peaks_batched(st, scratch, out, seeds, retain_used_views,
                                active_seed_count);
        return;
    }
    auto& ranked = scratch.ranked;
    auto& seed_views = scratch.seed_views;
    auto& seed_used = scratch.seed_used;
    auto& retained_seeds = scratch.retained_seeds;
    const std::vector<DnaLongSeedView>* exact_refine_seed_views =
        !scratch.pending_exact_refine_views.empty()
            ? &scratch.pending_exact_refine_views
            : (!scratch.seed_bundle.exact_refine_views.empty()
                   ? &scratch.seed_bundle.exact_refine_views
                   : nullptr);
    // The per-peak evidence interval bounds the tile mask's seeds only when this drain
    // walks every seed the mask enumerates.
    const bool evidence_interval_valid =
        exact_refine_seed_views == nullptr || retain_used_views;
    // The mask accepts postings within one diagonal bin; track the interval at
    // the wider of that and the vote radius so it bounds the mask's seeds.
    const int coarse_radius = std::max(1, st.vote_radius);
    while (!ranked.empty()) {
        if (static_cast<int>(out.size()) >= st.limit) break;
        // The heap order vote_rank built.
        std::pop_heap(ranked.begin(), ranked.end(), VoteHeapBelow{&st});
        const ChainWindowRankedBucket rb = ranked.back();
        ranked.pop_back();
        if (rb.support < st.min_support) continue;

        auto& ref_starts = scratch.ref_starts;
        ref_starts.clear();
        ref_starts.reserve(static_cast<size_t>(rb.support));
        int contributing_seed_count = 0;
        int contributing_score = 0;
        auto& occurrence_values = scratch.occurrence_values;
        occurrence_values.clear();
        occurrence_values.reserve(active_seed_count);
        int evidence_lo = std::numeric_limits<int>::max();
        int evidence_hi = std::numeric_limits<int>::min();
        auto visit_seed = [&](const QuerySeed& seed,
                              const KmerPostingView& v,
                              int /*seed_id*/) {
            if (!v.found() || v.count == 0) return;
            int chr_idx = chromosome_index_for_global_pos(
                st.chr_bounds, st.n_chr, v.positions[0]);
            if (chr_idx < 0 || chr_idx >= st.n_chr) return;
            uint64_t chr_lo = st.chr_bounds[static_cast<size_t>(chr_idx)];
            uint64_t chr_hi = st.chr_bounds[static_cast<size_t>(chr_idx + 1)];
            bool seed_contributed = false;
            int best_ref_delta = 0;
            const int seed_weight = 1;  // unweighted vote: one unit per distinct query seed
            // Unit weights guarantee the vote winner clears the ratio-admission floor, which
            // makes it catalogue candidate 0; weighted votes would break that.
            const int bin_center = rb.bin * st.W + st.W / 2;
            for (uint32_t i = 0; i < v.count; ++i) {
                const uint64_t g = static_cast<uint64_t>(v.positions[i]);
                advance_contig_cursor(st.chr_bounds, st.n_chr, g, chr_idx,
                                      chr_lo, chr_hi);
                if (g < chr_lo || g >= chr_hi || chr_idx != rb.chr) continue;
                // The accumulate stage's strand test, so every statistic here counts
                // compatible postings only; tested after the cursor advance.
                if (!packed_ref_orientation_compatible(
                        seed.z, v.positions.packed_at(i), st.is_rc))
                    continue;
                const int local = static_cast<int>(g - chr_lo);
                const int ref_start = local - seed.read_pos;
                const int bin_distance =
                    std::abs(vote_floor_div(ref_start, st.W) - rb.bin);
                if (bin_distance <= coarse_radius) {
                    evidence_lo = std::min(evidence_lo, seed.read_pos);
                    evidence_hi = std::max(evidence_hi, seed.read_pos);
                }
                if (bin_distance <= st.vote_radius) {
                    ref_starts.push_back(ref_start);
                    const int ref_delta = std::abs(ref_start - bin_center);
                    if (!seed_contributed || ref_delta < best_ref_delta) {
                        best_ref_delta = ref_delta;
                    }
                    seed_contributed = true;
                }
            }
            if (seed_contributed) {
                ++contributing_seed_count;
                contributing_score += seed_weight;
                occurrence_values.push_back(v.occurrence);
            }
        };
        if (exact_refine_seed_views) {
          for (size_t seed_i = 0; seed_i < exact_refine_seed_views->size();
               ++seed_i) {
            const auto& selected = (*exact_refine_seed_views)[seed_i];
            visit_seed(selected.seed, selected.view, static_cast<int>(seed_i));
          }
        } else if (retain_used_views) {
          for (size_t seed_i = 0; seed_i < retained_seeds.size(); ++seed_i) {
            const auto& retained = retained_seeds[seed_i];
            visit_seed(retained.seed, retained.view, static_cast<int>(seed_i));
          }
        } else {
          for (size_t seed_i = 0; seed_i < seeds.size(); ++seed_i) {
            if (seed_i >= seed_used.size() || !seed_used[seed_i])
              continue;
            visit_seed(seeds[seed_i], seed_views[seed_i],
                       static_cast<int>(seed_i));
          }
        }
        if (ref_starts.empty()) continue;
        uint32_t median_occurrence = 0;
        if (!occurrence_values.empty()) {
            const size_t median_idx = occurrence_values.size() / 2;
            std::nth_element(
                occurrence_values.begin(),
                occurrence_values.begin() + static_cast<ptrdiff_t>(median_idx),
                occurrence_values.end());
            median_occurrence = occurrence_values[median_idx];
        }
        std::nth_element(
            ref_starts.begin(), ref_starts.begin() + static_cast<ptrdiff_t>(ref_starts.size() / 2),
            ref_starts.end());
        int64_t ref_pos = ref_starts[ref_starts.size() / 2];
        const int64_t raw_ref_start = ref_pos;
        const int64_t chr_len =
            static_cast<int64_t>(st.chr_bounds[rb.chr + 1]) -
            static_cast<int64_t>(st.chr_bounds[rb.chr]);
        if (ref_pos < 0) ref_pos = 0;
        if (ref_pos > chr_len - static_cast<int64_t>(st.span)) {
            ref_pos = std::max<int64_t>(0, chr_len - static_cast<int64_t>(st.span));
        }

        // One peak per diagonal bin over the whole window.
        VotePeak cp;
        cp.read_lo = 0;
        cp.read_hi = st.span;
        cp.ref_pos = static_cast<int>(ref_pos);
        cp.chr = rb.chr;
        cp.is_rc = st.is_rc;
        // Support is the distinct-seed count recomputed over the folded neighbourhood by the
        // re-walk above.
        cp.support = contributing_seed_count;
        cp.center_support = std::min(rb.center, contributing_seed_count);
        cp.vote_score = contributing_score;
        cp.anchor.ref_start_bin = rb.bin;
        cp.anchor.ref_start_bin_width = st.W;
        cp.anchor.median_occurrence = median_occurrence;
        cp.raw_ref_start = raw_ref_start;
        if (evidence_interval_valid && evidence_lo <= evidence_hi) {
            cp.evidence_read_lo = evidence_lo;
            cp.evidence_read_hi = evidence_hi;
        }
        vote_append_peak_if_distinct(st, out, cp);
    }
    std::sort(out.begin(), out.end(), chain_peak_better);
}

}}}  // namespace fa::cpu::lr

/**** extract ****/

namespace fa { namespace cpu { namespace lr {

// Resolves the window's query seeds into scratch.seeds: from pre-resolved seed views, an
// explicit seed list, or a fresh extraction with occurrence-aware selection (else plain
// limiting). seed_view_override is in/out: a given override is consumed, or it is pointed
// at the fresh bundle's views. active_seed_count and direct_seed_view_input are outputs for
// the accumulate and refine stages.
static void vote_resolve_seeds(
        const LongReadSeedContext& ctx,
        const uint8_t* query_enc,
        int span,
        ChainWindowPeakScratch& scratch,
        ChainSeedLookupCache* lookup_cache,
        const std::vector<QuerySeed>* seed_override,
        bool seed_override_prelimited,
        const std::vector<DnaLongSeedView>*& seed_view_override,
        size_t& active_seed_count,
        bool& direct_seed_view_input) {
    auto& seeds = scratch.seeds;
    if (seed_view_override) {
        seeds.clear();
        direct_seed_view_input = true;
        active_seed_count = seed_view_override->size();
    } else if (seed_override) {
        seeds = *seed_override;
        if (!seed_override_prelimited)
            limit_query_seeds_with_lookup_cache(ctx, seeds, span, lookup_cache);
        active_seed_count = seeds.size();
    } else {
        DnaLongSeedBundle& bundle = scratch.seed_bundle;
        bool seeds_prelimited = false;
        seeds_prelimited =
            extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
                ctx, query_enc, span, lookup_cache, bundle, scratch);
        if (seeds_prelimited) {
            seed_view_override = bundle.views_ready
                ? &bundle.selected_views
                : nullptr;
            // With views ready, take the same direct seed-view path as a caller-supplied
            // override; otherwise copy bundle.seeds.
            const bool route_direct_seed_views = seed_view_override != nullptr;
            if (route_direct_seed_views) {
                seeds.clear();
                direct_seed_view_input = true;
            } else {
                seeds = bundle.seeds;
            }
        }
        if (!seeds_prelimited) {
            extract_chain_closed_syncmer_seeds_into(ctx, query_enc, span, seeds);
            limit_query_seeds_with_lookup_cache(ctx, seeds, span, lookup_cache);
        }
        // In direct-view mode scratch.seeds is empty; the live count is the view count.
        active_seed_count =
            direct_seed_view_input ? seed_view_override->size() : seeds.size();
    }
}

}}}  // namespace fa::cpu::lr

/**** extract_peaks ****/

namespace fa {
namespace cpu {
namespace lr {

const VotePeak *best_peak(const std::vector<VotePeak> &peaks) {
  const VotePeak *best = nullptr;
  for (const auto &cp : peaks) {
    if (!best || chain_peak_better(cp, *best))
      best = &cp;
  }
  return best;
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** vote ****/

namespace fa {
namespace cpu {
namespace lr {

std::vector<VotePeak> window_anchor_peaks(
    const LongReadSeedContext &ctx, const uint8_t *query_enc, int span,
    bool is_rc, ChainWindowPeakScratch &scratch,
    const WindowAnchorPeakParams &params) {
  std::vector<VotePeak> out;
  const int limit = std::max(1, ctx.chain_max_candidates_per_window);
  out.reserve(static_cast<size_t>(limit));
  scratch.clear_for_window();
  struct PendingExactRefineReset {
    std::vector<DnaLongSeedView>& views;
    ~PendingExactRefineReset() { views.clear(); }
  } pending_exact_refine_reset{scratch.pending_exact_refine_views};
  const bool has_resolved_seed_input =
      params.seed_view_override != nullptr || params.seed_override != nullptr;
  if (span < ctx.k || ctx.index->empty() ||
      (!query_enc && !has_resolved_seed_input))
    return out;

  auto &seeds = scratch.seeds;
  size_t active_seed_count = 0;
  bool direct_seed_view_input = false;
  const std::vector<DnaLongSeedView> *seed_view_override =
      params.seed_view_override;
  vote_resolve_seeds(ctx, query_enc, span, scratch, params.lookup_cache,
                     params.seed_override, params.seed_override_prelimited,
                     seed_view_override, active_seed_count,
                     direct_seed_view_input);
  if (active_seed_count == 0)
    return out;
  auto &seed_views = scratch.seed_views;
  auto &seed_used = scratch.seed_used;
  auto &retained_seeds = scratch.retained_seeds;
  const bool retain_used_views = direct_seed_view_input;
  if (!retain_used_views) {
    seed_views.reserve(active_seed_count);
    seed_used.reserve(active_seed_count);
  } else {
    retained_seeds.reserve(active_seed_count);
  }

  const uint64_t *chr_bounds = ctx.index->chrom_offsets_data();
  const int n_chr = static_cast<int>(ctx.chr_names->size());
  if (!chr_bounds || n_chr <= 0)
    return out;

  const int W = effective_vote_diag_bin_width(ctx, span);
  const int vote_radius = (span >= 128) ? 1 : 0;
  auto &buckets = scratch.buckets;
  {
    const size_t bucket_target = std::max<size_t>(8, active_seed_count * 2);
    if (buckets.bucket_count() < bucket_target) {
      buckets.reserve(bucket_target);
    }
  }
  const int min_support = std::max(1, ctx.min_support);
  const VoteWindowState st{ctx,
                           chr_bounds,
                           n_chr,
                           W,
                           vote_radius,
                           span,
                           is_rc,
                           limit,
                           min_support,
                           params.tie_seed};
  vote_accumulate(st, scratch, seeds, seed_view_override, retain_used_views,
                  params.lookup_cache);
  if (buckets.empty())
    return out;

  vote_rank(st, scratch);
  vote_emit_peaks(st, scratch, out, seeds, retain_used_views,
                  active_seed_count);
  return out;
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** vote_slope ****/

namespace fa { namespace cpu { namespace lr {
// The whole-read winner: the best peak of fwd, then rc, by
// chain_peak_better_seeded under the read's tie seed, the first on an exact
// tie. It is catalogue candidate 0.
static VotePeak* vote_slope_best_peak(std::vector<VotePeak>& fwd,
                                      std::vector<VotePeak>& rc,
                                      std::uint32_t tie_seed) {
  VotePeak* best = nullptr;
  for (std::vector<VotePeak>* lane : {&fwd, &rc})
    for (VotePeak& peak : *lane)
      if (best == nullptr || chain_peak_better_seeded(peak, *best, tie_seed))
        best = &peak;
  return best;
}

// Widens the winner's harvest window to its per-read line (vote/vote.h),
// fitted on the seeds its lane's exact refine walked. When the fit passes its
// gate, every chain pass harvests the winner over
// [min(raw, a) - pad, max(raw + L, a + L + stretch(L, b)) + pad] and keeps the
// line (VotePeak::line_a, line_b_q20). Nothing else about the winner changes,
// and every other peak keeps the plain window.
void vote_slope_widen_winner(const LongReadSeedContext& seed_ctx,
                             std::uint32_t tie_seed, int span,
                             const ChainAnchorScratch& scratch,
                             std::vector<VotePeak>& fwd,
                             std::vector<VotePeak>& rc) {
  VotePeak* winner = vote_slope_best_peak(fwd, rc, tie_seed);
  if (winner == nullptr)
    return;
  const int lane = winner->is_rc ? 1 : 0;
  const std::uint64_t* chr_bounds = seed_ctx.index->chrom_offsets_data();
  std::vector<VoteSlopePair> pairs;
  std::vector<std::int64_t> residuals;
  const std::vector<DnaLongSeedView>& refine = scratch.slope_refine_views[lane];
  const VoteSlopeFit fit =
      !refine.empty()
          ? vote_slope_fit_winner(refine, chr_bounds, *winner, span, pairs,
                                  residuals)
          : vote_slope_fit_winner(
                (lane ? scratch.rc : scratch.fwd).retained_seeds, chr_bounds,
                *winner, span, pairs, residuals);
  if (!fit.gate)
    return;
  // b1_q20 is clamped to [-6 %, +10 %], well inside int32.
  const std::int64_t line_end =
      fit.a + vote_slope_stretch(span, static_cast<std::int32_t>(fit.b1_q20));
  winner->harvest_below = static_cast<std::int32_t>(
      std::max<std::int64_t>(0, winner->raw_ref_start - fit.a));
  winner->harvest_above = static_cast<std::int32_t>(
      std::max<std::int64_t>(0, line_end - winner->raw_ref_start));
  winner->line_a = fit.a;
  winner->line_b_q20 = static_cast<std::int32_t>(fit.b1_q20);
  winner->line_gate = true;
}

}}}  // namespace fa::cpu::lr

/**** candidate_catalogue ****/

namespace fa::cpu::voting {

CandidateCatalogue build_candidate_catalogue(
    std::vector<CandidateInput> inputs, int max_per_lane,
    double admission_ratio, CandidateMaskSource* mask_source, int tile_count,
    bool unpartitioned) {
  const int lane_bound_ceiling =
      unpartitioned ? kAllChainsLaneBound : kMaxCatalogueLaneBound;
  if (max_per_lane < 1 || max_per_lane > lane_bound_ceiling) {
    throw std::invalid_argument(
        "candidate catalogue per-lane bound must be within [1," +
        std::to_string(lane_bound_ceiling) + "]");
  }
  if (admission_ratio < 0.0 || admission_ratio > 1.0)
    throw std::invalid_argument(
        "candidate catalogue admission ratio must be within [0,1]");
  if (mask_source) {
    for (const CandidateInput& input : inputs) {
      if (input.mask_slot < 0)
        throw std::invalid_argument(
            "candidate catalogue mask source requires a mask slot per input");
    }
  }

  // Catalogue order is the caller's ranking; inputs stay in place so admission decisions
  // and mask slots keep their identity.
  std::vector<std::size_t> order(inputs.size());
  for (std::size_t index = 0; index < order.size(); ++index)
    order[index] = index;
  std::stable_sort(order.begin(), order.end(),
                   [&inputs](std::size_t left, std::size_t right) {
    return std::tuple(inputs[left].catalogue_rank,
                      -inputs[left].vote_evidence, inputs[left].lane,
                      inputs[left].equivalence_key) <
           std::tuple(inputs[right].catalogue_rank,
                      -inputs[right].vote_evidence, inputs[right].lane,
                      inputs[right].equivalence_key);
  });

  const auto mask_of = [&](std::size_t index) -> const QueryTileMask& {
    return mask_source ? mask_source->mask(inputs[index].mask_slot)
                       : inputs[index].support;
  };

  // Ratio admission judges a candidate against the best vote among the candidates that
  // compete for its query interval (as minimap2's mask_level: their tile masks overlap by at
  // least half of the smaller one), not against the read's best. A chimeric segment
  // competes with nothing and is admitted on its own vote; a paralog of the primary is held
  // to the ratio. A candidate with no tiles competes with everything. Only strictly
  // stronger candidates can drop one, so the decisions are order-independent.
  std::vector<char> admitted(inputs.size(), 1);
  if (admission_ratio > 0.0 && !inputs.empty()) {
    int vote_max = 0;
    for (const CandidateInput& input : inputs)
      vote_max = std::max(vote_max, input.vote_evidence);
    const auto ratio_floor = [admission_ratio](int vote) {
      return static_cast<int>(
          std::ceil(admission_ratio * static_cast<double>(vote)));
    };
    const auto coarse_valid = [&inputs](std::size_t index) {
      return inputs[index].coarse_tile_lo <= inputs[index].coarse_tile_hi;
    };
    const auto coarse_range = [&inputs, tile_count](std::size_t index) {
      QueryTileMask range;
      const int lo = std::max(0, inputs[index].coarse_tile_lo);
      const int hi = std::min(tile_count - 1, inputs[index].coarse_tile_hi);
      for (int tile = lo; tile <= hi; ++tile) range.set(tile);
      return range;
    };
    // Judge in descending vote order, materializing the unconditional
    // admissions' masks first, so every mask that can decide a weaker
    // candidate is available to the containment shortcut when it fires.
    std::vector<std::size_t> by_vote(order);
    std::stable_sort(by_vote.begin(), by_vote.end(),
                     [&inputs](std::size_t left, std::size_t right) {
      return inputs[left].vote_evidence > inputs[right].vote_evidence;
    });
    std::vector<char> materialized(inputs.size(), 0);
    const auto materialize = [&](std::size_t index) -> const QueryTileMask& {
      materialized[index] = 1;
      return mask_of(index);
    };
    const int strong_floor = ratio_floor(vote_max);
    for (const std::size_t index : by_vote) {
      const CandidateInput& input = inputs[index];
      if (input.vote_evidence >= strong_floor) {
        // No competitor can raise the floor above the global best: admitted
        // with no competition test. The mask is needed downstream.
        materialize(index);
        continue;
      }
      bool dropped = false;
      // Containment shortcut: a coarse range inside a materialized stronger mask proves the
      // competition without this candidate's mask, which then overlaps it fully or is empty.
      if (coarse_valid(index)) {
        const QueryTileMask range = coarse_range(index);
        for (const std::size_t other_index : by_vote) {
          if (other_index == index) continue;
          if (input.vote_evidence >=
              ratio_floor(inputs[other_index].vote_evidence))
            break;
          if (!materialized[other_index]) continue;
          const QueryTileMask& other = mask_of(other_index);
          if (other.count() == 0) continue;
          if ((range & other) == range) {
            dropped = true;
            break;
          }
        }
      }
      if (!dropped) {
        const QueryTileMask& own = materialize(index);
        const int own_tiles = own.count();
        if (own_tiles == 0) {
          dropped = true;  // competes with everything, and vote_max is above
        } else {
          for (const std::size_t other_index : by_vote) {
            if (other_index == index) continue;
            if (input.vote_evidence >=
                ratio_floor(inputs[other_index].vote_evidence))
              break;
            // Disjoint coarse ranges bound disjoint masks: not competing.
            if (coarse_valid(index) && coarse_valid(other_index) &&
                (inputs[index].coarse_tile_hi <
                     inputs[other_index].coarse_tile_lo ||
                 inputs[other_index].coarse_tile_hi <
                     inputs[index].coarse_tile_lo))
              continue;
            const QueryTileMask& other = materialize(other_index);
            const int other_tiles = other.count();
            if (other_tiles == 0) continue;
            const int shared = (own & other).count();
            if (2 * shared >= std::min(own_tiles, other_tiles)) {
              dropped = true;
              break;
            }
          }
        }
      }
      admitted[index] = dropped ? 0 : 1;
    }
  }

  CandidateCatalogue result;
  std::set<std::uint64_t> seen;
  std::map<int, int> lane_counts;
  for (const std::size_t index : order) {
    const CandidateInput& input = inputs[index];
    if (!seen.insert(input.equivalence_key).second) {
      continue;
    }
    if (lane_counts[input.lane] >= max_per_lane) {
      continue;
    }
    if (!admitted[index]) {
      continue;
    }
    ++lane_counts[input.lane];
    QueryCandidate candidate;
    candidate.id = static_cast<CandidateId>(result.candidates.size());
    candidate.equivalence_key = input.equivalence_key;
    candidate.lane = input.lane;
    candidate.catalogue_rank = input.catalogue_rank;
    candidate.vote_evidence = input.vote_evidence;
    candidate.support = mask_of(index);
    result.candidates.push_back(candidate);
  }
  return result;
}

}  // namespace fa::cpu::voting
