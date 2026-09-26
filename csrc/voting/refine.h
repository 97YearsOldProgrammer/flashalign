// The vote's refine and emit stage (the heap drain): pops diagonals best first, re-walks
// each bucket's used seeds (the exact-refine pass) to recompute the median reference
// start and median occurrence, and emits one VotePeak per bucket.
#pragma once

#include "state.h"                 // VoteWindowState, vote_floor_div
#include "emit.h"                  // vote_append_peak_if_distinct
#include "../seeding/scratch.h"               // ChainWindowPeakScratch
#include "../seeding/types.h"                 // VotePeak, ChainWindowRankedBucket, comparators
#include "../index/seed.h"     // QuerySeed
#include "../index/format.h"   // KmerPostingView, chromosome_index_for_global_pos

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// Batched drain (ctx.vote_batched_refine): the same output as the sequential drain in
// vote_emit_peaks, but the exact-refine walk is paid once per chunk of drained buckets
// instead of once per bucket. Buckets are processed in heap order, so support counts,
// medians, the min_support skip and distinctness replay the sequential path, and a chunk
// never pops more buckets than emissions still needed.
inline void vote_emit_peaks_batched(const VoteWindowState& st,
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
    // The chunk's entries in (chr, bin) order, so a posting visits only the entries of its
    // chromosome within coarse_radius bins of it. Each entry still sees its postings in
    // order, so the emitted peaks match a full scan.
    std::vector<size_t> by_chr_bin;
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
            chunk_ref_starts[c].reserve(
                static_cast<size_t>(std::max(0, chunk[c].support)));
            chunk_occurrences[c].reserve(active_seed_count);
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
        auto visit_seed = [&](const QuerySeed& seed,
                              const KmerPostingView& v,
                              int /*seed_id*/) {
            if (!v.found() || v.count == 0) return;
            int chr_idx = chromosome_index_for_global_pos(
                st.chr_bounds, st.n_chr, v.positions[0]);
            if (chr_idx < 0 || chr_idx >= st.n_chr) return;
            uint64_t chr_lo = st.chr_bounds[static_cast<size_t>(chr_idx)];
            uint64_t chr_hi = st.chr_bounds[static_cast<size_t>(chr_idx + 1)];
            std::fill(seed_hit.begin(), seed_hit.end(), uint8_t{0});
            for (uint32_t i = 0; i < v.count; ++i) {
                const uint64_t g = static_cast<uint64_t>(v.positions[i]);
                advance_contig_cursor(st.chr_bounds, st.n_chr, g, chr_idx,
                                      chr_lo, chr_hi);
                if (g < chr_lo || g >= chr_hi) continue;
                // The accumulate stage's strand test, so every statistic here counts
                // compatible postings only; tested after the cursor advance.
                if (!packed_ref_orientation_compatible(
                        seed.z, v.positions.packed_at(i), st.is_rc))
                    continue;
                const int local = static_cast<int>(g - chr_lo);
                const int ref_start = local - seed.read_pos;
                const int bin = vote_floor_div(ref_start, st.W);
                // First entry of this chromosome at or past bin - coarse_radius.
                const int bin_lo = bin - coarse_radius;
                size_t lo = 0, hi = chunk_size;
                while (lo < hi) {
                    const size_t mid = lo + (hi - lo) / 2;
                    const ChainWindowRankedBucket& e = chunk[by_chr_bin[mid]];
                    if (e.chr < chr_idx || (e.chr == chr_idx && e.bin < bin_lo))
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                for (size_t o = lo; o < chunk_size; ++o) {
                    const size_t c = by_chr_bin[o];
                    if (chunk[c].chr != chr_idx) break;
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
                        seed_hit[c] = 1;
                    }
                }
            }
            for (size_t c = 0; c < chunk_size; ++c) {
                if (!seed_hit[c]) continue;
                // Unweighted vote: one unit per distinct query seed.
                ++chunk_seed_count[c];
                chunk_occurrences[c].push_back(v.occurrence);
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
inline void vote_emit_peaks(const VoteWindowState& st,
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
