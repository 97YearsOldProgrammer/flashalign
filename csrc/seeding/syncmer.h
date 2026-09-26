// Closed-syncmer seed extraction with occurrence-aware admission for long-read mapping, and
// the slicer that cuts a read's seed stream into windows without re-extracting.
#pragma once

#include "context.h"              // LongReadSeedContext, occurrence policy
#include "types.h"                // DnaLongSeedBundle, DnaLongSeedView
#include "scratch.h"              // ChainSeedLookupCache, ChainWindowPeakScratch
#include "../index/format.h"      // KmerPostingView
#include "../index/seed.h"        // QuerySeed, ClosedSyncmerConfig, extraction
#include "../core/hash.h"         // sketch_order_hash

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <vector>

namespace fa { namespace cpu { namespace lr {

inline bool chain_syncmer_occ_candidate_better(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b
) {
    if (!a.valid) return false;
    if (!b.valid) return true;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    return a.seed.read_pos < b.seed.read_pos;
}

inline bool chain_syncmer_occ_candidate_better_centred(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b,
    int tile_centre
) {
    if (!a.valid) return false;
    if (!b.valid) return true;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    const int a_distance = std::abs(a.seed.read_pos - tile_centre);
    const int b_distance = std::abs(b.seed.read_pos - tile_centre);
    if (a_distance != b_distance) return a_distance < b_distance;
    return a.seed.read_pos < b.seed.read_pos;
}

inline int syncmer_tile_centre(int tile, int n_kmers, int n_tiles) {
    if (n_tiles <= 1) return std::max(0, n_kmers - 1) / 2;
    const int64_t kmers = static_cast<int64_t>(std::max(1, n_kmers));
    const int64_t tiles = static_cast<int64_t>(n_tiles);
    const int64_t lo =
        (static_cast<int64_t>(tile) * kmers + tiles - 1) / tiles;
    const int64_t hi =
        (static_cast<int64_t>(tile + 1) * kmers + tiles - 1) / tiles - 1;
    return static_cast<int>((lo + std::max(lo, hi)) / 2);
}

inline constexpr int kVoteSeedOccurrenceStrata = 32;

inline int vote_seed_occurrence_stratum(uint32_t occurrence) {
    if (occurrence <= 1) return 0;
    const int floor_log2 =
        63 - __builtin_clzll(static_cast<unsigned long long>(occurrence));
    return std::min(floor_log2, kVoteSeedOccurrenceStrata - 1);
}

inline bool vote_seed_minimum_is_representative(
    const std::array<uint32_t, kVoteSeedOccurrenceStrata>& strata,
    uint32_t admitted,
    uint32_t minimum_occurrence
) {
    if (admitted == 0 || minimum_occurrence == 0) return false;
    const int minimum_stratum =
        vote_seed_occurrence_stratum(minimum_occurrence);
    uint32_t near = strata[static_cast<size_t>(minimum_stratum)];
    if (minimum_stratum + 1 < kVoteSeedOccurrenceStrata) {
        near += strata[static_cast<size_t>(minimum_stratum + 1)];
    }
    // The rare minimum is representative only when at least half of the tile shares its
    // floor(log2 occurrence) band or the next higher one.
    return static_cast<uint64_t>(near) * 2u >= admitted;
}

inline void extract_chain_closed_syncmer_seeds_into_downsample(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    int downsample,
    std::vector<QuerySeed>& out
) {
    ClosedSyncmerConfig cfg;
    cfg.k = ctx.k;
    cfg.s = ctx.syncmer_s;
    cfg.downsample = std::max(1, downsample);
    extract_closed_syncmer_query_seeds_into(query_enc, span, cfg, out);
}

inline void extract_chain_closed_syncmer_seeds_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    std::vector<QuerySeed>& out
) {
    extract_chain_closed_syncmer_seeds_into_downsample(
        ctx, query_enc, span, std::max(1, ctx.syncmer_downsample), out);
}

inline bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
) {
    out.clear();
    (void)scratch;
    const int k = ctx.k;
    const SeedIndex& index = *ctx.index;
    const LongOccPolicyConfig& cfg = ctx.occ_policy;
    if (!ctx.syncmer_occ_aware_active() || !lookup_cache ||
        input.empty() || !htable_valid_k(k) || span < k) {
        return false;
    }
    const int max_seeds = ctx.max_query_seeds_per_strand;
    const int n_kmers = std::max(1, span - k + 1);
    const int n_tiles = std::max(
        1,
        max_seeds > 0
            ? std::min(max_seeds, n_kmers)
            : n_kmers);
    auto make_candidate_from_view =
        [&](const QuerySeed& seed, const KmerPostingView& view) {
        ChainSyncmerOccCandidate cand;
        cand.seed = seed;
        cand.view = view;
        if (!view.found() || view.count == 0 || view.occurrence == 0) {
            return cand;
        }
        cand.occurrence = view.occurrence;
        cand.valid = seed_allowed_by_long_occ_policy(view, cfg);
        return cand;
    };

    // `seed_views`, when captured (ChainSeedLookupCache::warm), holds the cached view for
    // each position of `input`, read by position instead of by a probe.
    auto make_candidate = [&](const QuerySeed& seed, size_t position) {
        KmerPostingView view;
        if (seed_views != nullptr) {
            view = seed_views[position];
        } else {
            view = lookup_cache->lookup(index, seed.key);
        }
        return make_candidate_from_view(seed, view);
    };

    out.selected_views.reserve(static_cast<size_t>(n_tiles));
    out.exact_refine_views.reserve(static_cast<size_t>(3 * n_tiles));
    out.seeds.reserve(static_cast<size_t>(n_tiles));

    auto emit_tile = [&](const ChainSyncmerOccCandidate& first,
                         const ChainSyncmerOccCandidate& rarest,
                         const ChainSyncmerOccCandidate& centred,
                         const std::array<uint32_t,
                                          kVoteSeedOccurrenceStrata>& strata,
                         uint32_t admitted) {
      if (!first.valid)
        return;
      const ChainSyncmerOccCandidate& representative =
          vote_seed_minimum_is_representative(
              strata, admitted, rarest.occurrence)
              ? centred
              : first;
      out.selected_views.push_back(
          {representative.seed, representative.view});
      out.seeds.push_back(representative.seed);
      out.exact_refine_views.push_back({first.seed, first.view});
      if (rarest.valid && rarest.seed.read_pos != first.seed.read_pos) {
        out.exact_refine_views.push_back({rarest.seed, rarest.view});
      }
      if (representative.seed.read_pos != first.seed.read_pos &&
          representative.seed.read_pos != rarest.seed.read_pos) {
        out.exact_refine_views.push_back(
            {representative.seed, representative.view});
      }
    };

    // Streaming per-tile scan: each tile keeps its first, rarest and centred admitted seeds,
    // from which emit_tile picks the vote representative and the exact-refinement views.
    const int64_t tiles64 = static_cast<int64_t>(n_tiles);
    const int64_t kmers64 = static_cast<int64_t>(std::max(1, n_kmers));
    int tile = 0;
    ChainSyncmerOccCandidate first_admitted;
    ChainSyncmerOccCandidate rarest_admitted;
    ChainSyncmerOccCandidate centred_admitted;
    std::array<uint32_t, kVoteSeedOccurrenceStrata> tile_strata{};
    uint32_t tile_admitted = 0;
    int tile_centre = syncmer_tile_centre(0, n_kmers, n_tiles);
#ifndef NDEBUG
    int prev_pos = INT_MIN;
#endif
    for (size_t position = 0; position < input.size(); ++position) {
        const QuerySeed& seed = input[position];
#ifndef NDEBUG
        assert(seed.read_pos >= prev_pos &&
               "occ-stream selector requires read_pos-sorted input");
        prev_pos = seed.read_pos;
#endif
        const int cp =
            std::max(0, std::min(std::max(0, n_kmers - 1), seed.read_pos));
        const int prev_tile = tile;
        while (static_cast<int64_t>(cp) * tiles64 >=
               static_cast<int64_t>(tile + 1) * kmers64) {
            ++tile;
        }
        if (tile != prev_tile) {
          emit_tile(first_admitted, rarest_admitted, centred_admitted,
                    tile_strata, tile_admitted);
          first_admitted = ChainSyncmerOccCandidate{};
          rarest_admitted = ChainSyncmerOccCandidate{};
          centred_admitted = ChainSyncmerOccCandidate{};
          tile_strata.fill(0);
          tile_admitted = 0;
          tile_centre = syncmer_tile_centre(tile, n_kmers, n_tiles);
        }
        const ChainSyncmerOccCandidate cand = make_candidate(seed, position);
        if (cand.valid) {
          if (!first_admitted.valid)
            first_admitted = cand;
          if (chain_syncmer_occ_candidate_better(cand, rarest_admitted)) {
            rarest_admitted = cand;
          }
          if (chain_syncmer_occ_candidate_better_centred(
                  cand, centred_admitted, tile_centre)) {
            centred_admitted = cand;
          }
          ++tile_strata[static_cast<size_t>(
              vote_seed_occurrence_stratum(cand.occurrence))];
          ++tile_admitted;
        }
    }
    emit_tile(first_admitted, rarest_admitted, centred_admitted, tile_strata,
              tile_admitted);
    out.views_ready = true;
    return true;
}

inline bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch
) {
    out.clear();
    if (!query_enc || !htable_valid_k(ctx.k) || span < ctx.k) return false;
    scratch.extracted_seeds.clear();
    extract_chain_closed_syncmer_seeds_into(
        ctx, query_enc, span, scratch.extracted_seeds);
    return extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
        ctx, scratch.extracted_seeds, span, lookup_cache, out, scratch);
}

inline bool select_chain_closed_syncmer_occ_aware_seed_bundle_from_seeds_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
) {
    return extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
        ctx, input, span, lookup_cache, out, scratch, seed_views);
}

// Slices a read's seed stream to the window [start, end) and renumbers read_pos from the
// window start; key and z are copied unchanged. `k` sets the last valid read_pos (end - k).
// `downsample > 1` skips seeds whose sketch_order_hash(key) % downsample != 0. `seeds` must
// be sorted by read_pos, as the extractor emits them. Same seeds as re-extracting the window.
inline void rebase_syncmer_stream_window(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample,
    std::vector<QuerySeed>& out
) {
    out.clear();
    if (end - start < k) return;
    const int last_pos = end - k;
    auto it = std::lower_bound(
        seeds.begin(),
        seeds.end(),
        start,
        [](const QuerySeed& seed, int pos) {
            return seed.read_pos < pos;
        });
    for (; it != seeds.end() && it->read_pos <= last_pos; ++it) {
        if (downsample > 1 &&
            sketch_order_hash(it->key) %
                    static_cast<uint64_t>(downsample) != 0) {
            continue;
        }
        out.push_back(
            QuerySeed{it->key, it->read_pos - start, it->z});
    }
}

// True when rebase_syncmer_stream_window() with the same arguments would copy `seeds`
// unchanged (window from 0 holding every seed, no downsampling), so a caller can pass the
// stream through. O(1): the input is sorted, so its two ends decide.
inline bool syncmer_stream_window_is_identity(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample
) {
    if (downsample > 1 || start != 0 || end - start < k) return false;
    if (seeds.empty()) return true;
    return seeds.front().read_pos >= start &&
           seeds.back().read_pos <= end - k;
}

}}}  // namespace fa::cpu::lr
