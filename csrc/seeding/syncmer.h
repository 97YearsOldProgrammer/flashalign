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
#include <iterator>
#include <set>
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

// The vote's empty-tile rescue (options/dna_profile.h kDnaTileRescueOcc): per read and
// strand the rescued seeds' summed occurrence stays within kDnaTileRescueBudget, and reads
// shorter than kDnaTileRescueMinLen bp are not rescued.
inline constexpr uint64_t kDnaTileRescueBudget = 65536;
inline constexpr int kDnaTileRescueMinLen = 1000;

// One strand's rescue, walked in tile order: offer() every seed of the current tile, then
// close() the tile. A tile that admitted no seed votes with its rarest offered seed over the
// cap and at or under M, ranked as the tile's centred seed.
class VoteTileRescue {
public:
    VoteTileRescue(const LongReadSeedContext& ctx, int span)
        : max_occ_(span >= kDnaTileRescueMinLen
                       ? static_cast<uint32_t>(std::max(0, ctx.tile_rescue_occ))
                       : 0u) {}

    void offer(const ChainSyncmerOccCandidate& cand, int tile_centre) {
        // A seed with no postings keeps occurrence UINT32_MAX.
        if (cand.valid || cand.occurrence > max_occ_) return;
        ChainSyncmerOccCandidate over = cand;
        over.valid = true;
        if (chain_syncmer_occ_candidate_better_centred(over, best_, tile_centre))
            best_ = over;
    }

    // The closing tile's rescued seed, when it admitted none and the budget allows.
    bool close(uint32_t admitted, DnaLongSeedView& out) {
        const ChainSyncmerOccCandidate best = best_;
        best_ = ChainSyncmerOccCandidate{};
        if (admitted != 0 || !best.valid ||
            spent_ + best.occurrence > kDnaTileRescueBudget)
            return false;
        spent_ += best.occurrence;
        out = DnaLongSeedView{best.seed, best.view, true};
        return true;
    }

private:
    uint32_t max_occ_;
    uint64_t spent_ = 0;
    ChainSyncmerOccCandidate best_;
};

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

// One vote seed per tile of max_query_seeds_per_strand equal tiles (every
// seed when it is 0).
inline bool extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
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
        const uint32_t occurrence = vote_seed_occurrence(ctx, view);
        if (!view.found() || view.count == 0 || occurrence == 0) {
            return cand;
        }
        cand.occurrence = occurrence;
        cand.valid = occurrence_allowed_by_long_occ_policy(occurrence, cfg);
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

    // A tile that admits no seed may vote with a rescued one, also an exact-refine view.
    VoteTileRescue rescue(ctx, span);
    auto emit_rescue = [&](uint32_t admitted) {
      DnaLongSeedView rescued;
      if (!rescue.close(admitted, rescued))
        return;
      out.selected_views.push_back(rescued);
      out.seeds.push_back(rescued.seed);
      out.exact_refine_views.push_back(rescued);
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
          emit_rescue(tile_admitted);
          first_admitted = ChainSyncmerOccCandidate{};
          rarest_admitted = ChainSyncmerOccCandidate{};
          centred_admitted = ChainSyncmerOccCandidate{};
          tile_strata.fill(0);
          tile_admitted = 0;
          tile_centre = syncmer_tile_centre(tile, n_kmers, n_tiles);
        }
        const ChainSyncmerOccCandidate cand = make_candidate(seed, position);
        rescue.offer(cand, tile_centre);
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
    emit_rescue(tile_admitted);
    out.views_ready = true;
    return true;
}

// The vote seeds of one strand. With ctx.nested_vote_seeds and more than
// kVoteSeedNestBase seeds asked, the kVoteSeedNestBase selection is kept and
// filled from the full selection, up to its size: next is the seed farthest
// from every kept position, ties to the lower occurrence, then the lower
// position. A rescued seed that would take the rescued occurrence past
// kDnaTileRescueBudget is passed over, so the fill can stop short; when the
// kept selection is the larger, nothing is added. The exact-refine views are
// the union of both selections'.
inline bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
) {
    const bool ready =
        extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
            ctx, input, span, lookup_cache, out, scratch, seed_views);
    if (!ready || !ctx.nested_vote_seeds ||
        ctx.max_query_seeds_per_strand <= kVoteSeedNestBase)
        return ready;
    LongReadSeedContext base_ctx = ctx;
    base_ctx.max_query_seeds_per_strand = kVoteSeedNestBase;
    DnaLongSeedBundle base;
    extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
        base_ctx, input, span, lookup_cache, base, scratch, seed_views);

    const auto by_position = [](const DnaLongSeedView& a,
                                const DnaLongSeedView& b) {
        return a.seed.read_pos < b.seed.read_pos;
    };
    out.exact_refine_views.insert(out.exact_refine_views.end(),
                                  base.exact_refine_views.begin(),
                                  base.exact_refine_views.end());
    std::sort(out.exact_refine_views.begin(), out.exact_refine_views.end(),
              by_position);
    out.exact_refine_views.erase(
        std::unique(out.exact_refine_views.begin(),
                    out.exact_refine_views.end(),
                    [](const DnaLongSeedView& a, const DnaLongSeedView& b) {
                        return a.seed.read_pos == b.seed.read_pos;
                    }),
        out.exact_refine_views.end());

    const std::size_t target = out.selected_views.size();
    std::vector<DnaLongSeedView> selected = std::move(base.selected_views);
    std::set<int> kept;
    uint64_t rescued_occurrence = 0;
    for (const DnaLongSeedView& v : selected) {
        kept.insert(v.seed.read_pos);
        if (v.rescued) rescued_occurrence += vote_seed_occurrence(ctx, v.view);
    }
    const auto distance_to_kept = [&](int read_pos) {
        const auto hi = kept.lower_bound(read_pos);
        int distance = span;
        if (hi != kept.end()) distance = *hi - read_pos;
        if (hi != kept.begin())
            distance = std::min(distance, read_pos - *std::prev(hi));
        return distance;
    };
    // A max-heap; distances only shrink, so a popped seed whose distance has
    // shrunk since its push goes back with the new one.
    struct Fill {
        int distance;
        uint32_t occurrence;
        int read_pos;
        std::size_t index;
    };
    const auto after = [](const Fill& a, const Fill& b) {
        if (a.distance != b.distance) return a.distance < b.distance;
        if (a.occurrence != b.occurrence) return a.occurrence > b.occurrence;
        return a.read_pos > b.read_pos;
    };
    std::vector<Fill> heap;
    if (selected.size() < target) {
        for (std::size_t i = 0; i < out.selected_views.size(); ++i) {
            const DnaLongSeedView& v = out.selected_views[i];
            if (kept.count(v.seed.read_pos)) continue;
            heap.push_back({distance_to_kept(v.seed.read_pos),
                            vote_seed_occurrence(ctx, v.view),
                            v.seed.read_pos, i});
        }
        std::make_heap(heap.begin(), heap.end(), after);
    }
    while (selected.size() < target && !heap.empty()) {
        std::pop_heap(heap.begin(), heap.end(), after);
        Fill next = heap.back();
        heap.pop_back();
        const DnaLongSeedView& v = out.selected_views[next.index];
        if (v.rescued &&
            rescued_occurrence + vote_seed_occurrence(ctx, v.view) >
                kDnaTileRescueBudget)
            continue;
        const int distance = distance_to_kept(next.read_pos);
        if (distance != next.distance) {
            next.distance = distance;
            heap.push_back(next);
            std::push_heap(heap.begin(), heap.end(), after);
            continue;
        }
        selected.push_back(v);
        kept.insert(v.seed.read_pos);
        if (v.rescued) rescued_occurrence += vote_seed_occurrence(ctx, v.view);
    }
    std::sort(selected.begin(), selected.end(), by_position);
    out.selected_views = std::move(selected);
    out.seeds.clear();
    for (const DnaLongSeedView& v : out.selected_views)
        out.seeds.push_back(v.seed);
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
