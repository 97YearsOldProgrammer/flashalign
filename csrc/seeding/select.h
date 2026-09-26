// Query-seed occurrence evidence and low-occurrence seed selection.
#pragma once

#include "context.h"              // LongReadSeedContext
#include "scratch.h"              // ChainSeedLookupCache
#include "../index/index.h"
#include "../index/format.h"      // KmerPostingView
#include "../index/seed.h"        // QuerySeed

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace fa { namespace cpu {

struct QuerySeedOccurrenceEvidence {
    QuerySeed seed{};
    uint32_t occurrence = 0;     // Reference posting count for seed.key.
    float specificity = 0.0f;    // 1 / occurrence.
    bool found = false;          // Same as `retained`.
    bool present = false;        // Key is present in the index.
    bool retained = false;       // Retained postings are available for voting.
    int original_index = 0;
};

inline float query_seed_specificity(uint32_t occurrence) {
    return occurrence > 0 ? 1.0f / static_cast<float>(occurrence) : 0.0f;
}

inline QuerySeedOccurrenceEvidence query_seed_occurrence(
    const IndexView& index,
    const QuerySeed& seed,
    int original_index = 0
) {
    QuerySeedOccurrenceEvidence evidence;
    evidence.seed = seed;
    const KmerPostingView view = index.lookup(seed.key);
    evidence.present = view.found();
    evidence.retained = view.found() && view.count > 0;
    evidence.found = evidence.retained;
    evidence.occurrence = evidence.present ? view.occurrence : 0;
    evidence.specificity = query_seed_specificity(evidence.occurrence);
    evidence.original_index = original_index;
    return evidence;
}

inline bool better_query_seed_occurrence(
    const QuerySeedOccurrenceEvidence& a,
    const QuerySeedOccurrenceEvidence& b
) {
    if (a.found != b.found) return a.found && !b.found;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    if (a.seed.read_pos != b.seed.read_pos) return a.seed.read_pos < b.seed.read_pos;
    return a.original_index < b.original_index;
}

inline void restore_query_seed_read_order(std::vector<QuerySeed>& seeds) {
    std::sort(seeds.begin(), seeds.end(),
        [](const QuerySeed& a, const QuerySeed& b) {
            if (a.read_pos != b.read_pos) return a.read_pos < b.read_pos;
            return a.key < b.key;
        });
}

inline std::vector<QuerySeedOccurrenceEvidence> query_seed_occurrence_evidence(
    const std::vector<QuerySeed>& seeds,
    const IndexView& index
) {
    std::vector<QuerySeedOccurrenceEvidence> evidence;
    evidence.reserve(seeds.size());
    for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
        evidence.push_back(query_seed_occurrence(index, seeds[(size_t)i], i));
    }
    return evidence;
}

// Tile-spread pick: keep the lowest-occurrence seed in each of max_seeds read tiles, so kept
// seeds are rare and spread across the read; top up from a global low-occurrence sort if
// short, then restore read order.
inline void select_spread_from_evidence(
    std::vector<QuerySeed>& seeds,
    const std::vector<QuerySeedOccurrenceEvidence>& evidence,
    int read_len,
    int max_seeds
) {
    const int n_tiles = std::max(
        1,
        std::min(max_seeds, static_cast<int>(evidence.size())));
    std::vector<int> best_by_tile(static_cast<size_t>(n_tiles), -1);
    for (int i = 0; i < static_cast<int>(evidence.size()); ++i) {
        int tile = 0;
        if (read_len > 0) {
            tile = static_cast<int>(
                (static_cast<int64_t>(
                     std::max(0, evidence[static_cast<size_t>(i)].seed.read_pos)) *
                 static_cast<int64_t>(n_tiles)) /
                static_cast<int64_t>(read_len));
            tile = std::max(0, std::min(n_tiles - 1, tile));
        }
        const int prev = best_by_tile[static_cast<size_t>(tile)];
        if (prev < 0 ||
            better_query_seed_occurrence(
                evidence[static_cast<size_t>(i)],
                evidence[static_cast<size_t>(prev)])) {
            best_by_tile[static_cast<size_t>(tile)] = i;
        }
    }

    std::vector<char> selected(evidence.size(), 0);
    seeds.clear();
    seeds.reserve(static_cast<size_t>(max_seeds));
    for (int idx : best_by_tile) {
        if (idx < 0 || selected[static_cast<size_t>(idx)]) continue;
        selected[static_cast<size_t>(idx)] = 1;
        seeds.push_back(evidence[static_cast<size_t>(idx)].seed);
        if (static_cast<int>(seeds.size()) >= max_seeds) break;
    }

    if (static_cast<int>(seeds.size()) < max_seeds) {
        std::vector<int> order(evidence.size());
        for (int i = 0; i < static_cast<int>(order.size()); ++i) {
            order[static_cast<size_t>(i)] = i;
        }
        std::stable_sort(
            order.begin(),
            order.end(),
            [&](int a, int b) {
                return better_query_seed_occurrence(
                    evidence[static_cast<size_t>(a)],
                    evidence[static_cast<size_t>(b)]);
            });
        for (int idx : order) {
            if (selected[static_cast<size_t>(idx)]) continue;
            selected[static_cast<size_t>(idx)] = 1;
            seeds.push_back(evidence[static_cast<size_t>(idx)].seed);
            if (static_cast<int>(seeds.size()) >= max_seeds) break;
        }
    }
    restore_query_seed_read_order(seeds);
}

// Gathers evidence from the index, then spread-picks.
inline void select_low_occurrence_spread_query_seeds(
    std::vector<QuerySeed>& seeds,
    const IndexView& index,
    int read_len,
    int max_seeds
) {
    const auto evidence = query_seed_occurrence_evidence(seeds, index);
    select_spread_from_evidence(seeds, evidence, read_len, max_seeds);
}

}}  // namespace fa::cpu

namespace fa { namespace cpu { namespace lr {

inline QuerySeedOccurrenceEvidence chain_query_seed_occurrence_cached(
    const LongReadSeedContext& ctx,
    const QuerySeed& seed,
    int original_index,
    ChainSeedLookupCache* lookup_cache
) {
    if (!lookup_cache) {
        return query_seed_occurrence(*ctx.index, seed, original_index);
    }
    QuerySeedOccurrenceEvidence evidence;
    evidence.seed = seed;
    const KmerPostingView view = lookup_cache->lookup(*ctx.index, seed.key);
    evidence.present = view.found();
    evidence.retained = view.found() && view.count > 0;
    evidence.found = evidence.retained;
    evidence.occurrence = evidence.present ? view.occurrence : 0;
    evidence.specificity = query_seed_specificity(evidence.occurrence);
    evidence.original_index = original_index;
    return evidence;
}

inline std::vector<QuerySeedOccurrenceEvidence> chain_query_seed_occurrence_evidence_cached(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& seeds,
    ChainSeedLookupCache* lookup_cache
) {
    std::vector<QuerySeedOccurrenceEvidence> evidence;
    evidence.reserve(seeds.size());
    for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
        evidence.push_back(
            chain_query_seed_occurrence_cached(
                ctx, seeds[static_cast<size_t>(i)], i, lookup_cache));
    }
    return evidence;
}

inline void chain_select_low_occurrence_spread_query_seeds_cached(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len,
    int max_seeds,
    ChainSeedLookupCache* lookup_cache
) {
    const auto evidence = chain_query_seed_occurrence_evidence_cached(
        ctx, seeds, lookup_cache);
    select_spread_from_evidence(seeds, evidence, read_len, max_seeds);
}

// Query-seed limiting when occurrence-aware selection is off (RNA, or a window with no
// bundle): over budget, keep the lowest-occurrence seed per read strip. A no-op when
// uncapped (max_query_seeds <= 0, the RNA default).
inline void limit_query_seeds(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len
) {
    const int max_seeds = ctx.max_query_seeds_per_strand;
    if (max_seeds <= 0 || static_cast<int>(seeds.size()) <= max_seeds) return;
    select_low_occurrence_spread_query_seeds(
        seeds, *ctx.index, read_len, max_seeds);
}

inline void limit_query_seeds_with_lookup_cache(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len,
    ChainSeedLookupCache* lookup_cache
) {
    if (!lookup_cache) {
        limit_query_seeds(ctx, seeds, read_len);
        return;
    }
    const int max_seeds = ctx.max_query_seeds_per_strand;
    if (max_seeds <= 0 || static_cast<int>(seeds.size()) <= max_seeds) return;
    chain_select_low_occurrence_spread_query_seeds_cached(
        ctx, seeds, read_len, max_seeds, lookup_cache);
}

}}}  // namespace fa::cpu::lr
