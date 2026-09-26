#pragma once

#include "../seeding/context.h"               // LongReadSeedContext
#include "../seeding/scratch.h"               // ChainWindowPeakScratch, ChainSeedLookupCache, DnaLongSeedBundle
#include "../seeding/types.h"                 // DnaLongSeedView
#include "../seeding/syncmer.h"               // closed-syncmer extraction and selection
#include "../seeding/select.h"                // limit_query_seeds_with_lookup_cache
#include "../index/seed.h"   // QuerySeed

#include <algorithm>
#include <cstdint>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// Resolves the window's query seeds into scratch.seeds: from pre-resolved seed views, an
// explicit seed list, or a fresh extraction with occurrence-aware selection (else plain
// limiting). seed_view_override is in/out: a given override is consumed, or it is pointed
// at the fresh bundle's views. active_seed_count and direct_seed_view_input are outputs for
// the accumulate and refine stages.
inline void vote_resolve_seeds(
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
