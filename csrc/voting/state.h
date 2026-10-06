#pragma once

#include "../seeding/context.h"   // LongReadSeedContext, occurrence policy
#include "../seeding/types.h"     // ChainWindowRankedBucket, chain_window_ranked_better
#include "../index/format.h"    // KmerPostingView
#include "../index/seed.h"      // QuerySeed

#include <algorithm>
#include <cstdint>

namespace fa { namespace cpu { namespace lr {

inline bool chain_window_seed_occ_allowed(
    const LongReadSeedContext& ctx,
    const KmerPostingView& view
) {
    return seed_allowed_by_long_occ_policy(view, ctx.occ_policy);
}

// Per-window invariants shared by the vote stages (accumulate, rank, refine, emit): the
// context and the window geometry. Built once after seed extraction, never mutated.
struct VoteWindowState {
    const LongReadSeedContext& ctx;
    const uint64_t* chr_bounds;
    int n_chr;
    int W;                       // diagonal-bin width
    int vote_radius;             // neighbor-fold radius (span>=128 ? 1 : 0)
    int span;
    bool is_rc;
    int limit;                   // max emitted peaks per window
    int min_support;
    // The read's tie seed (tie_read_seed) for the drain heap's tie-break.
    std::uint32_t tie_seed;
};

// The drain heap's order: std::make_heap / std::pop_heap take "a ranks below b", so this is
// chain_window_ranked_better with its arguments swapped, under the read seed and strand.
struct VoteHeapBelow {
    const VoteWindowState* st;
    bool operator()(const ChainWindowRankedBucket& a,
                    const ChainWindowRankedBucket& b) const {
        return chain_window_ranked_better(b, a, st->tie_seed, st->is_rc);
    }
};

// Floor toward negative infinity: diagonal bin index for a (possibly negative)
// ref_start.
inline int vote_floor_div(int x, int w) {
    // floor(x / w) for every int x when w > 0, with no overflow at INT_MIN.
    int q = x / w;
    int r = x % w;
    return (r != 0 && ((r < 0) != (w < 0))) ? q - 1 : q;
}

// Pack (chr, bin) into one 64-bit bucket key.
inline int64_t vote_pack_key(int chr, int bin) {
    return (static_cast<int64_t>(chr) << 32)
         | static_cast<int64_t>(static_cast<uint32_t>(bin));
}

// Occurrence admission for one DNA vote seed, under the context's occurrence policy, on
// the occurrence the vote counts (vote_seed_occurrence).
inline bool vote_seed_occ_allowed(const LongReadSeedContext& ctx,
                                  const KmerPostingView& view) {
    if (ctx.self_contig < 0) return chain_window_seed_occ_allowed(ctx, view);
    return occurrence_allowed_by_long_occ_policy(
        vote_seed_occurrence(ctx, view), ctx.occ_policy);
}

inline bool vote_seed_occ_allowed(const VoteWindowState& st,
                                  const KmerPostingView& view) {
    return vote_seed_occ_allowed(st.ctx, view);
}

}}}  // namespace fa::cpu::lr
