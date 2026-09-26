#pragma once

#include "state.h"           // VoteWindowState, vote_pack_key
#include "../seeding/scratch.h"         // ChainWindowPeakScratch
#include "../seeding/types.h"           // ChainWindowRankedBucket, chain_window_ranked_better

#include <algorithm>
#include <cstdint>

namespace fa { namespace cpu { namespace lr {

// The vote's rank stage: folds each bucket of scratch.buckets with its vote_radius
// neighbours into scratch.ranked, drops diagonals below min_support, and heapifies by
// chain_window_ranked_better under the read's tie seed.
inline void vote_rank(const VoteWindowState& st,
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
