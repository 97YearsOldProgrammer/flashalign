#pragma once

#include "state.h"             // VoteWindowState
#include "../seeding/scratch.h"          // ChainWindowPeakScratch
#include "../seeding/types.h"            // VotePeak, ChainWindowRankedBucket

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// Appends cp unless an emitted peak on its contig lies within W * (vote_radius + 1)
// of its raw median diagonal (raw_ref_start, not the clamped ref_pos), so peaks that
// clamp to one start at a contig end stay distinct. Returns true iff appended.
inline bool vote_append_peak_if_distinct(const VoteWindowState& st,
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
