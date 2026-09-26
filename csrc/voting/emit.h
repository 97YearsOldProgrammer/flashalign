#pragma once

#include "state.h"             // VoteWindowState
#include "../seeding/scratch.h"          // ChainWindowPeakScratch
#include "../seeding/types.h"            // VotePeak, ChainWindowRankedBucket

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// Appends cp unless an emitted peak lies within W * (vote_radius + 1) of its diagonal.
// Returns true iff appended.
inline bool vote_append_peak_if_distinct(const VoteWindowState& st,
                                         std::vector<VotePeak>& out,
                                         VotePeak& cp) {
    bool distinct = true;
    for (const auto& prev : out) {
        if (prev.chr == cp.chr &&
            std::abs(prev.ref_pos - cp.ref_pos) <= st.W * (st.vote_radius + 1)) {
            distinct = false;
            break;
        }
    }
    if (distinct) out.push_back(cp);
    return distinct;
}

}}}  // namespace fa::cpu::lr
