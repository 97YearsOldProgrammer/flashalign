// The anchor-pool builder, internal to placement_chaining.cpp.
#pragma once

#include "placement_chaining.h"
#include "retained_seed_density.h"
#include "../chaining/anchor.h"
#include "../index/format.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {
namespace internal {

// Appends one anchor per posting in `interval` that passes the geometry tests
// and the orientation test for `reverse_lane`; no base is read.
// `chromosome_base` is the contig's offset in the flattened reference and
// `chromosome_length` its length, so postings off the contig are dropped. An
// anchor is flagged ANCHOR_TANDEM when its key has another posting on this
// contig within `tandem_window` bases (<= 0 disables). `skip_own_diagonal`
// drops the posting at the seed's own read position: the contig is the query
// read itself (DnaContext::self_contig) and the lane is forward. `query_length`
// bounds the query span.
void append_interval_anchors(
    const KmerPostingIntervalView& interval,
    const RetainedSeedRef& seed,
    std::uint64_t chromosome_base,
    int chromosome_length,
    int main_diagonal,
    int seed_length,
    int diagonal_band,
    int tandem_window,
    bool reverse_lane,
    bool skip_own_diagonal,
    int query_length,
    std::vector<chaining::Anchor>& anchors,
    DnaPlacementCandidateChain& record);

}  // namespace internal
}  // namespace fa::cpu::lr
