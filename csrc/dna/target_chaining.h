// The all-chains lane's chaining (the ava-ont and ava-hifi presets), one chain
// call per target of a read.
#pragma once

#include "placement_chaining.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

// The all-chains lane: the whole-query exact chains of the catalogue's
// candidates, parallel to family.candidates, with no screening pass and no
// partition. The candidates on one target, a contig and strand, are chained in
// one call, over the union of their pools, and take its chains best score
// first, in catalogue order, one each. A candidate that receives no chain
// keeps the refusing status, and one on a contig below DnaContext::dual_rank
// is not chained (NotSelected). Empty when the read's seeds cannot be indexed.
std::vector<DnaPlacementCandidateChain> build_dna_target_chains(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    const std::vector<ChainWindowRetainedSeed>* forward_seeds,
    const std::vector<ChainWindowRetainedSeed>* reverse_seeds,
    const std::vector<QuerySeed>* fine_forward_seeds,
    const std::vector<QuerySeed>* fine_reverse_seeds,
    ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots);

}  // namespace fa::cpu::lr
