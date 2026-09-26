// Writes the scored MAPQs onto the records of the emitted DNA family.
#pragma once

#include "result.h"
#include "placement_family_adapter.h"
#include "../voting/candidate_catalogue.h"

#include <vector>

namespace fa::cpu::lr {

// `supplementary_mapq` is parallel to realized.supplementary: an entry >= 0 is
// that record's own MAPQ, -1 inherits the primary's (records that own no
// block, such as terminal clips). Secondaries get 0, and an inversion middle
// gets min of its flanks' MAPQs, as minimap2's mm_set_inv_mapq.
void route_dna_mapq(const DnaPlacementFamily& catalogue,
                    ::fa::cpu::voting::CandidateId primary_candidate, int mapq,
                    const std::vector<int>& supplementary_mapq,
                    dna::Result& realized);

} // namespace fa::cpu::lr
