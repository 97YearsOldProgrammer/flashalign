#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include DNA results"
#endif

#include "../core/types.h"
#include "../voting/candidate_catalogue.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace dna {

struct Result : AlignResult {
  // Parallel to `supplementary`: the candidate whose block each record came
  // from, or kNullCandidate for a record that owns no block (terminal clip,
  // inversion middle). Read by the per-block MAPQ.
  std::vector<::fa::cpu::voting::CandidateId> supplementary_candidates;
};

inline void demote_unmapped(Result& result) {
  result.pos = -1;
  result.score = 0;
}

}  // namespace dna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
