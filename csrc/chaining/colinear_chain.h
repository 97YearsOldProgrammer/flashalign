// Colinear chaining as minimap2's mg_lchain_dp, for DNA and spliced (cDNA) anchors.
#pragma once

#include "colinear_params.h"
#include "result.h"

#include <cstddef>
#include <vector>

namespace fa {
namespace cpu {
namespace chaining {

ChainResult chain_colinear(std::vector<Anchor> anchors,
                           ColinearChainParams params);

ChainResult chain_spliced_colinear(std::vector<Anchor> anchors,
                                   SplicedChainParams params);

int best_chain_index(const ChainResult& result);

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
