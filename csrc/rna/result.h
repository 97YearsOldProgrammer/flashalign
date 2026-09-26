#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include RNA results"
#endif

#include "../core/types.h"

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct Result : AlignResult {};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
