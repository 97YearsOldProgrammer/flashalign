#pragma once

#include "splice_anchor.h"
#include "../anchoring/exact_anchor_path.h"

#include <string>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Builds the realization view of a selected anchor path; a non-empty `error` means the
// bundle was rejected.
SpliceAnchorView make_splice_anchor_view(
    const ExactAnchorPath& bundle, int index_k,
    std::string* error = nullptr);

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
