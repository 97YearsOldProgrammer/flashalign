// The option resolver shared by the CLI and the Python binding. Layers, lowest
// first: defaults, the named preset, the index's (k, s), explicit overrides,
// derived rules, then validation. req.preset is already the final preset name.
#pragma once

#include "types.h"

namespace fa {
namespace cpu {
namespace options {

// Throws std::invalid_argument on an unknown preset or occurrence policy, an
// explicit k or s that contradicts a prebuilt index, or an invalid value.
ResolvedMapOptions resolve_options(const ResolveRequest &req);

} // namespace options
} // namespace cpu
} // namespace fa
