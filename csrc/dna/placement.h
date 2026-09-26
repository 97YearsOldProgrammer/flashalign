#pragma once

#include "context.h"

#include <algorithm>

namespace fa::cpu::lr {

// Length of the whole-read vote window.
inline int effective_chain_segment_len(const DnaContext& context,
                                       int read_length) {
  return std::max(
      std::max(context.opts.k, 32), std::max(1, read_length));
}

}  // namespace fa::cpu::lr
