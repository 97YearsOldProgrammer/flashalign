// Per-worker mapping scratch. The batch driver owns one per worker thread and
// reuses it for every read of a mini-batch, like minimap2's mm_tbuf_t.
#pragma once

#include "../seeding/scratch.h" // ChainAnchorScratch

#include <cstddef>

namespace fa { namespace cpu { namespace lr {

struct DnaWorkerScratch {
  ChainAnchorScratch chain;

  // Clears per-read state, then trims retained capacity for the next read.
  void clear_for_read(std::size_t read_len) {
    chain.clear_for_read();
    chain.release_excess_for_read(read_len);
  }
};

}}} // namespace fa::cpu::lr
