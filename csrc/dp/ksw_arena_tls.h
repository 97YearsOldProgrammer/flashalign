// Per-thread ksw arena for the DNA DP controller, created on first use and reset once per
// kernel call, like minimap2's per-thread kalloc pool. ez->cigar grows by krealloc from NULL,
// so it stays on the system heap and the controller's free(ez.cigar) is correct.
// The RNA splice kernel owns its own arena.
#pragma once

#include "ksw_arena.h"

namespace fa {
namespace cpu {
namespace dp {
namespace detail {

inline ksw_arena_t* thread_arena() {
  struct Holder {
    ksw_arena_t* arena = ksw_arena_create();
    ~Holder() { ksw_arena_destroy(arena); }
    Holder() = default;
    Holder(const Holder&) = delete;
    Holder& operator=(const Holder&) = delete;
  };
  static thread_local Holder holder;
  return holder.arena; // may be NULL: the system allocator is used
}

} // namespace detail
} // namespace dp
} // namespace cpu
} // namespace fa
