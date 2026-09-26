// Open-addressed int64 -> Value map for the diagonal-bin vote, covering the subset of
// std::unordered_map it uses. Entries sit in a dense array in insertion order; the index
// holds generation-tagged entry numbers, so clear() is O(1). Iteration is in insertion order.
#pragma once

#include "hash.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

template <class Value> class FlatInt64Map {
public:
  struct value_type {
    int64_t first;
    Value second;
  };
  using iterator = value_type *;
  using const_iterator = const value_type *;

  FlatInt64Map() = default;

  iterator begin() { return entries_.data(); }
  iterator end() { return entries_.data() + entries_.size(); }
  const_iterator begin() const { return entries_.data(); }
  const_iterator end() const { return entries_.data() + entries_.size(); }

  size_t size() const { return entries_.size(); }
  bool empty() const { return entries_.empty(); }
  // Index capacity, as unordered_map::bucket_count.
  size_t bucket_count() const { return slots_.size(); }

  // Bumps the generation instead of clearing the index, except on wrap-around.
  void clear() {
    entries_.clear();
    if (++generation_ == 0) {
      std::fill(slots_.begin(), slots_.end(), uint64_t{0});
      generation_ = 1;
    }
  }

  // Room for n entries at load <= 0.5 without a rehash.
  void reserve(size_t n) {
    const size_t want = pow2_at_least(n == 0 ? 1 : n * 2);
    if (want > slots_.size())
      rebuild_index(want);
    entries_.reserve(n);
  }

  // Get-or-insert a value-initialized Value. The reference is valid until the next insert.
  Value &operator[](int64_t key) {
    ensure_index();
    size_t mask = slots_.size() - 1;
    size_t slot = hash(key) & mask;
    for (;;) {
      const uint64_t s = slots_[slot];
      if (!slot_live(s))
        break; // empty slot -> insert below
      if (entries_[slot_index(s)].first == key)
        return entries_[slot_index(s)].second;
      slot = (slot + 1) & mask;
    }
    if ((entries_.size() + 1) * 2 > slots_.size()) {
      rebuild_index(slots_.size() * 2);
      mask = slots_.size() - 1;
      slot = hash(key) & mask;
      while (slot_live(slots_[slot]))
        slot = (slot + 1) & mask;
    }
    const uint32_t new_idx = static_cast<uint32_t>(entries_.size());
    entries_.push_back(value_type{key, Value{}});
    slots_[slot] = pack(generation_, new_idx);
    return entries_.back().second;
  }

  iterator find(int64_t key) {
    const int64_t i = find_index(key);
    return i < 0 ? end() : entries_.data() + i;
  }
  const_iterator find(int64_t key) const {
    const int64_t i = find_index(key);
    return i < 0 ? end() : entries_.data() + i;
  }

  void swap(FlatInt64Map &o) {
    entries_.swap(o.entries_);
    slots_.swap(o.slots_);
    std::swap(generation_, o.generation_);
  }

private:
  static size_t pow2_at_least(size_t x) {
    size_t c = 1;
    while (c < x)
      c <<= 1;
    return c;
  }
  // Keys are (chr << 32) | bin. The default hash folds chr into the low bits, so contiguous
  // bins land in consecutive slots while equal bins on different chromosomes do not all
  // collide. FA_DNA_LONG_VOTE_HASH selects splitmix64 (0) or the identity (2) instead.
#ifndef FA_DNA_LONG_VOTE_HASH
#define FA_DNA_LONG_VOTE_HASH 1
#endif
  static uint64_t hash(int64_t key) {
#if FA_DNA_LONG_VOTE_HASH == 0
    return fa::cpu::detail::splitmix64(static_cast<uint64_t>(key));
#elif FA_DNA_LONG_VOTE_HASH == 2
    return static_cast<uint64_t>(key);
#else
    const uint64_t k = static_cast<uint64_t>(key);
    return k ^ (k >> 32);
#endif
  }
  static uint64_t pack(uint32_t gen, uint32_t idx) {
    return (static_cast<uint64_t>(gen) << 32) | idx;
  }
  bool slot_live(uint64_t s) const {
    return static_cast<uint32_t>(s >> 32) == generation_;
  }
  static uint32_t slot_index(uint64_t s) {
    return static_cast<uint32_t>(s & 0xFFFFFFFFu);
  }
  void ensure_index() {
    if (slots_.empty())
      rebuild_index(16);
  }

  int64_t find_index(int64_t key) const {
    if (slots_.empty())
      return -1;
    const size_t mask = slots_.size() - 1;
    size_t slot = hash(key) & mask;
    for (;;) {
      const uint64_t s = slots_[slot];
      if (!slot_live(s))
        return -1;
      if (entries_[slot_index(s)].first == key)
        return static_cast<int64_t>(slot_index(s));
      slot = (slot + 1) & mask;
    }
  }

  // Rebuilds the index at new_cap (a power of two, at least 16) from the live entries.
  void rebuild_index(size_t new_cap) {
    new_cap = pow2_at_least(new_cap < 16 ? 16 : new_cap);
    slots_.assign(new_cap, uint64_t{0});
    if (generation_ == 0)
      generation_ = 1; // 0 marks an empty slot
    const size_t mask = new_cap - 1;
    for (uint32_t i = 0; i < entries_.size(); ++i) {
      size_t slot = hash(entries_[i].first) & mask;
      while (slot_live(slots_[slot]))
        slot = (slot + 1) & mask;
      slots_[slot] = pack(generation_, i);
    }
  }

  std::vector<value_type> entries_; // dense, insertion order
  std::vector<uint64_t> slots_;     // gen<<32 | entry index; live iff gen match
  uint32_t generation_ = 1;
};

template <class Value>
inline void swap(FlatInt64Map<Value> &a, FlatInt64Map<Value> &b) {
  a.swap(b);
}

} // namespace lr
} // namespace cpu
} // namespace fa
