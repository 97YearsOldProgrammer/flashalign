// LSD radix sorts with the same output as std::stable_sort on the same key, without its
// temporary buffer or comparator calls. Sort composite keys least significant key first,
// one call per key.
#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace radix {

// The byte order that is the signed order.
inline std::uint32_t key_u32(std::int32_t value) {
  return static_cast<std::uint32_t>(value) ^ 0x80000000u;
}
inline std::uint32_t key_u32(std::uint32_t value) { return value; }

// Caller-owned scratch (ping-pong buffer and the four 256-bucket histograms), reused
// across calls.
template <class T>
struct Scratch {
  std::vector<T> tmp;
  std::vector<std::uint32_t> count;
};

// Stable ascending sort of `items` by the 32-bit key `key(item)` (int32 or uint32). A byte
// position shared by every key is skipped.
template <class T, class Key>
void stable_sort_u32(std::vector<T>& items, Scratch<T>& scratch, Key key) {
  const std::size_t n = items.size();
  if (n < 2)
    return;
  std::vector<T>& tmp = scratch.tmp;
  std::vector<std::uint32_t>& count = scratch.count;
  tmp.resize(n);
  count.assign(4 * 256, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::uint32_t k = key_u32(key(items[i]));
    ++count[k & 255u];
    ++count[256 + ((k >> 8) & 255u)];
    ++count[512 + ((k >> 16) & 255u)];
    ++count[768 + (k >> 24)];
  }
  std::vector<T>* src = &items;
  std::vector<T>* dst = &tmp;
  for (int pass = 0; pass < 4; ++pass) {
    std::uint32_t* bucket = count.data() + 256 * pass;
    bool trivial = false;
    for (int b = 0; b < 256; ++b) {
      if (bucket[b] == n) {
        trivial = true;
        break;
      }
    }
    if (trivial)
      continue;
    std::uint32_t sum = 0;
    for (int b = 0; b < 256; ++b) {
      const std::uint32_t c = bucket[b];
      bucket[b] = sum;
      sum += c;
    }
    const int shift = 8 * pass;
    const std::vector<T>& from = *src;
    std::vector<T>& to = *dst;
    for (std::size_t i = 0; i < n; ++i) {
      const std::uint32_t k = key_u32(key(from[i]));
      to[bucket[(k >> shift) & 255u]++] = from[i];
    }
    std::swap(src, dst);
  }
  if (src != &items)
    items.swap(tmp);
}

} // namespace radix
} // namespace cpu
} // namespace fa
