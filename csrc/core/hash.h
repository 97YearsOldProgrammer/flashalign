// Integer hashes shared by the index and seeding.
#pragma once

#include <cstdint>

namespace fa { namespace cpu {

// Seed-ordering hash for syncmer sketches.
inline uint64_t sketch_order_hash(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

namespace detail {

// Open-addressing probe hash for the k-mer table and the seeding lookup caches.
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Bijection on the low `twok` bits of a 2-bit-packed k-mer (twok = 2k): xorshifts and odd
// multipliers mod 2^twok. The sharded index takes the shard from the low bits of the mix and
// stores the rest as the tag, so (shard, tag) encodes the key losslessly. Mixing folds high
// bits into the shard bits, which keeps shards balanced despite GC-skewed trailing bases.
inline uint64_t shard_mix_2k(uint64_t key, int twok) {
    const uint64_t mask = (twok >= 64) ? ~0ULL : ((1ULL << twok) - 1ULL);
    const int s1 = (twok >= 3) ? twok / 3 : 1;
    const int s2 = (twok >= 2) ? twok / 2 : 1;
    uint64_t x = key & mask;
    x = (x ^ (x >> s1)) & mask;
    x = (x * 0x9E3779B97F4A7C15ULL) & mask;  // odd: invertible mod 2^twok
    x = (x ^ (x >> s2)) & mask;
    x = (x * 0xD6E8FEB86659FD93ULL) & mask;  // odd
    x = (x ^ (x >> s1)) & mask;
    return x;
}

}  // namespace detail

}}  // namespace fa::cpu
