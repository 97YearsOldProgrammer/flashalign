// minimap2's read-seeded tie-break. When two candidate loci tie on every evidence key, the
// DNA vote and the RNA catalogue order them by a hash seeded by the read, as minimap2 does,
// so on a diploid reference the haplotype whose contig sorts first does not win every tie.
// minimap2:
//
//   map.c mm_map_frag  hash  = X31(qname) (0 without a name);
//                      hash ^= Wang(qlen) + Wang(opt->seed);   opt->seed = 11
//                      hash  = Wang(hash)
//   hit.c mm_gen_regs  h = hash64((hash64(a.x) + hash64(a.y)) ^ hash)
//
// with a the chain's first anchor. tie_read_seed() is mm_map_frag's per-read seed and
// tie_locus_hash() mixes it with one locus key (strand, contig and position) where minimap2
// mixes an anchor. Comparators consult it only after every evidence key ties and fall back
// to contig and position on a 64-bit collision, so they stay strict weak orders. The
// result depends only on the read name, the read length and the locus.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace fa {
namespace cpu {
namespace lr {

// minimap2's opt->seed default (options.c mm_mapopt_init).
inline constexpr std::uint32_t kTieHashSeed = 11;

// khash.h __ac_X31_hash_string over the read name. Bytes are read unsigned, so the value
// is the same on every target (minimap2 reads plain char, which differs above 0x7f).
inline std::uint32_t tie_name_hash(std::string_view name) noexcept {
  if (name.empty())
    return 0;
  std::uint32_t hash = static_cast<unsigned char>(name[0]);
  for (std::size_t i = 1; i < name.size(); ++i)
    hash = (hash << 5) - hash + static_cast<unsigned char>(name[i]);
  return hash;
}

// khash.h __ac_Wang_hash.
inline std::uint32_t tie_wang_hash(std::uint32_t key) noexcept {
  key += ~(key << 15);
  key ^= (key >> 10);
  key += (key << 3);
  key ^= (key >> 6);
  key += ~(key << 11);
  key ^= (key >> 16);
  return key;
}

// hit.c hash64 (Thomas Wang's 64-bit integer mix).
inline std::uint64_t tie_hash64(std::uint64_t key) noexcept {
  key = (~key + (key << 21));
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8));
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4));
  key = key ^ key >> 28;
  key = (key + (key << 31));
  return key;
}

// mm_map_frag's per-read seed. `name_hash` is tie_name_hash(read name), or 0 without a
// name, as in minimap2; `read_len` is minimap2's qlen_sum for a single-segment read.
inline std::uint32_t tie_read_seed(std::uint32_t name_hash,
                                   int read_len) noexcept {
  std::uint32_t hash = name_hash;
  hash ^= tie_wang_hash(static_cast<std::uint32_t>(read_len)) +
          tie_wang_hash(kTieHashSeed);
  return tie_wang_hash(hash);
}

// One locus's tie rank under the read's seed; lower sorts first. The key packs strand
// (bit 63), contig (bits 32-62) and the position (low 32 bits: the DNA diagonal bin or the
// RNA window's reference start), so distinct loci give distinct keys; the mix is
// hash64(hash64(key) ^ seed).
inline std::uint64_t tie_locus_hash(std::uint32_t read_seed, int contig,
                                    bool reverse, std::int64_t bin) noexcept {
  const std::uint64_t key =
      (static_cast<std::uint64_t>(reverse ? 1u : 0u) << 63) |
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(contig) &
                                  0x7fffffffu)
       << 32) |
      static_cast<std::uint64_t>(static_cast<std::uint32_t>(bin));
  return tie_hash64(tie_hash64(key) ^ static_cast<std::uint64_t>(read_seed));
}

} // namespace lr
} // namespace cpu
} // namespace fa
