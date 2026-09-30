// The local chain that gates the late inversion probe of a seam or piece: the
// read's own fine seeds of the lane opposite the block, inside the fill's drop
// window, chained in that lane's frame. The seeds and their posting views are
// placement's (RetainedSeedDensity), so a window costs two lower bounds per
// seed and no key lookup. A posting counts only in the lane its orientation
// bit allows, and a key only within the whole-query harvest's occurrence gate.
// The chain gates the probe and nothing else.
#pragma once

#include "../index/index.h"

#include <cstdint>

namespace fa::cpu::lr {

class RetainedSeedDensity;

// The probe runs when the chain holds at least this many anchors.
inline constexpr int kDnaInvLocalMinAnchors = 3;
// Consecutive chained anchors differ by at most this many bases in diagonal,
inline constexpr int kDnaInvLocalDiagonal = 10;
// and by at most this many on either axis.
inline constexpr int kDnaInvLocalGap = 200;

// One fill of a block: the block's frame and the fill's origin in it.
struct DnaInvLocalChainSource {
  const RetainedSeedDensity* seeds = nullptr;
  const SeedIndex* index = nullptr;
  int chromosome = -1;
  // The block's orientation; the chain reads the other lane.
  bool block_reverse = false;
  int read_length = 0;
  int seed_length = 0;
  // The fill's first base in the block's oriented query and on the contig.
  int query_offset = 0;
  int target_offset = 0;
  // Keys over this genome-wide occurrence are skipped; 0 admits every key.
  std::uint32_t occurrence_cap = 0;
};

// Anchors in the longest chain of the opposite lane's seeds lying wholly inside
// the fill-local window [query_begin, query_end) x [target_begin, target_end).
// An anchor is (q', r): q' the seed's start in its lane, where the block's
// oriented window [a, b) is [L - b, L - a), and r its posting's start on the
// contig. Chained anchors have q' and r strictly increasing. `source` is a
// DnaInvLocalChainSource; the signature is dp/control.h DpInversionProbeGate.
int dna_inv_local_chain(const void* source, int query_begin, int query_end,
                        int target_begin, int target_end);

}  // namespace fa::cpu::lr
