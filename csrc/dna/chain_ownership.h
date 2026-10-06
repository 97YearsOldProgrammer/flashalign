// Which whole-query chain owns which stretch of the read, decided on the
// chains' anchors (chain_ownership.cpp). Reads anchors and scores only.
#pragma once

#include "../chaining/anchor.h"

#include <vector>

namespace fa::cpu::lr {

// minimap2's -n: a chain owns or attaches only with at least this many
// anchors, and a block keeps at least this many.
inline constexpr int kDnaOwnerMinAnchors = 3;

// One whole-query chain: a candidate's primary (path -1) or sibling j.
struct DnaOwnershipPath {
  int candidate = -1;
  int path = -1;
  int contig = -1;
  bool reverse = false;
  const std::vector<chaining::Anchor>* anchors = nullptr;
  int score = 0;
};

// One block of the read: item `item`'s anchors [a0, a1) in chain order, the
// bounds [forward_begin, forward_end) the blocks tile the read with, its score
// and anchor count, and minimap2's subsc and n_sub over the chains of the
// item's candidate attached to its owner.
struct DnaOwnershipBlock {
  int item = 0;
  int a0 = 0;
  int a1 = 0;
  int forward_begin = 0;
  int forward_end = 0;
  int score = 0;
  int kept = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct DnaOwnershipSelection {
  // The paths past the floors (at least kDnaOwnerMinAnchors anchors, score
  // >= min_chain_score), input order; DnaOwnershipBlock::item indexes it.
  std::vector<DnaOwnershipPath> items;
  // In query order; empty when no chain owns.
  std::vector<DnaOwnershipBlock> blocks;
};

// `paths` in catalogue order, each candidate's primary first, then its
// siblings by index.
DnaOwnershipSelection select_chain_owners(
    int read_length, int seed_length, int min_chain_score,
    const std::vector<DnaOwnershipPath>& paths);

}  // namespace fa::cpu::lr
