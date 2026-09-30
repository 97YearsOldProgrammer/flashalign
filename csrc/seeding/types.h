// Record types shared by long-read seeding and voting.
#pragma once

#include "../index/format.h" // KmerPostingView
#include "../index/seed.h"   // QuerySeed
#include "tie_hash.h"        // tie_locus_hash

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

// Vote table: 1 (default) is the flat open-addressed table, 0 a std::unordered_map.
#ifndef FA_DNA_LONG_VOTE_FLAT
#define FA_DNA_LONG_VOTE_FLAT 1
#endif

namespace fa {
namespace cpu {
namespace lr {

struct DnaLongSeedView {
  QuerySeed seed;
  KmerPostingView view;
  // An empty vote tile's seed over the cap, which the vote admits
  // (options/dna_profile.h kDnaTileRescueOcc).
  bool rescued = false;
};

struct DnaLongSeedBundle {
  std::vector<QuerySeed> seeds;
  std::vector<DnaLongSeedView> selected_views;
  // The vote uses the centred-rarest representative; exact candidate refinement uses the
  // deduplicated first, rarest and representative views from the same scan.
  std::vector<DnaLongSeedView> exact_refine_views;
  bool views_ready = false;

  void clear() {
    seeds.clear();
    selected_views.clear();
    exact_refine_views.clear();
    views_ready = false;
  }
};

struct ChainSyncmerOccCandidate {
  QuerySeed seed;
  KmerPostingView view;
  uint32_t occurrence = std::numeric_limits<uint32_t>::max();
  bool valid = false;
};

struct LongWindowAnchor {
  int ref_start_bin = 0;
  int ref_start_bin_width = 0;
  uint32_t median_occurrence = 0;
};

struct VotePeak {
  int read_lo = 0;
  int read_hi = 0;
  int chr = -1;
  int ref_pos = 0;
  bool is_rc = false;
  int support = 0;
  int center_support = 0;
  int vote_score = 0;
  LongWindowAnchor anchor;
  // How much further the chain's harvest window reaches below raw_ref_start - pad and
  // above raw_ref_start + L + pad, to cover the whole-read winner's per-read line
  // (vote_slope.h). Nonzero only on that winner, when its line passes the gate.
  std::int32_t harvest_below = 0;
  std::int32_t harvest_above = 0;
  // Median reference-start diagonal before whole-read projection clamps it to
  // chromosome bounds. Coarse block projection must not inherit that clamp.
  int64_t raw_ref_start = 0;
  // Oriented read range of the seeds with a posting on this peak's diagonal (bin distance
  // <= max(1, vote_radius)) in the exact-refine re-walk: a superset of the seeds the tile
  // support mask tests. Invalid (lo > hi) when the re-walk could not certify that.
  int evidence_read_lo = 0;
  int evidence_read_hi = -1;
};

static_assert(std::is_trivially_copyable_v<VotePeak>);

// The per-read line's stretch over q query bases, floor(q * b_q20 / 2^20) for a slope
// carried as b_q20 = round(b 2^20). Integer only, so no target's float contraction can
// move a window edge; both factors are 32-bit, so the product fits int64.
inline constexpr int kVoteSlopeFractionBits = 20;
inline int vote_slope_stretch(int q, std::int32_t b_q20) {
  const std::int64_t product = static_cast<std::int64_t>(q) * b_q20;
  const std::int64_t one = std::int64_t{1} << kVoteSlopeFractionBits;
  std::int64_t quotient = product / one;
  if (product % one != 0 && product < 0)
    --quotient;
  return static_cast<int>(quotient);
}

inline bool chain_peak_better(const VotePeak &a, const VotePeak &b) {
  if (a.vote_score != b.vote_score)
    return a.vote_score > b.vote_score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center_support != b.center_support) {
    return a.center_support > b.center_support;
  }
  if (a.read_lo != b.read_lo)
    return a.read_lo < b.read_lo;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.ref_pos < b.ref_pos;
}

// The DNA catalogue's order: chain_peak_better with its positional tail replaced by
// minimap2's read-seeded tie-break (tie_hash.h). A tie on vote, support, centre support and
// read_lo goes to the read-seeded hash of the peak's locus (contig, strand, diagonal bin),
// lower first; only a 64-bit collision falls back to contig and position.
inline bool chain_peak_better_seeded(const VotePeak &a, const VotePeak &b,
                                     std::uint32_t read_seed) {
  if (a.vote_score != b.vote_score)
    return a.vote_score > b.vote_score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center_support != b.center_support)
    return a.center_support > b.center_support;
  if (a.read_lo != b.read_lo)
    return a.read_lo < b.read_lo;
  const std::uint64_t a_hash =
      tie_locus_hash(read_seed, a.chr, a.is_rc, a.anchor.ref_start_bin);
  const std::uint64_t b_hash =
      tie_locus_hash(read_seed, b.chr, b.is_rc, b.anchor.ref_start_bin);
  if (a_hash != b_hash)
    return a_hash < b_hash;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.ref_pos < b.ref_pos;
}

struct ChainWindowBucketVote {
  int support = 0;
  int center = 0;
  int score = 0;
  // The last query-seed epoch that added support here, so each seed votes once per
  // diagonal. 0 means unstamped; the epoch counter restarts each window.
  uint32_t voter_epoch = 0;
};

struct ChainWindowRankedBucket {
  int chr = -1;
  int bin = 0;
  int score = 0;
  int support = 0;
  int center = 0;
};

struct ChainWindowRetainedSeed {
  QuerySeed seed;
  KmerPostingView view;
  bool rescued = false;  // DnaLongSeedView::rescued
};

// One posting captured in an accumulate pass: what vote_emit_peaks's re-walk would
// recompute, so a capturing path can replay it instead.
struct VoteHit {
  int seed_idx;
  int read_pos;
  int ref_start;
  uint32_t occurrence;
  uint32_t count;
};

// The vote drain order: score, support and centre decide; a tie goes to minimap2's
// read-seeded hash of the bucket's locus (contig, drain strand, bin), lower first, and a
// 64-bit collision to contig and bin. (chr, bin) is unique per bucket, so the order is
// strict and the drain depends only on the buckets and the read's seed.
inline bool chain_window_ranked_better(const ChainWindowRankedBucket &a,
                                       const ChainWindowRankedBucket &b,
                                       std::uint32_t read_seed, bool reverse) {
  if (a.score != b.score)
    return a.score > b.score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center != b.center)
    return a.center > b.center;
  const std::uint64_t a_hash = tie_locus_hash(read_seed, a.chr, reverse, a.bin);
  const std::uint64_t b_hash = tie_locus_hash(read_seed, b.chr, reverse, b.bin);
  if (a_hash != b_hash)
    return a_hash < b_hash;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.bin < b.bin;
}

} // namespace lr
} // namespace cpu
} // namespace fa
