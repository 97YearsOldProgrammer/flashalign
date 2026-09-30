// The DNA vote's accumulate stage: one distinct-query vote per seed into each reference
// start (diagonal) bin, into scratch.buckets, with the used-seed bookkeeping.
#pragma once

#include "../index/format.h" // KmerPostingView, chromosome_index_for_global_pos
#include "../index/seed.h"   // QuerySeed
#include "../seeding/scratch.h" // ChainWindowPeakScratch, ChainSeedLookupCache
#include "../seeding/types.h" // ChainWindowBucketVote, DnaLongSeedView
#include "state.h" // VoteWindowState, vote_pack_key, vote_floor_div

#include <cstdint>
#include <limits>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

// Accumulates distinct-query diagonal votes into scratch.buckets.
inline void
vote_accumulate(const VoteWindowState &st, ChainWindowPeakScratch &scratch,
                const std::vector<QuerySeed> &seeds,
                const std::vector<DnaLongSeedView> *seed_view_override,
                bool retain_used_views, ChainSeedLookupCache *lookup_cache) {
  const LongReadSeedContext &ctx = st.ctx;
  auto &seed_views = scratch.seed_views;
  auto &seed_used = scratch.seed_used;
  auto &retained_seeds = scratch.retained_seeds;
  auto &buckets = scratch.buckets;
  uint32_t seed_epoch = 0; // per-window, ++ once per voting seed (stamp)
  // A rescued seed (DnaLongSeedView::rescued) passes the occurrence cap.
  auto visit_seed_view_for_vote = [&](const QuerySeed &seed,
                                      const KmerPostingView &v, bool rescued) {
    if (!retain_used_views) {
      seed_views.push_back(v);
    }
    uint8_t used = 0;
    if (!v.found() || v.count == 0) {
      if (!retain_used_views)
        seed_used.push_back(used);
      return;
    }
    if (!rescued && !vote_seed_occ_allowed(st, v)) {
      // Rejected by the occurrence cap.
      if (!retain_used_views)
        seed_used.push_back(used);
      return;
    }
    used = 1;
    if (!retain_used_views) {
      seed_used.push_back(used);
    } else {
      retained_seeds.push_back({seed, v, rescued});
    }
    int chr_idx = chromosome_index_for_global_pos(st.chr_bounds, st.n_chr,
                                                  v.positions[0]);
    if (chr_idx < 0 || chr_idx >= st.n_chr)
      return;
    uint64_t chr_lo = st.chr_bounds[static_cast<size_t>(chr_idx)];
    uint64_t chr_hi = st.chr_bounds[static_cast<size_t>(chr_idx + 1)];
    const int seed_weight = 1; // unweighted vote: one unit per distinct query seed
    // Stamp the bucket with this seed's epoch so only the seed's first posting in a diagonal
    // adds support. Postings ascend, so a seed's hits in one bin are contiguous and the
    // bucket pointer is re-probed only when the bin changes.
    const uint32_t epoch = ++seed_epoch;
    int64_t cached_key = std::numeric_limits<int64_t>::min();
    ChainWindowBucketVote *cached = nullptr;
    for (uint32_t i = 0; i < v.count; ++i) {
      const uint64_t g = static_cast<uint64_t>(v.positions[i]);
      advance_contig_cursor(st.chr_bounds, st.n_chr, g, chr_idx, chr_lo,
                            chr_hi);
      if (g < chr_lo || g >= chr_hi)
        continue;
      // A posting votes only in the lane it is compatible with: the oriented query k-mer
      // must be the reference k-mer there. Tested after the cursor advance, which every
      // posting must make.
      if (!packed_ref_orientation_compatible(seed.z, v.positions.packed_at(i),
                                             st.is_rc))
        continue;
      const int local = static_cast<int>(g - chr_lo);
      const int ref_start = local - seed.read_pos;
      const int64_t bucket_key =
          vote_pack_key(chr_idx, vote_floor_div(ref_start, st.W));
      if (bucket_key != cached_key) {
        cached = &buckets[bucket_key];
        cached_key = bucket_key;
      }
      if (cached->voter_epoch != epoch) {
        cached->support++;
        cached->center++;
        cached->score += seed_weight;
        cached->voter_epoch = epoch;
      }
    }
  };
  if (seed_view_override) {
    // Views were resolved upstream, so this loop only accumulates.
    for (const auto &sv : *seed_view_override) {
      visit_seed_view_for_vote(sv.seed, sv.view, sv.rescued);
    }
  } else {
    for (const auto &seed : seeds) {
      const KmerPostingView v = lookup_cache
                                    ? lookup_cache->lookup(*ctx.index, seed.key)
                                    : ctx.index->lookup(seed.key);
      visit_seed_view_for_vote(seed, v, false);
    }
  }
}

} // namespace lr
} // namespace cpu
} // namespace fa
