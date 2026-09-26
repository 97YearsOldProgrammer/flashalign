// The per-window diagonal-bin vote (window_anchor_peaks): resolve the window's seeds,
// accumulate diagonal votes, rank, and drain the heap into VotePeak anchors.
#pragma once

#include "../index/format.h" // KmerPostingView, chromosome_index_for_global_pos
#include "../index/seed.h"   // QuerySeed
#include "../seeding/context.h" // LongReadSeedContext, occ policy
#include "../seeding/scratch.h" // ChainWindowPeakScratch, ChainSeedLookupCache
#include "../seeding/select.h"  // limit_query_seeds_with_lookup_cache
#include "../seeding/syncmer.h" // closed-syncmer extraction
#include "../seeding/types.h" // VotePeak, ChainWindow* records, comparators, DnaLongSeed*
#include "accumulate.h" // vote_accumulate
#include "extract.h" // vote_resolve_seeds
#include "extract_peaks.h" // best_peak
#include "rank.h" // vote_rank
#include "refine.h" // vote_emit_peaks
#include "state.h" // VoteWindowState

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

struct WindowAnchorPeakParams {
  ChainSeedLookupCache *lookup_cache = nullptr;
  const std::vector<QuerySeed> *seed_override = nullptr;
  bool seed_override_prelimited = true;
  const std::vector<DnaLongSeedView> *seed_view_override = nullptr;
  // The read's tie seed (tie_read_seed) for the drain heap's tie-break.
  std::uint32_t tie_seed = 0;
};

inline std::vector<VotePeak> window_anchor_peaks(
    const LongReadSeedContext &ctx, const uint8_t *query_enc, int span,
    bool is_rc, ChainWindowPeakScratch &scratch,
    const WindowAnchorPeakParams &params = {}) {
  std::vector<VotePeak> out;
  const int limit = std::max(1, ctx.chain_max_candidates_per_window);
  out.reserve(static_cast<size_t>(limit));
  scratch.clear_for_window();
  struct PendingExactRefineReset {
    std::vector<DnaLongSeedView>& views;
    ~PendingExactRefineReset() { views.clear(); }
  } pending_exact_refine_reset{scratch.pending_exact_refine_views};
  const bool has_resolved_seed_input =
      params.seed_view_override != nullptr || params.seed_override != nullptr;
  if (span < ctx.k || ctx.index->empty() ||
      (!query_enc && !has_resolved_seed_input))
    return out;

  auto &seeds = scratch.seeds;
  size_t active_seed_count = 0;
  bool direct_seed_view_input = false;
  const std::vector<DnaLongSeedView> *seed_view_override =
      params.seed_view_override;
  vote_resolve_seeds(ctx, query_enc, span, scratch, params.lookup_cache,
                     params.seed_override, params.seed_override_prelimited,
                     seed_view_override, active_seed_count,
                     direct_seed_view_input);
  if (active_seed_count == 0)
    return out;
  auto &seed_views = scratch.seed_views;
  auto &seed_used = scratch.seed_used;
  auto &retained_seeds = scratch.retained_seeds;
  const bool retain_used_views = direct_seed_view_input;
  if (!retain_used_views) {
    seed_views.reserve(active_seed_count);
    seed_used.reserve(active_seed_count);
  } else {
    retained_seeds.reserve(active_seed_count);
  }

  const uint64_t *chr_bounds = ctx.index->chrom_offsets_data();
  const int n_chr = static_cast<int>(ctx.chr_names->size());
  if (!chr_bounds || n_chr <= 0)
    return out;

  const int W = effective_vote_diag_bin_width(ctx, span);
  const int vote_radius = (span >= 128) ? 1 : 0;
  auto &buckets = scratch.buckets;
  {
    const size_t bucket_target = std::max<size_t>(8, active_seed_count * 2);
    if (buckets.bucket_count() < bucket_target) {
      buckets.reserve(bucket_target);
    }
  }
  const int min_support = std::max(1, ctx.min_support);
  const VoteWindowState st{ctx,
                           chr_bounds,
                           n_chr,
                           W,
                           vote_radius,
                           span,
                           is_rc,
                           limit,
                           min_support,
                           params.tie_seed};
  vote_accumulate(st, scratch, seeds, seed_view_override, retain_used_views,
                  params.lookup_cache);
  if (buckets.empty())
    return out;

  vote_rank(st, scratch);
  vote_emit_peaks(st, scratch, out, seeds, retain_used_views,
                  active_seed_count);
  return out;
}

} // namespace lr
} // namespace cpu
} // namespace fa
