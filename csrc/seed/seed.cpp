// Query-side seeds: closed-syncmer extraction under the occurrence cap, low-occurrence selection, tile rescue, the window slicer, minimap2's read-seeded tie hash.
#include "seed.h"
#include "scratch.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**** posting_density ****/

namespace fa::cpu {

bool posting_density_boundaries_valid(
    const std::uint64_t* offsets,
    std::uint64_t chromosome_count) noexcept {
  if (offsets == nullptr || chromosome_count == 0 || offsets[0] != 0)
    return false;
  for (std::uint64_t chromosome = 0; chromosome < chromosome_count;
       ++chromosome) {
    if (offsets[chromosome + 1] < offsets[chromosome])
      return false;
  }
  return true;
}

}  // namespace fa::cpu

/**** context ****/

namespace fa { namespace cpu { namespace lr {

static std::string long_occ_policy_token(std::string value) {
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char c) {
            const char lower = static_cast<char>(std::tolower(c));
            return lower == '-' ? '_' : lower;
        });
    return value;
}

std::string long_occ_policy_name(LongOccPolicy policy) {
    switch (policy) {
        case LongOccPolicy::Fixed:
            return "fixed";
        case LongOccPolicy::Platform:
            return "platform";
        case LongOccPolicy::Off:
        default:
            return "off";
    }
}

LongOccPolicy long_occ_policy_from_string(
    const std::string& value,
    LongOccPolicy fallback
) {
    const std::string mode = long_occ_policy_token(value);
    if (mode == "off" || mode == "none" || mode == "0") {
        return LongOccPolicy::Off;
    }
    if (mode == "fixed" || mode == "cap" || mode == "cap_fixed") {
        return LongOccPolicy::Fixed;
    }
    if (mode == "platform" || mode == "platform_fixed") {
        return LongOccPolicy::Platform;
    }
    return fallback;
}

bool long_occ_policy_string_valid(const std::string& value) {
    const std::string mode = long_occ_policy_token(value);
    return mode == "off" || mode == "none" || mode == "0" ||
           mode == "fixed" || mode == "cap" || mode == "cap_fixed" ||
           mode == "platform" || mode == "platform_fixed";
}

bool LongReadSeedContext::syncmer_occ_aware_active() const {
    return syncmer_occ_aware_enabled;
}

int effective_vote_diag_bin_width(
    const LongReadSeedContext& ctx, int read_span) {
    const int base_width = std::max(1, ctx.vote_diag_bin_width);
    if (!ctx.vote_diag_bin_width_adaptive || read_span <= 0)
        return base_width;
    // read_span > 0 makes the truncating division a floor; integer arithmetic keeps bin
    // boundaries identical on every target.
    const int drift = ctx.vote_diag_slope_den > 0
                          ? read_span / ctx.vote_diag_slope_den
                          : 0;
    return std::max(base_width, std::min(drift, ctx.vote_diag_width_max));
}

}}}  // namespace fa::cpu::lr

/**** syncmer ****/

namespace fa { namespace cpu { namespace lr {

bool chain_syncmer_occ_candidate_better(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b
) {
    if (!a.valid) return false;
    if (!b.valid) return true;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    return a.seed.read_pos < b.seed.read_pos;
}

bool chain_syncmer_occ_candidate_better_centred(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b,
    int tile_centre
) {
    if (!a.valid) return false;
    if (!b.valid) return true;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    const int a_distance = std::abs(a.seed.read_pos - tile_centre);
    const int b_distance = std::abs(b.seed.read_pos - tile_centre);
    if (a_distance != b_distance) return a_distance < b_distance;
    return a.seed.read_pos < b.seed.read_pos;
}

int syncmer_tile_centre(int tile, int n_kmers, int n_tiles) {
    if (n_tiles <= 1) return std::max(0, n_kmers - 1) / 2;
    const int64_t kmers = static_cast<int64_t>(std::max(1, n_kmers));
    const int64_t tiles = static_cast<int64_t>(n_tiles);
    const int64_t lo =
        (static_cast<int64_t>(tile) * kmers + tiles - 1) / tiles;
    const int64_t hi =
        (static_cast<int64_t>(tile + 1) * kmers + tiles - 1) / tiles - 1;
    return static_cast<int>((lo + std::max(lo, hi)) / 2);
}

int vote_seed_occurrence_stratum(uint32_t occurrence) {
    if (occurrence <= 1) return 0;
    const int floor_log2 =
        63 - __builtin_clzll(static_cast<unsigned long long>(occurrence));
    return std::min(floor_log2, kVoteSeedOccurrenceStrata - 1);
}

bool vote_seed_minimum_is_representative(
    const std::array<uint32_t, kVoteSeedOccurrenceStrata>& strata,
    uint32_t admitted,
    uint32_t minimum_occurrence
) {
    if (admitted == 0 || minimum_occurrence == 0) return false;
    const int minimum_stratum =
        vote_seed_occurrence_stratum(minimum_occurrence);
    uint32_t near = strata[static_cast<size_t>(minimum_stratum)];
    if (minimum_stratum + 1 < kVoteSeedOccurrenceStrata) {
        near += strata[static_cast<size_t>(minimum_stratum + 1)];
    }
    // The rare minimum is representative only when at least half of the tile shares its
    // floor(log2 occurrence) band or the next higher one.
    return static_cast<uint64_t>(near) * 2u >= admitted;
}

VoteTileRescue::VoteTileRescue(const LongReadSeedContext& ctx, int span)
    : max_occ_(span >= kDnaTileRescueMinLen
                   ? static_cast<uint32_t>(std::max(0, ctx.tile_rescue_occ))
                   : 0u) {}

void VoteTileRescue::offer(const ChainSyncmerOccCandidate& cand, int tile_centre) {
    // A seed with no postings keeps occurrence UINT32_MAX.
    if (cand.valid || cand.occurrence > max_occ_) return;
    ChainSyncmerOccCandidate over = cand;
    over.valid = true;
    if (chain_syncmer_occ_candidate_better_centred(over, best_, tile_centre))
        best_ = over;
}

bool VoteTileRescue::close(uint32_t admitted, DnaLongSeedView& out) {
    const ChainSyncmerOccCandidate best = best_;
    best_ = ChainSyncmerOccCandidate{};
    if (admitted != 0 || !best.valid ||
        spent_ + best.occurrence > kDnaTileRescueBudget)
        return false;
    spent_ += best.occurrence;
    out = DnaLongSeedView{best.seed, best.view, true};
    return true;
}

static void extract_chain_closed_syncmer_seeds_into_downsample(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    int downsample,
    std::vector<QuerySeed>& out
) {
    ClosedSyncmerConfig cfg;
    cfg.k = ctx.k;
    cfg.s = ctx.syncmer_s;
    cfg.downsample = std::max(1, downsample);
    extract_closed_syncmer_query_seeds_into(query_enc, span, cfg, out);
}

void extract_chain_closed_syncmer_seeds_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    std::vector<QuerySeed>& out
) {
    extract_chain_closed_syncmer_seeds_into_downsample(
        ctx, query_enc, span, std::max(1, ctx.syncmer_downsample), out);
}

// One vote seed per tile of max_query_seeds_per_strand equal tiles (every
// seed when it is 0).
static bool extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
) {
    out.clear();
    (void)scratch;
    const int k = ctx.k;
    const SeedIndex& index = *ctx.index;
    const LongOccPolicyConfig& cfg = ctx.occ_policy;
    if (!ctx.syncmer_occ_aware_active() || !lookup_cache ||
        input.empty() || !htable_valid_k(k) || span < k) {
        return false;
    }
    const int max_seeds = ctx.max_query_seeds_per_strand;
    const int n_kmers = std::max(1, span - k + 1);
    const int n_tiles = std::max(
        1,
        max_seeds > 0
            ? std::min(max_seeds, n_kmers)
            : n_kmers);
    auto make_candidate_from_view =
        [&](const QuerySeed& seed, const KmerPostingView& view) {
        ChainSyncmerOccCandidate cand;
        cand.seed = seed;
        cand.view = view;
        const uint32_t occurrence = vote_seed_occurrence(ctx, view);
        if (!view.found() || view.count == 0 || occurrence == 0) {
            return cand;
        }
        cand.occurrence = occurrence;
        cand.valid = occurrence_allowed_by_long_occ_policy(occurrence, cfg);
        return cand;
    };

    // `seed_views`, when captured (ChainSeedLookupCache::warm), holds the cached view for
    // each position of `input`, read by position instead of by a probe.
    auto make_candidate = [&](const QuerySeed& seed, size_t position) {
        KmerPostingView view;
        if (seed_views != nullptr) {
            view = seed_views[position];
        } else {
            view = lookup_cache->lookup(index, seed.key);
        }
        return make_candidate_from_view(seed, view);
    };

    out.selected_views.reserve(static_cast<size_t>(n_tiles));
    out.exact_refine_views.reserve(static_cast<size_t>(3 * n_tiles));
    out.seeds.reserve(static_cast<size_t>(n_tiles));

    auto emit_tile = [&](const ChainSyncmerOccCandidate& first,
                         const ChainSyncmerOccCandidate& rarest,
                         const ChainSyncmerOccCandidate& centred,
                         const std::array<uint32_t,
                                          kVoteSeedOccurrenceStrata>& strata,
                         uint32_t admitted) {
      if (!first.valid)
        return;
      const ChainSyncmerOccCandidate& representative =
          vote_seed_minimum_is_representative(
              strata, admitted, rarest.occurrence)
              ? centred
              : first;
      out.selected_views.push_back(
          {representative.seed, representative.view});
      out.seeds.push_back(representative.seed);
      out.exact_refine_views.push_back({first.seed, first.view});
      if (rarest.valid && rarest.seed.read_pos != first.seed.read_pos) {
        out.exact_refine_views.push_back({rarest.seed, rarest.view});
      }
      if (representative.seed.read_pos != first.seed.read_pos &&
          representative.seed.read_pos != rarest.seed.read_pos) {
        out.exact_refine_views.push_back(
            {representative.seed, representative.view});
      }
    };

    // A tile that admits no seed may vote with a rescued one, also an exact-refine view.
    VoteTileRescue rescue(ctx, span);
    auto emit_rescue = [&](uint32_t admitted) {
      DnaLongSeedView rescued;
      if (!rescue.close(admitted, rescued))
        return;
      out.selected_views.push_back(rescued);
      out.seeds.push_back(rescued.seed);
      out.exact_refine_views.push_back(rescued);
    };

    // Streaming per-tile scan: each tile keeps its first, rarest and centred admitted seeds,
    // from which emit_tile picks the vote representative and the exact-refinement views.
    const int64_t tiles64 = static_cast<int64_t>(n_tiles);
    const int64_t kmers64 = static_cast<int64_t>(std::max(1, n_kmers));
    int tile = 0;
    ChainSyncmerOccCandidate first_admitted;
    ChainSyncmerOccCandidate rarest_admitted;
    ChainSyncmerOccCandidate centred_admitted;
    std::array<uint32_t, kVoteSeedOccurrenceStrata> tile_strata{};
    uint32_t tile_admitted = 0;
    int tile_centre = syncmer_tile_centre(0, n_kmers, n_tiles);
#ifndef NDEBUG
    int prev_pos = INT_MIN;
#endif
    for (size_t position = 0; position < input.size(); ++position) {
        const QuerySeed& seed = input[position];
#ifndef NDEBUG
        assert(seed.read_pos >= prev_pos &&
               "occ-stream selector requires read_pos-sorted input");
        prev_pos = seed.read_pos;
#endif
        const int cp =
            std::max(0, std::min(std::max(0, n_kmers - 1), seed.read_pos));
        const int prev_tile = tile;
        while (static_cast<int64_t>(cp) * tiles64 >=
               static_cast<int64_t>(tile + 1) * kmers64) {
            ++tile;
        }
        if (tile != prev_tile) {
          emit_tile(first_admitted, rarest_admitted, centred_admitted,
                    tile_strata, tile_admitted);
          emit_rescue(tile_admitted);
          first_admitted = ChainSyncmerOccCandidate{};
          rarest_admitted = ChainSyncmerOccCandidate{};
          centred_admitted = ChainSyncmerOccCandidate{};
          tile_strata.fill(0);
          tile_admitted = 0;
          tile_centre = syncmer_tile_centre(tile, n_kmers, n_tiles);
        }
        const ChainSyncmerOccCandidate cand = make_candidate(seed, position);
        rescue.offer(cand, tile_centre);
        if (cand.valid) {
          if (!first_admitted.valid)
            first_admitted = cand;
          if (chain_syncmer_occ_candidate_better(cand, rarest_admitted)) {
            rarest_admitted = cand;
          }
          if (chain_syncmer_occ_candidate_better_centred(
                  cand, centred_admitted, tile_centre)) {
            centred_admitted = cand;
          }
          ++tile_strata[static_cast<size_t>(
              vote_seed_occurrence_stratum(cand.occurrence))];
          ++tile_admitted;
        }
    }
    emit_tile(first_admitted, rarest_admitted, centred_admitted, tile_strata,
              tile_admitted);
    emit_rescue(tile_admitted);
    out.views_ready = true;
    return true;
}

bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views
) {
    const bool ready =
        extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
            ctx, input, span, lookup_cache, out, scratch, seed_views);
    if (!ready || !ctx.nested_vote_seeds ||
        ctx.max_query_seeds_per_strand <= kVoteSeedNestBase)
        return ready;
    LongReadSeedContext base_ctx = ctx;
    base_ctx.max_query_seeds_per_strand = kVoteSeedNestBase;
    DnaLongSeedBundle base;
    extract_chain_closed_syncmer_occ_aware_seed_bundle_uniform_into(
        base_ctx, input, span, lookup_cache, base, scratch, seed_views);

    const auto by_position = [](const DnaLongSeedView& a,
                                const DnaLongSeedView& b) {
        return a.seed.read_pos < b.seed.read_pos;
    };
    out.exact_refine_views.insert(out.exact_refine_views.end(),
                                  base.exact_refine_views.begin(),
                                  base.exact_refine_views.end());
    std::sort(out.exact_refine_views.begin(), out.exact_refine_views.end(),
              by_position);
    out.exact_refine_views.erase(
        std::unique(out.exact_refine_views.begin(),
                    out.exact_refine_views.end(),
                    [](const DnaLongSeedView& a, const DnaLongSeedView& b) {
                        return a.seed.read_pos == b.seed.read_pos;
                    }),
        out.exact_refine_views.end());

    const std::size_t target = out.selected_views.size();
    std::vector<DnaLongSeedView> selected = std::move(base.selected_views);
    std::set<int> kept;
    uint64_t rescued_occurrence = 0;
    for (const DnaLongSeedView& v : selected) {
        kept.insert(v.seed.read_pos);
        if (v.rescued) rescued_occurrence += vote_seed_occurrence(ctx, v.view);
    }
    const auto distance_to_kept = [&](int read_pos) {
        const auto hi = kept.lower_bound(read_pos);
        int distance = span;
        if (hi != kept.end()) distance = *hi - read_pos;
        if (hi != kept.begin())
            distance = std::min(distance, read_pos - *std::prev(hi));
        return distance;
    };
    // A max-heap; distances only shrink, so a popped seed whose distance has
    // shrunk since its push goes back with the new one.
    struct Fill {
        int distance;
        uint32_t occurrence;
        int read_pos;
        std::size_t index;
    };
    const auto after = [](const Fill& a, const Fill& b) {
        if (a.distance != b.distance) return a.distance < b.distance;
        if (a.occurrence != b.occurrence) return a.occurrence > b.occurrence;
        return a.read_pos > b.read_pos;
    };
    std::vector<Fill> heap;
    if (selected.size() < target) {
        for (std::size_t i = 0; i < out.selected_views.size(); ++i) {
            const DnaLongSeedView& v = out.selected_views[i];
            if (kept.count(v.seed.read_pos)) continue;
            heap.push_back({distance_to_kept(v.seed.read_pos),
                            vote_seed_occurrence(ctx, v.view),
                            v.seed.read_pos, i});
        }
        std::make_heap(heap.begin(), heap.end(), after);
    }
    while (selected.size() < target && !heap.empty()) {
        std::pop_heap(heap.begin(), heap.end(), after);
        Fill next = heap.back();
        heap.pop_back();
        const DnaLongSeedView& v = out.selected_views[next.index];
        if (v.rescued &&
            rescued_occurrence + vote_seed_occurrence(ctx, v.view) >
                kDnaTileRescueBudget)
            continue;
        const int distance = distance_to_kept(next.read_pos);
        if (distance != next.distance) {
            next.distance = distance;
            heap.push_back(next);
            std::push_heap(heap.begin(), heap.end(), after);
            continue;
        }
        selected.push_back(v);
        kept.insert(v.seed.read_pos);
        if (v.rescued) rescued_occurrence += vote_seed_occurrence(ctx, v.view);
    }
    std::sort(selected.begin(), selected.end(), by_position);
    out.selected_views = std::move(selected);
    out.seeds.clear();
    for (const DnaLongSeedView& v : out.selected_views)
        out.seeds.push_back(v.seed);
    return true;
}

bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch
) {
    out.clear();
    if (!query_enc || !htable_valid_k(ctx.k) || span < ctx.k) return false;
    scratch.extracted_seeds.clear();
    extract_chain_closed_syncmer_seeds_into(
        ctx, query_enc, span, scratch.extracted_seeds);
    return extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
        ctx, scratch.extracted_seeds, span, lookup_cache, out, scratch);
}

bool select_chain_closed_syncmer_occ_aware_seed_bundle_from_seeds_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views
) {
    return extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
        ctx, input, span, lookup_cache, out, scratch, seed_views);
}

void rebase_syncmer_stream_window(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample,
    std::vector<QuerySeed>& out
) {
    out.clear();
    if (end - start < k) return;
    const int last_pos = end - k;
    auto it = std::lower_bound(
        seeds.begin(),
        seeds.end(),
        start,
        [](const QuerySeed& seed, int pos) {
            return seed.read_pos < pos;
        });
    for (; it != seeds.end() && it->read_pos <= last_pos; ++it) {
        if (downsample > 1 &&
            sketch_order_hash(it->key) %
                    static_cast<uint64_t>(downsample) != 0) {
            continue;
        }
        out.push_back(
            QuerySeed{it->key, it->read_pos - start, it->z});
    }
}

bool syncmer_stream_window_is_identity(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample
) {
    if (downsample > 1 || start != 0 || end - start < k) return false;
    if (seeds.empty()) return true;
    return seeds.front().read_pos >= start &&
           seeds.back().read_pos <= end - k;
}

}}}  // namespace fa::cpu::lr

/**** select ****/

namespace fa { namespace cpu {

static QuerySeedOccurrenceEvidence query_seed_occurrence(
    const IndexView& index,
    const QuerySeed& seed,
    int original_index = 0
) {
    QuerySeedOccurrenceEvidence evidence;
    evidence.seed = seed;
    const KmerPostingView view = index.lookup(seed.key);
    evidence.present = view.found();
    evidence.retained = view.found() && view.count > 0;
    evidence.found = evidence.retained;
    evidence.occurrence = evidence.present ? view.occurrence : 0;
    evidence.specificity = query_seed_specificity(evidence.occurrence);
    evidence.original_index = original_index;
    return evidence;
}

static bool better_query_seed_occurrence(
    const QuerySeedOccurrenceEvidence& a,
    const QuerySeedOccurrenceEvidence& b
) {
    if (a.found != b.found) return a.found && !b.found;
    if (a.occurrence != b.occurrence) return a.occurrence < b.occurrence;
    if (a.seed.read_pos != b.seed.read_pos) return a.seed.read_pos < b.seed.read_pos;
    return a.original_index < b.original_index;
}

static void restore_query_seed_read_order(std::vector<QuerySeed>& seeds) {
    std::sort(seeds.begin(), seeds.end(),
        [](const QuerySeed& a, const QuerySeed& b) {
            if (a.read_pos != b.read_pos) return a.read_pos < b.read_pos;
            return a.key < b.key;
        });
}

static std::vector<QuerySeedOccurrenceEvidence> query_seed_occurrence_evidence(
    const std::vector<QuerySeed>& seeds,
    const IndexView& index
) {
    std::vector<QuerySeedOccurrenceEvidence> evidence;
    evidence.reserve(seeds.size());
    for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
        evidence.push_back(query_seed_occurrence(index, seeds[(size_t)i], i));
    }
    return evidence;
}

void select_spread_from_evidence(
    std::vector<QuerySeed>& seeds,
    const std::vector<QuerySeedOccurrenceEvidence>& evidence,
    int read_len,
    int max_seeds
) {
    const int n_tiles = std::max(
        1,
        std::min(max_seeds, static_cast<int>(evidence.size())));
    std::vector<int> best_by_tile(static_cast<size_t>(n_tiles), -1);
    for (int i = 0; i < static_cast<int>(evidence.size()); ++i) {
        int tile = 0;
        if (read_len > 0) {
            tile = static_cast<int>(
                (static_cast<int64_t>(
                     std::max(0, evidence[static_cast<size_t>(i)].seed.read_pos)) *
                 static_cast<int64_t>(n_tiles)) /
                static_cast<int64_t>(read_len));
            tile = std::max(0, std::min(n_tiles - 1, tile));
        }
        const int prev = best_by_tile[static_cast<size_t>(tile)];
        if (prev < 0 ||
            better_query_seed_occurrence(
                evidence[static_cast<size_t>(i)],
                evidence[static_cast<size_t>(prev)])) {
            best_by_tile[static_cast<size_t>(tile)] = i;
        }
    }

    std::vector<char> selected(evidence.size(), 0);
    seeds.clear();
    seeds.reserve(static_cast<size_t>(max_seeds));
    for (int idx : best_by_tile) {
        if (idx < 0 || selected[static_cast<size_t>(idx)]) continue;
        selected[static_cast<size_t>(idx)] = 1;
        seeds.push_back(evidence[static_cast<size_t>(idx)].seed);
        if (static_cast<int>(seeds.size()) >= max_seeds) break;
    }

    if (static_cast<int>(seeds.size()) < max_seeds) {
        std::vector<int> order(evidence.size());
        for (int i = 0; i < static_cast<int>(order.size()); ++i) {
            order[static_cast<size_t>(i)] = i;
        }
        std::stable_sort(
            order.begin(),
            order.end(),
            [&](int a, int b) {
                return better_query_seed_occurrence(
                    evidence[static_cast<size_t>(a)],
                    evidence[static_cast<size_t>(b)]);
            });
        for (int idx : order) {
            if (selected[static_cast<size_t>(idx)]) continue;
            selected[static_cast<size_t>(idx)] = 1;
            seeds.push_back(evidence[static_cast<size_t>(idx)].seed);
            if (static_cast<int>(seeds.size()) >= max_seeds) break;
        }
    }
    restore_query_seed_read_order(seeds);
}

// Gathers evidence from the index, then spread-picks.
static void select_low_occurrence_spread_query_seeds(
    std::vector<QuerySeed>& seeds,
    const IndexView& index,
    int read_len,
    int max_seeds
) {
    const auto evidence = query_seed_occurrence_evidence(seeds, index);
    select_spread_from_evidence(seeds, evidence, read_len, max_seeds);
}

}}  // namespace fa::cpu

namespace fa { namespace cpu { namespace lr {

static QuerySeedOccurrenceEvidence chain_query_seed_occurrence_cached(
    const LongReadSeedContext& ctx,
    const QuerySeed& seed,
    int original_index,
    ChainSeedLookupCache* lookup_cache
) {
    if (!lookup_cache) {
        return query_seed_occurrence(*ctx.index, seed, original_index);
    }
    QuerySeedOccurrenceEvidence evidence;
    evidence.seed = seed;
    const KmerPostingView view = lookup_cache->lookup(*ctx.index, seed.key);
    evidence.present = view.found();
    evidence.retained = view.found() && view.count > 0;
    evidence.found = evidence.retained;
    evidence.occurrence = evidence.present ? view.occurrence : 0;
    evidence.specificity = query_seed_specificity(evidence.occurrence);
    evidence.original_index = original_index;
    return evidence;
}

static std::vector<QuerySeedOccurrenceEvidence> chain_query_seed_occurrence_evidence_cached(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& seeds,
    ChainSeedLookupCache* lookup_cache
) {
    std::vector<QuerySeedOccurrenceEvidence> evidence;
    evidence.reserve(seeds.size());
    for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
        evidence.push_back(
            chain_query_seed_occurrence_cached(
                ctx, seeds[static_cast<size_t>(i)], i, lookup_cache));
    }
    return evidence;
}

static void chain_select_low_occurrence_spread_query_seeds_cached(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len,
    int max_seeds,
    ChainSeedLookupCache* lookup_cache
) {
    const auto evidence = chain_query_seed_occurrence_evidence_cached(
        ctx, seeds, lookup_cache);
    select_spread_from_evidence(seeds, evidence, read_len, max_seeds);
}

// Query-seed limiting when occurrence-aware selection is off (RNA, or a window with no
// bundle): over budget, keep the lowest-occurrence seed per read strip. A no-op when
// uncapped (max_query_seeds <= 0, the RNA default).
static void limit_query_seeds(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len
) {
    const int max_seeds = ctx.max_query_seeds_per_strand;
    if (max_seeds <= 0 || static_cast<int>(seeds.size()) <= max_seeds) return;
    select_low_occurrence_spread_query_seeds(
        seeds, *ctx.index, read_len, max_seeds);
}

void limit_query_seeds_with_lookup_cache(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len,
    ChainSeedLookupCache* lookup_cache
) {
    if (!lookup_cache) {
        limit_query_seeds(ctx, seeds, read_len);
        return;
    }
    const int max_seeds = ctx.max_query_seeds_per_strand;
    if (max_seeds <= 0 || static_cast<int>(seeds.size()) <= max_seeds) return;
    chain_select_low_occurrence_spread_query_seeds_cached(
        ctx, seeds, read_len, max_seeds, lookup_cache);
}

}}}  // namespace fa::cpu::lr

/**** tie_hash ****/

namespace fa {
namespace cpu {
namespace lr {

std::uint32_t tie_name_hash(std::string_view name) noexcept {
  if (name.empty())
    return 0;
  std::uint32_t hash = static_cast<unsigned char>(name[0]);
  for (std::size_t i = 1; i < name.size(); ++i)
    hash = (hash << 5) - hash + static_cast<unsigned char>(name[i]);
  return hash;
}

// khash.h __ac_Wang_hash.
static std::uint32_t tie_wang_hash(std::uint32_t key) noexcept {
  key += ~(key << 15);
  key ^= (key >> 10);
  key += (key << 3);
  key ^= (key >> 6);
  key += ~(key << 11);
  key ^= (key >> 16);
  return key;
}

std::uint32_t tie_read_seed(std::uint32_t name_hash,
                            int read_len) noexcept {
  std::uint32_t hash = name_hash;
  hash ^= tie_wang_hash(static_cast<std::uint32_t>(read_len)) +
          tie_wang_hash(kTieHashSeed);
  return tie_wang_hash(hash);
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** query_filter ****/

namespace fa { namespace cpu { namespace lr {

const std::vector<QuerySeed>* drop_query_repeats(
    const std::vector<QuerySeed>& seeds, int cap, float frac,
    std::vector<std::uint64_t>& keys, std::vector<QuerySeed>& kept,
    std::size_t& dropped) {
  dropped = 0;
  const std::size_t n = seeds.size();
  if (cap <= 0 || frac <= 0.0f || n <= static_cast<std::size_t>(cap))
    return &seeds;
  keys.clear();
  keys.reserve(n);
  for (const QuerySeed& seed : seeds)
    keys.push_back(seed.key);
  std::sort(keys.begin(), keys.end());
  // minimap2's float test, cnt > n * q_occ_frac.
  const float frac_floor = static_cast<float>(n) * frac;
  // The repeated keys, compacted in order to the front of `keys`.
  std::size_t repeated = 0;
  for (std::size_t begin = 0, i = 1; i <= n; ++i) {
    if (i < n && keys[i] == keys[begin])
      continue;
    const std::size_t count = i - begin;
    if (count > static_cast<std::size_t>(cap) &&
        static_cast<float>(count) > frac_floor)
      keys[repeated++] = keys[begin];
    begin = i;
  }
  if (repeated == 0)
    return &seeds;
  keys.resize(repeated);
  kept.clear();
  kept.reserve(n);
  for (const QuerySeed& seed : seeds)
    if (!std::binary_search(keys.begin(), keys.end(), seed.key))
      kept.push_back(seed);
  dropped = n - kept.size();
  return &kept;
}

}}}  // namespace fa::cpu::lr

namespace fa { namespace cpu { namespace lr {

void seed_select(const std::vector<SeedSelectItem>& seeds, int query_length,
                 std::uint32_t cap, int dist,
                 std::vector<std::uint32_t>& kept) {
  if (dist <= 0)
    return;
  std::vector<std::size_t> run;
  const std::size_t n = seeds.size();
  // [st, i) is the streak of seeds over the cap that seed i (or the query
  // end) closes.
  std::size_t st = 0;
  for (std::size_t i = 0; i <= n; ++i) {
    if (i < n && seeds[i].occurrence > cap)
      continue;
    if (i > st) {
      const int ps = st == 0 ? 0 : seeds[st - 1].end;
      const int pe = i == n ? query_length : seeds[i].end;
      const int m = std::min(
          kSeedSelectMaxKept,
          static_cast<int>(static_cast<double>(pe - ps) / dist + .499));
      if (m > 0) {
        run.clear();
        for (std::size_t j = st; j < i; ++j)
          run.push_back(j);
        const std::size_t take =
            std::min(run.size(), static_cast<std::size_t>(m));
        // minimap2's heap keeps the earlier of two equal occurrences.
        std::partial_sort(run.begin(),
                          run.begin() + static_cast<std::ptrdiff_t>(take),
                          run.end(), [&](std::size_t a, std::size_t b) {
                            return seeds[a].occurrence != seeds[b].occurrence
                                       ? seeds[a].occurrence <
                                             seeds[b].occurrence
                                       : a < b;
                          });
        for (std::size_t j = 0; j < take; ++j)
          if (seeds[run[j]].occurrence <= kSeedSelectMaxOcc)
            kept.push_back(seeds[run[j]].id);
      }
    }
    st = i + 1;
  }
}

std::size_t mark_readmitted_seeds(const std::vector<QuerySeed>& seeds,
                                  const std::vector<KmerPostingView>& views,
                                  int read_len, int k, int cap, int dist,
                                  std::vector<std::uint8_t>& readmitted) {
  readmitted.assign(static_cast<std::size_t>(std::max(read_len, 0)), 0);
  if (cap <= 0 || dist <= 0)
    return 0;
  // minimap2's seed list: the seeds the index holds, in read order.
  std::vector<SeedSelectItem> items;
  std::vector<std::uint32_t> kept;
  const std::size_t n = std::min(seeds.size(), views.size());
  for (std::size_t i = 0; i < n; ++i)
    if (views[i].found() && views[i].count > 0)
      items.push_back(SeedSelectItem{seeds[i].read_pos + k,
                                     views[i].occurrence,
                                     static_cast<std::uint32_t>(i)});
  seed_select(items, read_len, static_cast<std::uint32_t>(cap), dist, kept);
  for (const std::uint32_t i : kept)
    readmitted[static_cast<std::size_t>(seeds[i].read_pos)] = 1;
  return kept.size();
}

}}}  // namespace fa::cpu::lr
