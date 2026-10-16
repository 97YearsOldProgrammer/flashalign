// The DNA conductor, DnaBackend::map_read: seeds, vote, catalogue, placement, ownership, projection or realization, MAPQ, one read at a time.
#include "dna.h"

#include "../core/core.h"
#include "../index/index.h"
#include "../index/seed.h"
#include "../seed/scratch.h"
#include "../seed/seed.h" // tie_read_seed (the vote's tie-break)
#include "../vote/vote.h" // vote_slope_widen_winner (the winner's window)

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace dna {

// Stands in for the reverse-complement query when no stage reads its bases
// (see reverse_query_bases_needed); consumers take the length from the read.
static const std::vector<uint8_t> kNoQueryBases;

} // namespace dna
} // namespace fa::cpu::lr

/**** unmapped ****/

namespace fa::cpu::lr {
namespace dna {

static void demote_unmapped(Result& result) {
  result.pos = -1;
  result.score = 0;
}

} // namespace dna
} // namespace fa::cpu::lr

/**** vote ****/

namespace fa::cpu::lr {

// Length of the whole-read vote window.
static int effective_chain_segment_len(const DnaContext& context,
                                       int read_length) {
  return std::max(
      std::max(context.opts.k, 32), std::max(1, read_length));
}
} // namespace fa::cpu::lr

/**** catalogue ****/

namespace fa::cpu::lr {
namespace dna {

// The all-chains lane (options/options.h all_chains). Every chain of
// `chains`, parallel to family.candidates, that reaches min_chain_score is
// projected whole as one record, best score first. A record on the contig and
// strand of a kept one, overlapping its reference span, is the same placement
// and is dropped. The best record is the read's head and the others its
// supplementary records, all MAPQ 0 and tp:A:S, as minimap2's -P prints them.
// A read left with none is unmapped.
static void emit_all_chains(
    const DnaContext& dctx, const DnaPlacementFamily& family,
    const std::vector<DnaPlacementCandidateChain>& chains, dna::Result& out) {
  std::vector<std::size_t> order;
  for (std::size_t index = 0; index < chains.size(); ++index)
    if (chains[index].status == DnaPlacementChainStatus::Accepted &&
        chains[index].chain_score >= dctx.opts.min_chain_score)
      order.push_back(index);
  // Catalogue order breaks a tie on score.
  std::stable_sort(order.begin(), order.end(),
                   [&chains](std::size_t left, std::size_t right) {
                     return chains[left].chain_score >
                            chains[right].chain_score;
                   });
  // The whole read as one block, so the projection trims nothing.
  ::fa::cpu::voting::QueryBlock whole;
  whole.query_tile_begin = 0;
  whole.query_tile_end = family.tile_count;
  whole.supporting_tiles = family.tile_count;
  std::vector<AlignResult> records;
  std::vector<std::size_t> kept;
  for (const std::size_t index : order) {
    const DnaPlacementCandidate& candidate = family.candidates[index];
    whole.candidate = candidate.id;
    AlignResult record;
    if (!dna_project_chained_block(dctx, family, chains[index], whole, record))
      continue;
    bool duplicate = false;
    for (std::size_t slot = 0; slot < kept.size() && !duplicate; ++slot) {
      const VotePeak& other = family.candidates[kept[slot]].peak;
      duplicate = other.chr == candidate.peak.chr &&
                  other.is_rc == candidate.peak.is_rc &&
                  records[slot].pos < record.target_end &&
                  record.pos < records[slot].target_end;
    }
    if (duplicate)
      continue;
    record.origin = AlignmentOrigin::DnaAllChains;
    record.chain_anchors = chains[index].chain_anchors;
    record.chain_score = chains[index].chain_score;
    records.push_back(std::move(record));
    kept.push_back(index);
  }
  if (records.empty())
    return;
  static_cast<AlignResult&>(out) = std::move(records.front());
  for (std::size_t slot = 1; slot < records.size(); ++slot) {
    out.supplementary.push_back(std::move(records[slot]));
    out.supplementary_candidates.push_back(family.candidates[kept[slot]].id);
  }
}

} // namespace dna
} // namespace fa::cpu::lr

/**** emission ****/

namespace fa::cpu::lr {
namespace dna {

// Map-only's floor. minimap2 keeps no chain under
// min_chain_score (mm_chain_dp) and does not test again what it splits off
// one (mm_split_reg), so each record is judged by the score of the block part
// it prints from. A failing head gives its place to the next record. `owners`
// and `parts` are parallel to `rest`, `head_part` is the head's part. Returns
// false when no record is left.
static bool keep_dna_min_chain_score(
    AlignResult& head, int head_chain_score, std::vector<AlignResult>& rest,
    std::vector<::fa::cpu::voting::CandidateId>& owners,
    const DnaPlacementFamily& family, int min_chain_score,
    std::vector<int>& parts, int& head_part) {
  std::size_t kept = 0;
  for (std::size_t index = 0; index < rest.size(); ++index) {
    if (family.block_parts[static_cast<std::size_t>(parts[index])].score <
        min_chain_score)
      continue;
    if (kept != index) {
      rest[kept] = std::move(rest[index]);
      owners[kept] = owners[index];
      parts[kept] = parts[index];
    }
    ++kept;
  }
  rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(kept), rest.end());
  owners.resize(kept);
  parts.resize(kept);
  if (head_chain_score >= min_chain_score)
    return true;
  if (rest.empty())
    return false;
  // `rest` aliases head.supplementary: move everything out first.
  std::vector<AlignResult> survivors = std::move(rest);
  std::vector<AlignResult> secondaries = std::move(head.secondary);
  AlignResult next = std::move(survivors.front());
  survivors.erase(survivors.begin());
  owners.erase(owners.begin());
  head_part = parts.front();
  parts.erase(parts.begin());
  // Read-level fields carry over; the record keeps its own MAPQ.
  next.read_len = head.read_len;
  next.median_occurrence = head.median_occurrence;
  head = std::move(next);
  head.supplementary = std::move(survivors);
  head.secondary = std::move(secondaries);
  return true;
}

// Applies it to a map-only read. Only a read left with no record is unmapped.
static void apply_dna_min_chain_score(
    const DnaPlacementFamily& family, int min_chain_score,
    dna::Result& realized) {
  if (keep_dna_min_chain_score(
          static_cast<AlignResult&>(realized),
          family.block_parts[static_cast<std::size_t>(realized.primary_part)]
              .score,
          realized.supplementary, realized.supplementary_candidates, family,
          min_chain_score, realized.supplementary_parts,
          realized.primary_part))
    return;
  demote_unmapped(realized);
  realized.supplementary.clear();
  realized.supplementary_candidates.clear();
  realized.supplementary_parts.clear();
  realized.secondary.clear();
}

} // namespace dna
} // namespace fa::cpu::lr

/**** map_read ****/

namespace fa::cpu::lr {
namespace dna {

AlignResult map_read(const DnaContext& base_dctx,
                     const LongReadSeedContext& seed_ctx,
                     DnaWorkerScratch& worker_scratch, const std::string& read,
                     const std::vector<uint8_t>& fwd_enc,
                     std::vector<uint8_t>& rc_enc,
                     const std::vector<QuerySeed>* shared_fwd_syncmer_seeds) {
  DnaContext dctx = base_dctx;
  auto effective_chain_segment_len = [&](int rl) {
    return ::fa::cpu::lr::effective_chain_segment_len(dctx, rl);
  };
  auto seed_window_estimate = [&]() {
    const int eff_s = effective_closed_syncmer_s(std::max(1, dctx.opts.k),
                                                 dctx.opts.chain_syncmer_s);
    return std::max(1, dctx.opts.k - eff_s + 1);
  };
  const int read_len = static_cast<int>(read.size());
  const int segment_len = effective_chain_segment_len(read_len);
  // The read's tie seed, from the name hash and the read length as minimap2;
  // every tie-break of the read uses it.
  dctx.vote_tie_seed = tie_read_seed(dctx.read_name_hash, read_len);
  if (dctx.opts.all_chains)
    dctx.self_contig = seed_ctx.self_contig;
  dna::Result out{};
  auto finish = [](dna::Result result) {
    return static_cast<AlignResult&&>(std::move(result));
  };
  out.read_len = read_len;
  if (read_len < std::max(2 * dctx.opts.k, segment_len)) {
    return finish(std::move(out));
  }

  ChainAnchorScratch& scratch = worker_scratch.chain;
  worker_scratch.clear_for_read(static_cast<size_t>(read_len));
  auto ensure_rc_enc = [&]() -> const std::vector<uint8_t>& {
    if (rc_enc.size() != fwd_enc.size())
      rc_enc = reverse_complement_encoded_u8(fwd_enc);
    return rc_enc;
  };
  // Only CIGAR realization reads reverse-strand query bases; map-only never
  // builds the buffer and passes kNoQueryBases instead.
  const bool reverse_query_bases_needed = dctx.opts.enable_full_read_cigar;
  auto reverse_query_stream = [&]() -> const std::vector<uint8_t>& {
    return reverse_query_bases_needed ? ensure_rc_enc() : kNoQueryBases;
  };
  auto& raw = scratch.raw;
  raw.reserve(2);
  ChainSeedLookupCache* lookup_cache = &scratch.lookup_cache;
  const size_t expected_keys =
      static_cast<size_t>(
          std::max(1, read_len / std::max(1, seed_window_estimate()))) *
          2u +
      8u + 64u;
  lookup_cache->reset(expected_keys);
  // Every consumer below, the vote, the tile rescue and the placement
  // harvest, reads the stream without the read's repeated keys, capped at the
  // pool gate's N. The all-chains lane keeps every seed.
  if (!dctx.opts.all_chains && shared_fwd_syncmer_seeds != nullptr) {
    std::size_t dropped = 0;
    shared_fwd_syncmer_seeds = drop_query_repeats(
        *shared_fwd_syncmer_seeds, dctx.opts.dna_pool_gate_occ,
        dctx.opts.dna_query_occ_frac, scratch.query_occ_keys,
        scratch.query_occ_seeds, dropped);
  }
  if (shared_fwd_syncmer_seeds == nullptr ||
      shared_fwd_syncmer_seeds->empty()) {
    return finish(std::move(out));
  }
  DnaQuerySeedPool query_seed_pool(seed_ctx, read_len, *lookup_cache);
  // Per-seed captured views and cache slots, read by position later instead
  // of probing the cache again.
  query_seed_pool.bind_capture_stores(
      &scratch.captured_fwd_views, &scratch.captured_fwd_slots,
      &scratch.captured_rc_views, &scratch.captured_rc_slots);
  query_seed_pool.bind(DnaQuerySeedStrand::Forward, *shared_fwd_syncmer_seeds);
  if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Forward))
    return finish(std::move(out));
  const std::vector<QuerySeed>* shared_rc_syncmer_seeds = nullptr;
  auto ensure_shared_rc_stream = [&]() -> const std::vector<QuerySeed>* {
    if (!scratch.full_syncmer_rc_seeds.empty()) {
      query_seed_pool.bind(DnaQuerySeedStrand::Reverse,
                           scratch.full_syncmer_rc_seeds);
      if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Reverse))
        return nullptr;
      return &scratch.full_syncmer_rc_seeds;
    }
    project_query_seeds_to_rc_coordinates(*shared_fwd_syncmer_seeds, read_len,
                                          dctx.opts.k,
                                          scratch.full_syncmer_rc_seeds);
    query_seed_pool.bind(DnaQuerySeedStrand::Reverse,
                         scratch.full_syncmer_rc_seeds);
    if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Reverse))
      return nullptr;
    return &scratch.full_syncmer_rc_seeds;
  };
  auto prepare_shared_window_views =
      [&](DnaQuerySeedStrand strand, int seed_start, int seed_end,
          int downsample, std::vector<QuerySeed>& window_seeds,
          const std::vector<DnaLongSeedView>*& exact_refine_views)
      -> const std::vector<DnaLongSeedView>* {
    scratch.shadow_rc.clear_for_window();
    DnaLongSeedBundle& bundle = scratch.shadow_rc.seed_bundle;
    const bool selected = query_seed_pool.prepare_representatives(
        strand, seed_start, seed_end, downsample, window_seeds, bundle,
        scratch.shadow_rc);
    if (!selected)
      return nullptr;
    exact_refine_views = bundle.exact_refine_views.empty()
                             ? nullptr
                             : &bundle.exact_refine_views;
    return bundle.views_ready ? &bundle.selected_views : nullptr;
  };
  // The vote runs over one whole-read window [0, read_len).
  if (read_len >= dctx.opts.k) {
    const int start = 0;
    const int end = read_len;
    const int span = read_len;
    const int actual_downsample =
        std::max(1, dctx.opts.chain_syncmer_downsample);
    const std::vector<DnaLongSeedView>* fwd_seed_views = nullptr;
    const std::vector<DnaLongSeedView>* fwd_exact_refine_views = nullptr;
    fwd_seed_views = prepare_shared_window_views(
        DnaQuerySeedStrand::Forward, start, end, actual_downsample,
        scratch.shared_window_fwd_seeds, fwd_exact_refine_views);
    if (!fwd_seed_views)
      return finish(std::move(out));
    scratch.fwd.pending_exact_refine_views =
        fwd_exact_refine_views ? *fwd_exact_refine_views
                               : std::vector<DnaLongSeedView>{};
    // Kept for vote_slope_widen_winner: the reverse prepare below clears the
    // bundle these views live in.
    scratch.slope_refine_views[0] = scratch.fwd.pending_exact_refine_views;
    WindowAnchorPeakParams fwd_peak_params;
    fwd_peak_params.lookup_cache = lookup_cache;
    fwd_peak_params.seed_view_override = fwd_seed_views;
    fwd_peak_params.tie_seed = dctx.vote_tie_seed;
    auto fwd =
        window_anchor_peaks(seed_ctx, fwd_enc.data() + start, span,
                            /*is_rc=*/false, scratch.fwd, fwd_peak_params);

    const int rc_start = read_len - end;
    const std::vector<DnaLongSeedView>* rc_seed_views = nullptr;
    const std::vector<DnaLongSeedView>* rc_exact_refine_views = nullptr;
    shared_rc_syncmer_seeds = ensure_shared_rc_stream();
    if (!shared_rc_syncmer_seeds)
      return finish(std::move(out));
    rc_seed_views = prepare_shared_window_views(
        DnaQuerySeedStrand::Reverse, rc_start, rc_start + span,
        actual_downsample, scratch.shared_window_rc_seeds,
        rc_exact_refine_views);
    if (!rc_seed_views)
      return finish(std::move(out));
    scratch.rc.pending_exact_refine_views =
        rc_exact_refine_views ? *rc_exact_refine_views
                              : std::vector<DnaLongSeedView>{};
    scratch.slope_refine_views[1] = scratch.rc.pending_exact_refine_views;
    WindowAnchorPeakParams rc_peak_params;
    rc_peak_params.lookup_cache = lookup_cache;
    rc_peak_params.seed_view_override = rc_seed_views;
    rc_peak_params.tie_seed = dctx.vote_tie_seed;
    auto rc = window_anchor_peaks(seed_ctx, nullptr, span, /*is_rc=*/true,
                                  scratch.rc, rc_peak_params);
    // The winner's harvest window and band, on its line.
    vote_slope_widen_winner(seed_ctx, dctx.vote_tie_seed, span, scratch, fwd,
                            rc);
    raw.insert(raw.end(), fwd.begin(), fwd.end());
    raw.insert(raw.end(), rc.begin(), rc.end());
  }
  if (raw.empty())
    return finish(std::move(out));
  // The placement family comes from the whole-read vote alone.
  std::vector<VotePeak> selected_chain;
  // Only whether a best peak exists matters here; the family's own order
  // decides catalogue candidate 0.
  const VotePeak* best = best_peak(raw);
  if (best)
    selected_chain.push_back(*best);
  DnaPlacementFamily placement_family;
  DnaPlacementChainingResult placement_chaining;
  bool placement_chaining_ran = false;
  if (best) {
    placement_family =
        build_dna_placement_family(dctx, fwd_enc, raw, scratch.fwd, scratch.rc);
    // The all-chains lane: the catalogue's chains are the records, and nothing
    // below runs.
    if (dctx.opts.all_chains) {
      if (placement_family.valid)
        emit_all_chains(
            dctx, placement_family,
            build_dna_target_chains(
                dctx, placement_family, fwd_enc, reverse_query_stream(),
                &scratch.fwd.retained_seeds, &scratch.rc.retained_seeds,
                shared_fwd_syncmer_seeds, shared_rc_syncmer_seeds,
                lookup_cache,
                query_seed_pool.captured_slots(DnaQuerySeedStrand::Forward),
                query_seed_pool.captured_slots(DnaQuerySeedStrand::Reverse)),
            out);
      return finish(std::move(out));
    }
    if (placement_family.valid) {
      const std::vector<uint8_t>& reverse_query = reverse_query_stream();
      placement_chaining = build_dna_placement_chains(
          dctx, placement_family, fwd_enc, reverse_query,
          &scratch.fwd.retained_seeds, &scratch.rc.retained_seeds,
          shared_fwd_syncmer_seeds, shared_rc_syncmer_seeds, lookup_cache,
          query_seed_pool.captured_slots(DnaQuerySeedStrand::Forward),
          query_seed_pool.captured_slots(DnaQuerySeedStrand::Reverse),
          worker_scratch.seed_density);
      placement_chaining_ran = true;
      dctx.inversion_gate_seeds = placement_chaining.inversion_gate_seeds;
      placement_family = placement_chaining.family;
    }
  }
  const auto& chain = selected_chain;
  if (chain.empty()) {
    return finish(std::move(out));
  }
  dna::Result realized;
  // Catalogue id 0 is the top-ranked vote candidate.
  ::fa::cpu::voting::CandidateId primary_candidate = 0;
  // The -c lane: the controller's records (realization.cpp) and what
  // the MAPQ reads for each primary.
  std::vector<DnaRegionPrice> region_prices;
  std::vector<DnaRealizedRival> region_rivals;
  if (dctx.opts.enable_full_read_cigar) {
    const std::vector<uint8_t>& reverse_query = ensure_rc_enc();
    DnaFamilyRealizationRequest request;
    request.family = &placement_family;
    request.placement = placement_chaining_ran ? &placement_chaining : nullptr;
    request.forward_query = &fwd_enc;
    request.reverse_query = &reverse_query;
    DnaRegionOutcome regions = realize_dna_regions(dctx, request);
    DnaFamilyRealizationOutcome& outcome = regions.family;
    primary_candidate = outcome.primary_candidate;
    realized = std::move(outcome.output);
    realized.read_len = read_len;
    region_prices = std::move(regions.prices);
    region_rivals = std::move(regions.rivals);
    if (!realized.mapped()) {
      if (regions.floored)
        demote_unmapped(realized);
    }
  } else {
    // Map-only: the placement family's projection replaces this carrier
    // wholesale, so no coarse vote geometry survives. A best peak without
    // vote support leaves the read unmapped without consulting the family.
    realized.read_len = read_len;
    realized.pos = 0;
    realized.score = std::max(0, chain.front().support);
  }
  std::vector<DnaJoinPiece> join_pieces;
  if (realized.mapped()) {
    if (!dctx.opts.enable_full_read_cigar) {
      // Placement chaining always ran here: a best peak exists past the
      // empty-chain return, so the catalogue is non-empty and the family is
      // valid by construction (catalogue.cpp).
      DnaFamilyProjectionResult projection =
          project_map_only_placement_family(dctx, placement_chaining, realized);
      primary_candidate = projection.primary_candidate;
      join_pieces = std::move(projection.pieces);
      realized = settle_map_only_projection(std::move(projection));
    }
    if (dctx.opts.enable_full_read_cigar) {
      price_dna_regions(dctx, placement_family, placement_chaining,
                        region_prices, region_rivals, read_len, realized);
    } else if (realized.mapped()) {
      // Map-only joins the block records that pass the bridge test, so the
      // MAPQ below scores a joined record on its owner's chain.
      dna_join_map_only_family(dctx, placement_chaining,
                               std::move(join_pieces), realized,
                               primary_candidate);
      price_dna_blocks(dctx, placement_family, placement_chaining,
                       placement_chaining_ran, primary_candidate, read_len,
                       realized);
      // The secondaries, after the MAPQ, which they do not enter.
      dna_project_map_only_secondaries(dctx, placement_chaining, realized);
    }
  }
  // The emitted family is final here, MAPQ included; the -c lane's emission
  // floor ran in its controller. One demotion follows on CIGAR output, and
  // one on map-only.
  //
  // (1) A mapped primary with an empty CIGAR is a realization refusal. Emit
  // it unmapped, as minimap2 does for a read it cannot align, rather than as
  // a whole-read match or all-soft-clip. Such a primary has no
  // supplementaries.
  if (dctx.opts.enable_full_read_cigar && realized.mapped() &&
      realized.cigar.empty()) {
    demote_unmapped(realized);
  }
  // (2) Map-only: no record under min_chain_score, judged on the block parts
  // the MAPQ reads (apply_dna_min_chain_score).
  if (!dctx.opts.enable_full_read_cigar && realized.mapped())
    apply_dna_min_chain_score(placement_chaining.family,
                              dctx.opts.min_chain_score, realized);
  // ms:i: the max-scoring segment of each emitted record's CIGAR, after the
  // demotion above (hence the second mapped() test). Map-only records have
  // no CIGAR and keep -1.
  if (dctx.opts.enable_full_read_cigar && realized.mapped()) {
    const std::vector<uint8_t>& reverse_query = ensure_rc_enc();
    auto stamp_dp_max_segment = [&](AlignResult& record) {
      record.dp_max_segment =
          dna_record_dp_max_segment(dctx, record, fwd_enc, reverse_query);
    };
    stamp_dp_max_segment(static_cast<AlignResult&>(realized));
    for (AlignResult& part : realized.supplementary)
      stamp_dp_max_segment(part);
    for (AlignResult& secondary : realized.secondary) {
      stamp_dp_max_segment(secondary);
      for (AlignResult& part : secondary.supplementary)
        stamp_dp_max_segment(part);
    }
  }
  return finish(std::move(realized));
}

} // namespace dna
} // namespace fa::cpu::lr
