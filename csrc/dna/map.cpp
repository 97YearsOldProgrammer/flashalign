// The DNA conductor, DnaBackend::map_read: seeds, vote, catalogue, placement, ownership, projection or realization, MAPQ, one read at a time.
#include "dna.h"

#include "../chain/chain.h" // dense_chain_scratch_bytes (kept_scratch)
#include "../core/core.h"
#include "../dp/dp.h"       // thread_arena_bytes (kept_scratch)
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

// The forward [begin, end) of each stretch of a query of `length`, `stretch`
// bp each; the last takes the remainder, and a remainder under stretch / 2
// joins the previous one.
static std::vector<std::pair<int, int>> query_stretches(int length,
                                                        int stretch) {
  int count = std::max(1, length / stretch);
  if (length > count * stretch &&
      2 * (length - count * stretch) >= stretch)
    ++count;
  std::vector<std::pair<int, int>> stretches;
  stretches.reserve(static_cast<std::size_t>(count));
  for (int at = 0; at < count; ++at)
    stretches.emplace_back(at * stretch,
                           at + 1 == count ? length : (at + 1) * stretch);
  return stretches;
}

// A stretch's vote ran on its window's rebased stream, `offset` the window's
// start on its lane's stream and `span` its length: its peaks and the seeds
// it retained move to the lane's oriented read coordinates. A diagonal
// d = r - q there is d + offset in the window's, so raw_ref_start loses the
// offset and the bin is taken again from it at the vote's width; read_lo /
// read_hi become the window's bounds, and ref_pos is clamped as the
// whole-read vote clamps it. The bin is floor(raw_ref_start / W) at offset 0
// too, not the whole-read vote's bucket bin. read_lo is the stretch's begin
// forward and L minus the window's end reverse, so within a stretch a forward
// and a reverse peak tied on vote_score, support and center_support go by
// read_lo in chain_peak_better_seeded before the locus hash.
static void stretch_to_contig(const DnaContext& dctx, int read_len,
                              int offset, int span,
                              std::vector<VotePeak>& peaks,
                              std::vector<ChainWindowRetainedSeed>& retained) {
  for (VotePeak& peak : peaks) {
    peak.raw_ref_start -= offset;
    const int width = std::max(1, peak.anchor.ref_start_bin_width);
    peak.anchor.ref_start_bin =
        vote_floor_div(static_cast<int>(peak.raw_ref_start), width);
    peak.read_lo = offset;
    peak.read_hi = offset + span;
    if (peak.evidence_read_lo <= peak.evidence_read_hi) {
      peak.evidence_read_lo += offset;
      peak.evidence_read_hi += offset;
    }
    const std::int64_t chromosome_length = dctx.ref.contig_length(peak.chr);
    peak.ref_pos = static_cast<int>(std::clamp<std::int64_t>(
        peak.raw_ref_start, 0,
        std::max<std::int64_t>(0, chromosome_length - read_len)));
  }
  for (ChainWindowRetainedSeed& seed : retained) seed.seed.read_pos += offset;
}

} // namespace dna
} // namespace fa::cpu::lr

/**** stage_release ****/
// A read's dead buffers, freed at their last reader once their group holds
// more than kKeptScratchReleaseBytes (dna.h kept_scratch).

namespace fa::cpu::lr {
namespace dna {
namespace {

template <class T>
void free_vector(std::vector<T>& buffer) {
  std::vector<T>().swap(buffer);
}

// Frees a read's buffers at their last reader once their group holds more
// than kKeptScratchReleaseBytes: the seed side, and the placement side (the
// chain scratch and the seed density). A group once due stays due for the
// read, and each of its points also frees what its earlier points free.
class StageRelease final : public DnaStageListener {
 public:
  StageRelease(const DnaContext& context, DnaWorkerScratch& worker,
               const std::vector<QuerySeed>* forward_stream)
      : worker_(worker),
        forward_stream_(context.stage_free_fwd_stream == forward_stream
                            ? context.stage_free_fwd_stream
                            : nullptr),
        // The inversion probe's gate reads the density during realization.
        keep_density_(context.opts.inversion_probe_local_gate &&
                      context.opts.enable_full_read_cigar) {}

  // After drop_query_repeats, which left `in_use` for the stages below.
  void after_query_repeats(const std::vector<QuerySeed>* in_use) {
    in_use_ = in_use;
    if (seed_due())
      free_query_repeats();
  }

  // After the read's last prepare_representatives, the captured views' last
  // reader: the whole-read vote's reverse prepare, or under the vote per
  // stretch the last stretch's.
  void after_last_prepare() {
    if (!seed_due())
      return;
    free_query_repeats();
    free_captured_views();
  }

  void at(DnaStagePoint point) override {
    switch (point) {
      case DnaStagePoint::Density:
        if (!seed_due())
          return;
        free_query_repeats();
        free_captured_views();
        free_density_inputs();
        return;
      case DnaStagePoint::Rivals:
        if (placement_due())
          free_chain_scratch();
        return;
      case DnaStagePoint::Retained:
        if (!placement_due())
          return;
        free_chain_scratch();
        if (!keep_density_)
          worker_.seed_density = RetainedSeedDensity{};
        return;
    }
  }

 private:
  std::size_t seed_bytes() const {
    const ChainAnchorScratch& chain = worker_.chain;
    std::size_t bytes =
        capacity_bytes(chain.query_occ_keys, chain.query_occ_seeds,
                       chain.captured_fwd_views, chain.captured_rc_views,
                       chain.captured_fwd_slots, chain.captured_rc_slots,
                       chain.full_syncmer_rc_seeds,
                       chain.shared_window_fwd_seeds,
                       chain.shared_window_rc_seeds) +
        chain.lookup_cache.heap_bytes() +
        worker_.seed_density.build_scratch_bytes();
    if (forward_stream_ != nullptr)
      bytes += capacity_bytes(*forward_stream_);
    return bytes;
  }

  std::size_t placement_bytes() const {
    return ::fa::cpu::chaining::dense_chain_scratch_bytes() +
           pass_slices_bytes() + pool_sort_scratch_bytes() +
           worker_.seed_density.heap_bytes();
  }

  bool seed_due() {
    if (!seed_fired_ && seed_bytes() > kKeptScratchReleaseBytes)
      seed_fired_ = true;
    return seed_fired_;
  }

  bool placement_due() {
    if (!placement_fired_ && placement_bytes() > kKeptScratchReleaseBytes)
      placement_fired_ = true;
    return placement_fired_;
  }

  // The keys are dead once the repeats are dropped; the kept seeds unless
  // they are the stream in use.
  void free_query_repeats() {
    ChainAnchorScratch& chain = worker_.chain;
    free_vector(chain.query_occ_keys);
    if (in_use_ != &chain.query_occ_seeds)
      free_vector(chain.query_occ_seeds);
  }

  void free_captured_views() {
    free_vector(worker_.chain.captured_fwd_views);
    free_vector(worker_.chain.captured_rc_views);
  }

  // What the density's build read last: the cache, the captured slots, the
  // seed streams and its own build scratch.
  void free_density_inputs() {
    ChainAnchorScratch& chain = worker_.chain;
    free_vector(chain.query_occ_seeds);
    free_vector(chain.captured_fwd_slots);
    free_vector(chain.captured_rc_slots);
    free_vector(chain.full_syncmer_rc_seeds);
    free_vector(chain.shared_window_fwd_seeds);
    free_vector(chain.shared_window_rc_seeds);
    worker_.seed_density.release_build_scratch();
    chain.lookup_cache = ChainSeedLookupCache{};
    if (forward_stream_ != nullptr)
      free_vector(*forward_stream_);
  }

  // Nothing chains after placement's groups (placement.cpp
  // chain_and_own_every_candidate).
  void free_chain_scratch() {
    ::fa::cpu::chaining::release_dense_chain_scratch();
    release_placement_scratch();
  }

  DnaWorkerScratch& worker_;
  std::vector<QuerySeed>* const forward_stream_;
  const std::vector<QuerySeed>* in_use_ = nullptr;
  const bool keep_density_;
  bool seed_fired_ = false;
  bool placement_fired_ = false;
};

} // namespace
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
  StageRelease stage_release(dctx, worker_scratch, shared_fwd_syncmer_seeds);
  dctx.stage_listener = &stage_release;
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
  stage_release.after_query_repeats(shared_fwd_syncmer_seeds);
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
  // minimap2's seed selection under -e, at the pool gate's N, over the stream
  // the filter left: the vote and every pool gate admit the re-admitted seeds
  // as rescued ones.
  if (dctx.opts.occ_dist > 0) {
    mark_readmitted_seeds(*shared_fwd_syncmer_seeds, scratch.captured_fwd_views,
                          read_len, dctx.opts.k, dctx.opts.dna_pool_gate_occ,
                          dctx.opts.occ_dist, scratch.readmitted);
    dctx.readmitted = &scratch.readmitted;
  }
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
  // `seed_start` is the window's start on the strand's stream, where its
  // rebased read positions begin.
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
    // A re-admitted representative is carried as a rescued one, which the
    // pool gates admit. The vote drew it from the seeds its own cap admits,
    // so the flag moves no vote.
    if (dctx.readmitted != nullptr)
      for (DnaLongSeedView& view : bundle.selected_views) {
        const int stream_position = seed_start + view.seed.read_pos;
        const int forward_position =
            strand == DnaQuerySeedStrand::Forward
                ? stream_position
                : read_len - dctx.opts.k - stream_position;
        if ((*dctx.readmitted)[static_cast<std::size_t>(forward_position)])
          view.rescued = true;
      }
    exact_refine_views = bundle.exact_refine_views.empty()
                             ? nullptr
                             : &bundle.exact_refine_views;
    return bundle.views_ready ? &bundle.selected_views : nullptr;
  };
  // The asm rows (options/options.h vote_stretch): the vote below, per
  // stretch and strand. A stretch's seeds are those starting in
  // its [begin, end), so its window reaches k - 1 bases further; the reverse
  // window mirrors the forward one. Each vote's peaks and retained seeds move
  // to contig coordinates (stretch_to_contig), each stretch's catalogue is
  // built from its own peaks and seeds (append_stretch_family), and each
  // lane's retained seeds become their union, in read order. No peak gets the
  // winner's line.
  const int stretch = dctx.opts.vote_stretch;
  DnaPlacementFamily stretch_family;
  if (stretch > 0 && read_len >= dctx.opts.k) {
    const int actual_downsample =
        std::max(1, dctx.opts.chain_syncmer_downsample);
    shared_rc_syncmer_seeds = ensure_shared_rc_stream();
    if (!shared_rc_syncmer_seeds)
      return finish(std::move(out));
    query_seed_pool.allow_windows();
    std::vector<ChainWindowRetainedSeed> retained[2];
    std::vector<VotePeak> stretch_peaks;
    const std::vector<std::pair<int, int>> stretches =
        query_stretches(read_len, stretch);
    for (const auto& [begin, end] : stretches) {
      const int window_end = std::min(read_len, end + dctx.opts.k - 1);
      const int span = window_end - begin;
      stretch_peaks.clear();
      for (int lane = 0; lane < 2; ++lane) {
        const bool reverse = lane == 1;
        const int offset = reverse ? read_len - window_end : begin;
        ChainWindowPeakScratch& lane_scratch =
            reverse ? scratch.rc : scratch.fwd;
        const std::vector<DnaLongSeedView>* refine_views = nullptr;
        const std::vector<DnaLongSeedView>* views =
            prepare_shared_window_views(
                reverse ? DnaQuerySeedStrand::Reverse
                        : DnaQuerySeedStrand::Forward,
                offset, offset + span, actual_downsample,
                reverse ? scratch.shared_window_rc_seeds
                        : scratch.shared_window_fwd_seeds,
                refine_views);
        // A window with no voter votes nothing.
        if (!views) {
          lane_scratch.retained_seeds.clear();
          continue;
        }
        lane_scratch.pending_exact_refine_views =
            refine_views ? *refine_views : std::vector<DnaLongSeedView>{};
        WindowAnchorPeakParams params;
        params.lookup_cache = lookup_cache;
        params.seed_view_override = views;
        params.tie_seed = dctx.vote_tie_seed;
        std::vector<VotePeak> peaks = window_anchor_peaks(
            seed_ctx, reverse ? nullptr : fwd_enc.data() + begin, span,
            reverse, lane_scratch, params);
        stretch_to_contig(dctx, read_len, offset, span, peaks,
                          lane_scratch.retained_seeds);
        stretch_peaks.insert(stretch_peaks.end(), peaks.begin(), peaks.end());
        retained[lane].insert(retained[lane].end(),
                              lane_scratch.retained_seeds.begin(),
                              lane_scratch.retained_seeds.end());
      }
      if (stretch_peaks.empty()) continue;
      raw.insert(raw.end(), stretch_peaks.begin(), stretch_peaks.end());
      // Its admission compares masks on the stretch's own tile grid.
      DnaPlacementFamily one =
          build_dna_placement_family(dctx, fwd_enc, stretch_peaks,
                                     scratch.fwd, scratch.rc, begin, span);
      // Its candidates take the contig's next ids: a contig they outgrow is
      // refused.
      const std::size_t candidates =
          stretch_family.candidates.size() + one.candidates.size();
      if (candidates > static_cast<std::size_t>(std::numeric_limits<
                           ::fa::cpu::voting::CandidateId>::max()))
        throw CandidateIdOverflow(dctx.read_name_hash, read_len, candidates,
                                  dctx.opts.catalogue_lane_bound);
      append_stretch_family(stretch_family, std::move(one));
    }
    // The reverse windows come in descending read order.
    std::stable_sort(retained[1].begin(), retained[1].end(),
                     [](const ChainWindowRetainedSeed& left,
                        const ChainWindowRetainedSeed& right) {
                       return left.seed.read_pos < right.seed.read_pos;
                     });
    scratch.fwd.retained_seeds.swap(retained[0]);
    scratch.rc.retained_seeds.swap(retained[1]);
    stage_release.after_last_prepare();
  } else if (read_len >= dctx.opts.k) {
    // The vote runs over one whole-read window [0, read_len).
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
    stage_release.after_last_prepare();
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
    // Under a vote per stretch the stretches' catalogues, built above.
    placement_family =
        stretch > 0 ? std::move(stretch_family)
                    : build_dna_placement_family(dctx, fwd_enc, raw,
                                                 scratch.fwd, scratch.rc);
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
      stage_release.at(DnaStagePoint::Retained);
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

/**** kept_scratch ****/
// dna.h kept_scratch: the bytes a mapping thread keeps across reads and their
// release after a read.

namespace fa::cpu::lr::dna {

std::vector<std::uint8_t>& ref_window_buffer() {
  static thread_local std::vector<std::uint8_t> buffer;
  return buffer;
}

namespace {

std::size_t kept_scratch_bytes(const DnaWorkerScratch& worker) {
  return worker.chain.heap_bytes() + worker.seed_density.heap_bytes() +
         ::fa::cpu::chaining::dense_chain_scratch_bytes() +
         pool_sort_scratch_bytes() + pass_slices_bytes() +
         ::fa::cpu::dp::thread_arena_bytes() + ava_read_targets_bytes() +
         ava_target_scratch_bytes() + capacity_bytes(ref_window_buffer());
}

void release_kept_scratch(DnaWorkerScratch& worker) {
  worker.chain = ChainAnchorScratch{};
  worker.seed_density = RetainedSeedDensity{};
  ::fa::cpu::chaining::release_dense_chain_scratch();
  release_placement_scratch();
  ::fa::cpu::dp::release_thread_arena();
  release_ava_scratch();
  ref_window_buffer() = std::vector<std::uint8_t>{};
}

} // namespace

void settle_kept_scratch(DnaWorkerScratch& worker) {
  if (kept_scratch_bytes(worker) > kKeptScratchReleaseBytes)
    release_kept_scratch(worker);
}

} // namespace fa::cpu::lr::dna
