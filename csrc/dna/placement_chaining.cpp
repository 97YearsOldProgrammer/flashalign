#include "placement_chaining.h"
#include "placement_chaining_internal.h"
#include "family_realization.h"  // plan_dna_family_join (the span rule)

#include "../chaining/colinear_chain.h"
#include "../chaining/dense_chain.h"
#include "../core/radix_sort.h"
#include "../chaining/partition.h"
#include "../dp/params.h"  // DpMapOpt::max_sw_mat (the span rule)
#include "../index/format.h"
#include "../seeding/tie_hash.h"  // tie_locus_hash (the residue clusters' tie-break)
#include "../split/query_geometry.h"
#include "../voting/query_tiles.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

constexpr int kCandidateReachDivisor = 4;
constexpr int kCandidateReachCap = 131072;

int forward_query_position(int oriented_seed_position, bool reverse,
                           int read_length, int seed_length) noexcept {
  return reverse ? read_length - seed_length - oriented_seed_position
                 : oriented_seed_position;
}

struct DeferredSlice {
  RetainedSeedRef seed;
  KmerPostingIntervalView interval;
};

::fa::cpu::voting::QueryTileMask valid_tiles(int read_length, int seed_length,
                                             int tile_count) {
  ::fa::cpu::voting::QueryTileMask result;
  for (int tile = 0; tile < tile_count; ++tile) {
    if (::fa::cpu::voting::query_tile_end(tile + 1, read_length, seed_length,
                                          tile_count) >
        ::fa::cpu::voting::query_tile_begin(tile, read_length, seed_length,
                                            tile_count)) {
      result.set(tile);
    }
  }
  return result;
}

int selected_tiles(
    const ::fa::cpu::voting::QueryPartitionPath& path,
    ::fa::cpu::voting::CandidateId candidate) {
  return static_cast<int>(std::count(
      path.assignment.begin(), path.assignment.end(), candidate));
}

::fa::cpu::voting::QueryPartitionResult repartition(
    const DnaContext& context, const DnaPlacementFamily& family) {
  ::fa::cpu::voting::QueryPartitionProblem problem;
  problem.tile_count = family.tile_count;
  problem.valid_tiles =
      valid_tiles(family.read_length, family.seed_length, family.tile_count);
  problem.parameters = context.opts.query_partition;
  problem.catalogue.candidates.reserve(family.candidates.size());
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    problem.catalogue.candidates.push_back(
        {candidate.id, candidate.equivalence_key, candidate.lane,
         candidate.catalogue_rank, candidate.vote_evidence,
         candidate.screening_chain_score, candidate.support});
  }
  return ::fa::cpu::voting::solve_query_partition(problem);
}

// Per-thread radix ping-pong buffer and histograms for the pool sort below.
::fa::cpu::radix::Scratch<chaining::Anchor>& pool_sort_scratch() {
  static thread_local ::fa::cpu::radix::Scratch<chaining::Anchor> scratch;
  return scratch;
}

// Per-thread slices of one chain_candidate pass, by density entry: an entry's
// interval is this pass's while its stamp equals `pass`.
struct PassSlices {
  std::vector<std::uint32_t> stamp;
  std::vector<KmerPostingIntervalView> interval;
  std::uint32_t pass = 0;
};
PassSlices& pass_slices(std::size_t entries) {
  static thread_local PassSlices slices;
  if (slices.stamp.size() < entries) {
    slices.stamp.resize(entries, 0);
    slices.interval.resize(entries);
  }
  if (++slices.pass == 0) {
    std::fill(slices.stamp.begin(), slices.stamp.end(), 0u);
    slices.pass = 1;
  }
  return slices;
}

// Sorts the pool by (r, q, span) ascending, flags descending, as a stable LSD
// radix sort. Flags descend so the std::unique that follows, which keeps the
// first of equal (r, q, span) anchors, keeps the tandem flag. A pass whose key
// is the same for every anchor would leave the order as it is and is skipped.
void sort_pool_anchors(std::vector<chaining::Anchor>& anchors) {
  if (anchors.size() < 2) return;
  ::fa::cpu::radix::Scratch<chaining::Anchor>& scratch = pool_sort_scratch();
  const chaining::Anchor first = anchors.front();
  std::uint32_t flags_differ = 0;
  std::int32_t span_differ = 0;
  for (const chaining::Anchor& anchor : anchors) {
    flags_differ |= anchor.flags ^ first.flags;
    span_differ |= anchor.span ^ first.span;
  }
  if (flags_differ != 0)
    ::fa::cpu::radix::stable_sort_u32(anchors, scratch,
                                      [](const chaining::Anchor& anchor) {
                                        return ~anchor.flags;
                                      });
  if (span_differ != 0)
    ::fa::cpu::radix::stable_sort_u32(
        anchors, scratch,
        [](const chaining::Anchor& anchor) { return anchor.span; });
  ::fa::cpu::radix::stable_sort_u32(
      anchors, scratch, [](const chaining::Anchor& anchor) { return anchor.q; });
  ::fa::cpu::radix::stable_sort_u32(
      anchors, scratch, [](const chaining::Anchor& anchor) { return anchor.r; });
}

}  // namespace

namespace internal {

void append_interval_anchors(
    const KmerPostingIntervalView& interval,
    const RetainedSeedRef& seed,
    std::uint64_t chromosome_base,
    int chromosome_length,
    int main_diagonal,
    int seed_length,
    int diagonal_band,
    int tandem_window,
    bool reverse_lane,
    bool skip_own_diagonal,
    int query_length,
    std::vector<chaining::Anchor>& anchors,
    DnaPlacementCandidateChain& record) {
  const std::uint64_t chromosome_end =
      chromosome_base + static_cast<std::uint64_t>(chromosome_length);
  record.interval_hits += interval.count;
  for (std::uint32_t posting = 0; posting < interval.count; ++posting) {
    const std::uint64_t global = interval.positions[posting];
    if (global < chromosome_base) {
      ++record.filtered_hits;
      continue;
    }
    const std::uint64_t local64 = global - chromosome_base;
    if (local64 > static_cast<std::uint64_t>(
                      std::numeric_limits<std::int32_t>::max())) {
      ++record.filtered_hits;
      continue;
    }
    const int reference_position = static_cast<int>(local64);
    const int query_position = seed.seed.read_pos;
    if (skip_own_diagonal && reference_position == query_position) {
      ++record.filtered_hits;
      continue;
    }
    const bool geometry_ok =
        query_position >= 0 &&
        query_position + seed_length <= query_length &&
        reference_position >= 0 &&
        reference_position + seed_length <= chromosome_length &&
        within_band(reference_position - query_position, main_diagonal,
                    diagonal_band);
    // The posting was found under this seed's canonical key, so the k-mers
    // match exactly when the orientation bits agree with the lane.
    const bool verified =
        geometry_ok && packed_ref_orientation_compatible(
                           seed.seed.z,
                           interval.positions.packed_at(posting), reverse_lane);
    if (!verified) {
      ++record.filtered_hits;
      continue;
    }
    const bool tandem = posting_is_tandem(
        posting, 0u, interval.count, global,
        [&interval](std::uint32_t at) -> std::uint64_t {
          return interval.positions[at];
        },
        chromosome_base, chromosome_end, tandem_window);
    const chaining::Anchor anchor{
        reference_position, query_position, seed_length,
        static_cast<std::int32_t>(anchors.size()),
        tandem ? static_cast<std::uint32_t>(chaining::ANCHOR_TANDEM) : 0u};
    anchors.push_back(anchor);
  }
}

}  // namespace internal

namespace {

void set_spans(DnaPlacementCandidateChain& record, bool reverse,
               int read_length) {
  if (record.primary.empty()) return;
  int query_begin = std::numeric_limits<int>::max();
  int query_end = std::numeric_limits<int>::min();
  int reference_begin = std::numeric_limits<int>::max();
  int reference_end = std::numeric_limits<int>::min();
  for (const chaining::Anchor& anchor : record.primary) {
    query_begin = std::min(query_begin, static_cast<int>(anchor.q));
    query_end = std::max(query_end, static_cast<int>(anchor.q_end()));
    reference_begin =
        std::min(reference_begin, static_cast<int>(anchor.r));
    reference_end =
        std::max(reference_end, static_cast<int>(anchor.r_end()));
  }
  record.oriented_query_begin = query_begin;
  record.oriented_query_end = query_end;
  record.forward_query_begin =
      reverse ? read_length - query_end : query_begin;
  record.forward_query_end =
      reverse ? read_length - query_begin : query_end;
  record.reference_begin = reference_begin;
  record.reference_end = reference_end;
}

enum class CandidateChainPass : std::uint8_t {
  BoundedScreening,
  WholeQueryExact,
};

}  // namespace

namespace internal {

void sort_unique_pool_anchors(std::vector<chaining::Anchor>& anchors) {
  sort_pool_anchors(anchors);
  anchors.erase(
      std::unique(
          anchors.begin(), anchors.end(),
          [](const chaining::Anchor& left,
             const chaining::Anchor& right) {
            return left.r == right.r && left.q == right.q &&
                   left.span == right.span;
          }),
      anchors.end());
}

void write_primary_chain(const DnaPlacementFamily& family, bool reverse,
                         const chaining::ChainResult& chained,
                         const chaining::Chain& chain,
                         DnaPlacementCandidateChain& record) {
  record.primary.clear();
  record.primary.reserve(chain.idx.size());
  for (const std::int32_t index : chain.idx)
    record.primary.push_back(chained.anchors[static_cast<std::size_t>(index)]);
  record.chain_score = chain.score;
  record.chain_anchors = static_cast<int>(record.primary.size());
  record.dense_support = {};
  for (const chaining::Anchor& anchor : record.primary)
    record.dense_support.set(
        dna_forward_query_tile(anchor.q, reverse, family.read_length,
                               family.seed_length, family.tile_count));
  set_spans(record, reverse, family.read_length);
}

}  // namespace internal

namespace {

// Writes into `to` what a chain_candidate call wrote into `from`, where the
// call returned `accepted`. The call leaves the record's candidate,
// sparse_support, screening_* and selected-tile counts alone, and its spans
// unless accepted.
void copy_chain_call(const DnaPlacementCandidateChain& from, bool accepted,
                     DnaPlacementCandidateChain& to) {
  DnaPlacementCandidateChain kept = std::move(to);
  to = from;
  to.candidate = kept.candidate;
  to.sparse_support = kept.sparse_support;
  to.screening_chain_score = kept.screening_chain_score;
  to.screening_chain_anchors = kept.screening_chain_anchors;
  to.screening_forward_query_begin = kept.screening_forward_query_begin;
  to.screening_forward_query_end = kept.screening_forward_query_end;
  to.initial_selected_tiles = kept.initial_selected_tiles;
  to.final_selected_tiles = kept.final_selected_tiles;
  if (!accepted) {
    to.oriented_query_begin = kept.oriented_query_begin;
    to.oriented_query_end = kept.oriented_query_end;
    to.forward_query_begin = kept.forward_query_begin;
    to.forward_query_end = kept.forward_query_end;
    to.reference_begin = kept.reference_begin;
    to.reference_end = kept.reference_end;
  }
}

bool chain_candidate(
    const DnaContext& context,
    const DnaPlacementFamily& family,
    const DnaPlacementCandidate& candidate,
    const std::vector<std::uint8_t>& query,
    RetainedSeedDensity& seed_index,
    DnaPlacementCandidateChain& record,
    CandidateChainPass pass) {
  // The whole-query pass defers every admissible slice and builds its anchors
  // after the scan. The screening pass defers nothing and uses chain_colinear.
  const bool whole_query_exact = pass != CandidateChainPass::BoundedScreening;
  record.dense_support = {};
  record.primary.clear();
  record.sibling_paths.clear();
  record.sibling_scores.clear();
  record.rival_sibling = -1;
  record.status = DnaPlacementChainStatus::NotSelected;
  record.interval_hits = 0;
  record.filtered_hits = 0;
  record.chain_score = 0;
  record.rival_chain_score = 0;
  record.overlapping_rivals = 0;
  record.chain_anchors = 0;
  record.exact = false;
  record.sparse_anchors = 0;
  record.deferred_anchors = 0;
  record.rescued_anchors = 0;
  record.dense_runs = 0;
  const bool reverse = candidate.peak.is_rc;
  const std::vector<RetainedSeedRef>& seeds =
      whole_query_exact
          ? (reverse ? seed_index.fine_reverse()
                     : seed_index.fine_forward())
          : (reverse ? seed_index.reverse() : seed_index.forward());
  const DnaPlacementChainStatus refusal = internal::contig_lane_refusal(
      context, candidate.peak.chr, seeds, query, family.read_length);
  if (refusal != DnaPlacementChainStatus::Accepted) {
    record.status = refusal;
    record.exact = whole_query_exact;
    return false;
  }

  const int seed_length = family.seed_length;
  const int chromosome_length =
      static_cast<int>(context.ref.contig_length(candidate.peak.chr));
  const std::int64_t expected = candidate.peak.raw_ref_start;
  const internal::HarvestWindow window = internal::harvest_window(
      context, candidate.peak, family.read_length, chromosome_length);
  if (window.high <= window.low) {
    record.status = DnaPlacementChainStatus::InvalidReference;
    record.exact = whole_query_exact;
    return false;
  }
  const std::uint64_t* offsets = context.ref.index->chrom_offsets_data();
  const std::uint64_t chromosome_base =
      offsets[static_cast<std::size_t>(candidate.peak.chr)];
  const internal::PoolGate gate =
      internal::pool_gate(context, whole_query_exact);
  const int diagonal_band =
      internal::pass_diagonal_band(context, whole_query_exact);
  const bool skip_own_diagonal = internal::skips_own_diagonal(
      context, candidate.peak.chr, reverse);

  std::vector<chaining::Anchor> sparse;
  std::vector<DeferredSlice> deferred;
  // The whole-query pass builds no anchors in the loop below; every
  // admissible slice is deferred instead.
  const bool defer_every_slice = whole_query_exact;
  sparse.reserve(defer_every_slice ? 0 : seeds.size());
  deferred.reserve(defer_every_slice ? seeds.size() : seeds.size() / 8 + 1);
  // Every seed resolves against the same (chr, low, high) window, so each
  // distinct entry is sliced once, in batches through slice_batch; seeds are
  // still processed in order.
  constexpr std::size_t kSliceBlock = RetainedSeedDensity::kSliceBatch;
  PassSlices& slices = pass_slices(seed_index.entry_count());
  std::uint32_t block_entries[kSliceBlock];
  KmerPostingIntervalView block_intervals[kSliceBlock];
  for (std::size_t next = 0; next < seeds.size();) {
    // Take seeds until kSliceBlock entries new to this pass are collected. A
    // skipped entry gets an interval that is not found, which the loop below
    // drops as it dropped the gated slice.
    const std::size_t first = next;
    std::size_t block_size = 0;
    for (; next < seeds.size(); ++next) {
      const std::uint32_t entry = seeds[next].entry;
      std::uint32_t& stamp = slices.stamp[entry];
      if (stamp == slices.pass) continue;
      if (gate.skips(seed_index, entry)) {
        stamp = slices.pass;
        slices.interval[entry] = {};
        continue;
      }
      if (block_size == kSliceBlock) break;
      stamp = slices.pass;
      block_entries[block_size++] = entry;
    }
    seed_index.slice_batch(
        *context.ref.index, block_entries, block_size, candidate.peak.chr,
        static_cast<std::uint32_t>(window.low),
        static_cast<std::uint32_t>(window.high), block_intervals);
    for (std::size_t lane = 0; lane < block_size; ++lane)
      slices.interval[block_entries[lane]] = block_intervals[lane];
    for (std::size_t at = first; at < next; ++at) {
      const RetainedSeedRef& seed = seeds[at];
      const KmerPostingIntervalView& interval = slices.interval[seed.entry];
      // A rescued seed passes the gate (RetainedSeedRef::rescued): the vote
      // seed on the screening pass, its fine twin on the whole-query pass. Its
      // key's interval is shared with the seeds that do not.
      const bool over_pool_gate =
          gate.drops(interval.global_count, seed.rescued);
      if (!interval.found() || interval.count == 0 || over_pool_gate) {
        if (interval.found() && over_pool_gate)
          record.filtered_hits += static_cast<int>(interval.count);
        continue;
      }
      // The screening pass builds its anchors directly; the whole-query pass
      // defers every slice and builds its anchors after the scan.
      if (whole_query_exact) {
        deferred.push_back({seed, interval});
        record.deferred_anchors += interval.count;
      } else {
        internal::append_interval_anchors(
            interval, seed, chromosome_base, chromosome_length,
            static_cast<int>(expected), seed_length, diagonal_band,
            context.opts.dna_tandem_window, reverse, skip_own_diagonal,
            family.read_length, sparse, record);
      }
    }
  }
  record.sparse_anchors = sparse.size();

  std::vector<chaining::Anchor> anchors = std::move(sparse);
  if (whole_query_exact) {
    // Restore every deferred slice, in order.
    for (std::size_t slot = 0; slot < deferred.size(); ++slot) {
      // Prefetch the next slice's first postings, which the binary search
      // that bounded it never touched.
      if (slot + 1 < deferred.size()) {
        const KmerPostingIntervalView& next = deferred[slot + 1].interval;
        __builtin_prefetch(next.positions.data(), 0, 1);
        if (next.count > 8)
          __builtin_prefetch(next.positions.data() + 8, 0, 1);
      }
      const DeferredSlice& item = deferred[slot];
      const std::size_t before = anchors.size();
      internal::append_interval_anchors(
          item.interval, item.seed, chromosome_base, chromosome_length,
          static_cast<int>(expected), seed_length, diagonal_band,
          context.opts.dna_tandem_window, reverse, skip_own_diagonal,
          family.read_length, anchors, record);
      record.rescued_anchors += anchors.size() - before;
    }
  }
  if (anchors.empty()) {
    record.status = whole_query_exact
                        ? (record.interval_hits == 0
                               ? DnaPlacementChainStatus::NoIntervalHits
                               : DnaPlacementChainStatus::NoChain)
                        : DnaPlacementChainStatus::NotSelected;
    record.exact = whole_query_exact;
    return false;
  }
  internal::sort_unique_pool_anchors(anchors);

  const chaining::ColinearChainParams chain_params =
      dna_candidate_chain_params(context, diagonal_band, seed_length,
                                 family.read_length);
  chaining::DenseChainStats dense_stats;
  chaining::DenseChainParams dense_params = dna_dense_chain_params(
      chain_params, seed_length, context.opts.dna_dense_diag_min_runs);
  // The whole-query pass runs the dense chain (the pool collapsed to exact
  // diagonal runs); the screening pass runs the plain colinear chain.
  const chaining::ChainResult chained =
      whole_query_exact
          ? chaining::chain_dense_colinear(std::move(anchors), dense_params,
                                           &dense_stats)
          : chaining::chain_colinear(std::move(anchors), chain_params);
  record.dense_runs = dense_stats.runs;
  const chaining::ChainPartition partition =
      chaining::partition_chains(chained, ::fa::cpu::split::kHalfFloor);
  if (partition.primary < 0) {
    record.status = DnaPlacementChainStatus::NoChain;
    record.exact = whole_query_exact;
    return false;
  }
  const chaining::Chain& primary =
      chained.chains[static_cast<std::size_t>(partition.primary)];
  internal::write_primary_chain(family, reverse, chained, primary, record);
  record.rival_chain_score = partition.f2;
  record.overlapping_rivals = partition.n_sub;
  if (whole_query_exact) {
    for (std::size_t which = 0; which < chained.chains.size(); ++which) {
      if (static_cast<int>(which) == partition.primary) continue;
      const chaining::Chain& sibling = chained.chains[which];
      if (sibling.idx.size() < 2) continue;
      std::vector<chaining::Anchor> path;
      path.reserve(sibling.idx.size());
      for (const std::int32_t index : sibling.idx)
        path.push_back(chained.anchors[static_cast<std::size_t>(index)]);
      // Record which sibling the partition's f2 came from, so the MAPQ can
      // realize it.
      if (static_cast<int>(which) == partition.f2_index)
        record.rival_sibling = static_cast<int>(record.sibling_paths.size());
      record.sibling_paths.push_back(std::move(path));
      record.sibling_scores.push_back(sibling.score);
    }
  }
  record.exact = whole_query_exact;
  record.status = whole_query_exact ? DnaPlacementChainStatus::Accepted
                                : DnaPlacementChainStatus::Sparse;
  return true;
}

// The span rule: which query tiles an accepted whole-query chain owns when
// the partition is re-solved (DnaTileOwnership). Under AnchorTiles an owner
// owns the tiles holding its anchors. Under Span it also owns every tile its
// chain traverses between its first and last anchor, as a minimap2 hit owns
// the query span it covers; otherwise a read crossing a high-copy stretch,
// with anchors in few tiles, hands the rest to another repeat copy. A
// traversed stretch stays with a rival on the owner's contig and strand in two
// cases:
//   1. Bridge. The rival is the same locus past a deletion, and the
//      realization would join it to the owner's block as one record with a
//      clean deletion. Each side where the rival's stretch meets an
//      uncontested owner tile is a seam, from the owner's nearest anchor to
//      the rival's first or last anchor inside the stretch, checked with
//      plan_dna_family_join. A stretch with at least one seam, all of them
//      eligible, keeps its rival; the bridge DP is not run.
//   2. Gap. The owner's chain crosses the stretch with a link whose reference
//      gap exceeds cigar_dp_max_gap, which one block could not fill. A long
//      query gap does not count: a block fills it as one insertion.
// The partition has no blocks yet, so a rival's block is approximated by a
// stretch of the owner's span without any uncontested owner anchor tile
// (anchor tiles the rival does not hold). A rival holds its vote tiles and its
// current chain's anchor tiles, never tiles widened by its own span.
// Owners take their supports after the sweep that chained them, so each
// sees every rival's post-sweep chain.

// One index per tile, -1 for none; the first 128 tiles allocate nothing.
struct TileAnchorIndices {
  std::array<int, ::fa::cpu::voting::kQueryTileCount> local;
  std::vector<int> extra;

  explicit TileAnchorIndices(int tile_count) {
    local.fill(-1);
    if (tile_count > static_cast<int>(local.size()))
      extra.assign(static_cast<std::size_t>(tile_count) - local.size(), -1);
  }
  int& operator[](std::size_t tile) {
    return tile < local.size() ? local[tile] : extra[tile - local.size()];
  }
  int operator[](std::size_t tile) const {
    return tile < local.size() ? local[tile] : extra[tile - local.size()];
  }
};

// Per forward-query tile, the anchors of one chain with the smallest and the
// largest oriented query start (indices into the chain, -1 for none).
struct TileAnchorEnds {
  TileAnchorIndices first;
  TileAnchorIndices last;
};

TileAnchorEnds tile_anchor_ends(const std::vector<chaining::Anchor>& anchors,
                                bool reverse, int read_length, int seed_length,
                                int tile_count) {
  TileAnchorEnds ends{TileAnchorIndices(tile_count),
                      TileAnchorIndices(tile_count)};
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    const int tile = dna_forward_query_tile(
        anchors[index].q, reverse, read_length, seed_length, tile_count);
    if (tile < 0 || tile >= tile_count) continue;
    const int at = static_cast<int>(index);
    const std::size_t tile_index = static_cast<std::size_t>(tile);
    if (ends.first[tile_index] < 0 ||
        anchors[index].q <
            anchors[static_cast<std::size_t>(ends.first[tile_index])].q)
      ends.first[tile_index] = at;
    if (ends.last[tile_index] < 0 ||
        anchors[index].q >
            anchors[static_cast<std::size_t>(ends.last[tile_index])].q)
      ends.last[tile_index] = at;
  }
  return ends;
}

// One side of a bridge seam: a forward-read range and the oriented-first and
// oriented-last anchors of the piece inside it.
struct SpanSeamPiece {
  int forward_begin = 0;
  int forward_end = 0;
  const chaining::Anchor* first = nullptr;
  const chaining::Anchor* last = nullptr;
};

// plan_dna_family_join over the seam from `earlier` to `later` (oriented
// query order), both on `chromosome` and strand `reverse`. The join rule
// takes its two blocks in forward-read order, so a reverse-strand seam is
// passed later-first; either way the gaps it measures are later.first minus
// earlier.last on both axes.
bool span_seam_bridges(const DnaContext& context, int chromosome, bool reverse,
                       const SpanSeamPiece& earlier,
                       const SpanSeamPiece& later) {
  const SpanSeamPiece& left = reverse ? later : earlier;
  const SpanSeamPiece& right = reverse ? earlier : later;
  const ::fa::cpu::DpMapOpt opt;  // bridge_geometry's matrix cap
  return plan_dna_family_join(DnaFamilyJoinGeometry{
             chromosome, chromosome, reverse, reverse, left.forward_begin,
             left.forward_end, right.forward_begin, right.forward_end,
             left.first->q, left.last->q_end(), right.first->q,
             right.last->q_end(), left.first->r, left.last->r_end(),
             right.first->r, right.last->r_end(),
             context.opts.cigar_dp_max_gap, context.opts.cigar_dp_max_gap,
             opt.max_sw_mat})
             .kind != DnaFamilyJoinKind::NotEligible;
}

// The tiles an accepted owner owns under DnaTileOwnership::Span: its anchor
// tiles, plus every tile its chain traverses that neither exception leaves
// with a same-locus rival.
::fa::cpu::voting::QueryTileMask span_owned_tiles(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& result,
    const DnaPlacementCandidate& owner,
    const DnaPlacementCandidateChain& chain) {
  namespace voting = ::fa::cpu::voting;
  const int tiles = family.tile_count;
  const voting::QueryTileMask& anchor_tiles = chain.dense_support;
  int first = -1;
  int last = -1;
  for (int tile = 0; tile < tiles; ++tile) {
    if (!anchor_tiles.test(tile)) continue;
    if (first < 0) first = tile;
    last = tile;
  }
  if (first < 0 || chain.primary.empty()) return anchor_tiles;
  const bool reverse = owner.peak.is_rc;
  const int chromosome = owner.peak.chr;
  const int read_length = family.read_length;
  const int seed_length = family.seed_length;
  const int max_gap = context.opts.cigar_dp_max_gap;
  const auto forward_begin = [&](int tile) {
    return tile >= tiles ? read_length
                          : voting::query_tile_begin(tile, read_length,
                                                     seed_length, tiles);
  };
  const auto forward_tile = [&](const chaining::Anchor& anchor) {
    return dna_forward_query_tile(anchor.q, reverse, read_length, seed_length,
                                  tiles);
  };
  const std::vector<chaining::Anchor>& own = chain.primary;
  const TileAnchorEnds own_ends =
      tile_anchor_ends(own, reverse, read_length, seed_length, tiles);
  const auto own_at = [&](int index) -> const chaining::Anchor& {
    return own[static_cast<std::size_t>(index)];
  };

  // The gap exception's links: consecutive anchors of the owner's chain
  // whose reference gap exceeds the join limit, as the forward tiles they
  // join (possibly one tile: a long deletion moves only the reference). The
  // chain must be in ascending query order.
  std::vector<std::pair<int, int>> long_links;
  bool ascending = true;
  for (std::size_t index = 1; index < own.size(); ++index) {
    const chaining::Anchor& from = own[index - 1];
    const chaining::Anchor& to = own[index];
    if (to.q < from.q) {
      ascending = false;
      break;
    }
    if (to.r - from.r_end() > max_gap) {
      const int a = forward_tile(from);
      const int b = forward_tile(to);
      long_links.emplace_back(std::min(a, b), std::max(a, b));
    }
  }
  if (!ascending) long_links.clear();

  voting::QueryTileMask keep;
  for (const DnaPlacementCandidate& rival : family.candidates) {
    if (rival.id == owner.id || rival.peak.chr != chromosome ||
        rival.peak.is_rc != reverse)
      continue;
    const DnaPlacementCandidateChain* rival_chain = result.find(rival.id);
    if (rival_chain == nullptr) continue;
    const voting::QueryTileMask held =
        rival_chain->sparse_support | rival_chain->dense_support;
    const auto at_stake = [&](int tile) {
      return !anchor_tiles.test(tile) && held.test(tile);
    };
    bool any = false;
    for (int tile = first; tile <= last && !any; ++tile) any = at_stake(tile);
    if (!any) continue;
    const std::vector<chaining::Anchor>& theirs = rival_chain->primary;
    const TileAnchorEnds their_ends =
        tile_anchor_ends(theirs, reverse, read_length, seed_length, tiles);
    const auto uncontested = [&](int tile) {
      return anchor_tiles.test(tile) && !held.test(tile);
    };
    for (int begin = first; begin <= last;) {
      if (uncontested(begin)) {
        ++begin;
        continue;
      }
      int end = begin;
      while (end < last && !uncontested(end + 1)) ++end;
      bool stake = false;
      int rival_first = -1;
      int rival_last = -1;
      for (int tile = begin; tile <= end; ++tile) {
        stake = stake || at_stake(tile);
        const int head = their_ends.first[static_cast<std::size_t>(tile)];
        const int tail = their_ends.last[static_cast<std::size_t>(tile)];
        if (head >= 0 &&
            (rival_first < 0 ||
             theirs[static_cast<std::size_t>(head)].q <
                 theirs[static_cast<std::size_t>(rival_first)].q))
          rival_first = head;
        if (tail >= 0 &&
            (rival_last < 0 ||
             theirs[static_cast<std::size_t>(tail)].q >
                 theirs[static_cast<std::size_t>(rival_last)].q))
          rival_last = tail;
      }
      if (!stake) {
        begin = end + 1;
        continue;
      }
      // An uncontested owner tile bounds the stretch on each side where it
      // does not reach the span's end.
      const bool left_seam = begin > first;
      const bool right_seam = end < last;
      bool bridges = rival_first >= 0 && (left_seam || right_seam);
      if (bridges) {
        const SpanSeamPiece stretch{
            forward_begin(begin), forward_begin(end + 1),
            &theirs[static_cast<std::size_t>(rival_first)],
            &theirs[static_cast<std::size_t>(rival_last)]};
        if (left_seam) {
          const int tile = begin - 1;
          const chaining::Anchor& near = own_at(
              reverse ? own_ends.first[static_cast<std::size_t>(tile)]
                      : own_ends.last[static_cast<std::size_t>(tile)]);
          const SpanSeamPiece piece{forward_begin(tile),
                                    forward_begin(tile + 1), &near, &near};
          bridges = reverse ? span_seam_bridges(context, chromosome, reverse,
                                                stretch, piece)
                            : span_seam_bridges(context, chromosome, reverse,
                                                piece, stretch);
        }
        if (bridges && right_seam) {
          const int tile = end + 1;
          const chaining::Anchor& near = own_at(
              reverse ? own_ends.last[static_cast<std::size_t>(tile)]
                      : own_ends.first[static_cast<std::size_t>(tile)]);
          const SpanSeamPiece piece{forward_begin(tile),
                                    forward_begin(tile + 1), &near, &near};
          bridges = reverse ? span_seam_bridges(context, chromosome, reverse,
                                                piece, stretch)
                            : span_seam_bridges(context, chromosome, reverse,
                                                stretch, piece);
        }
      }
      const bool gap = std::any_of(
          long_links.begin(), long_links.end(),
          [&](const std::pair<int, int>& link) {
            return link.first <= end && link.second >= begin;
          });
      if (bridges || gap) {
        for (int tile = begin; tile <= end; ++tile)
          if (at_stake(tile)) keep.set(tile);
      }
      begin = end + 1;
    }
  }

  voting::QueryTileMask owned = anchor_tiles;
  for (int tile = first; tile <= last; ++tile) {
    if (!anchor_tiles.test(tile) && !keep.test(tile))
      owned.set(tile);
  }
  return owned;
}

// An owner a sweep chained, settled after the sweep.
struct SweptOwner {
  DnaPlacementCandidate* candidate = nullptr;
  DnaPlacementCandidateChain* chain = nullptr;
  bool accepted = false;
};

// Stabilization: chain every selected owner that is not yet exact over the
// whole query, give each its support under `ownership`, re-solve the
// partition once, chain the owners the re-solve newly selected, and accept the
// family when every owner has an accepted whole-query chain. A family with an
// owner lacking one stays unaccepted, and the read is unmapped.
bool stabilize_selected_family(
    const DnaContext& context, DnaPlacementFamily& family,
    DnaPlacementChainingResult& result, RetainedSeedDensity& seed_index,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    DnaTileOwnership ownership) {
  bool resolved = false;
  for (;;) {
    bool changed = false;
    std::vector<SweptOwner> swept;
    for (const auto& block : family.partition.selected.blocks) {
      if (block.candidate == ::fa::cpu::voting::kNullCandidate) continue;
      auto record_it = std::find_if(
          result.candidates.begin(), result.candidates.end(),
          [&](const DnaPlacementCandidateChain& record) {
            return record.candidate == block.candidate;
          });
      auto candidate_it = std::find_if(
          family.candidates.begin(), family.candidates.end(),
          [&](const DnaPlacementCandidate& candidate) {
            return candidate.id == block.candidate;
          });
      if (record_it == result.candidates.end() ||
          candidate_it == family.candidates.end())
        return false;
      // An owner chained over the whole query is never chained again.
      if (record_it->exact) continue;
      const bool accepted = chain_candidate(
          context, family, *candidate_it,
          candidate_it->peak.is_rc ? reverse_query : forward_query,
          seed_index, *record_it, CandidateChainPass::WholeQueryExact);
      swept.push_back({&*candidate_it, &*record_it, accepted});
      changed = true;
    }
    // Supports are set after the sweep, since the span rule reads the
    // rivals' new chains: the vote tiles plus the owned tiles, or the vote
    // tiles alone when the chain was not accepted.
    for (const SweptOwner& owner : swept) {
      ::fa::cpu::voting::QueryTileMask restored;
      if (owner.accepted) {
        restored = ownership == DnaTileOwnership::Span
                       ? span_owned_tiles(context, family, result,
                                          *owner.candidate, *owner.chain)
                       : owner.chain->dense_support;
        if (restored.count() > owner.chain->dense_support.count())
          result.span_widened = true;
      }
      owner.candidate->support = owner.chain->sparse_support | restored;
    }
    if (changed && !resolved) {
      resolved = true;
      family.partition = repartition(context, family);
      continue;
    }
    result.accepted = true;
    for (const auto& block : family.partition.selected.blocks) {
      if (block.candidate == ::fa::cpu::voting::kNullCandidate) continue;
      const DnaPlacementCandidateChain* record = result.find(block.candidate);
      if (record == nullptr || !record->exact ||
          record->status != DnaPlacementChainStatus::Accepted)
        result.accepted = false;
    }
    return true;
  }
}

// Terminal-clip recovery, run after the family is committed, from the kept
// fine seeds and the posting views cached per read.

std::int64_t residue_floor_div(std::int64_t value, std::int64_t width) {
  const std::int64_t quotient = value / width;
  return (value % width != 0 && (value < 0) != (width < 0)) ? quotient - 1
                                                            : quotient;
}

// One query interval in the forward-query frame.
struct ResidueRun {
  int query_begin_bp = 0;
  int query_end_bp = 0;
};

// Distinct fine-seed keys inside the run whose cached view is under the
// occurrence cap, counted up to kDnaResidueMinClusterAnchors: the caller
// tests only that threshold. The index is never queried: recovery only rereads
// postings already fetched.
int residue_cached_supply(const DnaPlacementFamily& family,
                          const std::vector<RetainedSeedRef>& fine_forward,
                          const std::vector<RetainedSeedRef>& fine_reverse,
                          const ChainSeedLookupCache& lookup_cache,
                          std::uint32_t occurrence_cap,
                          const ResidueRun& run) {
  std::vector<std::uint64_t> keys;
  for (int strand = 0; strand < 2; ++strand) {
    const bool reverse = strand == 1;
    const std::vector<RetainedSeedRef>& seeds =
        reverse ? fine_reverse : fine_forward;
    for (const RetainedSeedRef& retained : seeds) {
      const QuerySeed& seed = retained.seed;
      const int forward_position = forward_query_position(
          seed.read_pos, reverse, family.read_length, family.seed_length);
      if (forward_position < run.query_begin_bp ||
          forward_position >= run.query_end_bp)
        continue;
      KmerPostingView view;
      if (!lookup_cache.find_cached_view(seed.key, view)) continue;
      if (!view.found() || view.occurrence > occurrence_cap) continue;
      if (std::find(keys.begin(), keys.end(), seed.key) != keys.end())
        continue;
      keys.push_back(seed.key);
      if (static_cast<int>(keys.size()) >= kDnaResidueMinClusterAnchors)
        return kDnaResidueMinClusterAnchors;
    }
  }
  return static_cast<int>(keys.size());
}

// Builds anchors as internal::append_interval_anchors does, except that the
// contig comes from the posting, there is no diagonal band (there is no prior
// diagonal) and no tandem flag (the view spans the whole reference, not one
// contig's window). `reverse_lane` is the strand of `seed`'s stream.
void append_residue_anchors(const DnaContext& context,
                            const KmerPostingView& view, const QuerySeed& seed,
                            int seed_length, bool reverse_lane,
                            int query_length,
                            std::vector<DnaResidueAnchor>& anchors) {
  const int query_position = seed.read_pos;
  if (query_position < 0 || query_position + seed_length > query_length)
    return;
  const int chromosome_count = std::min(
      context.ref.contig_count(),
      static_cast<int>(context.ref.index->chrom_count()));
  for (std::uint32_t posting = 0; posting < view.count; ++posting) {
    const PackedRefPos packed = view.positions.packed_at(posting);
    const int contig = static_cast<int>(packed_ref_contig(packed));
    if (contig < 0 || contig >= chromosome_count) continue;
    const std::uint64_t local = packed_ref_local(packed);
    if (local > static_cast<std::uint64_t>(
                    std::numeric_limits<std::int32_t>::max()))
      continue;
    const int reference_position = static_cast<int>(local);
    if (reference_position + seed_length >
        static_cast<int>(context.ref.contig_length(contig)))
      continue;
    const bool verified = packed_ref_orientation_compatible(
        seed.z, packed, reverse_lane);
    if (!verified) continue;
    anchors.push_back(
        {contig,
         chaining::Anchor{reference_position, query_position, seed_length,
                          static_cast<std::int32_t>(anchors.size()), 0u}});
  }
}

// Expands every cached, under-cap fine seed inside the run. Returns false when
// the posting budget is exceeded; the caller then discards both vectors.
bool collect_residue_anchors(const DnaContext& context,
                             const DnaPlacementFamily& family,
                             const std::vector<RetainedSeedRef>& fine_forward,
                             const std::vector<RetainedSeedRef>& fine_reverse,
                             const ChainSeedLookupCache& lookup_cache,
                             std::uint32_t occurrence_cap,
                             const ResidueRun& run,
                             std::vector<DnaResidueAnchor>& forward_anchors,
                             std::vector<DnaResidueAnchor>& reverse_anchors) {
  std::int64_t postings = 0;
  for (int strand = 0; strand < 2; ++strand) {
    const bool reverse = strand == 1;
    const std::vector<RetainedSeedRef>& seeds =
        reverse ? fine_reverse : fine_forward;
    std::vector<DnaResidueAnchor>& out =
        reverse ? reverse_anchors : forward_anchors;
    for (const RetainedSeedRef& retained : seeds) {
      const QuerySeed& seed = retained.seed;
      const int forward_position = forward_query_position(
          seed.read_pos, reverse, family.read_length, family.seed_length);
      if (forward_position < run.query_begin_bp ||
          forward_position >= run.query_end_bp)
        continue;
      KmerPostingView view;
      if (!lookup_cache.find_cached_view(seed.key, view)) continue;
      if (!view.found() || view.occurrence > occurrence_cap) continue;
      if (!dna_residue_posting_budget_allows(postings, view.count))
        return false;
      postings += view.count;
      append_residue_anchors(context, view, seed, family.seed_length, reverse,
                             family.read_length, out);
    }
  }
  return true;
}

// Forward-query bases a selected block owns, with the realizer's boundaries.
int residue_block_query_bp(const DnaPlacementFamily& family,
                           const ::fa::cpu::voting::QueryBlock& block) {
  const int begin = ::fa::cpu::voting::query_tile_begin(
      block.query_tile_begin, family.read_length, family.seed_length,
      family.tile_count);
  const int end = block.query_tile_end == family.tile_count
                      ? family.read_length
                      : ::fa::cpu::voting::query_tile_begin(
                            block.query_tile_end, family.read_length,
                            family.seed_length, family.tile_count);
  return std::max(0, end - begin);
}

int residue_owned_query_bp(const DnaPlacementFamily& family,
                           ::fa::cpu::voting::CandidateId candidate) {
  int total = 0;
  for (const auto& block : family.partition.selected.blocks) {
    if (block.candidate != candidate) continue;
    total += residue_block_query_bp(family, block);
  }
  return total;
}

// A partition-level stand-in for the record family's primary, which is the
// widest record: the candidate owning the most query bases.
struct ResiduePrimaryWitness {
  ::fa::cpu::voting::CandidateId dominant = ::fa::cpu::voting::kNullCandidate;
  int dominant_query_bp = 0;
};

ResiduePrimaryWitness residue_primary_witness(
    const DnaPlacementFamily& family) {
  ResiduePrimaryWitness witness;
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    const int owned = residue_owned_query_bp(family, candidate.id);
    // Ties keep the lower catalogue position.
    if (owned > witness.dominant_query_bp) {
      witness.dominant_query_bp = owned;
      witness.dominant = candidate.id;
    }
  }
  return witness;
}

using ResidueAdmission = DnaResidueAdmission;

// The best admissible cluster of one interval; at most one per interval.
bool best_residue_admission(const DnaContext& context,
                            const DnaPlacementFamily& family,
                            const std::vector<DnaResidueAnchor>& forward,
                            const std::vector<DnaResidueAnchor>& reverse,
                            ResidueAdmission& best,
                            const DnaResidueAdmissionBar& bar,
                            int max_clusters_per_strand) {
  bool found = false;
  for (int strand = 0; strand < 2; ++strand) {
    const bool is_reverse = strand == 1;
    // Ties between clusters keep the first found: forward strand first, then
    // the cluster order.
    const std::vector<DnaResidueCluster> clusters =
        dna_residue_diagonal_clusters(is_reverse ? reverse : forward,
                                      max_clusters_per_strand,
                                      kDnaResidueMinClusterAnchors,
                                      context.vote_tie_seed, is_reverse);
    for (const DnaResidueCluster& cluster : clusters) {
      const DnaResidueChainOutcome outcome = dna_residue_chain_cluster(
          context, cluster.anchors, family.seed_length, family.read_length,
          bar);
      if (!outcome.admitted) continue;
      const bool better =
          !found ||
          std::tie(outcome.chain_score, outcome.chain_anchors) >
              std::tie(best.outcome.chain_score, best.outcome.chain_anchors);
      if (!better) continue;
      best.cluster = cluster;
      best.reverse = is_reverse;
      best.outcome = outcome;
      found = true;
    }
  }
  return found;
}

// Chains the top kDnaMapqRivalChains catalogue rivals of the committed family
// over the whole query, into result.rival_exact; only the MAPQ reads them.
// Rivals own no selected block; they are ranked by vote, then catalogue rank.
// An existing whole-query chain (the retained alternative, or a candidate
// restored during stabilization and later deselected) is reused.
void chain_mapq_rivals(
    const DnaContext& context, const DnaPlacementFamily& family,
    DnaPlacementChainingResult& result, RetainedSeedDensity& seed_index,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query) {
  if (result.candidates.size() != family.candidates.size()) return;
  // A rival with little vote is not worth a whole-query chain: it still
  // enters the MAPQ through its vote factor. The strongest owner's vote
  // stands in for the primary's, which is decided later.
  int owner_vote = 0;
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    for (const auto& block : family.partition.selected.blocks) {
      if (block.candidate != candidate.id) continue;
      owner_vote = std::max(owner_vote, candidate.vote_evidence);
      break;
    }
  }
  std::vector<std::size_t> ranked;
  ranked.reserve(family.candidates.size());
  for (std::size_t index = 0; index < family.candidates.size(); ++index) {
    const DnaPlacementCandidate& candidate = family.candidates[index];
    if (static_cast<std::int64_t>(candidate.vote_evidence) *
            kDnaMapqRivalVoteDenominator <
        static_cast<std::int64_t>(owner_vote))
      continue;
    bool owns_selected_block = false;
    for (const auto& block : family.partition.selected.blocks) {
      if (block.candidate == candidate.id) {
        owns_selected_block = true;
        break;
      }
    }
    if (owns_selected_block) continue;
    ranked.push_back(index);
  }
  const std::size_t keep = std::min<std::size_t>(
      ranked.size(), static_cast<std::size_t>(kDnaMapqRivalChains));
  std::partial_sort(
      ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(keep),
      ranked.end(),
      [&family](std::size_t left, std::size_t right) {
        const DnaPlacementCandidate& a = family.candidates[left];
        const DnaPlacementCandidate& b = family.candidates[right];
        return std::tie(b.vote_evidence, a.catalogue_rank, a.id) <
               std::tie(a.vote_evidence, b.catalogue_rank, b.id);
      });
  ranked.resize(keep);
  // kDnaChainMapqStudyBlockRivals: on a split family each selected block
  // also gets its own top rivals, with the vote floor of that block's owner
  // and mask_level against the block's query span (using the candidate's
  // screening span). The union is chained below.
  if ((::fa::cpu::lr::dna_chain_mapq_study_bits(
           context.opts.chain_mapq_hifi_margin) &
       kDnaChainMapqStudyBlockRivals) != 0) {
    int selected_blocks = 0;
    for (const auto& block : family.partition.selected.blocks)
      if (block.candidate != ::fa::cpu::voting::kNullCandidate)
        ++selected_blocks;
    if (selected_blocks >= 2) {
      for (const auto& block : family.partition.selected.blocks) {
        if (block.candidate == ::fa::cpu::voting::kNullCandidate) continue;
        const DnaPlacementCandidate* owner = family.find(block.candidate);
        if (owner == nullptr) continue;
        const int block_begin = ::fa::cpu::voting::query_tile_begin(
            block.query_tile_begin, family.read_length, family.seed_length,
            family.tile_count);
        const int block_end = ::fa::cpu::voting::query_tile_end(
            block.query_tile_end, family.read_length, family.seed_length,
            family.tile_count);
        const int block_span = std::max(0, block_end - block_begin);
        std::vector<std::size_t> block_ranked;
        for (std::size_t index = 0; index < family.candidates.size();
             ++index) {
          const DnaPlacementCandidate& candidate = family.candidates[index];
          if (static_cast<std::int64_t>(candidate.vote_evidence) *
                  kDnaMapqRivalVoteDenominator <
              static_cast<std::int64_t>(owner->vote_evidence))
            continue;
          bool owns_selected_block = false;
          for (const auto& other : family.partition.selected.blocks) {
            if (other.candidate == candidate.id) {
              owns_selected_block = true;
              break;
            }
          }
          if (owns_selected_block) continue;
          const DnaPlacementCandidateChain& record = result.candidates[index];
          const int span = record.screening_forward_query_end -
                           record.screening_forward_query_begin;
          if (record.screening_forward_query_begin < 0 || span <= 0) continue;
          const int overlap =
              std::min(record.screening_forward_query_end, block_end) -
              std::max(record.screening_forward_query_begin, block_begin);
          if (overlap <= 0 || static_cast<std::int64_t>(overlap) * 2 <=
                                  std::min(span, block_span))
            continue;
          block_ranked.push_back(index);
        }
        const std::size_t block_keep = std::min<std::size_t>(
            block_ranked.size(),
            static_cast<std::size_t>(kDnaMapqRivalChains));
        std::partial_sort(
            block_ranked.begin(),
            block_ranked.begin() + static_cast<std::ptrdiff_t>(block_keep),
            block_ranked.end(),
            [&family](std::size_t left, std::size_t right) {
              const DnaPlacementCandidate& a = family.candidates[left];
              const DnaPlacementCandidate& b = family.candidates[right];
              return std::tie(b.vote_evidence, a.catalogue_rank, a.id) <
                     std::tie(a.vote_evidence, b.catalogue_rank, b.id);
            });
        block_ranked.resize(block_keep);
        for (const std::size_t index : block_ranked) {
          if (std::find(ranked.begin(), ranked.end(), index) == ranked.end())
            ranked.push_back(index);
        }
      }
    }
  }
  result.rival_exact.reserve(ranked.size());
  for (const std::size_t index : ranked) {
    const DnaPlacementCandidate& candidate = family.candidates[index];
    DnaRivalExactChain rival;
    rival.candidate = candidate.id;
    const bool is_retained_alternative =
        result.alternative.refusal == DnaAlternativeRefusal::None &&
        result.alternative.candidate == candidate.id;
    if (is_retained_alternative && result.alternative_exact.exact) {
      rival.chain = result.alternative_exact;
      rival.reused = true;
    } else if (result.candidates[index].exact) {
      rival.chain = result.candidates[index];
      rival.reused = true;
    } else {
      rival.chain.candidate = candidate.id;
      chain_candidate(context, family, candidate,
                      candidate.peak.is_rc ? reverse_query : forward_query,
                      seed_index, rival.chain,
                      CandidateChainPass::WholeQueryExact);
      ++result.mapq_rival_chains;
    }
    // Kept whatever the outcome; `status` says whether it chained.
    result.rival_exact.push_back(std::move(rival));
  }
}

// Restores an alternative's exact whole-query chain into `exact`. A candidate
// already holding an accepted whole-query chain would chain to the same
// result again, so that chain is copied. False unless the restore is
// accepted with a primary path.
bool restore_alternative_exact(const DnaContext& context,
                               const DnaPlacementFamily& family,
                               const DnaPlacementChainingResult& result,
                               RetainedSeedDensity& seed_index,
                               const std::vector<std::uint8_t>& forward_query,
                               const std::vector<std::uint8_t>& reverse_query,
                               ::fa::cpu::voting::CandidateId id,
                               DnaPlacementCandidateChain& exact) {
  const DnaPlacementCandidate* candidate = family.find(id);
  const DnaPlacementCandidateChain* attempted = result.find(id);
  exact.candidate = id;
  const bool chained_before =
      candidate != nullptr && attempted != nullptr && attempted->exact &&
      attempted->status == DnaPlacementChainStatus::Accepted;
  if (chained_before)
    copy_chain_call(*attempted, true, exact);
  const bool restored = chained_before ||
      (candidate != nullptr &&
       chain_candidate(
           context, family, *candidate,
           candidate->peak.is_rc ? reverse_query : forward_query, seed_index,
           exact, CandidateChainPass::WholeQueryExact));
  return restored && exact.status == DnaPlacementChainStatus::Accepted &&
         !exact.primary.empty();
}

}  // namespace

int dna_residue_cached_supply(
    const DnaPlacementFamily& family,
    const std::vector<RetainedSeedRef>& fine_forward,
    const std::vector<RetainedSeedRef>& fine_reverse,
    const ChainSeedLookupCache& lookup_cache, std::uint32_t occurrence_cap,
    int query_begin_bp, int query_end_bp) {
  ResidueRun run;
  run.query_begin_bp = query_begin_bp;
  run.query_end_bp = query_end_bp;
  return residue_cached_supply(family, fine_forward, fine_reverse,
                               lookup_cache, occurrence_cap, run);
}

bool dna_residue_collect_anchors(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<RetainedSeedRef>& fine_forward,
    const std::vector<RetainedSeedRef>& fine_reverse,
    const ChainSeedLookupCache& lookup_cache, std::uint32_t occurrence_cap,
    int query_begin_bp, int query_end_bp,
    std::vector<DnaResidueAnchor>& forward_anchors,
    std::vector<DnaResidueAnchor>& reverse_anchors) {
  ResidueRun run;
  run.query_begin_bp = query_begin_bp;
  run.query_end_bp = query_end_bp;
  return collect_residue_anchors(context, family, fine_forward, fine_reverse,
                                 lookup_cache, occurrence_cap, run,
                                 forward_anchors, reverse_anchors);
}

void dna_residue_collect_anchors_rarest_first(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<RetainedSeedRef>& fine_forward,
    const std::vector<RetainedSeedRef>& fine_reverse,
    const ChainSeedLookupCache& lookup_cache, std::uint32_t occurrence_cap,
    int query_begin_bp, int query_end_bp,
    std::vector<DnaResidueAnchor>& forward_anchors,
    std::vector<DnaResidueAnchor>& reverse_anchors) {
  struct Pending {
    std::uint32_t occurrence = 0;
    std::uint64_t key = 0;
    int read_pos = 0;
    bool reverse = false;
    QuerySeed seed;
    KmerPostingView view;
  };
  std::vector<Pending> pending;
  for (int strand = 0; strand < 2; ++strand) {
    const bool reverse = strand == 1;
    const std::vector<RetainedSeedRef>& seeds =
        reverse ? fine_reverse : fine_forward;
    for (const RetainedSeedRef& retained : seeds) {
      const QuerySeed& seed = retained.seed;
      const int forward_position = forward_query_position(
          seed.read_pos, reverse, family.read_length, family.seed_length);
      if (forward_position < query_begin_bp ||
          forward_position >= query_end_bp)
        continue;
      KmerPostingView view;
      if (!lookup_cache.find_cached_view(seed.key, view)) continue;
      if (!view.found() || view.occurrence > occurrence_cap) continue;
      pending.push_back({view.occurrence, seed.key, seed.read_pos, reverse,
                         seed, view});
    }
  }
  std::sort(pending.begin(), pending.end(),
            [](const Pending& left, const Pending& right) {
              return std::tie(left.occurrence, left.key, left.read_pos,
                              left.reverse) <
                     std::tie(right.occurrence, right.key, right.read_pos,
                              right.reverse);
            });
  std::int64_t postings = 0;
  for (const Pending& entry : pending) {
    if (!dna_residue_posting_budget_allows(postings, entry.view.count)) break;
    postings += entry.view.count;
    append_residue_anchors(
        context, entry.view, entry.seed, family.seed_length, entry.reverse,
        family.read_length,
        entry.reverse ? reverse_anchors : forward_anchors);
  }
}

bool dna_residue_best_admission(const DnaContext& context,
                                const DnaPlacementFamily& family,
                                const std::vector<DnaResidueAnchor>& forward,
                                const std::vector<DnaResidueAnchor>& reverse,
                                DnaResidueAdmission& best,
                                const DnaResidueAdmissionBar& bar,
                                int max_clusters_per_strand) {
  return best_residue_admission(context, family, forward, reverse, best, bar,
                                max_clusters_per_strand);
}

DnaResidueObservedDensity dna_residue_observed_density(
    const DnaPlacementChainingResult& placement) {
  DnaResidueObservedDensity observed;
  if (!placement.accepted || !placement.family.valid) return observed;
  const ResiduePrimaryWitness witness =
      residue_primary_witness(placement.family);
  if (witness.dominant == ::fa::cpu::voting::kNullCandidate) return observed;
  const DnaPlacementCandidateChain* winner = placement.find(witness.dominant);
  if (winner == nullptr || !winner->exact) return observed;
  const int span = winner->forward_query_end - winner->forward_query_begin;
  if (winner->chain_anchors <= 0 || span <= 0) return observed;
  observed.anchors = winner->chain_anchors;
  observed.query_span = span;
  return observed;
}

DnaResidueDetachedChain dna_residue_detached_chain(
    const DnaPlacementFamily& family, const DnaResidueCluster& cluster,
    bool reverse, const DnaResidueChainOutcome& outcome,
    ::fa::cpu::voting::CandidateId id) {
  DnaResidueDetachedChain pair;
  pair.candidate.id = id;
  pair.candidate.lane = reverse ? 1 : 0;
  pair.candidate.vote_evidence = outcome.chain_anchors;
  pair.candidate.screening_chain_score = outcome.chain_score;

  ::fa::cpu::voting::QueryTileMask support;
  int forward_begin = std::numeric_limits<int>::max();
  int forward_end = std::numeric_limits<int>::min();
  for (const chaining::Anchor& anchor : outcome.primary) {
    support.set(dna_forward_query_tile(anchor.q, reverse, family.read_length,
                                       family.seed_length, family.tile_count));
    const int position = forward_query_position(
        anchor.q, reverse, family.read_length, family.seed_length);
    forward_begin = std::min(forward_begin, position);
    forward_end = std::max(forward_end, position + family.seed_length);
  }
  pair.candidate.support = support;
  pair.candidate.peak.chr = cluster.contig;
  pair.candidate.peak.is_rc = reverse;
  pair.candidate.peak.raw_ref_start = cluster.peak_diagonal;
  pair.candidate.peak.ref_pos = static_cast<int>(
      std::max<std::int64_t>(0, cluster.peak_diagonal));
  pair.candidate.peak.read_lo = forward_begin;
  pair.candidate.peak.read_hi = forward_end;
  pair.candidate.peak.support = outcome.chain_anchors;
  pair.candidate.peak.center_support = outcome.chain_anchors;
  pair.candidate.peak.vote_score = outcome.chain_anchors;
  pair.candidate.peak.anchor.ref_start_bin =
      static_cast<int>(cluster.diagonal_bin);
  pair.candidate.peak.anchor.ref_start_bin_width = kDnaResidueDiagonalWidth;

  pair.chain.candidate = id;
  pair.chain.sparse_support = support;
  pair.chain.dense_support = support;
  pair.chain.primary = outcome.primary;
  pair.chain.status = DnaPlacementChainStatus::Accepted;
  pair.chain.chain_score = outcome.chain_score;
  pair.chain.chain_anchors = outcome.chain_anchors;
  pair.chain.screening_chain_score = outcome.chain_score;
  pair.chain.screening_chain_anchors = outcome.chain_anchors;
  pair.chain.exact = true;
  set_spans(pair.chain, reverse, family.read_length);
  pair.chain.screening_forward_query_begin = pair.chain.forward_query_begin;
  pair.chain.screening_forward_query_end = pair.chain.forward_query_end;
  return pair;
}

bool dna_residue_posting_budget_allows(std::int64_t used,
                                       std::uint32_t next) noexcept {
  return used >= 0 && used <= kDnaResidueMaxIntervalPostings &&
         static_cast<std::int64_t>(next) <=
             kDnaResidueMaxIntervalPostings - used;
}

std::vector<DnaResidueCluster> dna_residue_diagonal_clusters(
    const std::vector<DnaResidueAnchor>& anchors, int max_clusters,
    int min_anchors, std::uint32_t tie_seed, bool reverse) {
  std::vector<DnaResidueCluster> clusters;
  if (anchors.empty() || max_clusters <= 0) return clusters;

  struct Keyed {
    int contig;
    std::int64_t bin;
    std::int64_t diagonal;
    std::size_t index;
  };
  std::vector<Keyed> keyed;
  keyed.reserve(anchors.size());
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    const chaining::Anchor& anchor = anchors[index].anchor;
    const std::int64_t diagonal =
        static_cast<std::int64_t>(anchor.r) - static_cast<std::int64_t>(anchor.q);
    keyed.push_back({anchors[index].contig,
                     residue_floor_div(diagonal, kDnaResidueDiagonalWidth),
                     diagonal, index});
  }
  std::sort(keyed.begin(), keyed.end(),
            [](const Keyed& left, const Keyed& right) {
              return std::tie(left.contig, left.bin, left.diagonal,
                              left.index) <
                     std::tie(right.contig, right.bin, right.diagonal,
                              right.index);
            });

  struct BinRun {
    int contig;
    std::int64_t bin;
    std::size_t begin;
    std::size_t end;
  };
  std::vector<BinRun> bins;
  for (std::size_t index = 0; index < keyed.size();) {
    std::size_t end = index + 1;
    while (end < keyed.size() && keyed[end].contig == keyed[index].contig &&
           keyed[end].bin == keyed[index].bin)
      ++end;
    bins.push_back({keyed[index].contig, keyed[index].bin, index, end});
    index = end;
  }

  std::vector<std::size_t> ranked(bins.size());
  for (std::size_t index = 0; index < bins.size(); ++index) ranked[index] = index;
  std::sort(ranked.begin(), ranked.end(),
            [&bins, tie_seed, reverse](std::size_t left, std::size_t right) {
              const std::size_t left_count = bins[left].end - bins[left].begin;
              const std::size_t right_count =
                  bins[right].end - bins[right].begin;
              if (left_count != right_count) return left_count > right_count;
              // Ties go to the read-seeded hash of the bin's locus, then to
              // (contig, bin).
              const std::uint64_t left_hash = tie_locus_hash(
                  tie_seed, bins[left].contig, reverse, bins[left].bin);
              const std::uint64_t right_hash = tie_locus_hash(
                  tie_seed, bins[right].contig, reverse, bins[right].bin);
              if (left_hash != right_hash) return left_hash < right_hash;
              return std::tie(bins[left].contig, bins[left].bin) <
                     std::tie(bins[right].contig, bins[right].bin);
            });

  std::vector<std::size_t> taken;
  for (const std::size_t index : ranked) {
    if (static_cast<int>(clusters.size()) >= max_clusters) break;
    // Adjacent bins would share anchors, since a cluster absorbs its
    // neighbours.
    const bool adjacent = std::any_of(
        taken.begin(), taken.end(), [&](std::size_t other) {
          return bins[other].contig == bins[index].contig &&
                 std::llabs(bins[other].bin - bins[index].bin) <= 1;
        });
    if (adjacent) continue;
    DnaResidueCluster cluster;
    cluster.contig = bins[index].contig;
    cluster.diagonal_bin = bins[index].bin;
    std::vector<std::int64_t> diagonals;
    // Collect the bin and both neighbours.
    for (const BinRun& run : bins) {
      if (run.contig != bins[index].contig) continue;
      if (std::llabs(run.bin - bins[index].bin) > 1) continue;
      for (std::size_t at = run.begin; at < run.end; ++at) {
        cluster.anchors.push_back(anchors[keyed[at].index].anchor);
        diagonals.push_back(keyed[at].diagonal);
      }
    }
    if (static_cast<int>(cluster.anchors.size()) < min_anchors) continue;
    std::sort(diagonals.begin(), diagonals.end());
    cluster.peak_diagonal = diagonals[diagonals.size() / 2];
    taken.push_back(index);
    clusters.push_back(std::move(cluster));
  }
  return clusters;
}

DnaResidueChainOutcome dna_residue_chain_cluster(
    const DnaContext& context, std::vector<chaining::Anchor> anchors,
    int seed_length, int read_length, const DnaResidueAdmissionBar& bar) {
  DnaResidueChainOutcome outcome;
  if (anchors.empty()) return outcome;
  // The same ordering and deduplication as chain_candidate.
  internal::sort_unique_pool_anchors(anchors);
  const chaining::ColinearChainParams chain_params =
      dna_candidate_chain_params(context, context.opts.cigar_dp_bw_long,
                                 seed_length, read_length);
  // The dense chain, as in the whole-query pass. Clusters are small, so this
  // usually takes the linear scan.
  chaining::DenseChainStats dense_stats;
  const chaining::ChainResult chained = chaining::chain_dense_colinear(
      std::move(anchors),
      dna_dense_chain_params(chain_params, seed_length,
                             context.opts.dna_dense_diag_min_runs),
      &dense_stats);
  const chaining::ChainPartition partition =
      chaining::partition_chains(chained, ::fa::cpu::split::kHalfFloor);
  if (partition.primary < 0) return outcome;
  const chaining::Chain& primary =
      chained.chains[static_cast<std::size_t>(partition.primary)];
  outcome.chain_score = primary.score;
  outcome.chain_anchors = static_cast<int>(primary.idx.size());
  outcome.primary.reserve(primary.idx.size());
  int query_begin = std::numeric_limits<int>::max();
  int query_end = std::numeric_limits<int>::min();
  for (const std::int32_t index : primary.idx) {
    const chaining::Anchor& anchor =
        chained.anchors[static_cast<std::size_t>(index)];
    query_begin = std::min(query_begin, static_cast<int>(anchor.q));
    query_end = std::max(query_end, static_cast<int>(anchor.q_end()));
    outcome.primary.push_back(anchor);
  }
  // Oriented and forward query spans have the same length.
  outcome.chain_query_span = query_end - query_begin;
  // An absolute bar, compared with no other candidate. The density test, in
  // integers, admits a chain exactly at the bar and rejects spurious chains
  // in reads that do not contain the event.
  const int density_bar = bar.min_anchor_density_per_100bp < 0
                              ? context.opts
                                    .residue_min_anchor_density_per_100bp
                              : bar.min_anchor_density_per_100bp;
  const std::int64_t density_scaled =
      static_cast<std::int64_t>(outcome.chain_anchors) * 100;
  const std::int64_t density_required =
      static_cast<std::int64_t>(std::max(0, density_bar)) *
      static_cast<std::int64_t>(std::max(0, outcome.chain_query_span));
  outcome.admitted =
      outcome.chain_anchors >= bar.min_chain_anchors &&
      outcome.chain_score >=
          kDnaResidueScoreFloorMatches *
              std::min(std::max(seed_length, 0), 255) &&
      density_scaled >= density_required;
  return outcome;
}

chaining::ColinearChainParams dna_candidate_chain_params(
    const DnaContext& context, int band, int seed_length, int read_length) {
  chaining::ColinearChainParams params;
  params.bw = std::max(1, band);
  const int preset_reach = std::max(
      params.bw, context.opts.cigar_local_interval_anchor_chain_max_gap);
  const int read_reach =
      std::min(read_length / kCandidateReachDivisor, kCandidateReachCap);
  params.max_dist_x = params.max_dist_y = std::max(preset_reach, read_reach);
  params.bw = std::max(params.bw, read_reach);
  params.min_cnt = 1;
  params.min_sc = 1;
  params.chn_pen_gap = 0.008f * static_cast<float>(seed_length);
  params.chn_pen_skip = 0.0f;
  return params;
}

chaining::DenseChainParams
dna_dense_chain_params(const chaining::ColinearChainParams& params,
                       int seed_length, int diag_min_runs) {
  chaining::DenseChainParams dense;
  dense.span = seed_length;
  dense.max_dist_x = params.max_dist_x;
  dense.max_dist_y = params.max_dist_y;
  dense.bw = params.bw;
  dense.min_sc = params.min_sc;
  dense.min_cnt = params.min_cnt;
  dense.chn_pen_gap = params.chn_pen_gap;
  dense.chn_pen_skip = params.chn_pen_skip;
  // Pools at or above this run count use the diagonal-keyed search.
  dense.diag_min_runs = diag_min_runs;
  return dense;
}

const DnaPlacementCandidateChain* DnaPlacementChainingResult::find(
    ::fa::cpu::voting::CandidateId candidate) const noexcept {
  const auto found = std::find_if(
      candidates.begin(), candidates.end(),
      [candidate](const DnaPlacementCandidateChain& record) {
        return record.candidate == candidate;
      });
  return found == candidates.end() ? nullptr : &*found;
}

const DnaPlacementCandidateChain* DnaPlacementChainingResult::whole_query_chain(
    ::fa::cpu::voting::CandidateId candidate) const noexcept {
  for (const DnaRivalExactChain& chained : rival_exact)
    if (chained.candidate == candidate) return &chained.chain;
  if (alternative.refusal == DnaAlternativeRefusal::None &&
      alternative.candidate == candidate && alternative_exact.exact)
    return &alternative_exact;
  const DnaPlacementCandidateChain* record = find(candidate);
  return record != nullptr && record->exact ? record : nullptr;
}

DnaPlacementCandidateChain
dna_sibling_rival_chain(const DnaPlacementCandidateChain& winner, bool reverse,
                        int read_length, int seed_length, int tile_count) {
  DnaPlacementCandidateChain sibling;
  sibling.candidate = winner.candidate;
  if (winner.rival_sibling < 0 ||
      static_cast<std::size_t>(winner.rival_sibling) >=
          winner.sibling_paths.size() ||
      static_cast<std::size_t>(winner.rival_sibling) >=
          winner.sibling_scores.size())
    return sibling;
  const std::size_t which = static_cast<std::size_t>(winner.rival_sibling);
  sibling.primary = winner.sibling_paths[which];
  if (sibling.primary.empty()) return sibling;
  sibling.chain_score = winner.sibling_scores[which];
  sibling.chain_anchors = static_cast<int>(sibling.primary.size());
  for (const chaining::Anchor& anchor : sibling.primary)
    sibling.dense_support.set(dna_forward_query_tile(
        anchor.q, reverse, read_length, seed_length, tile_count));
  sibling.exact = true;
  sibling.status = DnaPlacementChainStatus::Accepted;
  set_spans(sibling, reverse, read_length);
  return sibling;
}

const DnaPlacementCandidateChain* dna_committed_winner_chain(
    const DnaPlacementChainingResult* stable,
    const DnaPlacementChainingResult* promoted_alternative,
    ::fa::cpu::voting::CandidateId primary_candidate,
    bool primary_is_alternative) noexcept {
  if (primary_is_alternative && promoted_alternative != nullptr &&
      promoted_alternative->family.original_candidate_id == primary_candidate)
    // The restricted family's only candidate has solver id 0.
    return promoted_alternative->find(0);
  return stable == nullptr ? nullptr : stable->find(primary_candidate);
}

DnaPlacementChainingResult build_dna_placement_chains(
    const DnaContext& context, DnaPlacementFamily family,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    const std::vector<ChainWindowRetainedSeed>* forward_seeds,
    const std::vector<ChainWindowRetainedSeed>* reverse_seeds,
    const std::vector<QuerySeed>* fine_forward_seeds,
    const std::vector<QuerySeed>* fine_reverse_seeds,
    ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots,
    DnaTileOwnership ownership) {
  DnaPlacementChainingResult result;
  // A deferred first partition is solved only where the family is returned
  // before the re-solve below. The initial_* fields then read an empty
  // partition, and nothing reads them.
  const auto settle_partition = [&context, &family] {
    if (!family.partition_deferred) return;
    family.partition = repartition(context, family);
    family.partition_deferred = false;
  };
  result.initial_score = family.partition.selected.score;
  result.initial_blocks = family.partition.selected.non_null_blocks;
  const std::vector<::fa::cpu::voting::CandidateId> initial_assignment =
      family.partition.selected.assignment;
  result.candidates.reserve(family.candidates.size());
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    DnaPlacementCandidateChain record;
    record.candidate = candidate.id;
    record.sparse_support = candidate.support;
    record.initial_selected_tiles =
        selected_tiles(family.partition.selected, candidate.id);
    result.candidates.push_back(std::move(record));
  }
  if (!family.valid || context.ref.index == nullptr ||
      context.ref.encoded == nullptr || fine_forward_seeds == nullptr ||
      fine_reverse_seeds == nullptr || lookup_cache == nullptr) {
    settle_partition();
    result.family = std::move(family);
    return result;
  }

  // Realization reads this density when it gates the late inversion probe
  // (inv_local_chain.h), so it then outlives placement.
  std::shared_ptr<RetainedSeedDensity> kept_seed_index;
  if (context.opts.inversion_probe_local_gate &&
      context.opts.enable_full_read_cigar)
    kept_seed_index = std::make_shared<RetainedSeedDensity>();
  RetainedSeedDensity local_seed_index;
  RetainedSeedDensity& seed_index =
      kept_seed_index ? *kept_seed_index : local_seed_index;
  if (!seed_index.build(
          *context.ref.index, forward_seeds, reverse_seeds, fine_forward_seeds,
          fine_reverse_seeds, lookup_cache, fine_forward_slots,
          fine_reverse_slots)) {
    settle_partition();
    result.family = std::move(family);
    return result;
  }
  result.inversion_gate_seeds = std::move(kept_seed_index);

  // The screening pass over the catalogue: one cheap chain per candidate from
  // the representative seeds, with no deferral.
  for (std::size_t index = 0; index < family.candidates.size(); ++index) {
    DnaPlacementCandidate& candidate = family.candidates[index];
    DnaPlacementCandidateChain& record = result.candidates[index];
    const bool sparse = chain_candidate(
        context, family, candidate,
        candidate.peak.is_rc ? reverse_query : forward_query,
        seed_index, record, CandidateChainPass::BoundedScreening);
    // Keep the screening result before the exact restore overwrites it.
    candidate.screening_chain_score = record.chain_score;
    record.screening_chain_score = record.chain_score;
    record.screening_chain_anchors = record.chain_anchors;
    record.screening_forward_query_begin = record.forward_query_begin;
    record.screening_forward_query_end = record.forward_query_end;
    if (sparse) {
      // Screening may add support but never removes vote support.
      candidate.support = candidate.support | record.dense_support;
    }
  }
  family.partition = repartition(context, family);
  family.partition_deferred = false;

  if (!stabilize_selected_family(context, family, result, seed_index,
                                 forward_query, reverse_query, ownership)) {
    result.family = std::move(family);
    return result;
  }

  // Keep the fine seeds for terminal-clip recovery, which runs after this
  // index is gone; the posting views stay in the per-read cache.
  result.residue_fine_forward = seed_index.fine_forward();
  result.residue_fine_reverse = seed_index.fine_reverse();
  for (DnaPlacementCandidateChain& record : result.candidates) {
    record.final_selected_tiles =
        selected_tiles(family.partition.selected, record.candidate);
  }
  result.alternative =
      select_dna_alternative_hypothesis(context, family, result);
  if (result.alternative.refusal == DnaAlternativeRefusal::None &&
      !restore_alternative_exact(context, family, result, seed_index,
                                 forward_query, reverse_query,
                                 result.alternative.candidate,
                                 result.alternative_exact)) {
    result.alternative.refusal = DnaAlternativeRefusal::ExactRestoreFailed;
    result.alternative.candidate = ::fa::cpu::voting::kNullCandidate;
    result.alternative_exact = {};
  }
  // After the alternative restore: the family is final, seed_index is alive
  // and the alternative's chain can be reused.
  chain_mapq_rivals(context, family, result, seed_index, forward_query,
                    reverse_query);
  // Ranks 2..n, after everything that decides the read, so that they
  // change none of it.
  if (context.opts.alternative_realize_max >= 2) {
    const std::vector<::fa::cpu::voting::CandidateId> ranked =
        rank_dna_alternative_hypotheses(
            context, family, result,
            static_cast<std::size_t>(context.opts.alternative_realize_max));
    for (std::size_t rank = 1; rank < ranked.size(); ++rank) {
      DnaRankedAlternative alternative;
      alternative.candidate = ranked[rank];
      if (restore_alternative_exact(context, family, result, seed_index,
                                    forward_query, reverse_query,
                                    alternative.candidate, alternative.exact))
        result.ranked_alternatives.push_back(std::move(alternative));
    }
  }
  result.selection_changed =
      initial_assignment != family.partition.selected.assignment;
  result.family = std::move(family);
  return result;
}

DnaPlacementChainingResult build_dna_alternative_placement(
    const DnaContext& context, const DnaPlacementFamily& stable_family,
    const DnaPlacementChainingResult& stable,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query) {
  if (stable.alternative.refusal != DnaAlternativeRefusal::None)
    return DnaPlacementChainingResult();
  return build_dna_rival_placement(context, stable_family,
                                   stable.alternative.candidate,
                                   stable.alternative_exact, forward_query,
                                   reverse_query);
}

DnaPlacementChainingResult build_dna_rival_placement(
    const DnaContext& context, const DnaPlacementFamily& stable_family,
    ::fa::cpu::voting::CandidateId original,
    const DnaPlacementCandidateChain& exact_chain,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query) {
  DnaPlacementChainingResult result;
  const DnaPlacementCandidate* source = stable_family.find(original);
  if (original == ::fa::cpu::voting::kNullCandidate || source == nullptr ||
      !stable_family.valid || exact_chain.primary.empty() ||
      context.ref.index == nullptr || context.ref.encoded == nullptr)
    return result;
  // The chain is given, so the seed index stays empty: nothing is rescanned.
  RetainedSeedDensity seed_index;

  DnaPlacementFamily family;
  family.read_length = stable_family.read_length;
  family.seed_length = stable_family.seed_length;
  family.tile_count = stable_family.tile_count;
  family.original_candidate_id = original;
  DnaPlacementCandidate only = *source;
  only.id = 0;
  only.catalogue_rank = 0;
  only.support = exact_chain.dense_support;
  if (only.peak.is_rc) family.reverse_candidates = 1;
  else family.forward_candidates = 1;
  family.candidates.push_back(std::move(only));
  family.valid = true;

  result.candidates.push_back(exact_chain);
  result.candidates.back().candidate = 0;
  family.partition = repartition(context, family);
  // The only candidate is already exact, so the ownership rule is unused.
  if (!stabilize_selected_family(context, family, result, seed_index,
                                 forward_query, reverse_query,
                                 DnaTileOwnership::Span))
    result.accepted = false;
  const int non_null_blocks = static_cast<int>(std::count_if(
      family.partition.selected.blocks.begin(),
      family.partition.selected.blocks.end(), [](const auto& block) {
        return block.candidate != ::fa::cpu::voting::kNullCandidate;
      }));
  if (non_null_blocks != 1) result.accepted = false;
  result.family = std::move(family);
  return result;
}

}  // namespace fa::cpu::lr
