#include "placement_chaining.h"
#include "placement_chaining_internal.h"
#include "chain_mapq.h" // kDnaChainMapqStudyBlockRivals
#include "chain_ownership.h"

#include "../chaining/colinear_chain.h"
#include "../chaining/dense_chain.h"
#include "../core/radix_sort.h"
#include "../chaining/partition.h"
#include "../index/format.h"
#include "../split/query_geometry.h"
#include "../voting/query_tiles.h"

#include <algorithm>
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
    const VotePeak* line_peak,
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
  // The winner's band follows its line at this seed.
  const std::int64_t line_diagonal =
      line_peak != nullptr
          ? line_peak->line_a +
                vote_slope_stretch(seed.seed.read_pos, line_peak->line_b_q20)
          : 0;
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
        (line_peak != nullptr
             ? within_band<std::int64_t>(
                   static_cast<std::int64_t>(reference_position) -
                       query_position,
                   line_diagonal, diagonal_band)
             : within_band(reference_position - query_position, main_diagonal,
                           diagonal_band));
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
  internal::HarvestWindow window = internal::harvest_window(
      context, candidate.peak, family.read_length, chromosome_length);
  // The whole-query pass of a candidate whose screening pass held a chain also
  // harvests [low - pad, high + L + pad] over that chain's lowest and highest
  // diagonal.
  if (whole_query_exact && candidate.screening_diagonals) {
    const int interval_pad =
        context.opts.cigar_local_interval_anchor_interval_pad;
    window.low = std::min(
        window.low,
        std::max<std::int64_t>(
            0, static_cast<std::int64_t>(candidate.screening_diagonal_low) -
                   interval_pad));
    window.high = std::max(
        window.high,
        std::min<std::int64_t>(
            chromosome_length,
            static_cast<std::int64_t>(candidate.screening_diagonal_high) +
                family.read_length + interval_pad));
  }
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
  const VotePeak* line_peak =
      candidate.peak.line_gate ? &candidate.peak : nullptr;
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
            static_cast<int>(expected), line_peak, seed_length, diagonal_band,
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
          static_cast<int>(expected), line_peak, seed_length, diagonal_band,
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

// Stabilization: chain every owner of the partition that is not yet exact
// over the whole query, and accept the family when every owner has an
// accepted whole-query chain. A family with an owner lacking one stays
// unaccepted, and the read is unmapped.
bool stabilize_selected_family(
    const DnaContext& context, DnaPlacementFamily& family,
    DnaPlacementChainingResult& result, RetainedSeedDensity& seed_index,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query) {
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
    chain_candidate(context, family, *candidate_it,
                    candidate_it->peak.is_rc ? reverse_query : forward_query,
                    seed_index, *record_it,
                    CandidateChainPass::WholeQueryExact);
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

// The blocks owned on the accepted whole-query chains' anchors
// (dna/chain_ownership.h), every candidate's primary and siblings, replace
// the tile owners: each block becomes one block of partition.selected, in
// query order, with its chain range and bounds in family.block_parts, and a
// tile goes to the block holding its middle. False, with the partition left
// as it is, when no chain qualifies to own. `kept` receives the selection's
// items and roles.
bool select_block_owners(const DnaContext& context, DnaPlacementFamily& family,
                         const DnaPlacementChainingResult& result,
                         DnaKeptSelection& kept) {
  namespace voting = ::fa::cpu::voting;
  std::vector<DnaOwnershipPath> paths;
  for (const DnaPlacementCandidateChain& record : result.candidates) {
    if (!record.exact || record.status != DnaPlacementChainStatus::Accepted)
      continue;
    const DnaPlacementCandidate* candidate = family.find(record.candidate);
    if (candidate == nullptr) continue;
    paths.push_back({record.candidate, -1, candidate->peak.chr,
                     candidate->peak.is_rc, &record.primary,
                     record.chain_score});
    for (std::size_t j = 0; j < record.sibling_paths.size(); ++j)
      paths.push_back({record.candidate, static_cast<int>(j),
                       candidate->peak.chr, candidate->peak.is_rc,
                       &record.sibling_paths[j],
                       j < record.sibling_scores.size()
                           ? record.sibling_scores[j]
                           : 0});
  }
  std::vector<DnaOwnershipRole> roles;
  const DnaOwnershipSelection selection =
      select_chain_owners(family.read_length, family.seed_length,
                          context.opts.min_chain_score, paths, roles);
  kept.items.clear();
  for (const DnaOwnershipPath& item : selection.items)
    kept.items.push_back(
        {item.candidate, item.path, item.contig, item.reverse, item.score});
  kept.roles = std::move(roles);
  if (selection.blocks.empty()) return false;
  const int read_length = family.read_length;
  const int seed_length = family.seed_length;
  const int tile_count = family.tile_count;
  std::vector<DnaBlockPart> parts;
  std::vector<voting::QueryBlock> blocks;
  parts.reserve(selection.blocks.size());
  blocks.reserve(selection.blocks.size());
  for (const DnaOwnershipBlock& owned : selection.blocks) {
    const DnaOwnershipPath& item =
        selection.items[static_cast<std::size_t>(owned.item)];
    const std::vector<chaining::Anchor>& anchors = *item.anchors;
    DnaBlockPart part;
    part.path = item.path;
    part.anchor_begin = owned.a0;
    part.anchor_end = owned.a1;
    part.forward_begin = owned.forward_begin;
    part.forward_end = owned.forward_end;
    part.score = owned.score;
    part.anchors = owned.kept;
    part.item_score = item.score;
    part.item_anchors = static_cast<int>(anchors.size());
    part.pool_subsc = owned.pool_subsc;
    part.pool_n_sub = owned.pool_n_sub;
    parts.push_back(part);
    voting::QueryBlock block;
    block.candidate = static_cast<voting::CandidateId>(item.candidate);
    block.query_tile_begin = voting::query_tile_for_position(
        owned.forward_begin, read_length, seed_length, tile_count);
    block.query_tile_end =
        voting::query_tile_for_position(owned.forward_end - 1, read_length,
                                        seed_length, tile_count) +
        1;
    // Anchor tiles are monotone in chain order.
    int supporting = 0;
    int previous = -1;
    for (int c = owned.a0; c < owned.a1; ++c) {
      const int tile = dna_forward_query_tile(
          anchors[static_cast<std::size_t>(c)].q, item.reverse, read_length,
          seed_length, tile_count);
      if (tile != previous) ++supporting;
      previous = tile;
    }
    block.supporting_tiles = std::max(1, supporting);
    blocks.push_back(block);
  }
  std::vector<voting::CandidateId> assignment(
      static_cast<std::size_t>(tile_count), voting::kNullCandidate);
  std::size_t holder = 0;
  for (int tile = 0; tile < tile_count; ++tile) {
    const int middle =
        (voting::query_tile_begin(tile, read_length, seed_length, tile_count) +
         voting::query_tile_begin(tile + 1, read_length, seed_length,
                                  tile_count)) /
        2;
    while (holder + 1 < parts.size() && middle >= parts[holder].forward_end)
      ++holder;
    assignment[static_cast<std::size_t>(tile)] = blocks[holder].candidate;
  }
  family.partition.selected.assignment = std::move(assignment);
  family.partition.selected.blocks = std::move(blocks);
  family.partition.selected.non_null_blocks = static_cast<int>(parts.size());
  family.block_parts = std::move(parts);
  return true;
}

// Chains the top kDnaMapqRivalChains catalogue rivals of the committed family
// over the whole query, into result.rival_exact; only the MAPQ reads them.
// Rivals own no selected block; they are ranked by vote, then catalogue rank.
// An existing whole-query chain (the alternative, or an owner of the
// tile partition that owns no block) is reused.
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

const std::vector<chaining::Anchor>*
DnaPlacementChainingResult::selection_anchors(
    const DnaSelectionItem& item) const noexcept {
  const DnaPlacementCandidateChain* record =
      find(static_cast<::fa::cpu::voting::CandidateId>(item.candidate));
  if (record == nullptr) return nullptr;
  if (item.path < 0) return &record->primary;
  return static_cast<std::size_t>(item.path) < record->sibling_paths.size()
             ? &record->sibling_paths[static_cast<std::size_t>(item.path)]
             : nullptr;
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
    const std::vector<std::uint32_t>* fine_reverse_slots) {
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
      int low = std::numeric_limits<int>::max();
      int high = std::numeric_limits<int>::min();
      for (const chaining::Anchor& anchor : record.primary) {
        low = std::min(low, anchor.r - anchor.q);
        high = std::max(high, anchor.r - anchor.q);
      }
      candidate.screening_diagonals = true;
      candidate.screening_diagonal_low = low;
      candidate.screening_diagonal_high = high;

      // Screening may add support but never removes vote support.
      candidate.support = candidate.support | record.dense_support;
    }
  }
  family.partition = repartition(context, family);
  family.partition_deferred = false;

  if (!stabilize_selected_family(context, family, result, seed_index,
                                 forward_query, reverse_query)) {
    result.family = std::move(family);
    return result;
  }
  // The blocks are owned on the dense chains' anchors. A family with owners
  // none of whose chains qualifies maps nothing.
  if (result.accepted &&
      std::any_of(family.partition.selected.blocks.begin(),
                  family.partition.selected.blocks.end(),
                  [](const auto& block) {
                    return block.candidate !=
                           ::fa::cpu::voting::kNullCandidate;
                  }) &&
      !select_block_owners(context, family, result, result.kept_selection)) {
    result.accepted = false;
    result.no_owner_chain = true;
  }
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
  // Every candidate the alternative ranking restores (the alternative, and
  // ranks 2..n under -N n) takes its whole-query chain into its own record,
  // as an owner of the partition does, and the ownership selection runs
  // again with those chains as items.
  std::vector<::fa::cpu::voting::CandidateId> entered;
  if (result.accepted &&
      result.alternative.refusal == DnaAlternativeRefusal::None) {
    const std::vector<::fa::cpu::voting::CandidateId> ranked =
        rank_dna_alternative_hypotheses(
            context, family, result,
            static_cast<std::size_t>(context.opts.alternative_realize_max));
    for (const ::fa::cpu::voting::CandidateId id : ranked) {
      const auto record = std::find_if(
          result.candidates.begin(), result.candidates.end(),
          [id](const DnaPlacementCandidateChain& chain) {
            return chain.candidate == id;
          });
      if (record == result.candidates.end() || record->exact) continue;
      DnaPlacementCandidateChain exact;
      if (id == result.alternative.candidate)
        exact = result.alternative_exact;
      else if (!restore_alternative_exact(context, family, result, seed_index,
                                          forward_query, reverse_query, id,
                                          exact))
        continue;
      copy_chain_call(exact, true, *record);
      entered.push_back(id);
    }
    if (!entered.empty()) {
      select_block_owners(context, family, result, result.kept_selection);
      for (DnaPlacementCandidateChain& record : result.candidates)
        record.final_selected_tiles =
            selected_tiles(family.partition.selected, record.candidate);
    }
  }
  // After the alternative restore: the family is final, seed_index is alive
  // and the alternative's chain can be reused.
  chain_mapq_rivals(context, family, result, seed_index, forward_query,
                    reverse_query);
  result.selection_changed =
      initial_assignment != family.partition.selected.assignment;
  result.family = std::move(family);
  return result;
}

}  // namespace fa::cpu::lr
