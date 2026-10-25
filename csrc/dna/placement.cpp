// Placement chaining from the catalogue: one whole-query dense pass per (contig, strand) of its candidates and the ownership of the read's blocks on those chains.
#include "dna.h"

#include "../chain/chain.h"
#include "../core/core.h"
#include "../index/format.h"
#include "../vote/vote.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

/**** placement_chaining_internal ****/

namespace fa::cpu::lr {
namespace internal {

int pass_diagonal_band(const DnaContext& context) {
  return std::max(1, context.opts.cigar_dp_bw_long);
}

HarvestWindow harvest_window(const DnaContext& context,
                             const VotePeak& peak, int read_length,
                             std::int64_t chromosome_length) {
  const int interval_pad =
      context.opts.cigar_local_interval_anchor_interval_pad;
  const std::int64_t expected = peak.raw_ref_start;
  return {std::max<std::int64_t>(
              0, expected - peak.harvest_below - interval_pad),
          std::min<std::int64_t>(chromosome_length,
                                 expected + read_length + peak.harvest_above +
                                     interval_pad)};
}

bool skips_own_diagonal(const DnaContext& context, int contig,
                        bool reverse) {
  return !reverse && contig == context.self_contig;
}

PoolGate pool_gate(const DnaContext& context) {
  const int pool_gate_occ = context.opts.dna_pool_gate_occ;
  PoolGate gate;
  gate.gates = pool_gate_occ > 0;
  gate.cap = static_cast<std::uint32_t>(std::max(1, pool_gate_occ));
  gate.skips_keys = gate.gates;
  return gate;
}

DnaPlacementChainStatus contig_lane_refusal(
    const DnaContext& context, int contig,
    const std::vector<RetainedSeedRef>& seeds,
    const std::vector<std::uint8_t>& query, int read_length) {
  if (seeds.empty()) return DnaPlacementChainStatus::MissingSeeds;
  if (context.ref.index == nullptr || contig < 0 ||
      contig >= context.ref.contig_count() ||
      contig >= static_cast<int>(context.ref.index->chrom_count()) ||
      (!query.empty() &&
       query.size() != static_cast<std::size_t>(read_length)) ||
      context.ref.index->chrom_offsets_data() == nullptr)
    return DnaPlacementChainStatus::InvalidReference;
  return DnaPlacementChainStatus::Accepted;
}

}  // namespace internal
}  // namespace fa::cpu::lr

/**** placement_chaining ****/

namespace fa::cpu::lr {
namespace {

struct DeferredSlice {
  RetainedSeedRef seed;
  KmerPostingIntervalView interval;
};

int selected_tiles(
    const ::fa::cpu::voting::QueryPartitionPath& path,
    ::fa::cpu::voting::CandidateId candidate) {
  return static_cast<int>(std::count(
      path.assignment.begin(), path.assignment.end(), candidate));
}

// Per-thread radix ping-pong buffer and histograms for the pool sort below.
::fa::cpu::radix::Scratch<chaining::Anchor>& pool_sort_scratch() {
  static thread_local ::fa::cpu::radix::Scratch<chaining::Anchor> scratch;
  return scratch;
}

// Per-thread slices of one merged interval of a group's harvest
// (harvest_group_anchors), by density entry: an entry's interval is this
// interval's while its stamp equals `pass`.
struct PassSlices {
  std::vector<std::uint32_t> stamp;
  std::vector<KmerPostingIntervalView> interval;
  std::uint32_t pass = 0;
};
PassSlices& thread_pass_slices() {
  static thread_local PassSlices slices;
  return slices;
}
PassSlices& pass_slices(std::size_t entries) {
  PassSlices& slices = thread_pass_slices();
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

// Clears what a whole-query pass writes into `record`.
void reset_pass_record(DnaPlacementCandidateChain& record) {
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
  record.deferred_anchors = 0;
  record.rescued_anchors = 0;
  record.dense_runs = 0;
}

// The candidate's whole-query window on its contig, of `chromosome_length`;
// empty when high <= low.
internal::HarvestWindow pass_window(const DnaContext& context,
                                    const DnaPlacementFamily& family,
                                    const DnaPlacementCandidate& candidate,
                                    int chromosome_length) {
  const std::int64_t expected = candidate.peak.raw_ref_start;
  internal::HarvestWindow window = internal::harvest_window(
      context, candidate.peak, family.read_length, chromosome_length);
  // The window reaches min(band, max(pad, L / 8)) further on both sides, and
  // the line-widened winner keeps its wider bound.
  {
    const int interval_pad =
        context.opts.cigar_local_interval_anchor_interval_pad;
    const std::int64_t widen = std::min<std::int64_t>(
        internal::pass_diagonal_band(context),
        std::max<std::int64_t>(interval_pad, family.read_length / 8));
    window.low = std::min(
        window.low,
        std::max<std::int64_t>(0, expected - widen - interval_pad));
    window.high = std::max(
        window.high,
        std::min<std::int64_t>(
            chromosome_length,
            expected + family.read_length + widen + interval_pad));
  }
  return window;
}

// A member of a group's pass: its candidate's index in the family and its
// pass window.
struct PassMember {
  std::size_t index = 0;
  internal::HarvestWindow window;
};

// Appends to `anchors` the harvest of a (contig, strand) group's `members`
// from `seeds`, its lane's fine seeds. The members' windows are merged into
// disjoint intervals, overlapping or touching ones joined, and each distinct
// seed entry is sliced once per interval through the pool gate; every
// admissible slice is deferred and its anchors are built after the scan. A
// member admits a posting inside its window that lies on the contig and the
// query, off the read's own diagonal, and whose orientation bit agrees with
// the lane, so mapping reads no base; no diagonal band is tested. A posting
// that some member admits gives one anchor, flagged ANCHOR_TANDEM when its
// seed is (RetainedSeedRef::tandem). Each member's interval_hits,
// filtered_hits, deferred_anchors and rescued_anchors are those of a harvest
// over its window alone.
void harvest_group_anchors(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<PassMember>& members,
    RetainedSeedDensity& seed_index, const std::vector<RetainedSeedRef>& seeds,
    DnaPlacementChainingResult& result,
    std::vector<chaining::Anchor>& anchors) {
  const VotePeak& lead = family.candidates[members.front().index].peak;
  const bool reverse = lead.is_rc;
  const int contig = lead.chr;
  const int seed_length = family.seed_length;
  const int query_length = family.read_length;
  const int chromosome_length =
      static_cast<int>(context.ref.contig_length(contig));
  const std::uint64_t chromosome_base =
      context.ref.index
          ->chrom_offsets_data()[static_cast<std::size_t>(contig)];
  const internal::PoolGate gate = internal::pool_gate(context);
  const bool skip_own_diagonal =
      internal::skips_own_diagonal(context, contig, reverse);

  // The members by window low, each window as contig offsets.
  struct Admitter {
    std::uint64_t low = 0;
    std::uint64_t high = 0;
    DnaPlacementCandidateChain* record = nullptr;
  };
  std::vector<Admitter> admitters;
  admitters.reserve(members.size());
  for (const PassMember& member : members)
    admitters.push_back({static_cast<std::uint64_t>(member.window.low),
                         static_cast<std::uint64_t>(member.window.high),
                         &result.candidates[member.index]});
  std::stable_sort(admitters.begin(), admitters.end(),
                   [](const Admitter& left, const Admitter& right) {
                     return left.low < right.low;
                   });
  // The merged intervals, each over admitters [first, last).
  struct Interval {
    std::uint64_t low = 0;
    std::uint64_t high = 0;
    std::size_t first = 0;
    std::size_t last = 0;
  };
  std::vector<Interval> intervals;
  for (std::size_t at = 0; at < admitters.size(); ++at) {
    if (intervals.empty() || admitters[at].low > intervals.back().high) {
      intervals.push_back(
          {admitters[at].low, admitters[at].high, at, at + 1});
    } else {
      intervals.back().high =
          std::max(intervals.back().high, admitters[at].high);
      intervals.back().last = at + 1;
    }
  }

  // No anchor is built in the scan; every admissible slice is deferred.
  std::vector<DeferredSlice> deferred;
  deferred.reserve(seeds.size());
  constexpr std::size_t kSliceBlock = RetainedSeedDensity::kSliceBatch;
  std::uint32_t block_entries[kSliceBlock];
  KmerPostingIntervalView block_intervals[kSliceBlock];
  for (const Interval& bounds : intervals) {
    Admitter* const first_admitter = admitters.data() + bounds.first;
    Admitter* const last_admitter = admitters.data() + bounds.last;
    // Every seed resolves against the same (chr, low, high) interval, so each
    // distinct entry is sliced once, in batches through slice_batch; seeds are
    // still processed in order.
    PassSlices& slices = pass_slices(seed_index.entry_count());
    deferred.clear();
    for (std::size_t next = 0; next < seeds.size();) {
      // Take seeds until kSliceBlock entries new to this interval are
      // collected. A skipped entry gets an interval that is not found, which
      // the loop below drops as it drops the gated slice.
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
      seed_index.slice_batch(*context.ref.index, block_entries, block_size,
                             contig, static_cast<std::uint32_t>(bounds.low),
                             static_cast<std::uint32_t>(bounds.high),
                             block_intervals);
      for (std::size_t lane = 0; lane < block_size; ++lane)
        slices.interval[block_entries[lane]] = block_intervals[lane];
      for (std::size_t at = first; at < next; ++at) {
        const RetainedSeedRef& seed = seeds[at];
        const KmerPostingIntervalView& interval = slices.interval[seed.entry];
        // A rescued seed's fine twin passes the gate
        // (RetainedSeedRef::rescued). Its key's interval is shared with the
        // seeds that do not.
        const bool over_pool_gate =
            gate.drops(interval.global_count, seed.rescued);
        if (!interval.found() || interval.count == 0) continue;
        if (over_pool_gate) {
          // Each member counts the gated postings inside its window.
          for (Admitter* admitter = first_admitter;
               admitter != last_admitter; ++admitter)
            admitter->record->filtered_hits += static_cast<int>(
                interval.positions.lower_bound(chromosome_base +
                                               admitter->high) -
                interval.positions.lower_bound(chromosome_base +
                                               admitter->low));
          continue;
        }
        deferred.push_back({seed, interval});
      }
    }

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
      const KmerPostingIntervalView& interval = deferred[slot].interval;
      const RetainedSeedRef& seed = deferred[slot].seed;
      const int query_position = seed.seed.read_pos;
      const bool query_ok = query_position >= 0 &&
                            query_position + seed_length <= query_length;
      for (std::uint32_t posting = 0; posting < interval.count; ++posting) {
        // Inside the interval, so inside some member's window and the contig.
        const std::uint64_t local =
            interval.positions[posting] - chromosome_base;
        const int reference_position = static_cast<int>(local);
        const bool geometry_ok =
            query_ok &&
            !(skip_own_diagonal && reference_position == query_position) &&
            reference_position + seed_length <= chromosome_length;
        // The orientation test, run once per posting: -1 until then.
        int verified = -1;
        bool admitted = false;
        for (Admitter* admitter = first_admitter; admitter != last_admitter;
             ++admitter) {
          // The windows ascend by low: no later one holds the posting.
          if (local < admitter->low) break;
          if (local >= admitter->high) continue;
          DnaPlacementCandidateChain& record = *admitter->record;
          ++record.interval_hits;
          ++record.deferred_anchors;
          if (!geometry_ok) {
            ++record.filtered_hits;
            continue;
          }
          if (verified < 0) {
            // The posting was found under this seed's canonical key, so the
            // k-mers match exactly when the orientation bits agree with the
            // lane.
            verified = packed_ref_orientation_compatible(
                           seed.seed.z, interval.positions.packed_at(posting),
                           reverse)
                           ? 1
                           : 0;
          }
          if (verified == 0) {
            ++record.filtered_hits;
            continue;
          }
          ++record.rescued_anchors;
          admitted = true;
        }
        if (!admitted) continue;
        const bool tandem = seed.tandem;
        anchors.push_back(
            {reference_position, query_position, seed_length,
             static_cast<std::int32_t>(anchors.size()),
             tandem ? static_cast<std::uint32_t>(chaining::ANCHOR_TANDEM)
                    : 0u});
      }
    }
  }
}

// One pass's chains over its pool: the chain call's result, its partition, and
// the call's statistics and parameters.
struct PassChains {
  chaining::ChainResult chained;
  chaining::ChainPartition partition;
  chaining::DenseChainStats stats;
  chaining::DenseChainParams params;
};

// The chain half of a pass: sorts and deduplicates the pool, chains it densely,
// and partitions the chains. Sets the record's run count.
void chain_pool(const DnaContext& context, const DnaPlacementFamily& family,
                std::vector<chaining::Anchor> anchors,
                DnaPlacementCandidateChain& record,
                PassChains& out) {
  internal::sort_unique_pool_anchors(anchors);

  const chaining::ColinearChainParams chain_params =
      dna_candidate_chain_params(context, internal::pass_diagonal_band(context),
                                 family.seed_length);
  chaining::DenseChainStats& dense_stats = out.stats;
  out.params = dna_dense_chain_params(chain_params, family.seed_length,
                                      context.opts.dna_dense_diag_min_runs);
  // The dense chain: the pool collapsed to exact diagonal runs.
  out.chained = chaining::chain_dense_colinear(std::move(anchors), out.params,
                                              &dense_stats);
  const chaining::ChainResult& chained = out.chained;
  record.dense_runs = dense_stats.runs;
  out.partition =
      chaining::partition_chains(chained, ::fa::cpu::split::kHalfFloor);
}

// The catalogue's whole-query passes, one per (contig, strand) group of its
// candidates, groups in catalogue order of their first member. The group is
// harvested once over its members' windows (harvest_group_anchors) into one
// pool, which is deduplicated, chained and partitioned once. The partition's
// primary and every other chain of at least two anchors go to one member:
// among the members whose window holds the chain's median anchor (else among
// all), the one whose expected diagonal is nearest that anchor's, ties to
// catalogue order. Each member's chains are written with its best as primary.
// The member holding the partition's primary keeps its f2 and n_sub; another
// takes its best other chain's score as f2 and n_sub 0. A member left without
// a chain has the status its own empty pass would have. Every record is
// exact. Each member keeps its harvest counters; the group's run count goes to
// its first member with a window. A one-member group is that candidate's own
// pass.
void chain_pooled_groups(const DnaContext& context,
                         const DnaPlacementFamily& family,
                         DnaPlacementChainingResult& result,
                         RetainedSeedDensity& seed_index,
                         const std::vector<std::uint8_t>& forward_query,
                         const std::vector<std::uint8_t>& reverse_query) {
  std::vector<std::int64_t> keys;
  std::vector<std::vector<std::size_t>> groups;
  for (std::size_t index = 0; index < family.candidates.size(); ++index) {
    const VotePeak& peak = family.candidates[index].peak;
    const std::int64_t key =
        static_cast<std::int64_t>(peak.chr) * 2 + (peak.is_rc ? 1 : 0);
    const auto found = std::find(keys.begin(), keys.end(), key);
    if (found == keys.end()) {
      keys.push_back(key);
      groups.emplace_back(1, index);
    } else {
      groups[static_cast<std::size_t>(found - keys.begin())].push_back(index);
    }
  }
  std::vector<PassMember> windowed;
  std::vector<std::vector<std::size_t>> owned;
  for (const std::vector<std::size_t>& group : groups) {
    const VotePeak& lead = family.candidates[group.front()].peak;
    const bool reverse = lead.is_rc;
    const std::vector<std::uint8_t>& query =
        reverse ? reverse_query : forward_query;
    for (const std::size_t index : group) {
      reset_pass_record(result.candidates[index]);
      result.candidates[index].exact = true;
    }
    // The group's pass.
    DnaPlacementCandidateChain pass;
    const std::vector<RetainedSeedRef>& seeds =
        reverse ? seed_index.fine_reverse() : seed_index.fine_forward();
    const DnaPlacementChainStatus refusal = internal::contig_lane_refusal(
        context, lead.chr, seeds, query, family.read_length);
    if (refusal != DnaPlacementChainStatus::Accepted) {
      for (const std::size_t index : group)
        result.candidates[index].status = refusal;
      continue;
    }
    const int chromosome_length =
        static_cast<int>(context.ref.contig_length(lead.chr));
    windowed.clear();
    for (const std::size_t index : group) {
      const internal::HarvestWindow member_window = pass_window(
          context, family, family.candidates[index], chromosome_length);
      if (member_window.high <= member_window.low) {
        result.candidates[index].status =
            DnaPlacementChainStatus::InvalidReference;
        continue;
      }
      windowed.push_back({index, member_window});
    }
    std::vector<chaining::Anchor> anchors;
    if (!windowed.empty())
      harvest_group_anchors(context, family, windowed, seed_index, seeds,
                            result, anchors);
    // A member's own pass with no chain: no interval hit, or no chain.
    const auto chainless = [&result](std::size_t index) {
      DnaPlacementCandidateChain& record = result.candidates[index];
      record.status = record.interval_hits == 0
                          ? DnaPlacementChainStatus::NoIntervalHits
                          : DnaPlacementChainStatus::NoChain;
    };
    if (anchors.empty()) {
      for (const PassMember& member : windowed) chainless(member.index);
      continue;
    }
    DnaPlacementCandidateChain& first_record =
        result.candidates[windowed.front().index];
    PassChains chains;
    chain_pool(context, family, std::move(anchors), pass, chains);
    first_record.dense_runs = pass.dense_runs;
    const chaining::ChainResult& chained = chains.chained;
    const chaining::ChainPartition& partition = chains.partition;

    owned.assign(windowed.size(), {});
    for (std::size_t which = 0;
         partition.primary >= 0 && which < chained.chains.size(); ++which) {
      const chaining::Chain& chain = chained.chains[which];
      if (static_cast<int>(which) != partition.primary && chain.idx.size() < 2)
        continue;
      const chaining::Anchor& median = chained.anchors[static_cast<std::size_t>(
          chain.idx[chain.idx.size() / 2])];
      const std::int64_t diagonal =
          static_cast<std::int64_t>(median.r) - median.q;
      std::size_t best = 0;
      std::pair<bool, std::int64_t> best_key;
      for (std::size_t at = 0; at < windowed.size(); ++at) {
        const internal::HarvestWindow& bounds = windowed[at].window;
        const VotePeak& peak = family.candidates[windowed[at].index].peak;
        const std::int64_t expected =
            peak.line_gate
                ? peak.line_a + vote_slope_stretch(median.q, peak.line_b_q20)
                : peak.raw_ref_start;
        const std::pair<bool, std::int64_t> key{
            median.r < bounds.low || median.r >= bounds.high,
            diagonal > expected ? diagonal - expected : expected - diagonal};
        if (at == 0 || key < best_key) {
          best = at;
          best_key = key;
        }
      }
      owned[best].push_back(which);
    }

    for (std::size_t at = 0; at < windowed.size(); ++at) {
      DnaPlacementCandidateChain& record =
          result.candidates[windowed[at].index];
      const std::vector<std::size_t>& mine = owned[at];
      if (mine.empty()) {
        chainless(windowed[at].index);
        continue;
      }
      const bool holds_primary =
          std::find(mine.begin(), mine.end(),
                    static_cast<std::size_t>(partition.primary)) != mine.end();
      std::size_t primary = static_cast<std::size_t>(partition.primary);
      if (!holds_primary) {
        primary = mine.front();
        for (const std::size_t which : mine)
          if (chained.chains[which].score > chained.chains[primary].score)
            primary = which;
      }
      internal::write_primary_chain(family, reverse, chained,
                                    chained.chains[primary], record);
      if (holds_primary) {
        record.rival_chain_score = partition.f2;
        record.overlapping_rivals = partition.n_sub;
      }
      for (const std::size_t which : mine) {
        if (which == primary) continue;
        const chaining::Chain& sibling = chained.chains[which];
        std::vector<chaining::Anchor> path;
        path.reserve(sibling.idx.size());
        for (const std::int32_t index : sibling.idx)
          path.push_back(chained.anchors[static_cast<std::size_t>(index)]);
        // The sibling the MAPQ can realize as the rival: the partition's f2
        // chain, or the member's best other chain.
        const bool rival =
            holds_primary ? static_cast<int>(which) == partition.f2_index
                          : (record.rival_sibling < 0 ||
                             sibling.score > record.rival_chain_score);
        if (rival) {
          if (!holds_primary) record.rival_chain_score = sibling.score;
          record.rival_sibling = static_cast<int>(record.sibling_paths.size());
        }
        record.sibling_paths.push_back(std::move(path));
        record.sibling_scores.push_back(sibling.score);
      }
      record.status = DnaPlacementChainStatus::Accepted;
    }
  }
}

// The blocks owned on the accepted whole-query chains' anchors
// (dna/dna.h), every candidate's primary and siblings: each block becomes one
// block of partition.selected, in query order, with its chain range and bounds
// in family.block_parts, and a tile goes to the block holding its middle.
// False, with the partition left as it is, when no chain qualifies to own.
// `kept` receives the selection's items and roles.
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

// Copies the whole-query chains of the top kDnaMapqRivalChains catalogue
// rivals of the committed family into result.rival_exact; only the MAPQ reads
// them. Rivals own no selected block; they are ranked by vote, then catalogue
// rank.
void chain_mapq_rivals(const DnaContext& context,
                       const DnaPlacementFamily& family,
                       DnaPlacementChainingResult& result) {
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
  // whole-query span). The union is copied below.
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
          const int span =
              record.forward_query_end - record.forward_query_begin;
          if (record.forward_query_begin < 0 || span <= 0) continue;
          const int overlap =
              std::min(record.forward_query_end, block_end) -
              std::max(record.forward_query_begin, block_begin);
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
    DnaRivalExactChain rival;
    rival.candidate = family.candidates[index].id;
    // Kept whatever the outcome; `status` says whether it chained.
    rival.chain = result.candidates[index];
    result.rival_exact.push_back(std::move(rival));
  }
}

// Placement: the catalogue's whole-query passes, one per
// (contig, strand) group of its candidates, then the ownership walk over every
// accepted chain's primary and siblings, which decides the blocks, then the
// MAPQ's rivals. The read is accepted when one chain is.
void chain_and_own_every_candidate(
    const DnaContext& context, DnaPlacementFamily& family,
    DnaPlacementChainingResult& result, RetainedSeedDensity& seed_index,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query) {
  result.strand_filter_seeds = &seed_index;
  chain_pooled_groups(context, family, result, seed_index, forward_query,
                      reverse_query);
  for (const DnaPlacementCandidateChain& record : result.candidates) {
    if (record.status == DnaPlacementChainStatus::Accepted)
      result.accepted = true;
  }
  if (result.accepted &&
      !select_block_owners(context, family, result, result.kept_selection)) {
    result.accepted = false;
    result.no_owner_chain = true;
  }
  for (DnaPlacementCandidateChain& record : result.candidates)
    record.final_selected_tiles =
        selected_tiles(family.partition.selected, record.candidate);
  chain_mapq_rivals(context, family, result);
}

}  // namespace

chaining::ColinearChainParams dna_candidate_chain_params(
    const DnaContext& context, int band, int seed_length) {
  chaining::ColinearChainParams params;
  params.bw = std::max(1, band);
  // minimap2's long join: the reach is at least the band.
  params.max_dist_x = params.max_dist_y =
      std::max(params.bw, context.opts.cigar_dp_max_gap);
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
    ::fa::cpu::voting::CandidateId candidate, bool realized) const noexcept {
  for (const DnaRivalExactChain& chained : rival_exact)
    if (chained.candidate == candidate) return &chained.chain;
  const DnaPlacementCandidateChain* record = find(candidate);
  if (record == nullptr || !record->exact) return nullptr;
  // Every chained record is exact; the MAPQ reads only the vote's picks
  // (rival_exact), the owners and the realized rival.
  if (!realized &&
      std::none_of(family.partition.selected.blocks.begin(),
                   family.partition.selected.blocks.end(),
                   [candidate](const auto& block) {
                     return block.candidate == candidate;
                   }))
    return nullptr;
  return record;
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
    const std::vector<std::uint32_t>* fine_reverse_slots,
    RetainedSeedDensity& seed_index) {
  DnaPlacementChainingResult result;
  result.candidates.reserve(family.candidates.size());
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    DnaPlacementCandidateChain record;
    record.candidate = candidate.id;
    record.sparse_support = candidate.support;
    result.candidates.push_back(std::move(record));
  }
  if (!family.valid || context.ref.index == nullptr ||
      context.ref.encoded == nullptr || fine_forward_seeds == nullptr ||
      fine_reverse_seeds == nullptr || lookup_cache == nullptr) {
    result.family = std::move(family);
    return result;
  }

  if (!seed_index.build(
          *context.ref.index, forward_seeds, reverse_seeds, fine_forward_seeds,
          fine_reverse_seeds, lookup_cache, fine_forward_slots,
          fine_reverse_slots, context.readmitted)) {
    result.family = std::move(family);
    return result;
  }
  // Realization reads the density when it gates the late inversion probe
  // (dna.h).
  if (context.opts.inversion_probe_local_gate &&
      context.opts.enable_full_read_cigar)
    result.inversion_gate_seeds = &seed_index;

  chain_and_own_every_candidate(context, family, result, seed_index,
                                forward_query, reverse_query);
  result.family = std::move(family);
  return result;
}

}  // namespace fa::cpu::lr

namespace fa::cpu::lr::dna {

std::size_t pool_sort_scratch_bytes() {
  const auto& scratch = pool_sort_scratch();
  return capacity_bytes(scratch.tmp, scratch.count);
}

std::size_t pass_slices_bytes() {
  const PassSlices& slices = thread_pass_slices();
  return capacity_bytes(slices.stamp, slices.interval);
}

void release_placement_scratch() {
  pool_sort_scratch() = {};
  thread_pass_slices() = PassSlices{};
}

}  // namespace fa::cpu::lr::dna
