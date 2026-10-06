#include "target_chaining.h"
#include "placement_chaining_internal.h"

#include "../chaining/dense_chain.h"
#include "../index/format.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

constexpr std::uint32_t kNoTarget = std::numeric_limits<std::uint32_t>::max();

// One seed's postings on one target's contig: [begin, end) of its key's posting
// list, the list's whole run on that contig. `seed` indexes the lane's fine
// seeds.
struct TargetHit {
  std::uint32_t seed = 0;
  std::uint32_t begin = 0;
  std::uint32_t end = 0;
};

// The hits of one target, in seed order.
struct TargetHits {
  const TargetHit* begin = nullptr;
  const TargetHit* end = nullptr;
};

// One key's postings on one candidate contig, [begin, end) of its list, and
// the contig's target per lane (kNoTarget where no candidate of the lane is).
struct KeyRun {
  std::uint32_t begin = 0;
  std::uint32_t end = 0;
  std::uint32_t target[2] = {kNoTarget, kNoTarget};
};

// The all-chains lane's postings of one read. A target is one (contig, lane)
// of the read's chained candidates, which are its members. Every key a
// target's pool can use is read once, from the read's lowest candidate contig
// to its highest, and each seed of the key gets one TargetHit per target of
// its lane that the key reaches. A target's hits are contiguous, in seed
// order. Per thread, reused across reads.
struct ReadTargets {
  // By contig * 2 + lane: the target, or kNoTarget.
  std::vector<std::uint32_t> target_of;
  // The target_of cells the last read set, by target.
  std::vector<std::uint32_t> target_cells;
  // By catalogue candidate: its target, or kNoTarget when it is not chained.
  std::vector<std::uint32_t> candidate_target;
  // The lanes that hold a target, and the lowest and highest target contig.
  bool lane_used[2] = {false, false};
  std::uint32_t first_contig = 0;
  std::uint32_t last_contig = 0;
  // The density entries to read, in order of first use.
  std::vector<std::uint32_t> entries;
  // By density entry: its runs, runs_begin kNoTarget when it is not read.
  std::vector<std::uint32_t> runs_begin;
  std::vector<std::uint32_t> runs_end;
  std::vector<KeyRun> runs;
  // Target t holds hits[hits_begin[t], hits_begin[t + 1]).
  std::vector<std::uint32_t> hits_begin;
  std::vector<std::uint32_t> cursor;
  std::vector<TargetHit> hits;
};

ReadTargets& read_targets() {
  static thread_local ReadTargets targets;
  return targets;
}

// Gives each chained candidate of the read its target, over the last read's,
// and sizes hits_begin by the targets.
void assign_targets(const DnaContext& context,
                    const DnaPlacementFamily& family, ReadTargets& out) {
  const std::size_t contigs = context.ref.index->chrom_count();
  for (const std::uint32_t cell : out.target_cells)
    out.target_of[cell] = kNoTarget;
  out.target_cells.clear();
  if (out.target_of.size() < 2 * contigs)
    out.target_of.resize(2 * contigs, kNoTarget);
  out.candidate_target.assign(family.candidates.size(), kNoTarget);
  out.hits.clear();
  out.lane_used[0] = false;
  out.lane_used[1] = false;
  out.first_contig = std::numeric_limits<std::uint32_t>::max();
  out.last_contig = 0;
  for (std::size_t index = 0; index < family.candidates.size(); ++index) {
    const VotePeak& peak = family.candidates[index].peak;
    if (peak.chr < 0 || static_cast<std::size_t>(peak.chr) >= contigs)
      continue;
    // --dual=no: a candidate on a contig ranked below the read is not chained,
    // so the pair prints from the read ranked first.
    if (context.dual_contig_rank != nullptr &&
        (*context.dual_contig_rank)[static_cast<std::size_t>(peak.chr)] <
            context.dual_rank)
      continue;
    const std::uint32_t contig = static_cast<std::uint32_t>(peak.chr);
    const int lane = peak.is_rc ? 1 : 0;
    const std::uint32_t cell = 2 * contig + static_cast<std::uint32_t>(lane);
    if (out.target_of[cell] == kNoTarget) {
      out.target_of[cell] = static_cast<std::uint32_t>(out.target_cells.size());
      out.target_cells.push_back(cell);
    }
    out.candidate_target[index] = out.target_of[cell];
    out.lane_used[lane] = true;
    out.first_contig = std::min(out.first_contig, contig);
    out.last_contig = std::max(out.last_contig, contig);
  }
  out.hits_begin.assign(out.target_cells.size() + 1, 0);
}

// Picks the keys the target pools use and splits each key's posting list into
// its runs on target contigs.
void split_key_runs(const DnaContext& context,
                    const RetainedSeedDensity& seed_index, ReadTargets& out) {
  // The keys the target pools use: every key of a lane with a candidate, but a
  // key over the pool gate is skipped as chain_candidate skips it.
  const internal::PoolGate gate =
      internal::pool_gate(context, /*whole_query=*/true);
  const std::vector<RetainedSeedRef>* lane_seeds[2] = {
      &seed_index.fine_forward(), &seed_index.fine_reverse()};
  out.runs_begin.assign(seed_index.entry_count(), kNoTarget);
  out.runs_end.assign(seed_index.entry_count(), 0);
  out.entries.clear();
  for (int lane = 0; lane < 2; ++lane) {
    if (!out.lane_used[lane]) continue;
    for (const RetainedSeedRef& seed : *lane_seeds[lane]) {
      const std::uint32_t entry = seed.entry;
      if (out.runs_begin[entry] != kNoTarget) continue;
      if (gate.skips(seed_index, entry)) continue;
      out.runs_begin[entry] = 0;
      out.entries.push_back(entry);
    }
  }

  // One pass over each list, split into its runs on candidate contigs. A
  // list ascends by contig, so it starts at the first candidate contig.
  out.runs.clear();
  const PackedRefPos first_posting = pack_ref_pos(out.first_contig, 0);
  constexpr std::size_t kPrefetchAhead = 4;
  for (std::size_t at = 0; at < out.entries.size(); ++at) {
    if (at + kPrefetchAhead < out.entries.size())
      __builtin_prefetch(
          seed_index.view(out.entries[at + kPrefetchAhead]).positions.data(),
          0, 1);
    const std::uint32_t entry = out.entries[at];
    const RefPosSpan& positions = seed_index.view(entry).positions;
    const PackedRefPos* postings = positions.data();
    const std::uint32_t count = positions.size();
    std::uint32_t posting =
        out.first_contig == 0
            ? 0
            : static_cast<std::uint32_t>(
                  std::lower_bound(postings, postings + count,
                                   first_posting) -
                  postings);
    out.runs_begin[entry] = static_cast<std::uint32_t>(out.runs.size());
    while (posting < count) {
      const std::uint32_t contig = packed_ref_contig(postings[posting]);
      if (contig > out.last_contig) break;
      std::uint32_t run_end = posting + 1;
      while (run_end < count &&
             packed_ref_contig(postings[run_end]) == contig)
        ++run_end;
      const std::uint32_t forward = out.target_of[2 * contig];
      const std::uint32_t reverse = out.target_of[2 * contig + 1];
      if (forward != kNoTarget || reverse != kNoTarget)
        out.runs.push_back({posting, run_end, {forward, reverse}});
      posting = run_end;
    }
    out.runs_end[entry] = static_cast<std::uint32_t>(out.runs.size());
  }
}

// Each lane's seeds in order, each with the runs of its key on its lane's
// targets: counted, then placed.
void place_target_hits(const RetainedSeedDensity& seed_index,
                       ReadTargets& out) {
  const std::size_t target_count = out.target_cells.size();
  const std::vector<RetainedSeedRef>* lane_seeds[2] = {
      &seed_index.fine_forward(), &seed_index.fine_reverse()};
  for (int place = 0; place < 2; ++place) {
    if (place == 1) {
      for (std::size_t target = 0; target < target_count; ++target)
        out.hits_begin[target + 1] += out.hits_begin[target];
      out.cursor.assign(out.hits_begin.begin(), out.hits_begin.end() - 1);
      out.hits.resize(out.hits_begin[target_count]);
    }
    for (int lane = 0; lane < 2; ++lane) {
      if (!out.lane_used[lane]) continue;
      const std::vector<RetainedSeedRef>& seeds = *lane_seeds[lane];
      for (std::size_t seed = 0; seed < seeds.size(); ++seed) {
        const std::uint32_t entry = seeds[seed].entry;
        if (out.runs_begin[entry] == kNoTarget) continue;
        for (std::uint32_t run = out.runs_begin[entry];
             run < out.runs_end[entry]; ++run) {
          const KeyRun& key_run = out.runs[run];
          const std::uint32_t target = key_run.target[lane];
          if (target == kNoTarget) continue;
          if (place == 0)
            ++out.hits_begin[target + 1];
          else
            out.hits[out.cursor[target]++] = {
                static_cast<std::uint32_t>(seed), key_run.begin, key_run.end};
        }
      }
    }
  }
}

void build_read_targets(const DnaContext& context,
                        const DnaPlacementFamily& family,
                        const RetainedSeedDensity& seed_index,
                        ReadTargets& out) {
  assign_targets(context, family, out);
  if (out.target_cells.empty()) return;
  split_key_runs(context, seed_index, out);
  place_target_hits(seed_index, out);
}

// Each target is chained once, over the union of its members' pools. A
// posting is an anchor when chain_candidate's whole-query pass would build it
// for at least one member, and it is flagged tandem when it is in any of those
// members' pools. A one-member target so builds chain_candidate's pool, in its
// order.

// A member's harvest window [low, high) on its contig (harvest_window) and its
// expected diagonal.
struct TargetMember {
  std::uint32_t candidate = 0;
  std::int64_t expected = 0;
  std::int64_t low = 0;
  std::int64_t high = 0;
};

struct TargetScratch {
  // By target: its members, in catalogue order.
  std::vector<std::uint32_t> member_begin;
  std::vector<std::uint32_t> member_cursor;
  std::vector<std::uint32_t> members;
  // One target's members with a window: in catalogue order, the ones without a
  // widened window by expected diagonal (and their diagonals), and the others.
  std::vector<TargetMember> windowed;
  std::vector<TargetMember> plain;
  std::vector<std::int64_t> plain_expected;
  std::vector<TargetMember> wide;
  std::vector<std::uint32_t> order;
};

TargetScratch& target_scratch() {
  static thread_local TargetScratch scratch;
  return scratch;
}

// Builds the union pool of one target into `anchors`; `scratch.windowed` is
// set.
void build_target_pool(const DnaContext& context,
                       const DnaPlacementFamily& family,
                       const RetainedSeedDensity& seed_index,
                       const std::vector<RetainedSeedRef>& seeds,
                       const TargetHits& hits, std::uint32_t contig,
                       bool reverse, int diagonal_band, TargetScratch& scratch,
                       std::vector<chaining::Anchor>& anchors) {
  const int read_length = family.read_length;
  const int seed_length = family.seed_length;
  const int chromosome_length =
      static_cast<int>(context.ref.contig_length(static_cast<int>(contig)));
  const int interval_pad =
      context.opts.cigar_local_interval_anchor_interval_pad;
  const int tandem_window = context.opts.dna_tandem_window;
  const internal::PoolGate gate =
      internal::pool_gate(context, /*whole_query=*/true);
  const bool skip_own_diagonal =
      internal::skips_own_diagonal(context, static_cast<int>(contig), reverse);

  scratch.plain.clear();
  scratch.wide.clear();
  std::int64_t target_low = std::numeric_limits<std::int64_t>::max();
  std::int64_t target_high = std::numeric_limits<std::int64_t>::min();
  for (const TargetMember& member : scratch.windowed) {
    const VotePeak& peak = family.candidates[member.candidate].peak;
    (peak.harvest_below == 0 && peak.harvest_above == 0 ? scratch.plain
                                                        : scratch.wide)
        .push_back(member);
    target_low = std::min(target_low, member.low);
    target_high = std::max(target_high, member.high);
  }
  std::stable_sort(scratch.plain.begin(), scratch.plain.end(),
                   [](const TargetMember& left, const TargetMember& right) {
                     return left.expected < right.expected;
                   });
  scratch.plain_expected.clear();
  for (const TargetMember& member : scratch.plain)
    scratch.plain_expected.push_back(member.expected);
  const std::int64_t* const expected_begin = scratch.plain_expected.data();
  const std::int64_t* const expected_end =
      expected_begin + scratch.plain_expected.size();

  const PackedRefPos window_low =
      pack_ref_pos(contig, static_cast<std::uint32_t>(target_low));
  const PackedRefPos window_high =
      pack_ref_pos(contig, static_cast<std::uint32_t>(target_high));
  for (const TargetHit* hit = hits.begin; hit != hits.end; ++hit) {
    const RetainedSeedRef& seed = seeds[hit->seed];
    const KmerPostingView& full = seed_index.view(seed.entry);
    if (gate.drops(full.occurrence, seed.rescued)) continue;
    const int query_position = seed.seed.read_pos;
    if (query_position < 0 || query_position + seed_length > read_length)
      continue;
    const PackedRefPos* postings = full.positions.data();
    const PackedRefPos* first =
        std::lower_bound(postings + hit->begin, postings + hit->end,
                         window_low);
    const PackedRefPos* last =
        std::lower_bound(first, postings + hit->end, window_high);
    for (const PackedRefPos* at = first; at != last; ++at) {
      const std::int64_t reference_position = packed_ref_local(*at);
      if (skip_own_diagonal && reference_position == query_position) continue;
      if (reference_position + seed_length > chromosome_length) continue;
      if (!packed_ref_orientation_compatible(seed.seed.z, *at, reverse))
        continue;
      const std::int64_t diagonal = reference_position - query_position;
      // The members admitting the posting, and the widest of their windows.
      std::int64_t admit_low = std::numeric_limits<std::int64_t>::max();
      std::int64_t admit_high = std::numeric_limits<std::int64_t>::min();
      // Without a widened window a member admits the posting exactly when
      // its expected diagonal is one of these starts.
      const internal::ExpectedStarts starts =
          internal::admitting_expected_starts(reference_position, diagonal,
                                              read_length, interval_pad,
                                              diagonal_band);
      if (starts.lowest <= starts.highest) {
        const std::int64_t* from =
            std::lower_bound(expected_begin, expected_end, starts.lowest);
        const std::int64_t* to =
            std::upper_bound(from, expected_end, starts.highest);
        if (from != to) {
          admit_low = scratch.plain[static_cast<std::size_t>(
                                        from - expected_begin)]
                          .low;
          admit_high = scratch.plain[static_cast<std::size_t>(
                                         to - expected_begin - 1)]
                           .high;
        }
      }
      for (const TargetMember& member : scratch.wide) {
        if (reference_position < member.low ||
            reference_position >= member.high ||
            !internal::within_band(diagonal, member.expected, diagonal_band))
          continue;
        admit_low = std::min(admit_low, member.low);
        admit_high = std::max(admit_high, member.high);
      }
      if (admit_low > admit_high) continue;
      // Over the union of the admitting members' windows.
      const bool tandem = internal::posting_is_tandem(
          at, postings + hit->begin, postings + hit->end, reference_position,
          [](const PackedRefPos* posting) -> std::int64_t {
            return packed_ref_local(*posting);
          },
          admit_low, admit_high, tandem_window);
      anchors.push_back(
          {static_cast<std::int32_t>(reference_position), query_position,
           seed_length, static_cast<std::int32_t>(anchors.size()),
           tandem ? static_cast<std::uint32_t>(chaining::ANCHOR_TANDEM) : 0u});
    }
  }
}

// Writes one chain of `chained` into a member's record as chain_candidate
// writes its accepted whole-query chain.
void accept_target_chain(const DnaPlacementFamily& family, bool reverse,
                         const chaining::ChainResult& chained,
                         const chaining::Chain& chain,
                         DnaPlacementCandidateChain& record) {
  internal::write_primary_chain(family, reverse, chained, chain, record);
  record.exact = true;
  record.status = DnaPlacementChainStatus::Accepted;
}

// Groups the read's chained candidates by target: target t's members are
// members[member_begin[t], member_begin[t + 1]), in catalogue order.
void group_members(const ReadTargets& targets, TargetScratch& scratch) {
  const std::size_t target_count = targets.target_cells.size();
  scratch.member_begin.assign(target_count + 1, 0);
  for (const std::uint32_t target : targets.candidate_target)
    if (target != kNoTarget) ++scratch.member_begin[target + 1];
  for (std::size_t target = 0; target < target_count; ++target)
    scratch.member_begin[target + 1] += scratch.member_begin[target];
  scratch.member_cursor.assign(scratch.member_begin.begin(),
                               scratch.member_begin.end() - 1);
  scratch.members.resize(scratch.member_begin[target_count]);
  for (std::size_t index = 0; index < targets.candidate_target.size();
       ++index) {
    const std::uint32_t target = targets.candidate_target[index];
    if (target != kNoTarget)
      scratch.members[scratch.member_cursor[target]++] =
          static_cast<std::uint32_t>(index);
  }
}

// Gives every member of `target`, on `contig`, the target's `refusal` when
// there is one, and otherwise its window: a member whose window is empty is
// refused, and the others go to scratch.windowed, in catalogue order, with
// NoChain until a chain reaches them.
void window_members(const DnaContext& context,
                    const DnaPlacementFamily& family, std::size_t target,
                    std::uint32_t contig, DnaPlacementChainStatus refusal,
                    TargetScratch& scratch,
                    std::vector<DnaPlacementCandidateChain>& chains) {
  const int read_length = family.read_length;
  const std::uint32_t* member_first =
      scratch.members.data() + scratch.member_begin[target];
  const std::uint32_t* member_last =
      scratch.members.data() + scratch.member_begin[target + 1];
  scratch.windowed.clear();
  for (const std::uint32_t* member = member_first; member != member_last;
       ++member) {
    DnaPlacementCandidateChain& record = chains[*member];
    record.exact = true;
    if (refusal != DnaPlacementChainStatus::Accepted) {
      record.status = refusal;
      continue;
    }
    const VotePeak& peak = family.candidates[*member].peak;
    const std::int64_t chromosome_length =
        context.ref.contig_length(static_cast<int>(contig));
    const internal::HarvestWindow window =
        internal::harvest_window(context, peak, read_length, chromosome_length);
    if (window.high <= window.low) {
      record.status = DnaPlacementChainStatus::InvalidReference;
      continue;
    }
    record.status = DnaPlacementChainStatus::NoChain;
    scratch.windowed.push_back(
        {*member, static_cast<std::int64_t>(
                      static_cast<int>(peak.raw_ref_start)),
         window.low, window.high});
  }
}

// Chains one target's pool with chain_candidate's whole-query parameters.
chaining::ChainResult chain_target_pool(const DnaContext& context,
                                        const DnaPlacementFamily& family,
                                        int diagonal_band,
                                        std::vector<chaining::Anchor> anchors) {
  const int read_length = family.read_length;
  const int seed_length = family.seed_length;
  const chaining::ColinearChainParams chain_params =
      dna_candidate_chain_params(context, diagonal_band, seed_length,
                                 read_length);
  return chaining::chain_dense_colinear(
      std::move(anchors),
      dna_dense_chain_params(chain_params, seed_length,
                             context.opts.dna_dense_diag_min_runs),
      nullptr);
}

// Hands one target's chains to its windowed members, best score first, in
// catalogue order, one each; chains past the member count are dropped.
void hand_out_chains(const DnaPlacementFamily& family, bool reverse,
                     const chaining::ChainResult& chained,
                     TargetScratch& scratch,
                     std::vector<DnaPlacementCandidateChain>& chains) {
  // Best score first; a tie keeps chain order, so the first is the chain
  // partition_chains makes primary.
  scratch.order.resize(chained.chains.size());
  for (std::size_t which = 0; which < scratch.order.size(); ++which)
    scratch.order[which] = static_cast<std::uint32_t>(which);
  std::stable_sort(scratch.order.begin(), scratch.order.end(),
                   [&chained](std::uint32_t left, std::uint32_t right) {
                     return chained.chains[left].score >
                            chained.chains[right].score;
                   });
  const std::size_t served =
      std::min(scratch.order.size(), scratch.windowed.size());
  for (std::size_t rank = 0; rank < served; ++rank)
    accept_target_chain(family, reverse, chained,
                        chained.chains[scratch.order[rank]],
                        chains[scratch.windowed[rank].candidate]);
}

// Chains every target of the read once, with chain_candidate's band and
// parameters. A target's chains, best score first, go to its members in
// catalogue order, at most one each; a member left without one keeps NoChain,
// and chains past the member count are dropped. `chains` is parallel to
// family.candidates.
void chain_read_targets(const DnaContext& context,
                        const DnaPlacementFamily& family,
                        const std::vector<std::uint8_t>& forward_query,
                        const std::vector<std::uint8_t>& reverse_query,
                        const RetainedSeedDensity& seed_index,
                        const ReadTargets& targets,
                        std::vector<DnaPlacementCandidateChain>& chains) {
  TargetScratch& scratch = target_scratch();
  group_members(targets, scratch);

  const std::size_t target_count = targets.target_cells.size();
  const int read_length = family.read_length;
  const int diagonal_band =
      internal::pass_diagonal_band(context, /*whole_query=*/true);
  for (std::size_t target = 0; target < target_count; ++target) {
    const std::uint32_t cell = targets.target_cells[target];
    const std::uint32_t contig = cell >> 1;
    const bool reverse = (cell & 1u) != 0;
    const std::vector<std::uint8_t>& query =
        reverse ? reverse_query : forward_query;
    const std::vector<RetainedSeedRef>& seeds =
        reverse ? seed_index.fine_reverse() : seed_index.fine_forward();
    // chain_candidate's refusals, which hold for every member alike.
    const DnaPlacementChainStatus refusal = internal::contig_lane_refusal(
        context, static_cast<int>(contig), seeds, query, read_length);
    window_members(context, family, target, contig, refusal, scratch, chains);
    if (scratch.windowed.empty()) continue;

    std::vector<chaining::Anchor> anchors;
    build_target_pool(context, family, seed_index, seeds,
                      {targets.hits.data() + targets.hits_begin[target],
                       targets.hits.data() + targets.hits_begin[target + 1]},
                      contig, reverse, diagonal_band, scratch, anchors);
    if (anchors.empty()) continue;
    internal::sort_unique_pool_anchors(anchors);
    const chaining::ChainResult chained =
        chain_target_pool(context, family, diagonal_band, std::move(anchors));
    hand_out_chains(family, reverse, chained, scratch, chains);
  }
}

}  // namespace

std::vector<DnaPlacementCandidateChain> build_dna_target_chains(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    const std::vector<ChainWindowRetainedSeed>* forward_seeds,
    const std::vector<ChainWindowRetainedSeed>* reverse_seeds,
    const std::vector<QuerySeed>* fine_forward_seeds,
    const std::vector<QuerySeed>* fine_reverse_seeds,
    ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots) {
  std::vector<DnaPlacementCandidateChain> chains;
  if (!family.valid || context.ref.index == nullptr ||
      fine_forward_seeds == nullptr || fine_reverse_seeds == nullptr ||
      lookup_cache == nullptr)
    return chains;
  RetainedSeedDensity seed_index;
  if (!seed_index.build(
          *context.ref.index, forward_seeds, reverse_seeds, fine_forward_seeds,
          fine_reverse_seeds, lookup_cache, fine_forward_slots,
          fine_reverse_slots))
    return chains;
  chains.resize(family.candidates.size());
  for (std::size_t index = 0; index < family.candidates.size(); ++index)
    chains[index].candidate = family.candidates[index].id;
  ReadTargets& targets = read_targets();
  build_read_targets(context, family, seed_index, targets);
  chain_read_targets(context, family, forward_query, reverse_query,
                     seed_index, targets, chains);
  return chains;
}

}  // namespace fa::cpu::lr
