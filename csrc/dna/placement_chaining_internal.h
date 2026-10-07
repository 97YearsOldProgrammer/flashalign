// The interface shared by the DNA chaining files, placement_chaining.cpp and
// target_chaining.cpp: the harvest rules, each stated once for
// chain_candidate's passes and the all-chains lane's target pools, and the
// pool and record steps both take.
#pragma once

#include "placement_chaining.h"
#include "retained_seed_density.h"
#include "../chaining/anchor.h"
#include "../chaining/result.h"
#include "../index/format.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace fa::cpu::lr {
namespace internal {

// The band of a chain pass, for both its harvest and its chain: -b's first
// value on the whole-query (dense) passes, its second on the screening pass.
inline int pass_diagonal_band(const DnaContext& context, bool whole_query) {
  return std::max(1, whole_query ? context.opts.chain_band
                                 : context.opts.screen_diag_band);
}

// A candidate's harvest window [low, high) on its contig: the read's span from
// its peak's expected start, padded on both sides by
// cigar_local_interval_anchor_interval_pad, widened by harvest_below /
// harvest_above and clamped to the contig. The widening reaches the whole-read
// winner's per-read line (vote_slope_widen_winner), which also centres that
// winner's band (append_interval_anchors), and is 0 on every other peak. Empty
// when high <= low.
struct HarvestWindow {
  std::int64_t low = 0;
  std::int64_t high = 0;
};
inline HarvestWindow harvest_window(const DnaContext& context,
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

// The band test: whether a posting on `diagonal` lies within `diagonal_band`
// of the candidate's expected diagonal, in the caller's integer type.
template <class Diagonal>
inline bool within_band(Diagonal diagonal, Diagonal expected,
                        int diagonal_band) {
  return std::abs(diagonal - expected) <= diagonal_band;
}

// The inverse of harvest_window and within_band for a peak with no widening:
// the expected starts [lowest, highest] whose window and band admit a posting
// at `reference_position` on `diagonal`, for the same read length, pad and
// band. Empty when lowest > highest.
struct ExpectedStarts {
  std::int64_t lowest = 0;
  std::int64_t highest = 0;
};
inline ExpectedStarts admitting_expected_starts(
    std::int64_t reference_position, std::int64_t diagonal, int read_length,
    int interval_pad, int diagonal_band) {
  return {std::max(reference_position - read_length - interval_pad + 1,
                   diagonal - diagonal_band),
          std::min(reference_position + interval_pad,
                   diagonal + diagonal_band)};
}

// Whether a pass over `contig` and the lane drops the posting at each seed's
// own read position: the read's own exact diagonal gives no anchor, so the
// contig is the read itself (DnaContext::self_contig) and the lane forward.
inline bool skips_own_diagonal(const DnaContext& context, int contig,
                               bool reverse) {
  return !reverse && contig == context.self_contig;
}

// Whether the posting at `at` of its key's sorted posting list [first, last),
// at `position`, has a neighbour in the list within `tandem_window` bases that
// lies in [low, high) (<= 0 disables). Such a seed pins one copy of a tandem
// array, chosen by accident, so its anchor is flagged ANCHOR_TANDEM and never
// becomes a realization corner. The raw list is used: a neighbour that fails
// the geometry tests still shows the key repeats nearby. `position_of(at)` is
// a posting's position in the coordinates of `position`, `low` and `high`.
template <class Position, class At, class PositionOf>
inline bool posting_is_tandem(At at, At first, At last, Position position,
                              const PositionOf& position_of, Position low,
                              Position high, int tandem_window) {
  if (tandem_window <= 0) return false;
  const Position window = static_cast<Position>(tandem_window);
  // Written to avoid unsigned wrap when `position` < `window`.
  if (at > first) {
    const Position previous = position_of(at - 1);
    if (previous >= low && previous + window >= position) return true;
  }
  if (at + 1 < last) {
    const Position next = position_of(at + 1);
    if (next < high && next <= position + window) return true;
  }
  return false;
}

// The harvest's occurrence gate. A seed whose key occurs more than `cap` times
// genome-wide gives no anchor unless it is rescued (RetainedSeedRef::rescued).
// The screening pass always gates, at the global cap; the whole-query pass
// gates at dna_pool_gate_occ (by default the vote's cap, or --max-chain-occ)
// and not at all when that is 0.
struct PoolGate {
  bool gates = false;
  std::uint32_t cap = 1;
  // A key over the cap is not sliced: nothing reads the postings the gate
  // drops. A key a rescued seed shares, which the gate admits, is sliced.
  bool skips_keys = false;

  // Whether the gate drops a seed whose key occurs `occurrence` times.
  bool drops(std::uint32_t occurrence, bool rescued) const {
    return gates && occurrence > cap && !rescued;
  }
  // Whether the key of the density entry is skipped unsliced.
  bool skips(const RetainedSeedDensity& seed_index,
             std::uint32_t entry) const {
    return skips_keys && seed_index.occurrence(entry) > cap &&
           !seed_index.holds_rescued(entry);
  }
};
inline PoolGate pool_gate(const DnaContext& context, bool whole_query) {
  const int pool_gate_occ = context.opts.dna_pool_gate_occ;
  PoolGate gate;
  gate.gates = !whole_query || pool_gate_occ > 0;
  gate.cap = static_cast<std::uint32_t>(
      std::max(1, whole_query ? pool_gate_occ
                              : context.opts.cigar_local_global_occ));
  gate.skips_keys = gate.gates;
  return gate;
}

// The refusal every candidate on `contig` in the lane of `seeds` meets alike,
// or Accepted when there is none: MissingSeeds when the lane has no seeds,
// InvalidReference when the index cannot serve the contig or `query` is not
// the read's length. Only metadata is needed: contig lengths come from the
// index offsets and the read length from the family, and a map-only run
// passes `query` empty.
inline DnaPlacementChainStatus contig_lane_refusal(
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

// Appends one anchor per posting in `interval` that passes the geometry tests
// and the orientation test for `reverse_lane`; no base is read.
// `chromosome_base` is the contig's offset in the flattened reference and
// `chromosome_length` its length, so postings off the contig are dropped. An
// anchor is flagged ANCHOR_TANDEM when posting_is_tandem holds over the
// contig. `skip_own_diagonal` is skips_own_diagonal's for the contig and lane.
// `query_length` bounds the query span.
// `line_peak`, when not null, is a peak whose line passed its gate: the band is
// then centred on that line at the seed's query position, not on
// `main_diagonal`.
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
    DnaPlacementCandidateChain& record);

// Sorts a pool by (r, q, span), flags descending, and keeps the first anchor of
// each (r, q, span), so the tandem flag survives.
void sort_unique_pool_anchors(std::vector<chaining::Anchor>& anchors);

// Writes `chain` of `chained` into `record` as its primary: the anchors, the
// score, the anchor count, the dense support and the spans.
void write_primary_chain(const DnaPlacementFamily& family, bool reverse,
                         const chaining::ChainResult& chained,
                         const chaining::Chain& chain,
                         DnaPlacementCandidateChain& record);

}  // namespace internal
}  // namespace fa::cpu::lr
