// Single-pass two-strand RNA vote capture. The reverse seed stream is the forward stream
// with projected query coordinates and the same canonical keys, so each forward seed's
// posting view is resolved once and reused by both strands. Postings append to a flat
// per-strand log while a reused counter tracks each bin's distinct-seed support, so bins
// below min_support are dropped before any per-bin container exists. Peaks are emitted
// in ascending bin key, forward strand first.
#pragma once

#include "../../index/format.h"    // KmerPostingView
#include "../../index/seed.h"      // QuerySeed
#include "../../seeding/context.h" // LongReadSeedContext + occ policy
#include "../../seeding/select.h" // QuerySeedOccurrenceEvidence, select_spread_from_evidence
#include "../../seeding/types.h" // VoteHit, SeedDrop
#include "types.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

// One captured posting with its diagonal-bin key and the bin's stable ordinal in the
// read's DiagonalBinCounter, so the drain reads bin support by index.
struct CapturedVote {
  std::int64_t key = 0;
  VoteHit hit{};
  std::uint32_t bin = 0;
};

// Reused open-addressing counter of distinct-seed support and postings per bin key. Keys
// are vote_pack_key values (>= 0), so -1 marks an empty slot; reset() once per read.
// Slots map key -> bin ordinal and the counts live in `bins_`. Ordinals are assigned in
// first-touch order and survive grow(), so add()'s return value can be stored.
class DiagonalBinCounter {
public:
  void reset(std::size_t expected_bins);
  // `first_posting_of_seed_in_bin` counts a new distinct seed; a seed's postings in one
  // bin are contiguous, so the caller passes key != previous key. Returns the key's
  // stable bin ordinal.
  std::uint32_t add(std::int64_t key, bool first_posting_of_seed_in_bin);
  std::uint32_t support(std::int64_t key) const;
  // support(key) by bin ordinal, with no hash probe.
  std::uint32_t support_at(std::uint32_t bin) const {
    return bins_[bin].support;
  }
  std::size_t bin_count() const { return bins_.size(); }
  std::uint64_t max_bin_postings() const { return max_bin_postings_; }

private:
  // Generation-stamped slots make reset() O(1): a slot is live only while its
  // generation matches the counter's.
  struct Slot {
    std::int64_t key = -1;
    std::uint32_t generation = 0;
    std::uint32_t bin = 0;
  };
  struct Bin {
    std::uint32_t support = 0;
    std::uint32_t postings = 0;
  };
  bool live(const Slot& slot) const { return slot.generation == generation_; }
  void grow();
  std::vector<Slot> slots_;
  std::vector<Bin> bins_;
  std::size_t mask_ = 0;
  std::uint32_t generation_ = 0;
  std::uint64_t max_bin_postings_ = 0;
};

struct FusedCaptureScratch {
  std::vector<CapturedVote> fwd_log;
  std::vector<CapturedVote> rc_log;
  DiagonalBinCounter fwd_bins;
  DiagonalBinCounter rc_bins;
  // One posting-list view per forward seed, resolved once and reused by the
  // per-strand occurrence evidence, the tile-spread selection, and the vote.
  std::vector<KmerPostingView> views;
  std::vector<::fa::cpu::QuerySeedOccurrenceEvidence> evidence;
  std::vector<QuerySeed> selected_fwd;
  std::vector<QuerySeed> selected_rc;
  // Aggregation scratch (survivor filter + per-bin peak math).
  std::vector<CapturedVote> survivors;
  std::vector<int> read_pos_scratch;
  std::vector<int> ref_start_scratch;
  std::vector<int> seed_id_scratch;
};

// Post-selection seed counts per strand.
struct FusedCaptureSummary {
  std::size_t fwd_voted_seeds = 0;
  std::size_t rc_voted_seeds = 0;
  // Strand-compatibility counts per lane: `inspected` postings reached the test,
  // `compatible` passed it and were logged, `rejected` failed it.
  std::uint64_t fwd_postings_inspected = 0;
  std::uint64_t fwd_postings_compatible = 0;
  std::uint64_t fwd_postings_rejected_strand = 0;
  std::uint64_t rc_postings_inspected = 0;
  std::uint64_t rc_postings_compatible = 0;
  std::uint64_t rc_postings_rejected_strand = 0;
};

// Captures votes over the forward canonical seed stream and its reverse projection with
// one posting lookup per forward seed: each strand's tile-spread selection
// (ctx.max_query_seeds_per_strand) is replayed from the cached evidence, then both
// strands vote under the per-seed gates (index hit, posting count, occurrence policy).
FusedCaptureSummary capture_votes_fused(
    const LongReadSeedContext& ctx, const std::vector<QuerySeed>& fwd_seeds,
    const std::vector<QuerySeed>& rc_seeds, int read_len,
    FusedCaptureScratch& scratch);

// Aggregates the captured logs into diagonal peaks: forward-strand bins in ascending key
// order, then reverse-strand bins.
std::vector<CoarseDiagonalPeak>
aggregate_fused_capture_to_peaks(FusedCaptureScratch& scratch, int k,
                                 int min_support);

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
