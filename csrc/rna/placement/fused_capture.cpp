#include "fused_capture.h"

#include "bin_peak_math.h"
#include "../../index/format.h" // chromosome_index_for_global_pos
#include "../../voting/state.h" // chain_window_seed_occ_allowed, vote_floor_div, vote_pack_key

#include <algorithm>
#include <climits>
#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

namespace {

inline std::uint64_t mix_bin_key(std::int64_t key) {
  std::uint64_t x = static_cast<std::uint64_t>(key);
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

// One-directional lower_bound cursor over the forward stream's read_pos. Each strand's
// selected seeds arrive in ascending read_pos, so the forward positions looked up ascend
// on the forward strand and descend on the reverse (fwd_pos = read_len - k - read_pos);
// a monotone cursor answers the binary search in O(selected + |fwd_seeds|).
class FwdIndexCursor {
public:
  FwdIndexCursor(const std::vector<QuerySeed>& fwd_seeds, bool descending)
      : seeds_(fwd_seeds), at_(descending ? fwd_seeds.size() : 0) {}

  std::size_t advance_to(int read_pos) {
    while (at_ < seeds_.size() && seeds_[at_].read_pos < read_pos)
      ++at_;
    return at_;
  }

  std::size_t retreat_to(int read_pos) {
    while (at_ > 0 && seeds_[at_ - 1].read_pos >= read_pos)
      --at_;
    return at_;
  }

private:
  const std::vector<QuerySeed>& seeds_;
  std::size_t at_;
};

} // namespace

void DiagonalBinCounter::reset(std::size_t expected_bins) {
  std::size_t capacity = 1024;
  const std::size_t target = expected_bins * 2 + 1;
  while (capacity < target)
    capacity <<= 1;
  if (slots_.size() < capacity) {
    slots_.assign(capacity, Slot{});
    generation_ = 0;
  }
  mask_ = slots_.size() - 1;
  bins_.clear();
  max_bin_postings_ = 0;
  if (++generation_ == 0) {
    // Generation wrap: one full clear per ~4 billion resets.
    std::fill(slots_.begin(), slots_.end(), Slot{});
    generation_ = 1;
  }
}

void DiagonalBinCounter::grow() {
  std::vector<Slot> old;
  old.swap(slots_);
  slots_.assign(old.size() * 2, Slot{});
  mask_ = slots_.size() - 1;
  for (const Slot& slot : old) {
    if (slot.generation != generation_ || slot.key < 0)
      continue;
    std::size_t idx = mix_bin_key(slot.key) & mask_;
    while (slots_[idx].generation == generation_)
      idx = (idx + 1) & mask_;
    slots_[idx] = slot;
  }
}

std::uint32_t DiagonalBinCounter::add(std::int64_t key,
                                      bool first_posting_of_seed_in_bin) {
  if ((bins_.size() + 1) * 10 >= slots_.size() * 7)
    grow();
  std::size_t idx = mix_bin_key(key) & mask_;
  while (live(slots_[idx]) && slots_[idx].key != key)
    idx = (idx + 1) & mask_;
  Slot& slot = slots_[idx];
  if (!live(slot)) {
    slot.key = key;
    slot.generation = generation_;
    slot.bin = static_cast<std::uint32_t>(bins_.size());
    bins_.push_back(Bin{});
  }
  Bin& bin = bins_[slot.bin];
  ++bin.postings;
  if (first_posting_of_seed_in_bin)
    ++bin.support;
  if (bin.postings > max_bin_postings_)
    max_bin_postings_ = bin.postings;
  return slot.bin;
}

std::uint32_t DiagonalBinCounter::support(std::int64_t key) const {
  std::size_t idx = mix_bin_key(key) & mask_;
  while (live(slots_[idx])) {
    if (slots_[idx].key == key)
      return bins_[slots_[idx].bin].support;
    idx = (idx + 1) & mask_;
  }
  return 0;
}

FusedCaptureSummary capture_votes_fused(
    const LongReadSeedContext& ctx, const std::vector<QuerySeed>& fwd_seeds,
    const std::vector<QuerySeed>& rc_seeds, int read_len,
    FusedCaptureScratch& scratch) {
  FusedCaptureSummary summary;
  scratch.fwd_log.clear();
  scratch.rc_log.clear();
  scratch.fwd_bins.reset(fwd_seeds.size());
  scratch.rc_bins.reset(rc_seeds.size());
  const int k = ctx.k;
  if (read_len < k || ctx.index == nullptr || ctx.index->empty() ||
      fwd_seeds.empty())
    return summary;
  const std::uint64_t* chr_bounds = ctx.index->chrom_offsets_data();
  const int n_chr = ctx.chr_names ? static_cast<int>(ctx.chr_names->size()) : 0;
  if (!chr_bounds || n_chr <= 0)
    return summary;
  const int W = effective_vote_diag_bin_width(ctx, read_len);
  // Strand-compatible vote: a posting of a selected seed's canonical key is logged on
  // lane h only when (z_i ^ z_ref(p)) == h (packed_ref_orientation_compatible), i.e. the
  // oriented query k-mer and the reference k-mer at the posting are the same bases.

  // One posting-list lookup per forward seed. The reverse projection keeps the canonical
  // key, so these views serve both strands.
  std::vector<KmerPostingView>& views = scratch.views;
  views.clear();
  views.reserve(fwd_seeds.size());
  for (const QuerySeed& seed : fwd_seeds)
    views.push_back(ctx.index->lookup(seed.key));

  // Per-strand low-occurrence tile-spread selection, replayed from the cached evidence
  // (limit_query_seeds semantics).
  const int max_seeds = ctx.max_query_seeds_per_strand;
  const auto select_strand = [&](const std::vector<QuerySeed>& stream,
                                 bool reversed_views,
                                 std::vector<QuerySeed>& selected) {
    if (max_seeds <= 0 || static_cast<int>(stream.size()) <= max_seeds) {
      selected = stream;
      return;
    }
    std::vector<::fa::cpu::QuerySeedOccurrenceEvidence>& evidence =
        scratch.evidence;
    evidence.clear();
    evidence.reserve(stream.size());
    const std::size_t n = stream.size();
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t view_index = reversed_views ? n - 1 - j : j;
      const KmerPostingView& view =
          view_index < views.size() ? views[view_index] : KmerPostingView{};
      ::fa::cpu::QuerySeedOccurrenceEvidence e;
      e.seed = stream[j];
      e.present = view.found();
      e.retained = view.found() && view.count > 0;
      e.found = e.retained;
      e.occurrence = e.present ? view.occurrence : 0;
      e.specificity = ::fa::cpu::query_seed_specificity(e.occurrence);
      e.original_index = static_cast<int>(j);
      evidence.push_back(e);
    }
    ::fa::cpu::select_spread_from_evidence(selected, evidence, read_len,
                                           max_seeds);
  };
  select_strand(fwd_seeds, /*reversed_views=*/false, scratch.selected_fwd);
  select_strand(rc_seeds, /*reversed_views=*/true, scratch.selected_rc);
  summary.fwd_voted_seeds = scratch.selected_fwd.size();
  summary.rc_voted_seeds = scratch.selected_rc.size();

  // Vote per strand from the cached views; no further index lookups.
  const auto vote_strand = [&](const std::vector<QuerySeed>& selected,
                               bool is_rc, std::vector<CapturedVote>& log,
                               DiagonalBinCounter& bins) {
    int sink_idx = 0;
    std::uint64_t inspected = 0;
    std::uint64_t compatible = 0;
    std::uint64_t rejected = 0;
    FwdIndexCursor cursor(fwd_seeds, /*descending=*/is_rc);
    for (const QuerySeed& seed : selected) {
      const int my_idx = sink_idx++;
      const int fwd_pos = is_rc ? read_len - k - seed.read_pos : seed.read_pos;
      const std::size_t view_index =
          is_rc ? cursor.retreat_to(fwd_pos) : cursor.advance_to(fwd_pos);
      if (view_index >= views.size() ||
          fwd_seeds[view_index].read_pos != fwd_pos ||
          fwd_seeds[view_index].key != seed.key)
        continue;
      const KmerPostingView& v = views[view_index];
      if (!v.found() || v.count == 0)
        continue;
      if (!chain_window_seed_occ_allowed(ctx, v))
        continue;
      int chr_idx =
          chromosome_index_for_global_pos(chr_bounds, n_chr, v.positions[0]);
      if (chr_idx < 0 || chr_idx >= n_chr)
        continue;
      std::uint64_t chr_lo = chr_bounds[static_cast<std::size_t>(chr_idx)];
      std::uint64_t chr_hi = chr_bounds[static_cast<std::size_t>(chr_idx + 1)];
      // A seed's postings ascend, so its postings within one diagonal bin are
      // contiguous: a key change marks the first posting of this seed in that
      // bin (the distinct-seed support signal).
      std::int64_t prev_key = INT64_MIN;
      for (std::uint32_t p = 0; p < v.count; ++p) {
        const std::uint64_t g = static_cast<std::uint64_t>(v.positions[p]);
        while (chr_idx + 1 < n_chr && g >= chr_hi) {
          ++chr_idx;
          chr_lo = chr_hi;
          chr_hi = chr_bounds[static_cast<std::size_t>(chr_idx + 1)];
        }
        if (g < chr_lo || g >= chr_hi)
          continue;
        // Strand-compatibility test, after the contig cursor advance (which must see
        // every posting) and the range check. `prev_key` is not updated on a rejected
        // posting, so it cannot suppress the bin's first compatible support.
        ++inspected;
        if (!packed_ref_orientation_compatible(
                seed.z, v.positions.packed_at(p), is_rc)) {
          ++rejected;
          continue;
        }
        ++compatible;
        const int local = static_cast<int>(g - chr_lo);
        const int ref_start = local - seed.read_pos;
        const std::int64_t key =
            vote_pack_key(chr_idx, vote_floor_div(ref_start, W));
        const std::uint32_t bin = bins.add(key, key != prev_key);
        prev_key = key;
        log.push_back(CapturedVote{
            key,
            VoteHit{my_idx, seed.read_pos, ref_start, v.occurrence, v.count},
            bin});
      }
    }
    if (is_rc) {
      summary.rc_postings_inspected = inspected;
      summary.rc_postings_compatible = compatible;
      summary.rc_postings_rejected_strand = rejected;
    } else {
      summary.fwd_postings_inspected = inspected;
      summary.fwd_postings_compatible = compatible;
      summary.fwd_postings_rejected_strand = rejected;
    }
  };
  vote_strand(scratch.selected_fwd, /*is_rc=*/false, scratch.fwd_log,
              scratch.fwd_bins);
  vote_strand(scratch.selected_rc, /*is_rc=*/true, scratch.rc_log,
              scratch.rc_bins);
  return summary;
}

std::vector<CoarseDiagonalPeak>
aggregate_fused_capture_to_peaks(FusedCaptureScratch& scratch, int k,
                                 int min_support) {
  std::vector<CoarseDiagonalPeak> peaks;
  const auto drain = [&](const std::vector<CapturedVote>& log,
                         const DiagonalBinCounter& bins, bool reverse) {
    std::vector<CapturedVote>& survivors = scratch.survivors;
    survivors.clear();
    for (const CapturedVote& cv : log) {
      if (bins.support_at(cv.bin) >= static_cast<std::uint32_t>(min_support))
        survivors.push_back(cv);
    }
    std::stable_sort(survivors.begin(), survivors.end(),
                     [](const CapturedVote& a, const CapturedVote& b) {
                       return a.key < b.key;
                     });
    // Reserve for the surviving bins only: the key-sorted survivors' runs are exactly
    // the bins that cleared min_support, an upper bound on this strand's pushes.
    std::size_t surviving_bins = 0;
    for (std::size_t s = 0; s < survivors.size(); ++s) {
      if (s == 0 || survivors[s].key != survivors[s - 1].key)
        ++surviving_bins;
    }
    peaks.reserve(peaks.size() + surviving_bins);
    std::size_t i = 0;
    while (i < survivors.size()) {
      std::size_t j = i;
      while (j < survivors.size() && survivors[j].key == survivors[i].key)
        ++j;
      CoarseDiagonalPeak peak;
      const CapturedVote* base = survivors.data() + i;
      if (build_bin_peak(
              survivors[i].key, j - i,
              [base](std::size_t idx) -> const VoteHit& {
                return base[idx].hit;
              },
              k, reverse, min_support, scratch.read_pos_scratch,
              scratch.ref_start_scratch, scratch.seed_id_scratch, peak)) {
        peaks.push_back(peak);
      }
      i = j;
    }
  };
  drain(scratch.fwd_log, scratch.fwd_bins, /*reverse=*/false);
  drain(scratch.rc_log, scratch.rc_bins, /*reverse=*/true);
  return peaks;
}

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
