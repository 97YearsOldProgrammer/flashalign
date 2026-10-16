// The RNA vote: single-pass two-strand capture, diagonal peaks, the coarse locus envelopes that form the catalogue, and the second coarse pass with its whole-read re-vote.
#include "rna.h"

#include "../index/format.h" // chromosome_index_for_global_pos, packed_ref_contig/local
#include "../seed/seed.h" // VoteHit
#include "../vote/vote.h" // chain_window_seed_occ_allowed, vote_floor_div, vote_pack_key

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

/**** fused_capture ****/

namespace fa::cpu::lr::rna::placement {

namespace {

// Inline: called per posting from DiagonalBinCounter::add (MurmurHash3's fmix64).
inline std::uint64_t mix_bin_key(std::int64_t key) {
  std::uint64_t x = static_cast<std::uint64_t>(key);
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

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

namespace {

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

// Per-strand low-occurrence tile-spread selection, replayed from the cached evidence
// (limit_query_seeds semantics).
void select_strand(FusedCaptureScratch& scratch, int read_len, int max_seeds,
                   const std::vector<QuerySeed>& stream, bool reversed_views,
                   std::vector<QuerySeed>& selected) {
  auto& views = scratch.views;
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
}

// capture_votes_fused's arguments and the values its prologue fixes, which both strands'
// vote reads under the same names.
struct CaptureCall {
  const LongReadSeedContext& ctx;
  const std::vector<QuerySeed>& fwd_seeds;
  int read_len;
  FusedCaptureSummary& summary; // each lane's posting counts
  int k;
  const std::uint64_t* chr_bounds;
  int n_chr;
  int W;
  std::vector<KmerPostingView>& views;
  // The seeds re-admitted past the cap, by forward read position, or null.
  const std::uint8_t* readmitted;
};

// Vote per strand from the cached views; no further index lookups.
void vote_strand(const CaptureCall& call,
                 const std::vector<QuerySeed>& selected, bool is_rc,
                 std::vector<CapturedVote>& log, DiagonalBinCounter& bins) {
  auto& fwd_seeds = call.fwd_seeds;
  auto& read_len = call.read_len;
  auto& k = call.k;
  auto& views = call.views;
  auto& ctx = call.ctx;
  auto& chr_bounds = call.chr_bounds;
  auto& n_chr = call.n_chr;
  auto& W = call.W;
  auto& summary = call.summary;
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
    if (!chain_window_seed_occ_allowed(ctx, v) &&
        !(call.readmitted != nullptr && call.readmitted[fwd_pos] != 0))
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
}

} // namespace

FusedCaptureSummary capture_votes_fused(
    const LongReadSeedContext& ctx, const std::vector<QuerySeed>& fwd_seeds,
    const std::vector<QuerySeed>& rc_seeds, int read_len, int max_seeds,
    FusedCaptureScratch& scratch,
    int occ_dist, std::vector<std::uint8_t>* readmitted) {
  FusedCaptureSummary summary;
  scratch.fwd_log.clear();
  scratch.rc_log.clear();
  if (readmitted != nullptr)
    readmitted->assign(static_cast<std::size_t>(std::max(read_len, 0)), 0);
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
  if (readmitted != nullptr)
    mark_readmitted_seeds(fwd_seeds, views, read_len, k,
                          effective_long_primary_occ_cap(ctx.occ_policy),
                          occ_dist, *readmitted);

  select_strand(scratch, read_len, max_seeds, fwd_seeds,
                /*reversed_views=*/false, scratch.selected_fwd);
  select_strand(scratch, read_len, max_seeds, rc_seeds,
                /*reversed_views=*/true, scratch.selected_rc);
  summary.fwd_voted_seeds = scratch.selected_fwd.size();
  summary.rc_voted_seeds = scratch.selected_rc.size();

  const CaptureCall call{ctx, fwd_seeds, read_len,
                         summary, k, chr_bounds, n_chr, W, views,
                         readmitted != nullptr ? readmitted->data() : nullptr};
  vote_strand(call, scratch.selected_fwd, /*is_rc=*/false, scratch.fwd_log,
              scratch.fwd_bins);
  vote_strand(call, scratch.selected_rc, /*is_rc=*/true, scratch.rc_log,
              scratch.rc_bins);
  return summary;
}

} // namespace fa::cpu::lr::rna::placement

/**** peaks ****/
// The capture log folded into coarse diagonal peaks: one per diagonal bin whose
// distinct-seed support clears min_support, with its seeds' query positions.

namespace fa::cpu::lr::rna::placement {

std::uint32_t union_seed_coverage(const std::vector<int>& sorted_read_pos,
                                  int k) {
  std::uint32_t covered = 0;
  int cur_lo = 0;
  int cur_hi = 0;
  bool open = false;
  for (int p : sorted_read_pos) {
    const int lo = p;
    const int hi = p + k;
    if (!open) {
      cur_lo = lo;
      cur_hi = hi;
      open = true;
    } else if (lo <= cur_hi) {
      cur_hi = std::max(cur_hi, hi);
    } else {
      covered += static_cast<std::uint32_t>(cur_hi - cur_lo);
      cur_lo = lo;
      cur_hi = hi;
    }
  }
  if (open)
    covered += static_cast<std::uint32_t>(cur_hi - cur_lo);
  return covered;
}

// Aggregates one diagonal bin's captured postings into a CoarseDiagonalPeak, templated
// over how a hit is fetched.
// Returns false when the bin is empty, carries a negative reference id, or
// fails the distinct-seed min_support gate. Scratch vectors are caller-owned
// and reused; their prior contents are discarded.
template <class HitAt>
static bool
build_bin_peak(std::int64_t key, std::size_t hit_count, HitAt&& hit_at, int k,
               bool reverse, int min_support, std::vector<int>& read_pos,
               std::vector<int>& ref_starts, std::vector<int>& seed_ids,
               CoarseDiagonalPeak& out) {
  if (hit_count == 0)
    return false;
  const int chr = static_cast<int>(key >> 32);
  if (chr < 0)
    return false;

  read_pos.clear();
  ref_starts.clear();
  seed_ids.clear();
  seed_ids.reserve(hit_count);
  std::int64_t ref_lo = INT64_MAX;
  std::int64_t ref_hi = INT64_MIN;
  int q_lo = INT32_MAX;
  int q_hi = INT32_MIN;
  for (std::size_t i = 0; i < hit_count; ++i) {
    const VoteHit& h = hit_at(i);
    read_pos.push_back(h.read_pos);
    ref_starts.push_back(h.ref_start);
    const std::int64_t genomic = static_cast<std::int64_t>(h.ref_start) +
                                 static_cast<std::int64_t>(h.read_pos);
    ref_lo = std::min(ref_lo, genomic);
    ref_hi = std::max(ref_hi, genomic + k);
    q_lo = std::min(q_lo, h.read_pos);
    q_hi = std::max(q_hi, h.read_pos + k);
    seed_ids.push_back(h.seed_idx);
  }
  std::sort(seed_ids.begin(), seed_ids.end());
  int distinct = 0;
  int last_seed = INT32_MIN;
  for (int s : seed_ids) {
    if (s != last_seed) {
      ++distinct;
      last_seed = s;
    }
  }
  if (distinct < min_support)
    return false;
  if (ref_lo < 0)
    ref_lo = 0;

  std::sort(read_pos.begin(), read_pos.end());
  std::nth_element(ref_starts.begin(),
                   ref_starts.begin() +
                       static_cast<long>(ref_starts.size() / 2),
                   ref_starts.end());
  const std::int64_t intercept = ref_starts[ref_starts.size() / 2];

  out.reference_id = chr;
  out.reverse = reverse;
  out.oriented_query_begin = static_cast<std::uint32_t>(std::max(0, q_lo));
  out.oriented_query_end = static_cast<std::uint32_t>(std::max(q_lo, q_hi));
  out.reference_begin = static_cast<std::uint64_t>(ref_lo);
  out.reference_end = static_cast<std::uint64_t>(std::max(ref_lo, ref_hi));
  out.median_diagonal = intercept;
  out.query_covered_bases = union_seed_coverage(read_pos, k);
  out.distinct_seed_support = static_cast<std::uint32_t>(distinct);
  return true;
}

// Drains one strand's log into peaks: one per bin that cleared min_support, in key order.
static void drain(FusedCaptureScratch& scratch, int k, int min_support,
                  std::vector<CoarseDiagonalPeak>& peaks,
                  const std::vector<CapturedVote>& log,
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
      // Sorted by build_bin_peak; one position is one seed.
      const std::vector<int>& positions = scratch.read_pos_scratch;
      for (std::size_t q = 0; q < positions.size(); ++q)
        if (q == 0 || positions[q] != positions[q - 1])
          scratch.peak_read_pos.push_back(positions[q]);
      scratch.peak_read_pos_offset.push_back(
          static_cast<std::uint32_t>(scratch.peak_read_pos.size()));
    }
    i = j;
  }
}

std::vector<CoarseDiagonalPeak>
aggregate_fused_capture_to_peaks(FusedCaptureScratch& scratch, int k,
                                 int min_support) {
  std::vector<CoarseDiagonalPeak> peaks;
  scratch.peak_read_pos.clear();
  scratch.peak_read_pos_offset.assign(1, 0);
  drain(scratch, k, min_support, peaks, scratch.fwd_log, scratch.fwd_bins,
        /*reverse=*/false);
  drain(scratch, k, min_support, peaks, scratch.rc_log, scratch.rc_bins,
        /*reverse=*/true);
  return peaks;
}

} // namespace fa::cpu::lr::rna::placement

/**** coarse_chain ****/
// Groups RNA diagonal peaks into candidate locus envelopes. An RNA read spans several
// exons, each on its own diagonal, so peaks are grouped rather than collapsed to one
// placement. This stage reads no reference bases. Envelopes are ranked by their seed
// evidence; the colinear score is computed for the kept loci only and is the coarse
// confidence the RNA MAPQ reads.

namespace fa::cpu::lr::rna::placement {

// Derives the internal score units from the seed span; only the intron bounds come from
// the user.
CoarseLocusOptions make_coarse_locus_options(uint32_t seed_span,
                                             uint32_t min_intron,
                                             uint32_t max_intron) {
  CoarseLocusOptions o;
  const uint32_t s = seed_span == 0 ? 15u : seed_span;
  o.diagonal_band = std::max<uint32_t>(16, s + s / 2); // ~1.5 seeds of drift
  o.max_query_gap = std::max<uint32_t>(64, s * 8);
  o.max_query_overlap = std::max<uint32_t>(16, s * 2);
  o.coalesce_min_piece_span = std::max<uint32_t>(16, 2 * s);
  o.min_intron = min_intron == 0 ? 20u : min_intron;
  o.max_intron = max_intron == 0 ? 200000u : max_intron;
  // Coalescing reach and ceiling in units of the maximum intron (see CoarseLocusOptions);
  // saturating, so an extreme -G cannot wrap them.
  const uint64_t mi = o.max_intron;
  o.coalesce_gap =
      static_cast<uint32_t>(std::min<uint64_t>(4u * mi, 0xffffffffull));
  o.max_coalesced_span =
      static_cast<uint32_t>(std::min<uint64_t>(16u * mi, 0xffffffffull));
  o.splice_open = static_cast<int32_t>(std::max<uint32_t>(4, s / 2));
  o.splice_log_scale = 2;
  o.query_gap_scale = 1;
  o.diagonal_drift_scale = 1;
  o.overlap_scale = 1;
  return o;
}

// Peak-to-peak transition scoring for the coarse colinear recurrence
// (best_colinear_score).
namespace detail {

static uint32_t ceil_log2_u64(uint64_t value) {
  uint32_t r = 0;
  uint64_t v = 1;
  while (v < value) {
    v <<= 1;
    ++r;
  }
  return r;
}

enum class TransitionKind : uint8_t { kSameExon, kSplice, kInvalid };

struct Transition {
  TransitionKind kind = TransitionKind::kInvalid;
  int64_t penalty = 0;
};

// Asymmetric transition score left -> right. The caller guarantees same contig and strand
// and reference-monotone order, so this scores geometry only.
static Transition score_transition(const CoarseDiagonalPeak& left,
                                   const CoarseDiagonalPeak& right,
                                   const CoarseLocusOptions& o) {
  Transition t;
  if (right.reference_begin <= left.reference_begin)
    return t; // non-monotone
  const int64_t dq = static_cast<int64_t>(right.oriented_query_begin) -
                     static_cast<int64_t>(left.oriented_query_end);
  const int64_t dr = static_cast<int64_t>(right.reference_begin) -
                     static_cast<int64_t>(left.reference_end);
  if (dr <= 0)
    return t;
  if (dq < -static_cast<int64_t>(o.max_query_overlap) ||
      dq > static_cast<int64_t>(o.max_query_gap)) {
    return t;
  }
  const int64_t overlap = std::max<int64_t>(0, -dq);
  const int64_t query_gap = std::max<int64_t>(0, dq);
  const int64_t diagonal_delta = std::llabs(dr - dq);
  if (diagonal_delta <= static_cast<int64_t>(o.diagonal_band)) {
    t.kind = TransitionKind::kSameExon;
    t.penalty = static_cast<int64_t>(o.diagonal_drift_scale) * diagonal_delta +
                static_cast<int64_t>(o.overlap_scale) * overlap;
    return t;
  }
  const int64_t intron = dr - dq;
  if (intron < static_cast<int64_t>(o.min_intron) ||
      intron > static_cast<int64_t>(o.max_intron)) {
    return t; // kInvalid: an over-long / impossible jump, never free-chained
  }
  // Sub-linear in intron length: a zero-cost ref jump over-chains noise, a
  // linear penalty erases real long introns.
  t.kind = TransitionKind::kSplice;
  t.penalty = static_cast<int64_t>(o.splice_open) +
              static_cast<int64_t>(o.splice_log_scale) *
                  static_cast<int64_t>(
                      ceil_log2_u64(static_cast<uint64_t>(intron) + 1)) +
              static_cast<int64_t>(o.query_gap_scale) * query_gap +
              static_cast<int64_t>(o.overlap_scale) * overlap;
  return t;
}

static int64_t incremental_query_reward(const CoarseDiagonalPeak& left,
                                        const CoarseDiagonalPeak& right) {
  const uint32_t start =
      std::max(left.oriented_query_end, right.oriented_query_begin);
  return right.oriented_query_end > start
             ? static_cast<int64_t>(right.oriented_query_end - start)
             : 0;
}

} // namespace detail

// A locus's rank in an exact tie: minimap2's read-seeded hash (tie_locus_hash) of its
// contig slot, strand and reference start, lower first. The rival lifecycle breaks ties
// by catalogue index, so this hash decides every RNA tie between windows equal on every
// vote key. The start enters by its low 32 bits, as the DNA bin does.
static std::uint64_t coarse_locus_tie_hash(std::uint32_t tie_seed,
                                           const CoarseLocus& locus) noexcept {
  return tie_locus_hash(tie_seed, locus.reference_id, locus.reverse,
                        static_cast<std::int64_t>(locus.reference_begin));
}

namespace {

// Union of the peaks' oriented query intervals, in bases. Only a tie-break: as a primary
// key it would reward a paralogous envelope for accumulating scattered peaks.
int64_t
union_query_coverage(const std::vector<CoarseDiagonalPeak>& peaks,
                     const std::vector<uint32_t>& indices,
                     std::vector<std::pair<uint32_t, uint32_t>>& spans) {
  spans.clear();
  spans.reserve(indices.size());
  for (uint32_t i : indices) {
    const CoarseDiagonalPeak& p = peaks[i];
    if (p.oriented_query_end > p.oriented_query_begin)
      spans.push_back({p.oriented_query_begin, p.oriented_query_end});
  }
  if (spans.empty())
    return 0;
  std::sort(spans.begin(), spans.end());
  int64_t total = 0;
  uint32_t lo = spans.front().first;
  uint32_t hi = spans.front().second;
  for (size_t i = 1; i < spans.size(); ++i) {
    if (spans[i].first > hi) {
      total += static_cast<int64_t>(hi) - lo;
      lo = spans[i].first;
      hi = spans[i].second;
    } else if (spans[i].second > hi) {
      hi = spans[i].second;
    }
  }
  return total + (static_cast<int64_t>(hi) - lo);
}

// How much of the read the envelope explains with seed matches: peaks in query order,
// each credited the smaller of its own seed coverage and the query it newly explains, so
// no query base counts twice. Unlike best_colinear_score's span increments, a sparse
// peak straddling a wide query range earns only the bases it matched. `indices` may be
// in any order; `scratch` is caller-owned.
int64_t seed_evidence_score(const std::vector<CoarseDiagonalPeak>& peaks,
                            const std::vector<uint32_t>& indices,
                            std::vector<uint32_t>& scratch) {
  scratch.assign(indices.begin(), indices.end());
  std::stable_sort(scratch.begin(), scratch.end(), [&](uint32_t a, uint32_t b) {
    const CoarseDiagonalPeak& pa = peaks[a];
    const CoarseDiagonalPeak& pb = peaks[b];
    if (pa.oriented_query_begin != pb.oriented_query_begin)
      return pa.oriented_query_begin < pb.oriented_query_begin;
    return pa.oriented_query_end < pb.oriented_query_end;
  });
  int64_t total = 0;
  int64_t reach = -1;
  for (uint32_t i : scratch) {
    const CoarseDiagonalPeak& p = peaks[i];
    const int64_t begin =
        std::max<int64_t>(reach, static_cast<int64_t>(p.oriented_query_begin));
    const int64_t novel = static_cast<int64_t>(p.oriented_query_end) - begin;
    if (novel <= 0)
      continue;
    total += std::min<int64_t>(p.query_covered_bases, novel);
    reach =
        std::max<int64_t>(reach, static_cast<int64_t>(p.oriented_query_end));
  }
  return total;
}

// Best colinear chain score over one complete envelope's peaks. It ranks envelopes
// (colinear evidence separates the true locus from a paralog of scattered peaks) but
// never bounds them. `indices` must be sorted by reference_begin.
int64_t best_colinear_score(const std::vector<CoarseDiagonalPeak>& peaks,
                            const std::vector<uint32_t>& indices,
                            const CoarseLocusOptions& o,
                            std::vector<int64_t>& score) {
  const size_t m = indices.size();
  if (m == 0)
    return 0;
  score.assign(m, 0);
  int64_t best = 0;
  for (size_t i = 0; i < m; ++i) {
    const CoarseDiagonalPeak& pi = peaks[indices[i]];
    score[i] = static_cast<int64_t>(pi.query_covered_bases);
    uint32_t scanned = 0;
    for (size_t j = i; j-- > 0 && scanned < o.max_predecessors;) {
      ++scanned;
      const CoarseDiagonalPeak& pj = peaks[indices[j]];
      const detail::Transition t = detail::score_transition(pj, pi, o);
      if (t.kind == detail::TransitionKind::kInvalid)
        continue;
      const int64_t candidate =
          score[j] + detail::incremental_query_reward(pj, pi) - t.penalty;
      if (candidate > score[i])
        score[i] = candidate;
    }
    if (score[i] > best)
      best = score[i];
  }
  return best;
}

// Reference-ordered count of transitions that look like an intron; context only, the
// realized junction count is the splice controller's.
uint32_t estimate_intron_count(const std::vector<CoarseDiagonalPeak>& peaks,
                               const std::vector<uint32_t>& indices,
                               const CoarseLocusOptions& o,
                               std::vector<uint32_t>& by_reference) {
  by_reference.assign(indices.begin(), indices.end());
  std::stable_sort(by_reference.begin(), by_reference.end(),
                   [&](uint32_t a, uint32_t b) {
                     return peaks[a].reference_begin < peaks[b].reference_begin;
                   });
  uint32_t introns = 0;
  for (size_t i = 1; i < by_reference.size(); ++i) {
    const CoarseDiagonalPeak& left = peaks[by_reference[i - 1]];
    const CoarseDiagonalPeak& right = peaks[by_reference[i]];
    if (right.reference_begin <= left.reference_end)
      continue;
    const int64_t dr = static_cast<int64_t>(right.reference_begin) -
                       static_cast<int64_t>(left.reference_end);
    const int64_t dq = static_cast<int64_t>(right.oriented_query_begin) -
                       static_cast<int64_t>(left.oriented_query_end);
    const int64_t intron = dr - dq;
    if (intron >= static_cast<int64_t>(o.min_intron) &&
        intron <= static_cast<int64_t>(o.max_intron)) {
      ++introns;
    }
  }
  return introns;
}

// Sets the forward-frame bases of [pos, pos + k) over the seeds of `peak_indices` that
// hold `from` in `cover` to `to`; returns how many it set.
int recolor_query_cover(const std::vector<CoarseDiagonalPeak>& peaks,
                        const std::vector<uint32_t>& peak_indices,
                        const CoarseLocusQueryCover& qc, uint8_t from,
                        uint8_t to, std::vector<uint8_t>& cover) {
  int changed = 0;
  for (uint32_t p : peak_indices) {
    const bool reverse = peaks[p].reverse;
    for (uint32_t s = (*qc.offset)[p]; s < (*qc.offset)[p + 1]; ++s) {
      const int pos = (*qc.read_pos)[s];
      const int begin = reverse ? qc.read_len - pos - qc.k : pos;
      const int hi = std::min(qc.read_len, begin + qc.k);
      for (int b = std::max(0, begin); b < hi; ++b)
        if (cover[static_cast<size_t>(b)] == from) {
          cover[static_cast<size_t>(b)] = to;
          ++changed;
        }
    }
  }
  return changed;
}

// Sorts the peaks by contig, strand and reference position and sweeps them into
// envelopes, one per run of same-contig, same-strand peaks within one legal intron.
void sweep_envelopes(const std::vector<CoarseDiagonalPeak>& peaks,
                     const CoarseLocusOptions& o, CoarseLocusScratch& scratch,
                     size_t n) {
  std::vector<uint32_t>& order = scratch.order;
  order.assign(n, 0);
  for (size_t i = 0; i < n; ++i)
    order[i] = static_cast<uint32_t>(i);
  std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
    const CoarseDiagonalPeak& pa = peaks[a];
    const CoarseDiagonalPeak& pb = peaks[b];
    if (pa.reference_id != pb.reference_id)
      return pa.reference_id < pb.reference_id;
    if (pa.reverse != pb.reverse)
      return pa.reverse < pb.reverse;
    if (pa.reference_begin != pb.reference_begin)
      return pa.reference_begin < pb.reference_begin;
    if (pa.reference_end != pb.reference_end)
      return pa.reference_end < pb.reference_end;
    return a < b; // original index: total order
  });

  // One linear sweep: a new envelope starts whenever the strand/contig changes
  // or the reference gap exceeds one legal intron.
  std::vector<CoarseLocus>& envelopes = scratch.envelopes;
  std::vector<uint64_t>& support = scratch.support;
  std::vector<int64_t>& coverage = scratch.coverage;
  std::vector<uint32_t>& query_order = scratch.query_order;
  envelopes.clear();
  support.clear();
  coverage.clear();
  size_t i = 0;
  while (i < n) {
    const CoarseDiagonalPeak& first = peaks[order[i]];
    if (first.reference_id < 0) {
      ++i;
      continue;
    }
    CoarseLocus locus;
    locus.reference_id = first.reference_id;
    locus.reverse = first.reverse;
    locus.reference_begin = first.reference_begin;
    locus.reference_end = first.reference_end;
    locus.oriented_query_begin = first.oriented_query_begin;
    locus.oriented_query_end = first.oriented_query_end;
    locus.peak_indices.push_back(order[i]);
    uint64_t locus_support = first.distinct_seed_support;
    size_t j = i + 1;
    for (; j < n; ++j) {
      const CoarseDiagonalPeak& p = peaks[order[j]];
      if (p.reference_id != locus.reference_id || p.reverse != locus.reverse) {
        break;
      }
      if (p.reference_begin >
          locus.reference_end + static_cast<uint64_t>(o.max_intron)) {
        break;
      }
      locus.reference_end = std::max(locus.reference_end, p.reference_end);
      locus.oriented_query_begin =
          std::min(locus.oriented_query_begin, p.oriented_query_begin);
      locus.oriented_query_end =
          std::max(locus.oriented_query_end, p.oriented_query_end);
      locus.peak_indices.push_back(order[j]);
      locus_support += p.distinct_seed_support;
    }
    // Rank by seed evidence; coverage and support are tie-breaks. `score` and
    // `intron_count` are not ranking terms and are computed for the kept loci only.
    locus.rank_score =
        seed_evidence_score(peaks, locus.peak_indices, query_order);
    locus.union_query_coverage_bases =
        union_query_coverage(peaks, locus.peak_indices, scratch.spans);
    coverage.push_back(locus.union_query_coverage_bases);
    envelopes.push_back(std::move(locus));
    support.push_back(locus_support);
    i = j;
  }
}

// Coalesces envelopes that are pieces of one transcript. The sweep splits
// an envelope when consecutive small exons each fall below min_support and the next
// peak lies beyond max_intron; the trailing piece would be spent as a rival and its
// junctions fall outside the harvest window. A pair merges only if it is colinear in
// the read (the later piece starts after the earlier one ends in oriented query
// coordinates); a paralog explains the same read interval and never merges. Each
// envelope is offered to every surviving envelope of its group, so a spurious upstream
// envelope neither blocks nor captures real pieces; the qualifying target reaching
// furthest right wins.
size_t coalesce_envelopes(const std::vector<CoarseDiagonalPeak>& peaks,
                          const CoarseLocusOptions& o,
                          CoarseLocusScratch& scratch) {
  auto& envelopes = scratch.envelopes;
  auto& support = scratch.support;
  auto& coverage = scratch.coverage;
  auto& query_order = scratch.query_order;
  size_t coalesced_count = 0;
  if (o.coalesce_gap != 0 && envelopes.size() > 1) {
    std::vector<uint8_t>& absorbed = scratch.coalesce_absorbed;
    absorbed.assign(envelopes.size(), 0); // every slot starts clean
    size_t write = 0;                     // next surviving slot
    size_t group_begin = 0; // first surviving slot of the current group
    // A target's evidence terms are re-derived once, when its group closes: both depend
    // on the whole peak set.
    auto close_group = [&]() {
      for (size_t t = group_begin; t < write; ++t) {
        if (!absorbed[t])
          continue;
        CoarseLocus& m = envelopes[t];
        m.rank_score = seed_evidence_score(peaks, m.peak_indices, query_order);
        m.union_query_coverage_bases =
            union_query_coverage(peaks, m.peak_indices, scratch.spans);
        coverage[t] = m.union_query_coverage_bases;
      }
    };
    for (size_t e = 0; e < envelopes.size(); ++e) {
      const CoarseLocus& c = envelopes[e];
      if (write > group_begin &&
          (envelopes[group_begin].reference_id != c.reference_id ||
           envelopes[group_begin].reverse != c.reverse)) {
        close_group();
        group_begin = write;
      }
      // The best earlier target: of this group's surviving envelopes that could
      // absorb `c`, the one whose reference_end is largest.
      size_t best = write; // == write means "none qualified"
      for (size_t t = group_begin; t < write; ++t) {
        const CoarseLocus& m = envelopes[t];
        // A group is strictly ordered and disjoint after the sweep; the test also keeps
        // the unsigned subtraction safe.
        if (c.reference_begin <= m.reference_end ||
            c.reference_begin - m.reference_end >
                static_cast<uint64_t>(o.coalesce_gap)) {
          continue;
        }
        // Both pieces must be exon-shaped (see coalesce_min_piece_span).
        if (m.oriented_query_end - m.oriented_query_begin <
                o.coalesce_min_piece_span ||
            c.oriented_query_end - c.oriented_query_begin <
                o.coalesce_min_piece_span) {
          continue;
        }
        if (static_cast<int64_t>(m.oriented_query_end) >
            static_cast<int64_t>(c.oriented_query_begin) +
                static_cast<int64_t>(o.max_query_overlap)) {
          continue; // not colinear in the read
        }
        const uint64_t span =
            std::max(m.reference_end, c.reference_end) - m.reference_begin;
        if (span > static_cast<uint64_t>(o.max_coalesced_span))
          continue;
        // Evidence gate: two scattered repeat peaks of a short junk read can land
        // colinear by chance, so one credible piece is required.
        if (std::max(support[t], support[e]) <
            static_cast<uint64_t>(o.coalesce_min_support)) {
          continue;
        }
        if (best == write || envelopes[best].reference_end < m.reference_end)
          best = t;
      }
      if (best != write) {
        CoarseLocus& m = envelopes[best];
        m.reference_end = std::max(m.reference_end, c.reference_end);
        m.oriented_query_begin =
            std::min(m.oriented_query_begin, c.oriented_query_begin);
        m.oriented_query_end =
            std::max(m.oriented_query_end, c.oriented_query_end);
        // Every peak of `c` begins past the target's reference end, so appending keeps
        // peak_indices reference-ordered, as best_colinear_score requires.
        m.peak_indices.insert(m.peak_indices.end(), c.peak_indices.begin(),
                              c.peak_indices.end());
        support[best] += support[e];
        absorbed[best] = 1;
        ++coalesced_count;
        continue;
      }
      // `e` survives as a target for what follows; compact it over the absorbed slots
      // together with its parallel entries.
      if (write != e) {
        envelopes[write] = std::move(envelopes[e]);
        support[write] = support[e];
        coverage[write] = coverage[e];
      }
      ++write;
    }
    close_group();
    envelopes.resize(write);
    support.resize(write);
    coverage.resize(write);
  }
  return coalesced_count;
}

// Ranks the envelopes by seed evidence, then coverage, support and peak count, then by
// minimap2's read-seeded locus hash.
void rank_envelopes(CoarseLocusScratch& scratch, uint32_t tie_seed) {
  auto& envelopes = scratch.envelopes;
  auto& coverage = scratch.coverage;
  auto& support = scratch.support;
  std::vector<uint32_t>& ranked = scratch.ranked;
  ranked.assign(envelopes.size(), 0);
  for (size_t e = 0; e < envelopes.size(); ++e)
    ranked[e] = static_cast<uint32_t>(e);
  std::stable_sort(ranked.begin(), ranked.end(), [&](uint32_t a, uint32_t b) {
    if (envelopes[a].rank_score != envelopes[b].rank_score)
      return envelopes[a].rank_score > envelopes[b].rank_score;
    if (coverage[a] != coverage[b])
      return coverage[a] > coverage[b];
    if (support[a] != support[b])
      return support[a] > support[b];
    if (envelopes[a].peak_indices.size() != envelopes[b].peak_indices.size())
      return envelopes[a].peak_indices.size() >
             envelopes[b].peak_indices.size();
    // Equal on every evidence key: minimap2's read-seeded locus hash decides rather
    // than contig order; contig slot and start settle only a 64-bit collision.
    const uint64_t a_hash = coarse_locus_tie_hash(tie_seed, envelopes[a]);
    const uint64_t b_hash = coarse_locus_tie_hash(tie_seed, envelopes[b]);
    if (a_hash != b_hash)
      return a_hash < b_hash;
    if (envelopes[a].reference_id != envelopes[b].reference_id)
      return envelopes[a].reference_id < envelopes[b].reference_id;
    return envelopes[a].reference_begin < envelopes[b].reference_begin;
  });
}

// Takes the ranked envelopes into `out` in rank order, one per distinct reference window,
// up to the cut, then into `overflow_out` or by query cover; finishes each kept locus.
void admit_envelopes(const std::vector<CoarseDiagonalPeak>& peaks,
                     const CoarseLocusOptions& o, CoarseLocusScratch& scratch,
                     std::vector<CoarseLocus>* overflow_out,
                     const CoarseLocusQueryCover* query_cover,
                     std::vector<CoarseLocus>& out) {
  auto& envelopes = scratch.envelopes;
  auto& ranked = scratch.ranked;
  const size_t keep = std::min<size_t>(o.max_locus_chains, ranked.size());
  // With overflow requested, up to `admit` loci are taken and everything past `keep`
  // goes to `overflow_out` unfinished.
  const size_t admit =
      overflow_out ? std::min<size_t>(keep + o.max_overflow_loci, ranked.size())
                   : keep;
  out.reserve(keep);
  // Past the cut, query_cover's admissions walk the rest of the ranking.
  const bool cover_admits = query_cover != nullptr && overflow_out == nullptr;
  std::vector<uint8_t>& cover = scratch.query_cover;
  bool cover_built = false;
  // Walk in rank order, keeping at most one envelope per distinct reference window: a
  // candidate overlapping a kept locus on the same contig, on either strand, names the
  // same place and is dropped, and the next ranked entry takes the slot. Ignoring strand
  // removes the winner's own opposite-strand mirror; same-strand envelopes are disjoint
  // by construction.
  for (size_t r = 0;
       r < ranked.size() &&
       (cover_admits ||
        out.size() + (overflow_out ? overflow_out->size() : 0) < admit);
       ++r) {
    const CoarseLocus& candidate = envelopes[ranked[r]];
    bool duplicate_window = false;
    for (const CoarseLocus& kept : out) {
      if (kept.reference_id != candidate.reference_id)
        continue;
      // Half-open intersection, strand ignored: the mirror is the same window.
      if (candidate.reference_begin < kept.reference_end &&
          kept.reference_begin < candidate.reference_end) {
        duplicate_window = true;
        break;
      }
    }
    // Overflow loci must be distinct places too.
    if (!duplicate_window && overflow_out) {
      for (const CoarseLocus& kept : *overflow_out) {
        if (kept.reference_id != candidate.reference_id)
          continue;
        if (candidate.reference_begin < kept.reference_end &&
            kept.reference_begin < candidate.reference_end) {
          duplicate_window = true;
          break;
        }
      }
    }
    if (duplicate_window)
      continue;
    if (cover_admits && out.size() >= keep) {
      if (!cover_built) {
        cover.assign(static_cast<size_t>(std::max(0, query_cover->read_len)), 0);
        for (const CoarseLocus& kept : out)
          recolor_query_cover(peaks, kept.peak_indices, *query_cover, 0, 1, cover);
        cover_built = true;
      }
      const std::vector<uint32_t>& candidate_peaks =
          envelopes[ranked[r]].peak_indices;
      const bool admitted =
          recolor_query_cover(peaks, candidate_peaks, *query_cover, 0, 2,
                              cover) >= query_cover->min_new_bases;
      recolor_query_cover(peaks, candidate_peaks, *query_cover, 2,
                          admitted ? 1 : 0, cover);
      if (!admitted)
        continue;
    }
    CoarseLocus locus = std::move(envelopes[ranked[r]]);
    if (out.size() >= keep && !cover_admits) {
      overflow_out->push_back(std::move(locus));
      continue;
    }
    // Non-ranking statistics, for the kept loci only (best_colinear_score is an
    // O(m * max_predecessors) DP). Both need the sweep's reference-ordered peak_indices,
    // so they run before the query-order sort. `score` is the coarse confidence the RNA
    // MAPQ reads.
    locus.score =
        best_colinear_score(peaks, locus.peak_indices, o, scratch.chain_score);
    locus.intron_count = estimate_intron_count(peaks, locus.peak_indices, o,
                                               scratch.by_reference);
    // Query order is what the skeleton and the downstream anchor priority read.
    std::stable_sort(locus.peak_indices.begin(), locus.peak_indices.end(),
                     [&](uint32_t a, uint32_t b) {
                       const CoarseDiagonalPeak& pa = peaks[a];
                       const CoarseDiagonalPeak& pb = peaks[b];
                       if (pa.oriented_query_begin != pb.oriented_query_begin)
                         return pa.oriented_query_begin <
                                pb.oriented_query_begin;
                       if (pa.reference_begin != pb.reference_begin)
                         return pa.reference_begin < pb.reference_begin;
                       return a < b;
                     });
    out.push_back(std::move(locus));
  }
}

} // namespace

std::vector<CoarseLocus>
select_coarse_loci(const std::vector<CoarseDiagonalPeak>& peaks,
                   const CoarseLocusOptions& o, CoarseLocusScratch& scratch,
                   uint32_t tie_seed, size_t* coalesced,
                   std::vector<CoarseLocus>* overflow_out,
                   const CoarseLocusQueryCover* query_cover) {
  std::vector<CoarseLocus> out;
  if (overflow_out)
    overflow_out->clear();
  // Zeroed up front so the early returns report no coalescing.
  if (coalesced)
    *coalesced = 0;
  const size_t n = peaks.size();
  if (n == 0)
    return out;

  sweep_envelopes(peaks, o, scratch, n);
  auto& envelopes = scratch.envelopes;
  if (envelopes.empty())
    return out;

  size_t coalesced_count = coalesce_envelopes(peaks, o, scratch);
  if (coalesced)
    *coalesced = coalesced_count;

  rank_envelopes(scratch, tie_seed);

  admit_envelopes(peaks, o, scratch, overflow_out, query_cover, out);
  return out;
}

} // namespace fa::cpu::lr::rna::placement

/**** second_pass ****/

namespace fa::cpu::lr::rna {

// Re-votes every extracted seed whose span lies inside one of `intervals`, selected or
// not, with the fused capture's admission (a held, non-empty posting list within the
// occurrence policy) and bin key, on both strand frames. `views` is the capture's
// posting view per forward seed, so nothing probes the index. The posting budget is
// checked before any posting is walked; over it the read is declined. Returns the peaks
// drained at `min_support`; `capture` is cleared and used as scratch. A seed marked in
// `readmitted` (by forward read position, or null) passes the occurrence policy.
static std::vector<placement::CoarseDiagonalPeak>
rna_revote_intervals(const LongReadSeedContext& seed_ctx,
                     const std::vector<QuerySeed>& forward_seeds,
                     const std::vector<KmerPostingView>& views,
                     const std::vector<RnaQuerySpan>& intervals, int read_len,
                     int k, int min_support, int n_chr,
                     const std::uint8_t* readmitted,
                     placement::FusedCaptureScratch& capture,
                     RnaRevoteStats& stats) {
  std::vector<placement::CoarseDiagonalPeak> peaks;
  const auto inside = [&](int forward_position) {
    for (const RnaQuerySpan& interval : intervals)
      if (forward_position >= interval.begin &&
          forward_position + k <= interval.end)
        return true;
    return false;
  };
  const auto occ_allowed = [&](std::size_t i) {
    return chain_window_seed_occ_allowed(seed_ctx, views[i]) ||
           (readmitted != nullptr &&
            readmitted[forward_seeds[i].read_pos] != 0);
  };
  const std::size_t seed_count = std::min(forward_seeds.size(), views.size());
  for (std::size_t i = 0; i < seed_count; ++i) {
    if (!inside(forward_seeds[i].read_pos))
      continue;
    const KmerPostingView& view = views[i];
    if (!view.found() || view.count == 0 || !occ_allowed(i))
      continue;
    stats.postings += static_cast<std::int64_t>(view.count);
    ++stats.seeds_revoted;
  }
  if (stats.postings > kRnaExplainMaxPostings) {
    stats.over_budget = true;
    stats.postings = 0;
    stats.seeds_revoted = 0;
    return peaks;
  }
  if (stats.seeds_revoted == 0)
    return peaks;
  capture.fwd_log.clear();
  capture.rc_log.clear();
  capture.fwd_bins.reset(static_cast<std::size_t>(stats.seeds_revoted));
  capture.rc_bins.reset(static_cast<std::size_t>(stats.seeds_revoted));
  const int W = effective_vote_diag_bin_width(seed_ctx, read_len);
  // Strand-compatible re-vote, the capture's rule (vote_strand): each posting is
  // read once and emits a forward and a reverse record, each only on the lane the
  // posting is compatible with.
  for (std::size_t i = 0; i < seed_count; ++i) {
    const int forward_position = forward_seeds[i].read_pos;
    if (!inside(forward_position))
      continue;
    const KmerPostingView& view = views[i];
    if (!view.found() || view.count == 0 || !occ_allowed(i))
      continue;
    const int reverse_position = read_len - k - forward_position;
    const int seed_index = static_cast<int>(i);
    // A seed's postings ascend, so its postings within one diagonal bin are
    // contiguous: a key change marks the first posting of this seed in that
    // bin (the distinct-seed support signal), as in the fused capture.
    std::int64_t previous_forward = INT64_MIN;
    std::int64_t previous_reverse = INT64_MIN;
    for (std::uint32_t posting = 0; posting < view.count; ++posting) {
      const PackedRefPos packed = view.positions.packed_at(posting);
      const int contig = static_cast<int>(packed_ref_contig(packed));
      if (contig < 0 || contig >= n_chr)
        continue;
      const std::uint32_t local = packed_ref_local(packed);
      if (local > static_cast<std::uint32_t>(INT_MAX))
        continue;
      // One posting, two candidate records; exactly one lane is compatible. The
      // per-lane "first posting of this seed in this bin" signals update only on their
      // own lane's emission, so a declined record cannot suppress the first compatible
      // support.
      const bool emit_forward = packed_ref_orientation_compatible(
          forward_seeds[i].z, packed, /*reverse_lane=*/false);
      const bool emit_reverse = packed_ref_orientation_compatible(
          forward_seeds[i].z, packed, /*reverse_lane=*/true);
      ++stats.postings_inspected;
      if (!emit_forward)
        ++stats.rejected_forward;
      if (!emit_reverse)
        ++stats.rejected_reverse;
      if (emit_forward) {
        const int forward_start = static_cast<int>(local) - forward_position;
        const std::int64_t forward_key =
            vote_pack_key(contig, vote_floor_div(forward_start, W));
        const std::uint32_t forward_bin =
            capture.fwd_bins.add(forward_key, forward_key != previous_forward);
        previous_forward = forward_key;
        capture.fwd_log.push_back(placement::CapturedVote{
            forward_key,
            VoteHit{seed_index, forward_position, forward_start,
                    view.occurrence, view.count},
            forward_bin});
      }
      if (emit_reverse) {
        const int reverse_start = static_cast<int>(local) - reverse_position;
        const std::int64_t reverse_key =
            vote_pack_key(contig, vote_floor_div(reverse_start, W));
        const std::uint32_t reverse_bin =
            capture.rc_bins.add(reverse_key, reverse_key != previous_reverse);
        previous_reverse = reverse_key;
        capture.rc_log.push_back(placement::CapturedVote{
            reverse_key,
            VoteHit{seed_index, reverse_position, reverse_start,
                    view.occurrence, view.count},
            reverse_bin});
      }
    }
  }
  stats.emitted_forward = static_cast<std::int64_t>(capture.fwd_log.size());
  stats.emitted_reverse = static_cast<std::int64_t>(capture.rc_log.size());
  peaks = placement::aggregate_fused_capture_to_peaks(capture, k, min_support);
  return peaks;
}

// The second coarse pass. The budgeted vote needs two distinct seeds in one diagonal
// bin, so a spliced read with one exact hit per exon may have no peak at its true
// locus, leaving an empty catalogue or only a decoy that Stage 4 refuses. Such a read
// gets one more pass: every extracted seed of the whole read re-voted from the
// capture's posting views (within kRnaExplainMaxPostings), drained at
// kRnaCoprimaryMinAnchors (a re-vote sees every seed, so two-seed agreement is a
// coincidence), and the catalogue rebuilt by select_coarse_loci with the same options.
bool build_second_pass_catalogue(const ReadFacts& facts, ReadState& rs) {
  const auto& read_len = facts.read_len;
  const auto& seed_ctx = facts.seed_ctx;
  const auto& shared_fwd_syncmer_seeds = facts.shared_fwd_syncmer_seeds;
  auto& capture = facts.capture;
  const auto& k = facts.k;
  const auto& reference_count = facts.reference_count;
  auto& worker_scratch = facts.worker_scratch;
  const auto& copt = facts.copt;
  const auto& tie_seed = facts.tie_seed;
  auto& peaks = rs.peaks;
  auto& loci = rs.loci;
  auto& coarse_pass = rs.coarse_pass;
  // Every extracted seed of the whole read, both strand frames: -u is a
  // realization-time hypothesis and constrains no vote.
  const std::vector<rna::RnaQuerySpan> whole_read{
      rna::RnaQuerySpan{0, read_len}};
  rna::RnaRevoteStats revote_stats;
  std::vector<rna::placement::CoarseDiagonalPeak> pass_peaks =
      rna::rna_revote_intervals(
          seed_ctx, *shared_fwd_syncmer_seeds, capture.views, whole_read,
          read_len, k, rna::kRnaCoprimaryMinAnchors, reference_count,
          facts.readmitted, worker_scratch.explain.revote_capture,
          revote_stats);
  if (revote_stats.over_budget)
    return false;
  const rna::placement::CoarseLocusQueryCover pass_query_cover{
      &worker_scratch.explain.revote_capture.peak_read_pos_offset,
      &worker_scratch.explain.revote_capture.peak_read_pos, read_len, k,
      rna::kRnaChimeraMinQueryBases};
  std::vector<rna::placement::CoarseLocus> pass_loci =
      rna::placement::select_coarse_loci(
          pass_peaks, copt, worker_scratch.coarse, tie_seed,
          /*coalesced=*/nullptr, /*overflow_out=*/nullptr, &pass_query_cover);
  if (pass_loci.empty())
    return false;
  peaks = std::move(pass_peaks);
  loci = std::move(pass_loci);
  coarse_pass = 2;
  return true;
}

} // namespace fa::cpu::lr::rna
