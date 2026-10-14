// Seed, posting and vote-hit record types, the resolved seeding parameters, and the query-side seed selection's declarations.
#pragma once

#include "../core/core.h"
#include "../index/format.h"
#include "../index/index.h"
#include "../index/seed.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

/**** tie_hash ****/
// minimap2's read-seeded tie-break. When two candidate loci tie on every evidence key, the
// DNA vote and the RNA catalogue order them by a hash seeded by the read, as minimap2 does,
// so on a diploid reference the haplotype whose contig sorts first does not win every tie.
// minimap2:
//
//   map.c mm_map_frag  hash  = X31(qname) (0 without a name);
//                      hash ^= Wang(qlen) + Wang(opt->seed);   opt->seed = 11
//                      hash  = Wang(hash)
//   hit.c mm_gen_regs  h = hash64((hash64(a.x) + hash64(a.y)) ^ hash)
//
// with a the chain's first anchor. tie_read_seed() is mm_map_frag's per-read seed and
// tie_locus_hash() mixes it with one locus key (strand, contig and position) where minimap2
// mixes an anchor. Comparators consult it only after every evidence key ties and fall back
// to contig and position on a 64-bit collision, so they stay strict weak orders. The
// result depends only on the read name, the read length and the locus.

namespace fa {
namespace cpu {
namespace lr {

// minimap2's opt->seed default (options.c mm_mapopt_init).
inline constexpr std::uint32_t kTieHashSeed = 11;

// khash.h __ac_X31_hash_string over the read name. Bytes are read unsigned, so the value
// is the same on every target (minimap2 reads plain char, which differs above 0x7f).
std::uint32_t tie_name_hash(std::string_view name) noexcept;

// hit.c hash64 (Thomas Wang's 64-bit integer mix).
// Inline: called per comparison from tie_locus_hash, below.
inline std::uint64_t tie_hash64(std::uint64_t key) noexcept {
  key = (~key + (key << 21));
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8));
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4));
  key = key ^ key >> 28;
  key = (key + (key << 31));
  return key;
}

// mm_map_frag's per-read seed. `name_hash` is tie_name_hash(read name), or 0 without a
// name, as in minimap2; `read_len` is minimap2's qlen_sum for a single-segment read.
std::uint32_t tie_read_seed(std::uint32_t name_hash,
                            int read_len) noexcept;

// One locus's tie rank under the read's seed; lower sorts first. The key packs strand
// (bit 63), contig (bits 32-62) and the position (low 32 bits: the DNA diagonal bin or the
// RNA window's reference start), so distinct loci give distinct keys; the mix is
// hash64(hash64(key) ^ seed).
// Inline: called per bucket comparison from the vote's drain heap (vote/vote.h).
inline std::uint64_t tie_locus_hash(std::uint32_t read_seed, int contig,
                                    bool reverse, std::int64_t bin) noexcept {
  const std::uint64_t key =
      (static_cast<std::uint64_t>(reverse ? 1u : 0u) << 63) |
      (static_cast<std::uint64_t>(static_cast<std::uint32_t>(contig) &
                                  0x7fffffffu)
       << 32) |
      static_cast<std::uint64_t>(static_cast<std::uint32_t>(bin));
  return tie_hash64(tie_hash64(key) ^ static_cast<std::uint64_t>(read_seed));
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** types ****/
// Record types shared by long-read seeding and voting.

namespace fa {
namespace cpu {
namespace lr {

struct DnaLongSeedView {
  QuerySeed seed;
  KmerPostingView view;
  // An empty vote tile's seed over the cap, which the vote admits
  // (options/options.h kDnaTileRescueOcc).
  bool rescued = false;
};

struct DnaLongSeedBundle {
  std::vector<QuerySeed> seeds;
  std::vector<DnaLongSeedView> selected_views;
  // The vote uses the centred-rarest representative; exact candidate refinement uses the
  // deduplicated first, rarest and representative views from the same scan.
  std::vector<DnaLongSeedView> exact_refine_views;
  bool views_ready = false;

  void clear() {
    seeds.clear();
    selected_views.clear();
    exact_refine_views.clear();
    views_ready = false;
  }
};

struct ChainSyncmerOccCandidate {
  QuerySeed seed;
  KmerPostingView view;
  uint32_t occurrence = std::numeric_limits<uint32_t>::max();
  bool valid = false;
};

struct LongWindowAnchor {
  int ref_start_bin = 0;
  int ref_start_bin_width = 0;
  uint32_t median_occurrence = 0;
};

struct VotePeak {
  int read_lo = 0;
  int read_hi = 0;
  int chr = -1;
  int ref_pos = 0;
  bool is_rc = false;
  int support = 0;
  int center_support = 0;
  int vote_score = 0;
  LongWindowAnchor anchor;
  // The whole-read winner's per-read line d(q) = line_a + stretch(q, line_b_q20)
  // (vote/vote.h), set only on that winner, when its line passes the gate (line_gate).
  // Its chain harvest window then reaches harvest_below further below raw_ref_start - pad
  // and harvest_above further above raw_ref_start + L + pad, and its chain band follows
  // the line.
  std::int32_t harvest_below = 0;
  std::int32_t harvest_above = 0;
  std::int64_t line_a = 0;
  std::int32_t line_b_q20 = 0;
  bool line_gate = false;
  // Median reference-start diagonal before whole-read projection clamps it to
  // chromosome bounds. Coarse block projection must not inherit that clamp.
  int64_t raw_ref_start = 0;
  // Oriented read range of the seeds with a posting on this peak's diagonal (bin distance
  // <= max(1, vote_radius)) in the exact-refine re-walk: a superset of the seeds the tile
  // support mask tests. Invalid (lo > hi) when the re-walk could not certify that.
  int evidence_read_lo = 0;
  int evidence_read_hi = -1;
};

static_assert(std::is_trivially_copyable_v<VotePeak>);

// The per-read line's stretch over q query bases, floor(q * b_q20 / 2^20) for a slope
// carried as b_q20 = round(b 2^20). Integer only, so no target's float contraction can
// move a window edge; both factors are 32-bit, so the product fits int64.
inline constexpr int kVoteSlopeFractionBits = 20;
inline int vote_slope_stretch(int q, std::int32_t b_q20) {
  const std::int64_t product = static_cast<std::int64_t>(q) * b_q20;
  const std::int64_t one = std::int64_t{1} << kVoteSlopeFractionBits;
  std::int64_t quotient = product / one;
  if (product % one != 0 && product < 0)
    --quotient;
  return static_cast<int>(quotient);
}

inline bool chain_peak_better(const VotePeak &a, const VotePeak &b) {
  if (a.vote_score != b.vote_score)
    return a.vote_score > b.vote_score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center_support != b.center_support) {
    return a.center_support > b.center_support;
  }
  if (a.read_lo != b.read_lo)
    return a.read_lo < b.read_lo;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.ref_pos < b.ref_pos;
}

// The DNA catalogue's order: chain_peak_better with its positional tail replaced by
// minimap2's read-seeded tie-break (above). A tie on vote, support, centre support and
// read_lo goes to the read-seeded hash of the peak's locus (contig, strand, diagonal bin),
// lower first; only a 64-bit collision falls back to contig and position.
inline bool chain_peak_better_seeded(const VotePeak &a, const VotePeak &b,
                                     std::uint32_t read_seed) {
  if (a.vote_score != b.vote_score)
    return a.vote_score > b.vote_score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center_support != b.center_support)
    return a.center_support > b.center_support;
  if (a.read_lo != b.read_lo)
    return a.read_lo < b.read_lo;
  const std::uint64_t a_hash =
      tie_locus_hash(read_seed, a.chr, a.is_rc, a.anchor.ref_start_bin);
  const std::uint64_t b_hash =
      tie_locus_hash(read_seed, b.chr, b.is_rc, b.anchor.ref_start_bin);
  if (a_hash != b_hash)
    return a_hash < b_hash;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.ref_pos < b.ref_pos;
}

struct ChainWindowBucketVote {
  int support = 0;
  int center = 0;
  int score = 0;
  // The last query-seed epoch that added support here, so each seed votes once per
  // diagonal. 0 means unstamped; the epoch counter restarts each window.
  uint32_t voter_epoch = 0;
};

struct ChainWindowRankedBucket {
  int chr = -1;
  int bin = 0;
  int score = 0;
  int support = 0;
  int center = 0;
};

struct ChainWindowRetainedSeed {
  QuerySeed seed;
  KmerPostingView view;
  bool rescued = false;  // DnaLongSeedView::rescued
};

// One posting captured in an accumulate pass: what vote_emit_peaks's re-walk would
// recompute, so a capturing path can replay it instead.
struct VoteHit {
  int seed_idx;
  int read_pos;
  int ref_start;
  uint32_t occurrence;
  uint32_t count;
};

// The vote drain order: score, support and centre decide; a tie goes to minimap2's
// read-seeded hash of the bucket's locus (contig, drain strand, bin), lower first, and a
// 64-bit collision to contig and bin. (chr, bin) is unique per bucket, so the order is
// strict and the drain depends only on the buckets and the read's seed.
inline bool chain_window_ranked_better(const ChainWindowRankedBucket &a,
                                       const ChainWindowRankedBucket &b,
                                       std::uint32_t read_seed, bool reverse) {
  if (a.score != b.score)
    return a.score > b.score;
  if (a.support != b.support)
    return a.support > b.support;
  if (a.center != b.center)
    return a.center > b.center;
  const std::uint64_t a_hash = tie_locus_hash(read_seed, a.chr, reverse, a.bin);
  const std::uint64_t b_hash = tie_locus_hash(read_seed, b.chr, reverse, b.bin);
  if (a_hash != b_hash)
    return a_hash < b_hash;
  if (a.chr != b.chr)
    return a.chr < b.chr;
  return a.bin < b.bin;
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** posting_density ****/

namespace fa::cpu {

// Validate the chromosome directory used to interpret global posting
// coordinates. Empty chromosomes are permitted; the global coordinate space
// must start at zero and remain nondecreasing.
bool posting_density_boundaries_valid(
    const std::uint64_t* chromosome_offsets,
    std::uint64_t chromosome_count) noexcept;

}  // namespace fa::cpu

/**** context ****/
// Resolved long-read seeding parameters (LongReadSeedContext) and the vote occurrence-cap
// policy.

namespace fa { namespace cpu { namespace lr {

enum class LongOccPolicy {
    Off,
    Fixed,
    Platform,
};

struct LongOccPolicyConfig {
    LongOccPolicy policy = LongOccPolicy::Off;
    int primary_occ_cap = 0;
    int platform_occ_cap = 12;
};

// The DNA vote occurrence cap. Under LongOccPolicy::Platform, the presets' policy, it is
//     max(platform_occ_cap, mm_idx_cal_max_occ(attached index, kLongOccQuantileF))
// limited by the preset's ceiling, if any. The quantile only raises the preset floor, so a
// repeat-rich reference gets a higher cap; f is chosen so a human reference resolves to the
// floor. --max-vote-occ selects LongOccPolicy::Fixed (or Off at 0), which skips the rule
// and its index scan.
inline constexpr double kLongOccQuantileF = 1.81e-4;

// How the cap was resolved, filled at index-attach time for reporting only. `resolved` is
// false until an index is attached.
struct LongOccIndexResolution {
    bool resolved = false;
    double quantile_f = 0.0;
    std::uint64_t distinct_keys = 0;
    std::uint64_t total_postings = 0;
    std::int64_t raw_cap = 0;   // faix_cal_max_occ, before the floor
    int floor_cap = 0;          // the preset/caller cap the quantile floors to
    int ceiling_cap = 0;        // the ceiling in force, 0 = none
    int resolved_cap = 0;       // min(ceiling_cap, max(floor_cap, raw_cap))
};

std::string long_occ_policy_name(LongOccPolicy policy);

LongOccPolicy long_occ_policy_from_string(
    const std::string& value,
    LongOccPolicy fallback
);

bool long_occ_policy_string_valid(const std::string& value);

// Inline: called per seed from vote/vote.h, through the two admission tests below.
inline int effective_long_primary_occ_cap(
    const LongOccPolicyConfig& cfg
) {
    switch (cfg.policy) {
        case LongOccPolicy::Fixed:
            return std::max(0, cfg.primary_occ_cap);
        case LongOccPolicy::Platform:
            // The rule resolves into platform_occ_cap; before an index is attached this
            // is the preset floor, never 0 (no cap).
            return std::max(0, cfg.platform_occ_cap);
        case LongOccPolicy::Off:
        default:
            return 0;
    }
}

// Inline: called per seed from vote/vote.h (chain_window_seed_occ_allowed).
inline bool seed_allowed_by_long_occ_policy(
    const KmerPostingView& view,
    const LongOccPolicyConfig& cfg
) {
    const int cap = effective_long_primary_occ_cap(cfg);
    if (cap <= 0) return true;
    if (!view.found() || view.count == 0) return false;
    return view.occurrence > 0 &&
           view.occurrence <= static_cast<uint32_t>(cap);
}

// The same test for a found seed, on a given occurrence.
// Inline: called per seed from vote/vote.h (vote_seed_occ_allowed).
inline bool occurrence_allowed_by_long_occ_policy(
    uint32_t occurrence,
    const LongOccPolicyConfig& cfg
) {
    const int cap = effective_long_primary_occ_cap(cfg);
    if (cap <= 0) return true;
    return occurrence > 0 && occurrence <= static_cast<uint32_t>(cap);
}

// Seeds per strand that a larger DNA vote selection keeps (below).
inline constexpr int kVoteSeedNestBase = 128;

struct LongReadSeedContext {
    int k = 0;
    const SeedIndex* index = nullptr;

    // Closed-syncmer extraction.
    int syncmer_s = 0;                   // ClosedSyncmerConfig::s
    int syncmer_downsample = 1;          // ClosedSyncmerConfig::downsample

    // Occurrence-aware seed admission. When enabled (DNA), the kernel keeps one
    // lowest-occ seed per read-position strip (occ <= cap). When disabled (RNA),
    // the window votes every extracted seed (capped only by max_query_seeds).
    int max_query_seeds_per_strand = 0;
    // DNA: above kVoteSeedNestBase seeds per strand the selection keeps the
    // kVoteSeedNestBase selection and fills from the larger one (below).
    bool nested_vote_seeds = false;
    bool syncmer_occ_aware_enabled = false;  // DNA: true; RNA: false
    LongOccPolicyConfig occ_policy;

    // Positive preset (or explicit override) reference-start diagonal width.
    int vote_diag_bin_width = 64;
    // Adaptive width (DNA presets), in integer arithmetic:
    //   W(L) = clamp(L / vote_diag_slope_den, vote_diag_bin_width, vote_diag_width_max)
    // a length term for indel drift, floored at the base width. A denominator <= 0 drops the
    // length term. False for RNA and for a fixed --dw override.
    bool vote_diag_bin_width_adaptive = false;
    int vote_diag_slope_den = 128;
    int vote_diag_width_max = 2048;
    int min_support = 3;
    int chain_max_candidates_per_window = 0;
    // DNA ratio admission, and the all-chains lane under either admission rule: drain the
    // vote heap in batches, so the exact refine walk is paid per chunk of buckets rather
    // than per bucket. Same output as the sequential drain.
    bool vote_batched_refine = false;
    // DNA: a vote tile that admits no seed votes with its rarest found seed at or under
    // this occurrence (options/options.h kDnaTileRescueOcc); 0 = none, as on RNA.
    int tile_rescue_occ = 0;
    // DNA, options/options.h skip_self: the contig that is the query read itself,
    // found by its exact name; -1 when the read is not in the index or the option is off.
    int self_contig = -1;
    const std::vector<std::string>* chr_names = nullptr;

    bool syncmer_occ_aware_active() const;
};

// The occurrence the DNA vote counts for a seed: the index's, less the key's postings on
// the query read's own contig (ctx.self_contig), so 0 for a key that occurs only in the
// read itself.
// Inline: called per seed from vote/vote.h (vote_seed_occ_allowed).
inline uint32_t vote_seed_occurrence(const LongReadSeedContext& ctx,
                                     const KmerPostingView& view) {
    if (ctx.self_contig < 0 || !view.found()) return view.occurrence;
    const uint32_t self = static_cast<uint32_t>(ctx.self_contig);
    const uint32_t lo =
        view.positions.lower_bound_packed(pack_ref_pos(self, 0));
    const uint32_t hi =
        view.positions.lower_bound_packed(pack_ref_pos(self + 1, 0));
    return view.occurrence - (hi - lo);
}

int effective_vote_diag_bin_width(
    const LongReadSeedContext& ctx, int read_span);

}}}  // namespace fa::cpu::lr

/**** syncmer ****/
// Closed-syncmer seed extraction with occurrence-aware admission for long-read mapping, and
// the slicer that cuts a read's seed stream into windows without re-extracting.

namespace fa { namespace cpu { namespace lr {

struct ChainSeedLookupCache;    // scratch.h
struct ChainWindowPeakScratch;  // scratch.h

bool chain_syncmer_occ_candidate_better(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b
);

bool chain_syncmer_occ_candidate_better_centred(
    const ChainSyncmerOccCandidate& a,
    const ChainSyncmerOccCandidate& b,
    int tile_centre
);

int syncmer_tile_centre(int tile, int n_kmers, int n_tiles);

inline constexpr int kVoteSeedOccurrenceStrata = 32;

int vote_seed_occurrence_stratum(uint32_t occurrence);

bool vote_seed_minimum_is_representative(
    const std::array<uint32_t, kVoteSeedOccurrenceStrata>& strata,
    uint32_t admitted,
    uint32_t minimum_occurrence
);

// The vote's empty-tile rescue (options/options.h kDnaTileRescueOcc): per read and
// strand the rescued seeds' summed occurrence stays within kDnaTileRescueBudget, and reads
// shorter than kDnaTileRescueMinLen bp are not rescued.
inline constexpr uint64_t kDnaTileRescueBudget = 65536;
inline constexpr int kDnaTileRescueMinLen = 1000;

// One strand's rescue, walked in tile order: offer() every seed of the current tile, then
// close() the tile. A tile that admitted no seed votes with its rarest offered seed over the
// cap and at or under M, ranked as the tile's centred seed.
class VoteTileRescue {
public:
    VoteTileRescue(const LongReadSeedContext& ctx, int span);

    void offer(const ChainSyncmerOccCandidate& cand, int tile_centre);

    // The closing tile's rescued seed, when it admitted none and the budget allows.
    bool close(uint32_t admitted, DnaLongSeedView& out);

private:
    uint32_t max_occ_;
    uint64_t spent_ = 0;
    ChainSyncmerOccCandidate best_;
};

void extract_chain_closed_syncmer_seeds_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    std::vector<QuerySeed>& out
);

// The vote seeds of one strand. With ctx.nested_vote_seeds and more than
// kVoteSeedNestBase seeds asked, the kVoteSeedNestBase selection is kept and
// filled from the full selection, up to its size: next is the seed farthest
// from every kept position, ties to the lower occurrence, then the lower
// position. A rescued seed that would take the rescued occurrence past
// kDnaTileRescueBudget is passed over, so the fill can stop short; when the
// kept selection is the larger, nothing is added. The exact-refine views are
// the union of both selections'.
bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
);

bool extract_chain_closed_syncmer_occ_aware_seed_bundle_into(
    const LongReadSeedContext& ctx,
    const uint8_t* query_enc,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch
);

bool select_chain_closed_syncmer_occ_aware_seed_bundle_from_seeds_into(
    const LongReadSeedContext& ctx,
    const std::vector<QuerySeed>& input,
    int span,
    ChainSeedLookupCache* lookup_cache,
    DnaLongSeedBundle& out,
    ChainWindowPeakScratch& scratch,
    const KmerPostingView* seed_views = nullptr
);

// Slices a read's seed stream to the window [start, end) and renumbers read_pos from the
// window start; key and z are copied unchanged. `k` sets the last valid read_pos (end - k).
// `downsample > 1` skips seeds whose sketch_order_hash(key) % downsample != 0. `seeds` must
// be sorted by read_pos, as the extractor emits them. Same seeds as re-extracting the window.
void rebase_syncmer_stream_window(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample,
    std::vector<QuerySeed>& out
);

// True when rebase_syncmer_stream_window() with the same arguments would copy `seeds`
// unchanged (window from 0 holding every seed, no downsampling), so a caller can pass the
// stream through. O(1): the input is sorted, so its two ends decide.
bool syncmer_stream_window_is_identity(
    const std::vector<QuerySeed>& seeds,
    int start,
    int end,
    int k,
    int downsample
);

}}}  // namespace fa::cpu::lr

/**** select ****/
// Query-seed occurrence evidence and low-occurrence seed selection.

namespace fa { namespace cpu {

struct QuerySeedOccurrenceEvidence {
    QuerySeed seed{};
    uint32_t occurrence = 0;     // Reference posting count for seed.key.
    float specificity = 0.0f;    // 1 / occurrence.
    bool found = false;          // Same as `retained`.
    bool present = false;        // Key is present in the index.
    bool retained = false;       // Retained postings are available for voting.
    int original_index = 0;
};

// Inline: called per seed from rna/vote.cpp.
inline float query_seed_specificity(uint32_t occurrence) {
    return occurrence > 0 ? 1.0f / static_cast<float>(occurrence) : 0.0f;
}

// Tile-spread pick: keep the lowest-occurrence seed in each of max_seeds read tiles, so kept
// seeds are rare and spread across the read; top up from a global low-occurrence sort if
// short, then restore read order.
void select_spread_from_evidence(
    std::vector<QuerySeed>& seeds,
    const std::vector<QuerySeedOccurrenceEvidence>& evidence,
    int read_len,
    int max_seeds
);

}}  // namespace fa::cpu

namespace fa { namespace cpu { namespace lr {

void limit_query_seeds_with_lookup_cache(
    const LongReadSeedContext& ctx,
    std::vector<QuerySeed>& seeds,
    int read_len,
    ChainSeedLookupCache* lookup_cache
);

}}}  // namespace fa::cpu::lr
