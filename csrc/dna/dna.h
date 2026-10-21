// The DNA lane's context, per-worker scratch, result, block and gap types, and the declarations of its stages.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include DNA context"
#endif

#include "../chain/chain.h"
#include "../core/core.h" // output::CigarReplayRequest (cs/MD)
#include "../dp/dp.h" // ::fa::cpu::GotohCigarResult
#include "../index/format.h"
#include "../index/index.h" // engine::ReferenceContext
#include "../index/seed.h" // QuerySeed
#include "../seed/scratch.h" // ChainAnchorScratch
#include "../seed/seed.h" // LongReadSeedContext
#include "../vote/vote.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

/**** context ****/
// What the DNA stages depend on: non-owning reference views and a flat copy
// of the resolved options.

namespace fa {
namespace cpu {
namespace lr {

// Projection and realization refuse a block with fewer supporting tiles.
inline constexpr int kDnaPostCommitRecordMinBlockTiles = 1;

// minimap2's default opt->min_chain_score (-m), which map-ont and map-hifi
// keep.
inline constexpr int kDnaMinChainScore = 40;

class RetainedSeedDensity;

// The options the DNA stages read, copied from the typed options at the call
// boundary. Defaults are the shipped settings.
struct ResolvedDnaOptions {
  // Local interval-anchor harvest.
  int cigar_local_interval_anchor_interval_pad = 512;
  // Run count at or above which the dense chain uses the diagonal-keyed
  // predecessor search instead of the linear scan. Both are exact; 0 always
  // uses the diagonal search, -1 never does.
  int dna_dense_diag_min_runs = 2048;
  // Occurrence gate of the whole-query exact pass (--max-chain-occ, by
  // default the vote's resolved cap); 0 = off.
  int dna_pool_gate_occ = 0;
  // The within-read repeat filter's fraction (--q-occ-frac); 0 = off.
  float dna_query_occ_frac = 0.0f;
  // minimap2's occ_dist (-e), the seed selection's streak spacing; 0 = off.
  int occ_dist = 0;

  // Realization (CIGAR output).
  int k = 21; // seed length
  int cigar_dp_ambi = 1;
  int cigar_dp_bw = 500;
  int cigar_dp_bw_long = 20000;
  int cigar_dp_gap_extend1 = 2;
  int cigar_dp_gap_extend2 = 1;
  int cigar_dp_gap_open1 = 4;
  int cigar_dp_gap_open2 = 24;
  int cigar_dp_match = 2;
  int cigar_dp_max_gap = 5000;
  // minimap2 -s: the DP half of the per-record emission floor and the
  // local-inversion segment gate. 0 disables both.
  int cigar_dp_min_dp_max = 80;
  // minimap2 opt->zdrop_inv, for the local-inversion probe (mm_test_zdrop).
  int cigar_dp_inversion_zdrop = 200;
  // minimap2 opt->pri_ratio (-p): the credibility ratio of the -c lane's
  // secondary records (realization.cpp).
  double pri_ratio = 0.8;
  // minimap2 opt->min_chain_score (-m): the emission floor, the MAPQ's subsc
  // floor and the minimum size of an inversion middle.
  int min_chain_score = kDnaMinChainScore;
  int cigar_dp_split_min_anchors = 3;  // minimap2 opt->min_cnt
  // minimap2 opt->min_ksw_len: a stretch piece ends once both spans reach it.
  int cigar_dp_min_ksw_len = 200;
  int cigar_dp_mismatch = 4;
  int cigar_dp_tail_end_bonus = -1;
  int cigar_dp_tail_zdrop = 400;       // -z
  // The gap-fill row, -A -B -O -E -z --score-N; the cigar_dp_* row above is
  // the preset's end row, which prices every path. fill_dp_min_dp_max is -S
  // in the fill row's units (realization.cpp).
  int fill_dp_match = 4;
  int fill_dp_mismatch = 8;
  int fill_dp_ambi = 2;
  int fill_dp_gap_open1 = 8;
  int fill_dp_gap_extend1 = 4;
  int fill_dp_gap_open2 = 48;
  int fill_dp_gap_extend2 = 1;
  int fill_dp_tail_zdrop = 800;
  int fill_dp_inversion_zdrop = 200;
  int fill_dp_min_dp_max = 160;
  bool enable_full_read_cigar = true;
  // Optional cs:Z / MD:Z output (minimap2 --cs / --MD); empty by default.
  ::fa::cpu::output::CigarReplayRequest cigar_replay_request;
  // Chain MAPQ, HiFi presets only: on the DP branch the margin dp_max -
  // dp_max2 over the dp2 owner replaces the ratio form.
  bool chain_mapq_hifi_margin = false;
  // HiFi presets only: a seam's or piece's late inversion probe runs only
  // over a local chain of the read's opposite-lane seeds in its drop window
  // (dna.h).
  bool inversion_probe_local_gate = false;
  // -N, minimap2's best_n: the cap on secondary alignments in both lanes
  // (realization.cpp retain_regions).
  int secondary_max = 5;
  // lr and lr:hq: a secondary-role region ends at its first split, with no
  // continuation (realization.cpp plan_block_split).
  bool rival_split_ends = false;
  // Chain MAPQ own-locus rules, both on: a chained shadow's vote counts for
  // the winner, and the chained shadow window grows with read length. No
  // option sets them.
  bool dna_chain_mapq_shadow_vote_credit = true;
  bool dna_chain_mapq_chain_shadow_read_slack = true;
  // Vote seeding.
  int chain_syncmer_s = 9;
  int chain_syncmer_downsample = 1;
  // Admit vote peaks by their ratio to the read's best vote; 0 admits by
  // count. See catalogue_lane_bound.
  double vote_admission_ratio = 0.0;
  // The query tiles per read, the preset's (options/options.h query_tiles).
  int query_tiles = ::fa::cpu::voting::kQueryTileCount;
  // The all-chains lane (options/options.h all_chains).
  bool all_chains = false;
  // The catalogue's candidates per strand (options/options.h
  // dna_chain_max_candidates).
  int catalogue_lane_bound = ::fa::cpu::voting::kCatalogueLaneBound;
};

// Non-owning views, valid for the duration of one stage call.
struct DnaContext {
  engine::ReferenceContext ref;
  ResolvedDnaOptions opts;
  // Hash of the read name (tie_name_hash); 0 without a name, as minimap2
  // without a qname.
  std::uint32_t read_name_hash = 0;
  // tie_read_seed(read_name_hash, read length), set by map_read before the
  // vote. Every tie-break of the read uses this value.
  std::uint32_t vote_tie_seed = 0;
  // The all-chains lane with options/options.h skip_self: the contig that
  // is the query read itself, set by map_read. A chain takes no anchor on the
  // read's own exact diagonal there, as minimap2's -D. -1 otherwise.
  int self_contig = -1;
  // The all-chains lane under --dual=no (options/options.h dual): the
  // read's rank in the pair order and each contig's rank in it, set by the
  // engine (engine/engine.cpp bind_dual_order). A candidate on a contig ranked
  // below the read is not chained. 0 and null otherwise.
  int dual_rank = 0;
  const std::vector<int>* dual_contig_rank = nullptr;
  // The read's seed density for the inversion probe's local chain, kept by
  // placement (DnaPlacementChainingResult::inversion_gate_seeds) and set by
  // map_read for every realization of the read; null when not kept.
  const RetainedSeedDensity* inversion_gate_seeds = nullptr;
  // The read's seeds re-admitted past the pool gate's N (seed/seed.h
  // mark_readmitted_seeds), one byte per forward read position, set by
  // map_read when opts.occ_dist > 0, null otherwise. The vote and every pool
  // gate admit them as rescued seeds (RetainedSeedRef::rescued).
  const std::vector<std::uint8_t>* readmitted = nullptr;
};

} // namespace lr
} // namespace cpu
} // namespace fa

/**** result ****/

namespace fa {
namespace cpu {
namespace lr {
namespace dna {

struct Result : AlignResult {
  // Parallel to `supplementary`: the candidate whose block each record came
  // from, or kNullCandidate for a record that owns no block (an inversion
  // middle). Read by the per-block MAPQ.
  std::vector<::fa::cpu::voting::CandidateId> supplementary_candidates;
  // The family.block_parts entry the head prints from, -1 when the family
  // has no parts; then supplementary_parts is parallel to `supplementary` (-1
  // for a record from no part), else empty.
  int primary_part = -1;
  std::vector<int> supplementary_parts;
};

}  // namespace dna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** region_gap ****/
// A gap inside a verified region: the bases between one anchor's k-mer and
// the next one's, which passed the ungapped certificate. Offsets are relative
// to the region start and apply to both axes, since a verified region lies on
// one diagonal.

namespace fa::cpu::lr {

struct RegionGap {
  int begin = 0;
  int end = 0;
};

} // namespace fa::cpu::lr

/**** retained_seed_density ****/
// Read-local posting views of the retained and fine seeds, deduplicated by
// canonical key while each strand keeps its own query positions. A candidate
// slices a view to its reference window with two lower bounds and no new key
// lookup.

namespace fa::cpu::lr {

struct RetainedSeedDensityEntry {
  std::uint64_t key = 0;
  KmerPostingView view;
  // Some RetainedSeedRef of this key is rescued.
  bool rescued = false;
  bool retained = false;
};

struct RetainedSeedRef {
  QuerySeed seed;
  std::uint32_t entry = 0;
  // A vote seed the empty-tile rescue admitted past the cap
  // (options/options.h kDnaTileRescueOcc), or its fine twin, which
  // placement's harvest gate admits too.
  bool rescued = false;
  // minimap2's seed tandem: a seed whose query neighbour has the same key (seed.c), on every hit.
  bool tandem = false;
};

class RetainedSeedDensity {
 public:
  // A rescued vote seed's fine twin is the fine seed of its strand with its
  // read_pos and key: both come from the strand's one seed stream. A fine seed
  // is tandem by its neighbours in its strand's whole fine stream, a vote seed
  // as its fine twin. A fine seed marked in `readmitted` (DnaContext::
  // readmitted, at its forward read position) is rescued too.
  bool build(const SeedIndex& index,
             const std::vector<ChainWindowRetainedSeed>* forward,
             const std::vector<ChainWindowRetainedSeed>* reverse,
             const std::vector<QuerySeed>* fine_forward = nullptr,
             const std::vector<QuerySeed>* fine_reverse = nullptr,
             const ChainSeedLookupCache* lookup_cache = nullptr,
             const std::vector<std::uint32_t>* fine_forward_slots = nullptr,
             const std::vector<std::uint32_t>* fine_reverse_slots = nullptr,
             const std::vector<std::uint8_t>* readmitted = nullptr);

  const std::vector<RetainedSeedRef>& fine_forward() const {
    return fine_forward_;
  }
  const std::vector<RetainedSeedRef>& fine_reverse() const {
    return fine_reverse_;
  }

  // Frees buffers a much longer read left; build() reuses the rest.
  void release_excess_for_read(std::size_t read_len);

  // Entries use capture IDs when available; absent views have no refs.
  std::size_t entry_count() const { return entries_.size(); }

  // The genome-wide occurrence of an entry's key: the global_count of every
  // interval slice() finds for it.
  std::uint32_t occurrence(std::uint32_t entry_index) const {
    return entries_[entry_index].view.occurrence;
  }

  // Whether a rescued seed refers to the entry.
  bool holds_rescued(std::uint32_t entry_index) const {
    return entries_[entry_index].rescued;
  }

  // The entry's whole posting list.
  const KmerPostingView& view(std::uint32_t entry_index) const {
    return entries_[entry_index].view;
  }

  // Intervals resolved per lock-step round of slice_batch, which accepts any
  // n and chunks internally.
  static constexpr std::size_t kSliceBatch = 16;

  KmerPostingIntervalView slice(const SeedIndex& index,
                                std::uint32_t entry_index,
                                int chromosome,
                                std::uint32_t low,
                                std::uint32_t high) const;

  // Resolves n intervals for the same (chromosome, low, high) window; result i
  // equals slice(index, entry_indices[i], chromosome, low, high). The binary
  // searches of up to kSliceBatch entries are interleaved so their cache
  // misses overlap.
  void slice_batch(const SeedIndex& index,
                   const std::uint32_t* entry_indices,
                   std::size_t n,
                   int chromosome,
                   std::uint32_t low,
                   std::uint32_t high,
                   KmerPostingIntervalView* out);

 private:
  std::vector<RetainedSeedDensityEntry> entries_;
  // The forward fine seeds' entries by position, -1 where skipped.
  std::vector<std::int64_t> fine_entry_;
  std::vector<RetainedSeedRef> forward_;
  std::vector<RetainedSeedRef> reverse_;
  std::vector<RetainedSeedRef> fine_forward_;
  std::vector<RetainedSeedRef> fine_reverse_;
  // The chromosome directory build() validated; slicing with the same
  // directory skips re-validation.
  const std::uint64_t* validated_offsets_ = nullptr;
  std::uint64_t validated_chrom_count_ = 0;
};

}  // namespace fa::cpu::lr

/**** query_seed_pool ****/
// Read-local access to the query seeds of both strands. Non-owning: the
// syncmer streams stay with the caller and the posting views in the
// ChainSeedLookupCache.

namespace fa {
namespace cpu {
namespace lr {
namespace dna {

enum class DnaQuerySeedStrand : unsigned char { Forward, Reverse };

class DnaQuerySeedPool {
public:
  DnaQuerySeedPool(const LongReadSeedContext& context, int read_length,
                   ChainSeedLookupCache& lookup_cache)
      : context_(context), read_length_(read_length),
        lookup_cache_(lookup_cache) {}

  void bind(DnaQuerySeedStrand strand,
            const std::vector<QuerySeed>& full_stream);

  // Optional caller-owned stores for the capture's views and cache slots, one
  // entry per seed of the bound stream. Without them the capture only warms
  // the cache.
  void bind_capture_stores(std::vector<KmerPostingView>* forward_views,
                           std::vector<std::uint32_t>* forward_slots,
                           std::vector<KmerPostingView>* reverse_views,
                           std::vector<std::uint32_t>* reverse_slots);

  // The captured views / slots of a strand, or null when the strand was not
  // captured into a store.
  const std::vector<KmerPostingView>* captured_views(
      DnaQuerySeedStrand strand) const;
  const std::vector<std::uint32_t>* captured_slots(
      DnaQuerySeedStrand strand) const;

  bool bound(DnaQuerySeedStrand strand) const;

  // Resolves every query key into the read-local cache. Later exact chains
  // read only these cached views.
  bool capture_posting_views(DnaQuerySeedStrand strand);

  // Selects the representative vote seeds of the whole read, once per strand;
  // a second call returns false.
  bool prepare_representatives(
      DnaQuerySeedStrand strand, int start, int end, int downsample,
      std::vector<QuerySeed>& rebased_seeds, DnaLongSeedBundle& bundle,
      ChainWindowPeakScratch& selection_scratch);

private:
  const std::vector<QuerySeed>* stream(DnaQuerySeedStrand strand) const {
    return strand == DnaQuerySeedStrand::Forward ? forward_ : reverse_;
  }

  const std::vector<QuerySeed>*& stream_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_ : reverse_;
  }

  bool& prepared_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_prepared_
                                                 : reverse_prepared_;
  }

  bool& captured_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_captured_
                                                 : reverse_captured_;
  }

  template <class T>
  static void mirror(const std::vector<T>& from, std::vector<T>& to) {
    to.resize(from.size());
    for (std::size_t i = 0; i < from.size(); ++i)
      to[i] = from[from.size() - 1 - i];
  }

  const LongReadSeedContext& context_;
  int read_length_ = 0;
  ChainSeedLookupCache& lookup_cache_;
  const std::vector<QuerySeed>* forward_ = nullptr;
  const std::vector<QuerySeed>* reverse_ = nullptr;
  bool forward_prepared_ = false;
  bool reverse_prepared_ = false;
  bool forward_captured_ = false;
  bool reverse_captured_ = false;
  std::vector<KmerPostingView>* forward_views_ = nullptr;
  std::vector<std::uint32_t>* forward_slots_ = nullptr;
  std::vector<KmerPostingView>* reverse_views_ = nullptr;
  std::vector<std::uint32_t>* reverse_slots_ = nullptr;
};

} // namespace dna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** worker_scratch ****/
// Per-worker mapping scratch. The batch driver owns one per worker thread and
// reuses it for every read of a mini-batch, like minimap2's mm_tbuf_t.

namespace fa { namespace cpu { namespace lr {

struct DnaWorkerScratch {
  ChainAnchorScratch chain;
  // The read's seed density, built by placement and read by realization
  // until the worker's next read.
  RetainedSeedDensity seed_density;

  // Clears per-read state, then trims retained capacity for the next read.
  void clear_for_read(std::size_t read_len) {
    chain.clear_for_read();
    chain.release_excess_for_read(read_len);
    seed_density.release_excess_for_read(read_len);
  }
};

}}} // namespace fa::cpu::lr

/**** placement_family_adapter ****/

namespace fa::cpu::lr {

struct DnaPlacementCandidate {
  ::fa::cpu::voting::CandidateId id =
      ::fa::cpu::voting::kNullCandidate;
  VotePeak peak;
  ::fa::cpu::voting::QueryTileMask support;
  std::uint64_t equivalence_key = 0;
  int lane = 0;
  int catalogue_rank = 0;
  int vote_evidence = 0;
};

// The chain a block prints from and its bounds in forward-read bp
// (dna/dna.h). path -1 is the candidate's primary, j its
// sibling_paths[j]; [anchor_begin, anchor_end) in chain order. score and
// anchors are the block's own; item_score and item_anchors its whole chain's;
// pool_subsc and pool_n_sub minimap2's subsc and n_sub over the chains of the
// same candidate attached to the block's owner.
struct DnaBlockPart {
  int path = -1;
  int anchor_begin = 0;
  int anchor_end = 0;
  int forward_begin = 0;
  int forward_end = 0;
  int score = 0;
  int anchors = 0;
  int item_score = 0;
  int item_anchors = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct DnaPlacementFamily {
  int read_length = 0;
  int seed_length = 0;
  // The read's query tiles; masks, blocks and assignments use them.
  int tile_count = ::fa::cpu::voting::kQueryTileCount;
  int forward_candidates = 0;
  int reverse_candidates = 0;
  std::uint64_t exact_posting_tests = 0;
  std::vector<DnaPlacementCandidate> candidates;
  ::fa::cpu::voting::QueryPartitionResult partition;
  // Parallel to partition.selected.blocks once the owners are decided on the
  // dense chains' anchors; empty before.
  std::vector<DnaBlockPart> block_parts;
  bool valid = false;

  const DnaPlacementCandidate* find(
      ::fa::cpu::voting::CandidateId id) const noexcept;
};

// Maps an oriented seed position to its forward-read query tile of
// tile_count. The generic catalogue and partition code only sees forward-read
// tiles.
int dna_forward_query_tile(
    int oriented_seed_position, bool reverse, int read_length, int seed_length,
    int tile_count) noexcept;

// Builds the candidate catalogue from the whole-read vote peaks. A candidate's
// support is the set of forward query tiles holding a retained seed with a
// posting compatible with its geometry.
DnaPlacementFamily build_dna_placement_family(
    const DnaContext& context, const std::vector<std::uint8_t>& forward_query,
    const std::vector<VotePeak>& whole_read_peaks,
    const ChainWindowPeakScratch& forward_scratch,
    const ChainWindowPeakScratch& reverse_scratch);

}  // namespace fa::cpu::lr

/**** chain_ownership ****/
// Which whole-query chain owns which stretch of the read, decided on the
// chains' anchors (ownership.cpp). Reads anchors and scores only.

namespace fa::cpu::lr {

// minimap2's -n: a chain owns or attaches only with at least this many
// anchors, and a block keeps at least this many.
inline constexpr int kDnaOwnerMinAnchors = 3;

// One whole-query chain: a candidate's primary (path -1) or sibling j.
struct DnaOwnershipPath {
  int candidate = -1;
  int path = -1;
  int contig = -1;
  bool reverse = false;
  const std::vector<chaining::Anchor>* anchors = nullptr;
  int score = 0;
};

// One block of the read: item `item`'s anchors [a0, a1) in chain order, the
// bounds [forward_begin, forward_end) the blocks tile the read with, its score
// and anchor count, and minimap2's subsc and n_sub over the chains of the
// item's candidate attached to its owner.
struct DnaOwnershipBlock {
  int item = 0;
  int a0 = 0;
  int a1 = 0;
  int forward_begin = 0;
  int forward_end = 0;
  int score = 0;
  int kept = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

// The parent walk's outcome in its order: a chain that became an owner
// (owner_item == item) or was attached to one. An owner's row carries
// minimap2's subsc and n_sub over the chains attached to it, and the same two
// over those of its own candidate, as they stand at the end; an attachment's
// are 0. Items index DnaOwnershipSelection::items.
struct DnaOwnershipRole {
  int item = 0;
  int owner_item = 0;
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct DnaOwnershipSelection {
  // The paths past the floors (at least kDnaOwnerMinAnchors anchors, score
  // >= min_chain_score), input order; DnaOwnershipBlock::item indexes it.
  std::vector<DnaOwnershipPath> items;
  // In query order; empty when no chain owns.
  std::vector<DnaOwnershipBlock> blocks;
};

// `paths` in catalogue order, each candidate's primary first, then its
// siblings by index. `roles` receives the parent walk's roles.
DnaOwnershipSelection select_chain_owners(
    int read_length, int seed_length, int min_chain_score,
    const std::vector<DnaOwnershipPath>& paths,
    std::vector<DnaOwnershipRole>& roles);

}  // namespace fa::cpu::lr

/**** placement_chaining ****/
// Placement chaining. The whole-read vote catalogue is the only source of
// candidates; this stage chains every candidate over the whole query and
// decides the read's blocks on those chains' anchors. Both map-only
// projection and CIGAR realization consume the resulting family and anchor
// paths. No DP runs here.

namespace fa::cpu::lr {

enum class DnaPlacementChainStatus : std::uint8_t {
  NotSelected = 0,
  Accepted = 2,
  InvalidReference,
  MissingSeeds,
  NoIntervalHits,
  NoChain,
};

struct DnaPlacementCandidateChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  ::fa::cpu::voting::QueryTileMask sparse_support;
  ::fa::cpu::voting::QueryTileMask dense_support;
  std::vector<chaining::Anchor> primary;
  // Sibling paths from the same exact DP, items of the ownership selection
  // beside the primary.
  std::vector<std::vector<chaining::Anchor>> sibling_paths;
  // Parallel to sibling_paths: each sibling's chain score. rival_sibling is
  // the index of the sibling rival_chain_score came from
  // (ChainPartition::f2_index, or the best sibling of a record that does not
  // hold its group's primary; placement.cpp chain_pooled_groups), or -1 when
  // there is none or it was not materialized. Whole-query passes only; read by
  // the ownership selection (select_block_owners).
  std::vector<int> sibling_scores;
  int rival_sibling = -1;
  DnaPlacementChainStatus status = DnaPlacementChainStatus::NotSelected;
  std::int64_t interval_hits = 0;
  // Undercounts: the harvest gate skips slicing the seeds it drops.
  int filtered_hits = 0;
  int chain_score = 0;
  int rival_chain_score = 0;
  int overlapping_rivals = 0;
  int chain_anchors = 0;
  int oriented_query_begin = -1;
  int oriented_query_end = -1;
  int forward_query_begin = -1;
  int forward_query_end = -1;
  int reference_begin = -1;
  int reference_end = -1;
  int final_selected_tiles = 0;
  bool exact = false;
  // Harvest counters. deferred_anchors: postings of the slices deferred as
  // metadata, within this candidate's window. rescued_anchors: those it
  // admits, the anchors a harvest of its window alone would build.
  std::uint64_t deferred_anchors = 0;
  std::uint64_t rescued_anchors = 0;
  // Runs produced by the whole-query pass's collapse; on the
  // group's first member with a window only.
  std::int64_t dense_runs = 0;
};

// Catalogue rivals chained over the whole query for the MAPQ, per ranking.
// mm_set_mapq2 reads one best rival and a count; two let the best rival by
// vote lose to the runner-up on chain score. Under
// kDnaChainMapqStudyBlockRivals each block also ranks its own.
inline constexpr int kDnaMapqRivalChains = 2;
// Only rivals with at least 1/kDnaMapqRivalVoteDenominator of the strongest
// owner's vote are taken; the rest still enter the MAPQ through their vote.
inline constexpr int kDnaMapqRivalVoteDenominator = 4;

// A catalogue rival chained over the whole query for the MAPQ; a group's
// chain partition only sees rivals inside its members' windows. Evidence only: kept
// out of `candidates`, which is parallel to family.candidates and drives
// selection.
struct DnaRivalExactChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  DnaPlacementCandidateChain chain;
};

// One chain of the final ownership selection, kept for the -c lane's
// controller (realization.cpp): a candidate's primary (path -1) or
// sibling path.
struct DnaSelectionItem {
  int candidate = -1;
  int path = -1;
  int contig = -1;
  bool reverse = false;
  int score = 0;
};

// The final selection's chains and the parent walk's roles over them, in its
// order (dna/dna.h); role items index `items`.
struct DnaKeptSelection {
  std::vector<DnaSelectionItem> items;
  std::vector<DnaOwnershipRole> roles;
};

struct DnaPlacementChainingResult {
  DnaPlacementFamily family;
  std::vector<DnaPlacementCandidateChain> candidates;
  // The MAPQ's cross-locus rivals; empty when none is worth chaining.
  std::vector<DnaRivalExactChain> rival_exact;
  // The read's seed density when realization gates the late inversion probe
  // on it (dna.h); null otherwise. It is the caller's and lives
  // until the caller's next build.
  const RetainedSeedDensity* inversion_gate_seeds = nullptr;
  // The same density, for the -c controller's divergence filter
  // (dna/realization.cpp retain_regions).
  const RetainedSeedDensity* strand_filter_seeds = nullptr;
  bool accepted = false;
  // Not accepted because no accepted whole-query chain qualifies to own a
  // block (dna/dna.h).
  bool no_owner_chain = false;
  // The final ownership selection.
  DnaKeptSelection kept_selection;

  const DnaPlacementCandidateChain*
  find(::fa::cpu::voting::CandidateId candidate) const noexcept;

  // A candidate's whole-query exact chain: its MAPQ rival chain (any status),
  // else its own record, read only for a block owner or a `realized` rival
  // (the CIGAR lane's dp_max2 owner). Null when there is none. Callers that
  // need an accepted chain check the status.
  const DnaPlacementCandidateChain*
  whole_query_chain(::fa::cpu::voting::CandidateId candidate,
                    bool realized = false) const noexcept;

  // A kept selection item's anchors, in chain order; null when gone.
  const std::vector<chaining::Anchor>*
  selection_anchors(const DnaSelectionItem& item) const noexcept;
};

// Chain parameters for the candidate paths under diagonal band `band`: reach
// max(band, cigar_dp_max_gap).
chaining::ColinearChainParams
dna_candidate_chain_params(const DnaContext& context, int band,
                           int seed_length);

// The same parameters for the dense run chain: the scoring is unchanged, and
// max_iter and max_skip have no counterpart. The diagonal-keyed search runs
// when diag_min_runs >= 0 && runs >= diag_min_runs, the linear scan
// otherwise; both are exact.
chaining::DenseChainParams
dna_dense_chain_params(const chaining::ColinearChainParams& params,
                       int seed_length, int diag_min_runs);

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
    RetainedSeedDensity& seed_index);

} // namespace fa::cpu::lr

/**** placement_chaining_internal ****/
// The interface shared by the DNA chaining files, placement.cpp and
// ava.cpp: the harvest rules, each stated once for
// placement's whole-query passes and the all-chains lane's target pools, and
// the pool and record steps both take.

namespace fa::cpu::lr {
namespace internal {

// The band of a whole-query chain pass, for its chain and its window's
// widening, and of the ava lane's harvest: -r's second value.
int pass_diagonal_band(const DnaContext& context);

// A candidate's harvest window [low, high) on its contig: the read's span from
// its peak's expected start, padded on both sides by
// cigar_local_interval_anchor_interval_pad, widened by harvest_below /
// harvest_above and clamped to the contig. The widening reaches the whole-read
// winner's per-read line (vote_slope_widen_winner) and is 0 on every other
// peak. Empty when high <= low.
struct HarvestWindow {
  std::int64_t low = 0;
  std::int64_t high = 0;
};
HarvestWindow harvest_window(const DnaContext& context,
                             const VotePeak& peak, int read_length,
                             std::int64_t chromosome_length);

// The band test: whether a posting on `diagonal` lies within `diagonal_band`
// of the candidate's expected diagonal, in the caller's integer type.
template <class Diagonal>
inline bool within_band(Diagonal diagonal, Diagonal expected,
                        int diagonal_band) {
  return std::abs(diagonal - expected) <= diagonal_band;
}

// Whether a pass over `contig` and the lane drops the posting at each seed's
// own read position: the read's own exact diagonal gives no anchor, so the
// contig is the read itself (DnaContext::self_contig) and the lane forward.
bool skips_own_diagonal(const DnaContext& context, int contig,
                        bool reverse);
// The harvest's occurrence gate. A seed whose key occurs more than `cap` times
// genome-wide gives no anchor unless it is rescued (RetainedSeedRef::rescued).
// The cap is dna_pool_gate_occ (by default the vote's cap, or
// --max-chain-occ); there is no gate when that is 0.
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
PoolGate pool_gate(const DnaContext& context);

// The refusal every candidate on `contig` in the lane of `seeds` meets alike,
// or Accepted when there is none: MissingSeeds when the lane has no seeds,
// InvalidReference when the index cannot serve the contig or `query` is not
// the read's length. Only metadata is needed: contig lengths come from the
// index offsets and the read length from the family, and a map-only run
// passes `query` empty.
DnaPlacementChainStatus contig_lane_refusal(
    const DnaContext& context, int contig,
    const std::vector<RetainedSeedRef>& seeds,
    const std::vector<std::uint8_t>& query, int read_length);

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

/**** target_chaining ****/
// The all-chains lane's chaining (the ava-ont and ava-hifi presets), one chain
// call per target of a read.

namespace fa::cpu::lr {

// The all-chains lane: the whole-query exact chains of the catalogue's
// candidates, parallel to family.candidates, with no partition. The
// candidates on one target, a contig and strand, are chained in one call,
// over the union of their pools, and take its chains best score
// first, in catalogue order, one each. A candidate that receives no chain
// keeps the refusing status, and one on a contig ranked below the read
// (DnaContext::dual_rank) is not chained (NotSelected). Empty when the read's
// seeds cannot be indexed.
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
    const std::vector<std::uint32_t>* fine_reverse_slots);

}  // namespace fa::cpu::lr

/**** record_family ****/

namespace fa::cpu::lr {

struct DnaSegmentRecord {
  AlignResult alignment;
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  // The family.block_parts entry, -1 for none.
  int part = -1;
};

struct DnaRecordFamily {
  AlignResult primary;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  int primary_part = -1;
  // Each supplementary keeps the candidate whose block it came from, because
  // MAPQ is scored per block owner.
  std::vector<DnaSegmentRecord> supplementary;
  bool valid = false;
};

} // namespace fa::cpu::lr

/**** family_projection ****/

namespace fa::cpu::lr {

struct DnaPlacementChainingResult;
struct DnaPlacementCandidateChain;

// Projects the part of an accepted chain inside `block` (forward-query tiles)
// to a record without CIGAR. Returns false on invalid geometry.
bool dna_project_chained_block(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementCandidateChain& chain,
    const ::fa::cpu::voting::QueryBlock& block, AlignResult& output);

enum class DnaProjectionFallbackReason : std::uint8_t {
  None,
  ChainingUnavailableOrRefused,
  NoOwnerChain,
  NoProjectedRecords,
  MissingOrInvalidSelectedChain,
  InvalidBlockGeometry,
  OverlappingProjectedQueryBlocks,
  TransactionValidationFailure,
  InvalidPlacementFamily,
};

const char* dna_projection_fallback_name(
    DnaProjectionFallbackReason reason) noexcept;

// A projected block record with what the join below reads: its owner, contig
// and strand, the first and last anchors it was projected from (oriented
// order), and the fuzzy walk over those anchors before the record's trims.
struct DnaJoinPiece {
  AlignResult record;
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  int chromosome = -1;
  bool reverse = false;
  chaining::Anchor first;
  chaining::Anchor last;
  int walk_matches = 0;
  int walk_length = 0;
  // The family.block_parts entry, -1 for none.
  int part = -1;
};

struct DnaFamilyProjectionResult {
  dna::Result output;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  bool committed = false;
  DnaProjectionFallbackReason fallback =
      DnaProjectionFallbackReason::None;
  // Every projected block record, in partition order.
  std::vector<DnaJoinPiece> pieces;

  bool fell_back() const noexcept {
    return fallback != DnaProjectionFallbackReason::None;
  }
};

// Returns the committed projection, or the incumbent with score 0 and the
// refusal reason recorded.
dna::Result settle_map_only_projection(DnaFamilyProjectionResult result);

DnaFamilyProjectionResult project_map_only_placement_family(
    const DnaContext& context,
    const DnaPlacementChainingResult& placement,
    const dna::Result& incumbent);

// Joins consecutive block records of the committed map-only family, in
// forward-query order, into one record wherever they pass the bridge test (no
// duplicate seam anchor, then plan_dna_family_join) with the centers of the
// records' corner anchors as the cuts. A record that owns no block never
// joins. The joined record spans both, its columns 10 and 11 are the walk over
// both anchor runs, and it is owned by the piece of the higher whole chain
// score, the leftmost on a tie. No DP and no reference base. `realized` and `primary_candidate`
// change only when two records join.
void dna_join_map_only_family(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    std::vector<DnaJoinPiece> pieces, dna::Result& realized,
    ::fa::cpu::voting::CandidateId& primary_candidate);

// Each attached chain dna_retained_secondary_items keeps, appended to
// realized.secondary as one record projected from its whole path, MAPQ 0,
// cm:i and s1:i its path's. A path that does not project is left out.
void dna_project_map_only_secondaries(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    dna::Result& realized);

}  // namespace fa::cpu::lr

/**** ordered_anchor_path ****/
// The anchor path of one block and the pre-DP control minimap2's mm_align1
// applies to it: k-mer center adjustment, terminal DP windows, bad-seed
// filtering and the verified realization geometry. No DP runs here; the
// geometry planner reads bases only to test gaps against the certificate.
//
// Anchor coordinates are exact k-mer starts. minimap2 stores a minimizer end
// and mm_adjust_minier() subtracts k/2 to reach the center;
// adjusted_endpoint() applies the equivalent start-to-center shift.

namespace fa {
namespace cpu {
namespace lr {
namespace ordered_anchor {

enum class CoordinateConvention : std::uint8_t {
  ExactKmerStart,
};

// A window [data, data + count) of an anchor array the caller owns, as
// minimap2's (r->as, r->cnt) into a[]. Indices are window-local.
struct AnchorWindow {
  ::fa::cpu::chaining::Anchor* data = nullptr;
  std::size_t count = 0;

  std::size_t size() const { return count; }
  bool empty() const { return count == 0; }
  ::fa::cpu::chaining::Anchor& operator[](std::size_t i) const {
    return data[i];
  }
  ::fa::cpu::chaining::Anchor& front() const { return data[0]; }
  ::fa::cpu::chaining::Anchor& back() const { return data[count - 1]; }
};

struct OrderedAnchorPath {
  int reference_id = -1;
  bool reverse_complemented = false;
  int query_length = 0;
  int target_length = 0;
  int query_bound_begin = 0;
  int query_bound_end = 0;
  int target_bound_begin = 0;
  int target_bound_end = 0;
  int minimizer_k = 0;
  CoordinateConvention coordinate_convention =
      CoordinateConvention::ExactKmerStart;
  bool homopolymer_compressed = false;

  // Coordinates are in the oriented query and the forward reference.
  // [realization_begin, realization_end) is the subrange mm_fix_bad_ends
  // keeps; the anchors outside it still bound the terminal windows.
  // Normalization flags the anchors in place.
  AnchorWindow selected;
  std::size_t realization_begin = 0;
  std::size_t realization_end = 0;
};

enum class ValidationCode : std::uint8_t {
  Valid,
  EmptySelectedPath,
  InvalidReference,
  InvalidBounds,
  UnsupportedCoordinateConvention,
  UnsupportedHomopolymerCompression,
  InvalidAnchorSpan,
  InconsistentMinimizerSpan,
  AnchorOutOfBounds,
  NonMonotoneSelectedPath,
  InvalidAnchorFlags,
  InvalidSelectedSubrange,
  InvalidControl,
};

struct ValidationResult {
  ValidationCode code = ValidationCode::Valid;
  std::size_t anchor_index = 0;

  explicit operator bool() const { return code == ValidationCode::Valid; }
};

ValidationResult validate(const OrderedAnchorPath& path);

struct AdjustedEndpoint {
  int query = 0;
  int target = 0;
};

// mm_adjust_minier for exact k-mer starts, without HPC:
// start + (k - 1) - floor(k / 2).
AdjustedEndpoint adjusted_endpoint(const OrderedAnchorPath& path,
                                   const ::fa::cpu::chaining::Anchor& anchor);

struct NormalizationControl {
  int bandwidth = 0;
  int min_chain_score = 0;
  int max_gap = 0;
  bool trim_left = true;
  bool trim_right = true;
};

// Pre-DP normalization as minimap2's mm_fix_bad_ends, mm_filter_bad_seeds and
// mm_filter_bad_seeds_alt, in place: the trimmed range is written to `path`
// and the bad-seed filters, which always run, flag anchors ANCHOR_IGNORE /
// ANCHOR_LONG_JOIN in its window. End trimming is enabled per side. `path`
// must pass validate(); a window after a cut of a valid one does, since the
// filters add only known flags. Only the control is checked here.
ValidationResult normalize_for_realization(OrderedAnchorPath& path,
                                           const NormalizationControl& control);

// The verified realization geometry (see align.cpp). A block is laid
// out left to right as steps of three kinds.
enum class GeometryStepKind : std::uint8_t {
  // A verified region, emitted ungapped from its first anchor's center to its
  // last's, with no DP.
  Verified,
  // The fill between two verified regions with no retained anchor between
  // them.
  Seam,
  // One fill of a stretch: unverified material between verified regions that
  // are not neighbours, before the first or after the last. minimap2's loop
  // cuts a stretch into pieces.
  Piece,
};

struct GeometryRegion {
  std::size_t first = 0; // selected index of its first retained anchor
  std::size_t last = 0;  // selected index of its last retained anchor
  int anchors = 0;       // retained anchors in it
  bool verified = false;
  // Gaps that failed the certificate, counting all of them.
  int failing_gaps = 0;
};

// One executor call, in path order. Consecutive steps abut and tile
// [center(first retained anchor), center(last retained anchor)].
struct GeometryStep {
  GeometryStepKind kind = GeometryStepKind::Piece;
  // The right corner anchor, where a Z-drop split cuts (minimap2's `i`).
  std::size_t selected_index = 0;
  int query_begin = 0; // [center(left corner), center(right corner))
  int query_end = 0;
  int target_begin = 0;
  int target_end = 0;
  // Seam / Piece: the right corner carries ANCHOR_LONG_JOIN.
  bool long_join = false;
  // Seam / Piece: bw_long_eff(), or max(q_span, r_span) on a long join.
  // Verified: bw_long_eff(), recorded; no DP runs.
  int selected_band = -1;
  // Verified: its gaps, [gap_begin, gap_begin + gap_count) in
  // GeometryPlan::gaps, relative to the step start on both axes.
  std::size_t gap_begin = 0;
  std::size_t gap_count = 0;
  // Piece: its index in GeometryPlan::stretches and its ordinal within the
  // stretch. -1 elsewhere.
  int stretch = -1;
  int piece = -1;

  // Counters only. Verified: the region's retained anchors. Piece: the
  // retained anchors it crosses. Seam: 0.
  int anchor_count = 0;
  // Seam / Piece: max(next.q - prev.q_end, next.r - prev.r_end) between the
  // corner anchors, negative when they overlap. 0 on a verified region.
  int anchor_free = 0;
};

// One stretch of unverified material. Counters only.
struct GeometryStretch {
  std::size_t first_step = 0; // its first piece in GeometryPlan::steps
  int pieces = 0;
  int query_bases = 0;    // center(end corner) - center(start corner)
  int lone_anchors = 0;   // regions of one retained anchor inside it
  int failed_regions = 0; // regions of >= 2 anchors with a failing gap
  int failed_anchors = 0; // the retained anchors of those regions
  int failed_bases = 0;   // their query span, center to center
  int failing_gaps = 0;   // the gaps in them that failed the certificate
};

struct GeometryPlan {
  int initial_query_cursor = 0; // center(first retained anchor)
  int initial_target_cursor = 0;
  int final_query_cursor = 0; // center(last retained anchor)
  int final_target_cursor = 0;
  std::vector<GeometryStep> steps;
  std::vector<GeometryStretch> stretches;
  // The gaps of every verified region.
  std::vector<::fa::cpu::lr::RegionGap> gaps;
};

struct GeometryControl {
  int normal_long_band = 0; // DpMapOpt::bw_long_eff()
  // minimap2 opt->min_ksw_len: inside a stretch a piece ends at a retained
  // anchor once both spans since the previous corner reach it.
  int min_ksw_len = 0;
};

// Inputs of the ungapped certificate, borrowed for the call. The ungapped
// alignment of an L-base gap is optimal when its score S0 reaches
// (L - 1) * a - 2 * g(1): a gapped alignment of two equal-length windows keeps
// at most L - 1 match columns and pays two one-base gaps. L < 2 always passes.
// This allows up to 3 mismatches per gap with map-hifi scoring, 2 with map-ont.
struct GapCertificate {
  const std::uint8_t* query = nullptr;  // the whole oriented query (anchor.q)
  const std::uint8_t* target = nullptr; // the whole chromosome (anchor.r)
  // ksw2_simple_mat(a, b, sc_ambi): N scores sc_ambi, even against N.
  std::array<std::int8_t, 25> matrix{};
  int match = 0;        // a
  int one_base_gap = 0; // g(1) = min(q + e, q2 + e2)
};

// Plans the verified geometry of one normalized block, left to right:
//   1. interior anchors flagged ANCHOR_IGNORE or ANCHOR_TANDEM are skipped, as
//      in minimap2; the first and last anchors never are;
//   2. the retained anchors are cut into regions, and each region of two or
//      more anchors is tested gap by gap, a gap being the bases between
//      consecutive k-mers, [p.q_end(), n.q); bases of skipped anchors inside a
//      gap are tested like any other;
//   3. a verified region becomes a Verified step, neighbouring verified
//      regions are joined by a Seam, and everything else is a stretch cut into
//      Pieces by minimap2's loop.
// Regions are laid out one at a time, only as far as the steps asked for, so
// a caller that stops at a Z-drop, as minimap2's mm_align1 loop breaks there,
// plans nothing past it. `path` must be valid normalize_for_realization
// output; nothing is re-checked.
class GeometryPlanner {
public:
  GeometryPlanner(const OrderedAnchorPath& path, const GeometryControl& control,
                  const GapCertificate& certificate, GeometryPlan& plan);
  // True once plan.steps holds step `s`; false when the block ends before it.
  bool reach(std::size_t s);

private:
  bool skipped(std::size_t i) const;
  GeometryStep fill(GeometryStepKind kind, std::size_t left,
                    std::size_t right) const;
  void emit_stretch(std::size_t start, std::size_t end,
                    GeometryStretch stretch);
  // Lays out the next region; false once the block is laid out.
  bool advance();

  const OrderedAnchorPath& path_;
  const GeometryControl control_;
  const GapCertificate& certificate_;
  GeometryPlan& plan_;
  std::size_t begin_ = 0;
  std::size_t last_ = 0;
  std::size_t next_ = 0;  // the next anchor the region scan reads
  GeometryRegion open_;   // the region the scan is growing
  bool done_ = false;
  // The block's first anchor or the last anchor of the verified region laid
  // out last; `pending` collects the unverified regions since then, which
  // belong to the next stretch.
  std::size_t corner_ = 0;
  bool previous_verified_ = false;
  GeometryStretch pending_;
};

struct TerminalWindowControl {
  int max_gap = 0;
  int match = 0;
  int gap_open = 0;
  int gap_extend = 0;
  // minimap2's cap at nearby seeds (mm_align1): where set, the left window
  // starts at or after (left_cap_query, left_cap_target) and the right one
  // ends at or before (right_cap_query, right_cap_target).
  bool left_cap = false;
  int left_cap_query = 0;
  int left_cap_target = 0;
  bool right_cap = false;
  int right_cap_query = 0;
  int right_cap_target = 0;
};

struct TerminalWindowPlan {
  ValidationResult validation;
  AdjustedEndpoint first_adjusted;
  AdjustedEndpoint last_adjusted;

  int query_start = 0;  // qs0
  int target_start = 0; // rs0
  int query_end = 0;    // qe0
  int target_end = 0;   // re0

  bool has_left_packet() const {
    return query_start < first_adjusted.query &&
           target_start < first_adjusted.target;
  }
  bool has_right_packet() const {
    return last_adjusted.query < query_end && last_adjusted.target < target_end;
  }

  explicit operator bool() const {
    return validation.code == ValidationCode::Valid;
  }
};

// `path` must be valid normalize_for_realization output; only the control is
// checked.
TerminalWindowPlan plan_terminal_windows(const OrderedAnchorPath& path,
                                         const TerminalWindowControl& control);

} // namespace ordered_anchor
} // namespace lr
} // namespace cpu
} // namespace fa

/**** realization_controller ****/
// The realization controller: validates one DP packet, picks an executor
// (ksw2, or base comparison for a verified region) and returns a packet-local
// result. Slicing, concatenation, clipping and commit are the caller's.

namespace fa {
namespace cpu {
namespace lr {
namespace realization {

enum class RealizationRole : std::uint8_t {
  DnaLeftExtension,
  DnaInternalFill,
  DnaRightExtension,
  DnaSupplementaryFill,
  DnaLocalInversionMiddle,
  RnaExonFill,
};

enum class EndpointObjective : std::uint8_t {
  OpenEndedLeft,
  PinnedGlobal,
  OpenEndedRight,
};

enum class SequenceOrientation : std::uint8_t {
  Forward,
  ReverseComplement,
};

struct SequenceSlice {
  const std::uint8_t* data = nullptr;
  int length = 0;
  int original_begin = 0;
  int original_end = 0;
  SequenceOrientation orientation = SequenceOrientation::Forward;
  bool reversed_for_executor = false;
};

struct AffineGapRow {
  int open = 0;
  int extend = 0;
};

enum class ScoreModel : std::uint8_t {
  SingleAffine,
  DualAffine,
  Splice,
};

struct ScoringContract {
  int match = 0;
  int mismatch = 0;
  int ambiguity = 0;
  std::array<std::int8_t, 25> substitution_matrix{};
  AffineGapRow gap1;
  AffineGapRow gap2;
  ScoreModel model = ScoreModel::DualAffine;
};

// dp/dp.h DpInversionProbeControl's gate: the late probe runs only where
// count(source, window) is at least min_count; source is borrowed for the call.
struct InversionProbeGate {
  int (*count)(const void* source, int query_begin, int query_end,
               int target_begin, int target_end) = nullptr;
  const void* source = nullptr;
  int min_count = 0;
};

struct RealizationRequest {
  RealizationRole role = RealizationRole::DnaInternalFill;
  EndpointObjective objective = EndpointObjective::PinnedGlobal;
  SequenceSlice query;
  SequenceSlice target;
  ScoringContract scoring;

  int caller_requested_band = -1;
  int configured_band = -1;
  int configured_long_band = -1;
  int selected_band = -1;
  int zdrop = -1;
  int inversion_zdrop = -1;
  bool inversion_probe_enabled = false;
  int inversion_max_gap = 0;
  int inversion_min_chain_score = 0;
  int inversion_min_dp_max = 0;
  InversionProbeGate inversion_probe_gate;
  // A gap fill takes dp_fill_gap's exact first pass.
  bool exact_first_pass = false;
  int end_bonus = -1;
  std::int64_t matrix_cell_cap = 0;
  bool long_join = false;
  std::uint32_t ksw_flags = 0;

  // A verified region: one diagonal, every base under an exact-match anchor
  // or in one of the certified gaps below. Runs on ExecutorId::Verified.
  bool verified = false;
  // Region-relative certified gaps, borrowed for the call. Empty when the
  // region's anchors overlap end to end.
  const ::fa::cpu::lr::RegionGap* gaps = nullptr;
  int gap_count = 0;
};

enum class ControllerDecision : std::uint8_t {
  FireKsw2,
  FireExecutor,
  PureInsertion,
  PureDeletion,
  RefuseMatrixCap,
  RejectInvalid,
  RejectInvalidResult,
  UnsupportedExecutor,
};

enum class OutcomeKind : std::uint8_t {
  Aligned,
  Zdropped,
  SplitRequested,
  PureInsertion,
  PureDeletion,
  RefusedMatrixCap,
  InvalidPacket,
  InvalidResult,
  UnsupportedExecutor,
};

enum class ExecutorId : std::uint8_t {
  None,
  Ksw2,
  // Base comparison over a verified region; no DP.
  Verified,
};

// The ksw2 kernels dp/dp.h can run (DpKernelId), plus None.
enum class KernelId : std::uint8_t {
  None,
  Extz2,
  Extd2,
};

enum class SupportReason : std::uint8_t {
  Supported,
  InvalidPacket,
  UnsupportedObjective,
  UnsupportedOrientation,
  UnsupportedScoreModel,
  UnsupportedSubstitutionMatrix,
  UnsupportedFlags,
  UnsupportedBand,
  ExecutorUnavailable,
  // A verified region handed to ksw2, or anything else handed to Verified.
  UnsupportedPacketClass,
};

enum class SelectionReason : std::uint8_t {
  ProductionDefault,
  RequestedExecutor,
  UnsupportedExecutorFallback,
};

struct RawAlignmentResult {
  bool has_cigar = false;
  bool zdropped = false;
  bool reach_end = false;
  int max_score = 0;
  int max_query = -1;
  int max_target = -1;
  int query_end_score = 0;
  int query_end_target = -1;
  int target_end_score = 0;
  int target_end_query = -1;
  int corner_score = 0;
  std::vector<std::uint32_t> packed_cigar;
};

struct KernelAttempt {
  ExecutorId executor = ExecutorId::None;
  KernelId kernel = KernelId::None;
  std::uint32_t flags = 0;
  int band = -1;
  int zdrop = -1;
  bool approximate = false;
  bool exact = false;
  // raw.packed_cigar is empty even when raw.has_cigar: an attempt keeps only
  // metadata and scores.
  // (dp/dp.h notify_attempt): release consumers get only metadata and scores.
  RawAlignmentResult raw;
};

enum class InversionProbeStatus : std::uint8_t {
  Disabled,
  DropBelowInversionThreshold,
  InvalidBounds,
  WeakReverseScore,
  Accepted,
  Gated,
};

struct ZdropTestEvidence {
  bool tested = false;
  int code = 0;
  int max_drop = 0;
  int target_begin = -1;
  int target_end = -1;
  int query_begin = -1;
  int query_end = -1;
  InversionProbeStatus inversion_status = InversionProbeStatus::Disabled;
  bool reverse_probe_attempted = false;
  int reverse_score = 0;
  int reverse_query_end = -1;
  int reverse_target_end = -1;
};

// Work of one verified region: bases scored and gaps emitted. It runs no
// kernel, so it records no attempt.
struct VerifiedWork {
  std::int64_t bases = 0;
  int gaps = 0;
};

struct ExecutionTrace {
  static constexpr std::size_t kMaxAttempts = 4;
  std::array<KernelAttempt, kMaxAttempts> attempts{};
  std::size_t attempt_count = 0;
  bool overflow = false;
  ZdropTestEvidence zdrop_test;
  bool exact_retry = false;
  int exact_retry_zdrop = -1;
  VerifiedWork verified;

  void push(KernelAttempt attempt);
};

// The traceback is `raw.packed_cigar` (BAM-packed, ops past N folded to M);
// `alignment.cigar` stays empty and `alignment` supplies only the accounting
// (ref_offset, ref_consumed, matches, score, zdropped, max_q, max_t).
struct RealizationOutcome {
  OutcomeKind kind = OutcomeKind::InvalidPacket;
  ControllerDecision decision = ControllerDecision::RejectInvalid;
  ExecutorId requested_executor = ExecutorId::None;
  ExecutorId selected_executor = ExecutorId::None;
  ExecutorId committed_executor = ExecutorId::None;
  SelectionReason selection_reason = SelectionReason::ProductionDefault;
  SupportReason support_reason = SupportReason::InvalidPacket;
  SupportReason fallback_reason = SupportReason::Supported;
  RawAlignmentResult raw;
  ::fa::cpu::GotohCigarResult alignment{};
  ExecutionTrace trace;
};

struct ControllerObserver {
  using Callback = void (*)(void* user_data, const RealizationRequest& request,
                            const RealizationOutcome& outcome);

  Callback callback = nullptr;
  void* user_data = nullptr;
};

struct ExecutionContext {
  ControllerObserver observer;
};

struct ExecutorResult {
  RawAlignmentResult raw;
  ::fa::cpu::GotohCigarResult alignment{};
};

using SupportsFn = bool (*)(const RealizationRequest& request,
                            SupportReason& reason);
using ExecuteFn = ExecutorResult (*)(const RealizationRequest& request,
                                     ExecutionContext& context,
                                     ExecutionTrace& trace);

struct ExecutorEntry {
  ExecutorId id = ExecutorId::None;
  const char* name = "none";
  SupportsFn supports = nullptr;
  ExecuteFn execute = nullptr;
};

struct ExecutorTable {
  const ExecutorEntry* entries = nullptr;
  std::size_t size = 0;
};

struct ControllerPolicy {
  ExecutorId selected_executor = ExecutorId::Ksw2;
};

const char* role_name(RealizationRole role);

SequenceSlice make_query_slice(const std::uint8_t* data, int length,
                               int strand_begin, int read_length,
                               bool reverse_complemented);
SequenceSlice make_target_slice(const std::uint8_t* data, int length,
                                int reference_begin);

RealizationOutcome run(const RealizationRequest& request,
                       const ControllerPolicy& policy,
                       const ExecutorTable& executors,
                       ExecutionContext& context);

} // namespace realization
} // namespace lr
} // namespace cpu
} // namespace fa

/**** dp_runner ****/
// The banded DP entry point of DNA realization.

namespace fa {
namespace cpu {
namespace lr {

// Dual-affine scoring and bands, built once by the caller.
struct DpScoringParams {
  int match = 2;     // minimap2 -A
  int mismatch = 4;  // minimap2 -B, positive
  int ambi = 1;      // minimap2 --score-N
  int gap_open1 = 4; // minimap2 -O
  int gap_extend1 = 2;  // minimap2 -E
  int gap_open2 = 24;   // minimap2 -O2
  int gap_extend2 = 1;  // minimap2 -E2
  int zdrop = 400; // minimap2 -z
  int tail_end_bonus =
      -1;       // minimap2 opt->end_bonus
  int bw = 500; // minimap2 -r, first value
  int bw_long = 20000; // minimap2 -r, second value
  // Local-inversion parameters, read only by DNA internal and supplementary
  // fills.
  int inversion_zdrop = 200;    // minimap2 zdrop_inv
  int inversion_max_gap = 5000; // the preset's max_gap
  int inversion_min_chain_score = 40;
  int inversion_min_dp_max = 80;
};

// The certified gaps of a verified region, region-relative. Passed only for
// a verified region; null for every other packet.
struct VerifiedRegionInputs {
  const ::fa::cpu::lr::RegionGap *gaps = nullptr;
  int gap_count = 0;
};

// Runs one realization packet through the controller. Slicing the packet and
// committing the resulting CIGAR are left to the caller. `inversion_probe_gate`,
// when set, gates the late probe of a fill whose probe is enabled.
// `exact_first_pass` asks a gap fill for dp_fill_gap's exact first pass.
realization::RealizationOutcome run_dna_long_realization(
    realization::RealizationRole role, const DpScoringParams &dp,
    realization::SequenceSlice query, realization::SequenceSlice target,
    int caller_requested_band, int interior_zdrop = -1, bool long_join = false,
    realization::ExecutionContext *execution_context = nullptr,
    bool inversion_probe_enabled = false,
    const VerifiedRegionInputs *verified_region = nullptr,
    const realization::InversionProbeGate *inversion_probe_gate = nullptr,
    bool exact_first_pass = false);

} // namespace lr
} // namespace cpu
} // namespace fa

/**** postdp_scoring ****/
// minimap2's DP terms of one emitted record, read from its CIGAR: the
// maximum-scoring segment under the log-gap cost (mm_update_extra), which is
// the ms:i tag and the -c lane's dp_max, the event identity
// (mm_event_identity) and dp_max recomputed under b2 (mm_update_dp_max).

namespace fa::cpu::lr {

// The ms:i value of one record: the maximum-scoring segment of its CIGAR under
// the log-gap cost, as minimap2 computes ms. Returns -1 (no tag) for an
// unmapped record, an empty CIGAR or an unknown contig.
int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc);

// The terms of minimap2's mm_update_dp_max (align.c) for one record, under
// the same pricing as dna_record_dp_max_segment: its event identity
// (mm_event_identity) and its dp_max recomputed under b2 (mm_recal_max_dp, a
// gap costing b2 + log2(1 + len), floored at 0). False, and -1, for a record
// that pricing cannot read.
bool dna_record_event_identity(const DnaContext& context,
                               const AlignResult& record,
                               const std::vector<std::uint8_t>& fwd,
                               const std::vector<std::uint8_t>& rc,
                               double& identity);
int dna_record_recal_dp_max(const DnaContext& context,
                            const AlignResult& record,
                            const std::vector<std::uint8_t>& fwd,
                            const std::vector<std::uint8_t>& rc, double b2);
// mm_update_dp_max's b2 from the best record's event identity.
double dna_rank_b2(double identity, int match, int mismatch) noexcept;

} // namespace fa::cpu::lr

/**** geometry_work ****/
// Counters for the steps of verified-geometry realization. DP cells are not
// counted here; they go to ksw2_attempts / estimated_cells.

namespace fa::cpu::lr {

struct DnaGeometryWork {
  int verified_regions = 0;        // verified regions emitted, no DP
  std::int64_t verified_bases = 0; // their bases
  int seams = 0;                   // fills between two verified neighbours
  int stretches = 0;               // stretches whose first piece fired
  int stretch_pieces = 0;          // pieces fired, one gap fill each
  std::int64_t stretch_bases = 0;  // those stretches' query bases

  // Regions of two or more anchors that went to a stretch because a gap failed
  // the certificate: their anchors, query bases and failing gaps.
  int failed_regions = 0;
  int failed_region_anchors = 0;
  std::int64_t failed_region_bases = 0;
  int failing_gaps = 0;

  // Seams and pieces whose right corner is a long join, filled at
  // minimap2's max(q_span, r_span) band.
  int join_calls = 0;

  void add(const DnaGeometryWork& other);
};

} // namespace fa::cpu::lr

/**** family_realization ****/

namespace fa::cpu::lr {

struct DnaPlacementChainingResult;

enum class DnaFamilyFailure : std::uint8_t {
  None,
  InvalidInput,
  PlacementChainRefused,
  NoOwnerChain,
  NoSelectedFamily,
  NoAnchors,
  InvalidOrderedPath,
  InteriorRefused,
  TerminalRefused,
  FamilyValidation,
};

enum class DnaFamilyJoinKind : std::uint8_t {
  NotEligible,
  ZeroSeam,
  PureInsertion,
  PureDeletion,
  TwoAxis,
};

struct DnaFamilyJoinGeometry {
  int left_reference = -1;
  int right_reference = -1;
  bool left_reverse = false;
  bool right_reverse = false;
  int left_forward_begin = 0;
  int left_forward_end = 0;
  int right_forward_begin = 0;
  int right_forward_end = 0;
  int left_query_begin_cut = 0;
  int left_query_end_cut = 0;
  int right_query_begin_cut = 0;
  int right_query_end_cut = 0;
  int left_target_begin_cut = 0;
  int left_target_end_cut = 0;
  int right_target_begin_cut = 0;
  int right_target_end_cut = 0;
  int max_query_gap = 0;
  int max_reference_gap = 0;
  std::int64_t matrix_cell_cap = 0;
};

struct DnaFamilyJoinPlan {
  DnaFamilyJoinKind kind = DnaFamilyJoinKind::NotEligible;
  int query_gap = -1;
  int reference_gap = -1;
};

DnaFamilyJoinPlan plan_dna_family_join(
    const DnaFamilyJoinGeometry& geometry) noexcept;

// Whether the two anchors that meet at a family seam are the same anchor.
// Reverse-strand blocks run backwards in the oriented query, so their seam is
// left.begin / right.end rather than left.end / right.begin.
bool dna_family_seam_has_duplicate_anchor(
    const ordered_anchor::OrderedAnchorPath& left,
    const ordered_anchor::OrderedAnchorPath& right,
    bool reverse) noexcept;

struct DnaFamilyRealizationRequest {
  const DnaPlacementFamily* family = nullptr;
  const DnaPlacementChainingResult* placement = nullptr;
  const std::vector<std::uint8_t>* forward_query = nullptr;
  const std::vector<std::uint8_t>* reverse_query = nullptr;
};

// An interior fill takes dp_fill_gap's exact first pass when the read's last
// interior fills hold enough exact retries that did not confirm the Z-drop.
// The window, in the read's last interior fills (fewer than 32).
inline constexpr int kDnaExactFirstWindow = 16;
// The unconfirmed retries in the window that make the first pass exact.
inline constexpr int kDnaExactFirstUnconfirmed = 2;

// The -c lane's realization of one read (realization.cpp) and its work.
struct DnaFamilyRealizationOutcome {
  dna::Result output;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  DnaFamilyFailure failure = DnaFamilyFailure::InvalidInput;
  int failed_block = -1;
  int block_count = 0;
  // Realized segments. Differs from block_count when a Z-drop split a region
  // into further segments, or a segment was cut at its start and realized
  // nothing.
  int segment_count = 0;
  // Interior fills truncated at a Z-drop, and how many of those continued in
  // a new segment.
  int zdrop_truncations = 0;
  int zdrop_splits = 0;
  // Local inversions (minimap2 mm_test_zdrop, mm_split_reg, mm_align1_inv).
  // Probes: reverse scans run over an interior fill's largest drop. Splits:
  // probes that found an inversion. Middles: inversion middles attempted.
  // Records: middles emitted as records.
  int inversion_probes = 0;
  int inversion_splits = 0;
  int inversion_middles = 0;
  int inversion_records = 0;
  // The read's interior gap fills in realization order, the latest in bit 0:
  // set where the fill's exact retry returned a fill that did not Z-drop
  // (kDnaExactFirstWindow).
  std::uint32_t unconfirmed_window = 0;
  int unit_count = 0;
  int ksw2_attempts = 0;
  std::int64_t estimated_cells = 0;
  DnaGeometryWork geometry;
};

const char* dna_family_failure_name(DnaFamilyFailure failure) noexcept;

}  // namespace fa::cpu::lr

/**** region_realization ****/
// The -c lane's controller in minimap2's form, from the final ownership
// selection to records. Every owner and the attached
// chains mm_select_sub keeps are regions; each is realized from its whole
// anchor vector as mm_align1 realizes one (both ends trimmed, terminal windows
// to the read ends capped at nearby seeds, a Z-drop splitting off the rest as
// a region of its own), the records are filtered (mm_filter_regs,
// mm_update_dp_max) and roles are assigned again on the aligned spans
// (mm_set_parent, mm_select_sub). Records may overlap.

namespace fa::cpu::lr {

// What the MAPQ reads for one primary record: the SAM primary first, then each
// supplementary in order.
struct DnaRegionPrice {
  // The catalogue candidate of its chain; -1 for an inversion middle.
  int candidate = -1;
  bool inversion = false;
  // minimap2's score and cnt (the region's share of its chain after Z-drop
  // splits, and its anchors), and score0 with the whole chain's anchors.
  int score = 0;
  int cnt = 0;
  int score0 = 0;
  int anchors0 = 0;
  // The whole chain's forward-query and reference spans.
  int chain_q_begin = -1;
  int chain_q_end = -1;
  int chain_ref_begin = -1;
  int chain_ref_end = -1;
  int dp_max = 0;
  int dp_max2 = 0;
  // The catalogue candidate of the record dp_max2 came from, -1 for none.
  int dp_max2_candidate = -1;
  // dp_max and dp_max2 before the rank rescale: the sweep scores (ms:i) of
  // the record and of the record dp_max2 came from.
  int dp_max0 = 0;
  int dp_max2_0 = 0;
  // mm_set_parent's subsc and n_sub, and the same two over the secondaries
  // from its own candidate's chains.
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
  // The entry of DnaRegionOutcome::rivals that dp_max2 came from, -1 for
  // none.
  int dp_max2_rival = -1;
};

// A record the controller realized that mm_set_parent made a secondary,
// whether select_records kept it or not, with its region's whole chain as the
// MAPQ reads a rival's.
struct DnaRealizedRival {
  int candidate = -1;
  // Its chain's path in the candidate (DnaSelectionItem::path).
  int path = -1;
  int region = -1;
  // The chain's score and anchors (DnaRegionPrice score0 / anchors0) and its
  // forward-query and reference spans.
  int score0 = 0;
  int anchors0 = 0;
  int chain_q_begin = -1;
  int chain_q_end = -1;
  int chain_ref_begin = -1;
  int chain_ref_end = -1;
};

struct DnaRegionOutcome {
  // `output` holds the SAM primary, the other primaries as supplementaries
  // (supplementary_candidates parallel, no parts) and each secondary as a
  // hypothesis of its own. When no record survives it is unmapped, with
  // `failure` and `failed_block` naming the first region refused, or `floored`
  // set when records were realized and the filter dropped them all.
  // block_count counts the regions, segment_count the segments realized and
  // unit_count the records before the filter.
  DnaFamilyRealizationOutcome family;
  std::vector<DnaRegionPrice> prices;
  // Every realized secondary.
  std::vector<DnaRealizedRival> rivals;
  bool floored = false;
};

DnaRegionOutcome realize_dna_regions(const DnaContext& context,
                                     const DnaFamilyRealizationRequest& request);

// The attached chains of placement.kept_selection that the -c lane's
// selection keeps (realization.cpp retain_regions), as item indices in the
// parent walk's order. Reads no reference base.
std::vector<int> dna_retained_secondary_items(
    const DnaContext& context, const DnaPlacementChainingResult& placement);

}  // namespace fa::cpu::lr

/**** inv_local_chain ****/
// The local chain that gates the late inversion probe of a seam or piece: the
// read's own fine seeds of the lane opposite the block, inside the fill's drop
// window, chained in that lane's frame. The seeds and their posting views are
// placement's (RetainedSeedDensity), so a window costs two lower bounds per
// seed and no key lookup. A posting counts only in the lane its orientation
// bit allows, and a key only within the whole-query harvest's occurrence gate.
// The chain gates the probe and nothing else.

namespace fa::cpu::lr {

class RetainedSeedDensity;

// The probe runs when the chain holds at least this many anchors.
inline constexpr int kDnaInvLocalMinAnchors = 3;
// Consecutive chained anchors differ by at most this many bases in diagonal,
inline constexpr int kDnaInvLocalDiagonal = 10;
// and by at most this many on either axis.
inline constexpr int kDnaInvLocalGap = 200;

// One fill of a block: the block's frame and the fill's origin in it.
struct DnaInvLocalChainSource {
  const RetainedSeedDensity* seeds = nullptr;
  const SeedIndex* index = nullptr;
  int chromosome = -1;
  // The block's orientation; the chain reads the other lane.
  bool block_reverse = false;
  int read_length = 0;
  int seed_length = 0;
  // The fill's first base in the block's oriented query and on the contig.
  int query_offset = 0;
  int target_offset = 0;
  // Keys over this genome-wide occurrence are skipped; 0 admits every key.
  std::uint32_t occurrence_cap = 0;
};

}  // namespace fa::cpu::lr

/**** chain_mapq ****/
// The specification of the DNA MAPQ is the comment under mapq.cpp's chain_mapq banner.

namespace fa::cpu::lr {

// minimap2 constants.
inline constexpr double kDnaChainMapqCoef = 40.0; // q_coef
// minimap2's mask_level: a rival overlapping the committed hypothesis over
// more than half the shorter query span competes with it; less overlap means
// another segment of the same read.
inline constexpr double kDnaChainMapqMaskLevel = 0.5;
// The vote floor under which pen_cm is discounted linearly, and the vote cap
// q <= (int)(4.343 * w / 3).
inline constexpr double kDnaChainMapqVoteFloor = 40.0;
inline constexpr double kDnaChainMapqVoteCapCoef = 4.343;
inline constexpr double kDnaChainMapqVoteCapDiv = 3.0;
// The vote cap applies only while the winner has fewer votes than this; above
// it the chain's own score and anchors carry the absolute evidence.
inline constexpr int kDnaChainMapqVoteCapMinVote = 10;
// Divergence contrast (dna.h): a record whose event divergence
// (PAF de:f) exceeds the family minimum by z > kDnaChainMapqDivergentRecordZ
// is likely on the wrong repeat copy and has its MAPQ capped at
// kDnaChainMapqDivergentRecordMapq. Only records that clear the emission
// floor take part, and the primary is judged like any other record. The event
// denominator is block_len + ambiguities - gap_bases + gap_opens.
inline constexpr double kDnaChainMapqDivergentRecordZ = 8.0;
inline constexpr int kDnaChainMapqDivergentRecordMapq = 3;
// HiFi margin rule (R1-R6 above).
inline constexpr int kDnaChainMapqHifiRuleVersion = 2;
// R2: every competing rival other than the seat must be strictly below this
// strength.
inline constexpr double kDnaChainMapqHifiRivalStrengthMax = 1.0;
// R3: the cap for a margin above one substitution (A + B) and below a
// two-base gap (O1 + 2 E1), which one gap event (the HiFi homopolymer error)
// explains.
inline constexpr int kDnaChainMapqHifiSingleEventMapq = 9;
// R6: the split-read cap; (a) the rival strength that triggers it; (b) the
// minimum event count (mismatches + indel events) and mismatch share. (b)
// applies with or without a seat.
inline constexpr int kDnaChainMapqHifiSplitCapMapq = 3;
inline constexpr double kDnaChainMapqHifiSplitRivalStrength = 0.4;
inline constexpr int kDnaChainMapqHifiSubstitutionEvents = 20;
inline constexpr double kDnaChainMapqHifiSubstitutionShare = 0.45;
// Shadow rule, chained rivals: a chained rival is the winner's own locus when
// its chain overlaps the winner's reference span over more than half the
// shorter span, lies within that span extended by the slack on each side, and
// starts and ends within the slack of the winner's diagonals. A repeat copy
// one period away sits on a diagonal of its own. The default slack is the
// interval pad of the whole-query scan window.
inline constexpr int kDnaChainMapqDefaultShadowSlackBp = 512;
// With chain_shadow_read_slack the slack is max(shadow_slack_bp, floor,
// read_len / fraction): a rival chain harvests its own anchors and can reach a
// few kb past the winner's span at an SV flank. The fraction is wider than the
// peak window below because the chain's containment must also absorb that
// overhang.
inline constexpr int kDnaChainMapqChainShadowFloorBp = 2000;
inline constexpr int kDnaChainMapqChainShadowFraction = 5;
// Shadow rule, un-chained rivals: a vote peak on the winner's contig and strand
// within max(floor, read_len / fraction) of the winner's diagonal, with no
// more vote than the winner, is the winner's own locus.
inline constexpr int kDnaChainMapqPeakShadowFloorBp = 2000;
inline constexpr int kDnaChainMapqPeakShadowFraction = 10;
// The rule bits (header comment).
inline constexpr int kDnaChainMapqStudyRelativeShareCap = 2;
inline constexpr int kDnaChainMapqStudyRecordShadow = 8;
inline constexpr int kDnaChainMapqStudyFamilyVouch = 64;
inline constexpr int kDnaChainMapqStudyBlockRivals = 256;

// The bits in force under the HiFi presets.
inline constexpr int kDnaChainMapqStudy =
    kDnaChainMapqStudyRelativeShareCap | kDnaChainMapqStudyRecordShadow |
    kDnaChainMapqStudyFamilyVouch | kDnaChainMapqStudyBlockRivals;
static_assert(kDnaChainMapqStudy == 330, "the HiFi MAPQ rule set");

// The rule bits for a run: kDnaChainMapqStudy under the HiFi presets, else 0.
inline constexpr int dna_chain_mapq_study_bits(bool hifi_margin_rule) noexcept {
  return hifi_margin_rule ? kDnaChainMapqStudy : 0;
}

// One catalogue rival of the committed hypothesis. A chained rival carries its
// whole-query chain; an un-chained one only its vote evidence and query
// span, enough for the vote factor but not for f2 / n_sub.
struct DnaChainMapqRival {
  // candidate and status are not read by the formula.
  int candidate = -1;
  int status = 0;
  int chain_score = 0;
  int chain_anchors = 0;
  int q_begin = -1;
  int q_end = -1;
  int vote_evidence = 0;
  bool chained = false;
  // The evidence's dp2 is this rival's DP score, as minimap2 attributes
  // dp_max2 to one rival; n_sub then counts it once under either rule.
  bool owns_dp2 = false;
  // Locus, for the shadow rule: the vote peak's contig and strand, the
  // chain's contig-local half-open reference span (chained rivals only), and
  // the vote peak's reference-start diagonal.
  int contig = -1;
  bool reverse = false;
  int ref_begin = -1;
  int ref_end = -1;
  std::int64_t peak_diagonal = 0;
  // Bit 8 / 64: the committed alignment maps both ends of this rival's chain
  // to within the chained shadow slack of its reference span, so the rival is
  // part of this placement. The formula treats it as a shadow.
  bool record_shadow = false;
  // Which record vouched: 0 none, 1 the block's own record (bit 8), 2 another
  // record or the family jointly (bit 64). Not read by the formula.
  int vouch = 0;
};

struct DnaChainMapqEvidence {
  // minimap2's min_chain_score (-m): a rival scoring under it counts as it.
  int min_chain_score = kDnaMinChainScore;
  int f1 = 0;  // committed whole-query chain score      (mm2 score0 / score)
  int cnt = 0; // committed whole-query chain anchors    (mm2 cnt)
  // The committed candidate's sibling f2 from its own chain partition, a free
  // lower bound on subsc.
  int sib_f2 = 0;
  int winner_q_begin = -1;
  int winner_q_end = -1;
  int winner_vote = 0;
  int read_len = 0;
  // The winner's locus, for the shadow rule; see DnaChainMapqRival.
  int winner_contig = -1;
  bool winner_reverse = false;
  int winner_ref_begin = -1;
  int winner_ref_end = -1;
  std::int64_t winner_peak_diagonal = 0;
  // The chained shadow slack: the whole-query scan window's interval pad.
  int shadow_slack_bp = kDnaChainMapqDefaultShadowSlackBp;
  std::vector<DnaChainMapqRival> rivals;
  // CIGAR output only; 0 otherwise. After the controller's rank rescale
  // (minimap2's dp_max): the pair the ratio reads.
  double dp1 = 0.0;
  double dp2 = 0.0;
  // dp1 and dp2 before the controller's rank rescale (minimap2's dp_max0):
  // the pair the margin terms read. Equal to dp1 and dp2 when the rescale
  // did not fire.
  double dp1_raw = 0.0;
  double dp2_raw = 0.0;
  // mlen/blen of the realized primary; 1.0 when there is no base-level
  // accounting (map-only), which makes the identity factor inert.
  double identity = 1.0;
  int match_sc = 1; // minimap2 opt->a
  // minimap2's opt->a * 2 + opt->b: the dp2 owner within this of dp1 counts
  // toward n_sub (mm_set_parent).
  int sub_diff = 0;
  // mm_set_parent's `ri->cnt >= rp->cnt` clause: a chained rival with at
  // least the winner's anchors counts toward n_sub.
  bool n_sub_cnt_clause = true;
  // Shadow refinements, set from ResolvedDnaOptions (on by default) and
  // off on a bare struct.
  //   shadow_vote_credit: a chained shadow's vote counts for the winner, in
  //   the absolute term and the strength denominators, since one locus's
  //   vote can be split between twin candidates.
  //   chain_shadow_read_slack: see dna_chain_mapq_chain_shadow_slack.
  bool shadow_vote_credit = false;
  bool chain_shadow_read_slack = false;
  // HiFi margin rule (header comment), HiFi presets only: the margin
  // dp1_raw - dp2_raw is read against the scoring row below.
  bool hifi_margin_rule = false;
  int substitution_cost = 0; // A + B
  int gap_open1 = 0;         // O1
  int gap_extend1 = 0;       // E1
  // R6: the family's block-owning records, and this record's mismatches (-1
  // without accounting) and indel events (-1 without a CIGAR).
  int family_blocks = 1;
  int record_mismatches = -1;
  int record_ins_events = -1;
  int record_del_events = -1;
  // The rule bits (dna_chain_mapq_study_bits); the formula reads only bit 2.
  // For it: this record's aligned bases, and the (mismatches, aligned bases)
  // of the family's cleanest other accounted record, -1 when there is none.
  int study_bits = 0;
  int record_aligned_bases = -1;
  int family_min_mismatches = -1;
  int family_min_aligned = -1;
  // The record prints from a block part (DnaBlockPart). f1 and cnt are then
  // its whole chain's and sib_f2 the part's same-candidate subsc; pen_cm and
  // the map-only logarithm read part_score and part_anchors (minimap2's score
  // and cnt against score0), and n_sub starts at the part's pool_n_sub.
  bool part_priced = false;
  int part_score = 0;
  int part_anchors = 0;
  int part_n_sub = 0;
};

// The formula's verdict on one rival, for the caller.
struct DnaChainMapqRivalVerdict {
  // minimap2's mask_level test on the query spans.
  bool admissible = false;
  // The rival is the winner's own locus.
  bool shadow = false;
  // Not read by the formula: the half-overlap test without containment.
  bool shadow_half_only = false;
  // Chained exact tie at another locus.
  bool tie = false;
  // The rival's contribution to x (0 when it contributes nothing).
  double strength = 0.0;
  // The rival counted toward n_sub. The standalone verdict applies only the
  // anchor rule; dna_chain_mapq adds minimap2's DP rule for the dp2 owner.
  bool counted = false;
};

// Every intermediate the formula forms. Filled only when the caller passes a
// non-null pointer.
struct DnaChainMapqBreakdown {
  int rivals_total = 0;
  int rivals_admissible = 0; // minimap2's mask_level test, shadows included
  int rivals_shadow = 0;     // admissible rivals the shadow rule excluded
  int rivals_tie = 0;
  int f2 = 0;
  int subsc = 0;
  int n_sub = 0;
  double vote_ratio = 0.0; // largest vote_i / w over the competing rivals
  double x_chain = 0.0;    // minimap2's subsc / f1 over the competing rivals
  double x_used = 0.0;     // the x the emitted branch consumed
  double pen_cm = 0.0;     // after every discount, as consumed
  int cap = 60;            // the absolute vote cap as applied
  bool dp_branch = false;
  // HiFi margin rule: the largest strength among competing rivals other than
  // the seat, whether the rule decided the MAPQ, and its Phred before the
  // tail (-1 when it did not).
  double x_other = 0.0;
  bool hifi_margin_applied = false;
  int hifi_margin_mapq = -1;
  // The seat (0 none, 1 the dp2 owner) and its margin on the rescaled pair
  // (R1 reads the pair before the rescale), the pen the margin rule uses,
  // and the split-read cap applied (0 none, 1 rival strength, 2 substitution
  // share, 3 both).
  int seat_kind = 0;
  double seat_m_eff = 0.0;
  double hifi_pen = 0.0;
  int split_cap = 0;
};

// The chained shadow slack in bp: max(shadow_slack_bp,
// kDnaChainMapqChainShadowFloorBp, read_len / kDnaChainMapqChainShadowFraction)
// with chain_shadow_read_slack, else shadow_slack_bp.
std::int64_t dna_chain_mapq_chain_shadow_slack(
    const DnaChainMapqEvidence& evidence) noexcept;

// The chain part of the shadow rule as geometry: the winner's chain against a
// rival's, each a contig-local reference span, a forward-query span and a
// strand, with `slack` bp of tolerance. The caller checks contig and strand.
// Without query spans on_diagonal stays false; without reference spans
// nothing is set. `shadow` is half && contained && on_diagonal.
struct DnaChainMapqChainShadow {
  bool half = false;
  bool contained = false;
  bool on_diagonal = false;
  bool shadow = false;
};

DnaChainMapqChainShadow dna_chain_mapq_chain_shadow(
    int winner_ref_begin, int winner_ref_end, int winner_q_begin,
    int winner_q_end, bool winner_reverse, int rival_ref_begin,
    int rival_ref_end, int rival_q_begin, int rival_q_end, bool rival_reverse,
    int read_len, std::int64_t slack) noexcept;

// The per-rival verdict, without the parts that depend on the whole
// evidence: the DP factor of `strength`, the DP rule of `counted` and the
// shadow vote credit, which dna_chain_mapq adds.
DnaChainMapqRivalVerdict
dna_chain_mapq_rival_verdict(const DnaChainMapqEvidence& evidence,
                             const DnaChainMapqRival& rival) noexcept;

// minimap2's MAPQ for the committed hypothesis, in [0,60]. Returns 0 when
// there is no committed whole-query chain (f1 <= 0). `verdicts`, when
// non-null, receives one entry per evidence.rivals element, in order.
int dna_chain_mapq(const DnaChainMapqEvidence& evidence,
                   DnaChainMapqBreakdown* breakdown = nullptr,
                   std::vector<DnaChainMapqRivalVerdict>* verdicts = nullptr);

} // namespace fa::cpu::lr

/**** record_divergence ****/
// The divergence contrast between the records of one read. All records of a
// read share one error rate, so a record whose divergence sits far above the
// family's minimum is likely aligned to the wrong copy of a repeat.

namespace fa::cpu::lr {
namespace dna {

// `accounted` and `divergence` (gap-compressed, as PAF de:f) are inputs: the
// caller sets `accounted` for a record with base-level accounting that clears
// the emission floor and is not an inversion middle. `z` and `divergent` are
// outputs.
struct DnaRecordDivergence {
  bool accounted = false;
  bool divergent = false;
  double divergence = 0.0;
  double z = 0.0;
};

// The family: index 0 is the primary, index i + 1 the i-th supplementary.
struct DnaFamilyDivergence {
  int accounted = 0;
  std::vector<DnaRecordDivergence> records;
};

// Fills `family.accounted` and each accounted record's `z` and `divergent`.
// `denominators` holds each record's event denominator L, in record order.
// Against the accounted record with the lowest divergence, record i scores the
// pooled two-sample rate difference
//   z = (de_i - de_min) / sqrt(p (1 - p) (1 / L_i + 1 / L_min)),
//   p = (de_i L_i + de_min L_min) / (L_i + L_min),
// and is divergent above kDnaChainMapqDivergentRecordZ. With fewer than two
// accounted records nothing is flagged.
void dna_divergence_contrast(DnaFamilyDivergence& family,
                             const std::vector<double>& denominators);

} // namespace dna
} // namespace fa::cpu::lr

/**** pricing ****/

namespace fa::cpu::lr {
namespace dna {

// The MAPQ of the -c lane's records, each priced on its own region; `prices`
// is parallel to the realized family (index 0 the primary) and `rivals` holds
// the realized secondaries (DnaRegionOutcome::rivals).
void price_dna_regions(const DnaContext& dctx,
                       const DnaPlacementFamily& family,
                       const DnaPlacementChainingResult& placement,
                       const std::vector<DnaRegionPrice>& prices,
                       const std::vector<DnaRealizedRival>& rivals,
                       int read_len, dna::Result& realized);

// The MAPQ of map-only's records, each priced on the block part it prints from.
void price_dna_blocks(const DnaContext& dctx,
                      const DnaPlacementFamily& placement_family,
                      const DnaPlacementChainingResult& placement_chaining,
                      bool placement_chaining_ran,
                      ::fa::cpu::voting::CandidateId primary_candidate,
                      int read_len, dna::Result& realized);

} // namespace dna
} // namespace fa::cpu::lr

/**** backend ****/
// The DNA mapping backend. Internal header, not public API.

namespace fa {
namespace cpu {
namespace lr {

namespace dna {

// Maps one read: whole-read vote, bounded-candidate chaining, projection
// (map-only) or realization (CIGAR), then MAPQ and output assembly.
AlignResult map_read(const DnaContext &dctx, const LongReadSeedContext &seed_ctx,
                  DnaWorkerScratch &worker_scratch, const std::string &read,
                  const std::vector<uint8_t> &fwd_enc,
                  std::vector<uint8_t> &rc_enc,
                  const std::vector<QuerySeed> *shared_fwd_syncmer_seeds);

// Stateless backend tag for the templated engine; the RNA counterpart is
// rna::RnaBackend.
struct DnaBackend {
  using WorkerScratch = DnaWorkerScratch;
  using Context = DnaContext;

  static AlignResult map_read(const DnaContext &dctx,
                           const LongReadSeedContext &seed_ctx,
                           DnaWorkerScratch &worker_scratch,
                           const std::string &read,
                           const std::vector<uint8_t> &fwd_enc,
                           std::vector<uint8_t> &rc_enc,
                           const std::vector<QuerySeed> *shared_fwd_syncmer_seeds) {
    return ::fa::cpu::lr::dna::map_read(dctx, seed_ctx, worker_scratch, read,
                                        fwd_enc, rc_enc,
                                        shared_fwd_syncmer_seeds);
  }
};

} // namespace dna
} // namespace lr
} // namespace cpu
} // namespace fa
