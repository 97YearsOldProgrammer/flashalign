// What the DNA stages depend on: non-owning reference views and a flat copy
// of the resolved options.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include DNA context"
#endif

#include "../core/cigar.h"              // output::CigarReplayRequest (cs/MD)
#include "../index/reference_context.h" // engine::ReferenceContext
#include "../voting/query_partition.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

// Projection and realization refuse a block with fewer supporting tiles; the
// partition's own blocks need two.
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
  int cigar_local_interval_anchor_chain_max_gap = 5000;
  // The chain passes' diagonal bands (-b): the dense (whole-query) passes'
  // and the screening pass's.
  int chain_band = 20000;
  int screen_diag_band = 20000;
  // Global occurrence cap: the screening pass's pool gate. Set to the vote's
  // resolved cap (INT_MAX when it has none); 200 is the vote's cap on a human
  // index.
  int cigar_local_global_occ = 200;
  // Run count at or above which the dense chain uses the diagonal-keyed
  // predecessor search instead of the linear scan. Both are exact; 0 always
  // uses the diagonal search, -1 never does.
  int dna_dense_diag_min_runs = 2048;
  // Occurrence gate of the whole-query exact pass (--max-chain-occ, by
  // default the vote's resolved cap); 0 = off.
  int dna_pool_gate_occ = 0;

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
  // minimap2 opt->pri_ratio (-p): the credibility ratio of the alternative
  // (alternative_hypothesis.h) and of the -c lane's secondary records
  // (region_realization.h).
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
  // in the fill row's units (family_realization.cpp).
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
  // (inv_local_chain.h).
  bool inversion_probe_local_gate = false;
  // The alternatives whose whole-query chains enter the ownership selection,
  // best first (placement_chaining.cpp).
  int alternative_realize_max = 1;
  // -N, minimap2's best_n: the -c lane's cap on secondary alignments
  // (region_realization.cpp).
  int secondary_max = 1;
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
  ::fa::cpu::voting::QueryPartitionParameters query_partition;
  // The query partition's tiles per read, kMinQueryTiles..kMaxQueryTiles.
  int query_tiles = ::fa::cpu::voting::kQueryTileCount;
  // The all-chains lane (options/dna_profile.h all_chains).
  bool all_chains = false;
  // The catalogue's candidates per strand (options/dna_profile.h
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
  // The all-chains lane with options/dna_profile.h skip_self: the contig that
  // is the query read itself, set by map_read. A chain takes no anchor on the
  // read's own exact diagonal there, as minimap2's -D. -1 otherwise.
  int self_contig = -1;
  // The all-chains lane under --dual=no (options/dna_profile.h dual): the
  // read's rank in the pair order and each contig's rank in it, set by the
  // engine (engine/aligner.h bind_dual_order). A candidate on a contig ranked
  // below the read is not chained. 0 and null otherwise.
  int dual_rank = 0;
  const std::vector<int>* dual_contig_rank = nullptr;
  // The read's seed density for the inversion probe's local chain, kept by
  // placement (DnaPlacementChainingResult::inversion_gate_seeds) and set by
  // map_read for every realization of the read; null when not kept.
  const RetainedSeedDensity* inversion_gate_seeds = nullptr;
};

} // namespace lr
} // namespace cpu
} // namespace fa
