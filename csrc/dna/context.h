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

namespace fa {
namespace cpu {
namespace lr {

// A post-commit terminal-clip record may own a single tile once it clears its
// absolute anchor and span floors; stable partition blocks need two.
inline constexpr int kDnaPostCommitRecordMinBlockTiles = 1;

// minimap2's opt->min_chain_score (-m); map-ont and map-hifi keep the default.
inline constexpr int kDnaMinChainScore = 40;

// The options the DNA stages read, copied from the typed options at the call
// boundary. Defaults are the shipped settings.
struct ResolvedDnaOptions {
  // Local interval-anchor harvest.
  int cigar_local_interval_anchor_interval_pad = 512;
  int cigar_local_interval_anchor_chain_max_gap = 5000;
  int cigar_local_diag_band = 20000;
  // Global occurrence cap: the screening pass's pool gate, the admissibility
  // bar inside an over-share query tile, residue recovery and terminal-clip
  // recovery. Set to the vote's resolved cap (INT_MAX when it has none); 200
  // is the vote's cap on a human index.
  int cigar_local_global_occ = 200;
  // Run count at or above which the dense chain uses the diagonal-keyed
  // predecessor search instead of the linear scan. Both are exact; 0 always
  // uses the diagonal search, -1 never does.
  int dna_dense_diag_min_runs = 2048;
  // Tandem half-window W in bp: an anchor whose key has another posting on
  // the same contig within W is flagged ANCHOR_TANDEM and is never a
  // realization corner. 0 = off.
  int dna_tandem_window = 1000;
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
  // minimap2 opt->zdrop_inv and opt->min_chain_score, for the local-inversion
  // probe (mm_test_zdrop) and the minimum size of an inversion middle.
  int cigar_dp_inversion_zdrop = 200;
  int cigar_dp_inversion_min_chain_score = kDnaMinChainScore;
  int cigar_dp_split_min_anchors = 3;  // minimap2 opt->min_cnt
  // minimap2 opt->min_ksw_len: a stretch piece ends once both spans reach it.
  int cigar_dp_min_ksw_len = 200;
  int cigar_dp_mismatch = 4;
  int cigar_dp_tail_end_bonus = -1;
  int cigar_dp_tail_zdrop = 400;       // -z
  float cigar_band_frac = 0.10f;
  bool enable_full_read_cigar = true;
  // Optional cs:Z / MD:Z output (minimap2 --cs / --MD); empty by default.
  ::fa::cpu::output::CigarReplayRequest cigar_replay_request;
  // Preset bounds for residue recovery from cached fine-seed postings.
  int residue_recovery_anchor_floor = 800;
  int residue_min_interval_bp = 200;
  int residue_min_anchor_density_per_100bp = 9;
  // Post-DP rescoring (postdp_scoring.h) supplies the MAPQ's dp1 / dp2. HiFi
  // presets only.
  bool postdp_rescoring = false;
  // Chain MAPQ, HiFi presets only: on the DP branch the raw ksw2 margin
  // replaces the ratio form.
  bool chain_mapq_hifi_margin = false;
  // Off by default: nominate a second locus from the read's unclaimed
  // terminal query.
  bool clip_nominate = false;
  // Chain MAPQ own-locus rules, both on: a chained shadow's vote counts for
  // the winner, and the chained shadow window grows with read length. No
  // option sets them.
  bool dna_chain_mapq_shadow_vote_credit = true;
  bool dna_chain_mapq_chain_shadow_read_slack = true;
  // Vote seeding.
  int chain_syncmer_s = 9;
  int chain_syncmer_downsample = 1;
  // Admit vote peaks by their ratio to the read's best vote; 0 admits by
  // count. See kCatalogueLaneBound.
  double vote_admission_ratio = 0.0;
  ::fa::cpu::voting::QueryPartitionParameters query_partition;
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
};

} // namespace lr
} // namespace cpu
} // namespace fa
