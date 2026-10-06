// DNA mapping options. Seeding (k, s) and the options shared with RNA are in
// options/common_profile.h.
#pragma once

#include "../seeding/context.h"            // LongOccPolicy
#include "../voting/query_partition.h"

#include <algorithm>
#include <limits>

namespace fa {
namespace cpu {
namespace lr {

// Vote admission ratio installed by the DNA presets: a whole-read vote peak is
// nominated when its vote is at least this fraction of the read's best.
inline constexpr double kDnaProductionVoteAdmissionRatio = 0.25;

// Collapsed-run count at or above which the dense chain's predecessor search
// is the exact branch-and-bound search instead of the linear scan
// (chaining::kDenseDiagMinRuns, restated to keep chaining out of this header).
// Both return the recurrence's optimum; branch-and-bound is cheaper on large
// pools. 0 sends every pool to branch-and-bound, -1 every pool to the scan.
inline constexpr int kDnaDenseDiagMinRuns = 2048;

// Tandem release half-window W, in bp. A DNA anchor whose seed key has another
// posting on the same chromosome within +-W is flagged ANCHOR_TANDEM
// (dna/placement_chaining.cpp append_interval_anchors), and the realizer skips
// it as a corner, as minimap2 skips MM_SEED_TANDEM seeds. Otherwise an anchor
// on another copy of a tandem array would pin the DP to the wrong copy.
// 0 disables the rule.
inline constexpr int kDnaTandemWindow = 1000;

// minimap2's min_ksw_len, in bases: the piece length of its gap-filling loop
// over unverified material, which dna/ordered_anchor_path.h
// plan_verified_geometry reproduces. A piece ends at the first retained anchor
// where both the query and the reference span since the previous corner reach
// it.
inline constexpr int kDnaMinKswLen = 200;

// The pool occurrence gate (--max-chain-occ). At N > 0 the whole-query pass
// drops every seed whose exact reference occurrence
// (KmerPostingIntervalView::global_count) exceeds N, the same test the
// screening pass applies. Like minimap2's cap before chaining, it matters most
// on ONT, whose errors break exact runs so that an ungated pool is mostly
// short, noisy runs. The gate and the run collapse are complementary: the
// high-occurrence postings it removes are the ones that scatter across
// diagonals and collapse poorly. 0 is the ungated pass. Nothing else bounds
// the pool, so at 0 nothing does.
inline constexpr int kDnaPoolGateOcc = 0;

// Every DNA occurrence threshold is one number, N: the vote's seed occurrence
// cap as resolved against the index (seeding/context.h, engine/aligner.h
// resolve_index_occ_cap_into):
//   lr     N = max(200, mm_idx_cal_max_occ(index, 1.81e-4))
//   lr:hq  N = min(500, max(200, mm_idx_cal_max_occ(index, 1.81e-4)))
// where 500 is minimap2's map-hifi max_mid_occ. dna_chain_occ_thresholds()
// hands N to the pool gate and to the global cap read by the screening pass.
// --max-vote-occ INT sets N for the vote and the chain alike, and 0 removes
// every cap; --max-chain-occ INT replaces N on the pool gate only.
//
// kDnaPoolGateVoteCap is the dna_pool_gate_occ value that means N. It is
// outside --max-chain-occ's domain.
inline constexpr int kDnaPoolGateVoteCap = -1;
// Ceiling of the vote's cap rule (long_occ_ceiling); 0 is none.
inline constexpr int kDnaOccCeiling = 0;
inline constexpr int kDnaHiFiOccCeiling = 500;

// The vote's empty-tile rescue, M. A vote tile (seeding/syncmer.h) whose found
// seeds are all over N votes with the rarest of them at or under M, ranked as a
// tile ranks its centred seed (occurrence, distance to the tile centre,
// read_pos), and that seed is also an exact-refine view. Per read and strand,
// in tile order, a rescue that would take the rescued occurrence past
// kDnaTileRescueBudget is skipped; reads shorter than kDnaTileRescueMinLen bp
// are not rescued (both constants in seeding/syncmer.h). Placement's screening
// pass admits a rescued seed past its gate, and the whole-query pool admits
// exactly the seeds the vote rescued (dna/retained_seed_density.h); the
// inversion probe keeps its gate. With no cap, or N >= M, nothing is rescued.
// lr's M is minimap2's max_max_occ, a constant there too; the RNA presets
// install 0.
inline constexpr int kDnaTileRescueOcc = 4095;
inline constexpr int kDnaHiFiTileRescueOcc = 1024;

// minimap2's seed-repetitiveness discount (uniq_ratio from rep_len, hit.c) is
// deliberately not applied: pen_s1 is the bare term (dna/chain_mapq.cpp).

struct DnaLongOptions {
  // --vote-ratio R: nominate every vote peak with vote >= R * best, with the
  // lane bound (dna_chain_max_candidates) as a cost bound only; 0 admits by
  // count. DNA presets install kDnaProductionVoteAdmissionRatio.
  double vote_admission_ratio = 0.0;
  // Objective of the query partition (--tile-score), per tile.
  ::fa::cpu::voting::QueryPartitionParameters query_partition;
  // The query partition's tiles per read (--tiles),
  // kMinQueryTiles..kMaxQueryTiles. DNA presets only; the splice presets
  // partition into kQueryTileCount.
  int query_tiles = ::fa::cpu::voting::kQueryTileCount;
  // A read that is itself in the index, found by its exact name, is left out
  // of its own vote: a seed's occurrence does not count the key's postings on
  // the read's own contig (seeding/context.h vote_seed_occurrence), and the
  // read's own exact diagonal casts no vote, as minimap2's -D. In the
  // all-chains lane the chain takes no anchor on that diagonal either
  // (dna/context.h DnaContext::self_contig); everything else after the vote
  // still reads the read's own postings.
  bool skip_self = false;
  // The all-chains lane, minimap2's -P for map-only PAF: the catalogue
  // candidates on each contig and strand share one whole-query chain call and
  // take its chains best first, one each, and each chain at -m or above is
  // printed once, MAPQ 0 and tp:A:S (dna/backend.cpp emit_all_chains). No
  // tile mask, partition, -p, -N or MAPQ runs, so its lane bound (max_cands)
  // may reach kAllChainsLaneBound. Map-only.
  bool all_chains = false;
  // In the all-chains lane, whether a pair prints from both of its reads.
  // When false (--dual=no) the reads are ordered by length, then name, and a
  // candidate on a contig ordered before the read is not chained, so a pair
  // prints once, from its shorter read (engine/aligner.h bind_dual_order,
  // dna/context.h DnaContext::dual_rank).
  bool dual = true;
  // --max-cands: the lane bound, 1..kMaxCatalogueLaneBound, or
  // 1..kAllChainsLaneBound in the all-chains lane; 0 leaves it to the
  // admission rule (dna_chain_max_candidates).
  int max_cands = 0;
  int vote_diag_bin_width = 64;
  // DNA presets widen the vote with read length; --dw turns this off.
  bool vote_diag_bin_width_adaptive = false;
  //   W(L) = clamp(L / vote_diag_slope_den, vote_diag_bin_width,
  //                vote_diag_width_max)
  // Defaults are lr's. A denominator of 0 pins W at the base width.
  int vote_diag_slope_den = 128;
  int vote_diag_width_max = 2048;
  int chain_max_candidates_per_window = 2;
  int long_primary_occ_cap = 0;
  // Under LongOccPolicy::Platform, the floor of the vote's cap: the engine
  // raises it to the index's mid-occurrence quantile when that is larger
  // (seeding/context.h). --max-vote-occ selects Fixed or Off instead.
  int long_occ_cap = 12;
  // Ceiling of the same rule, 0 = none: the cap is
  // min(ceiling, max(floor, quantile)), and a ceiling below the floor does
  // not apply, as minimap2's mm_mapopt_update clamps mid_occ.
  int long_occ_ceiling = kDnaOccCeiling;
  LongOccPolicy long_occ_policy = LongOccPolicy::Platform;
  // How the engine resolved the cap against the index; reported only.
  LongOccIndexResolution long_occ_index_resolution;
  // DP scoring as minimap2's; defaults are map-ont's.
  int cigar_dp_match = 2;
  int cigar_dp_mismatch = 4;
  int cigar_dp_ambi = 1;
  int cigar_dp_gap_open1 = 4;
  int cigar_dp_gap_extend1 = 2;
  int cigar_dp_gap_open2 = 24;
  int cigar_dp_gap_extend2 = 1;
  int cigar_dp_tail_zdrop = 400; // minimap2 -z
  int cigar_dp_tail_end_bonus = -1;
  // Bandwidths (minimap2 -r) before minimap2's 1.5x inflation, which
  // dp/params.h bw_eff() and bw_long_eff() apply. On a DNA preset
  // cigar_dp_bw_long, uninflated, is also the dense chains' diagonal band.
  // On a splice preset cigar_dp_bw is only the fine chain's band and
  // cigar_dp_bw_long only bounds it.
  int cigar_dp_bw = 500;
  int cigar_dp_bw_long = 20000;
  int cigar_dp_max_gap = 5000;
  // minimap2 -s: the DP half of the DNA per-record emission floor
  // (mm_filter_regs) and the inversion-segment gate. 0 disables.
  int cigar_dp_min_dp_max = 80;
  // minimap2's zdrop_inv; on a DNA preset the end row's, fixed.
  int cigar_dp_inversion_zdrop = 200;
  // DNA presets: the gap-fill row, which -A -B -O -E -z --score-N set; the
  // cigar_dp_* row above is then the preset's end row, which prices every
  // path. fill_dp_min_dp_max is -S scaled by fill_dp_match / cigar_dp_match.
  // Defaults are lr's.
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
  // minimap2's -p (pri_ratio): the alternative and the -c lane's secondary
  // records are credible at p times the owner's chain score, or within 2k of
  // it.
  double pri_ratio = 0.8;
  // minimap2's -m (min_chain_score): the DNA emission floor, the MAPQ's subsc
  // floor and the minimum size of an inversion middle. DNA presets only.
  int min_chain_score = 40;
  // As minimap2, a Z-dropped region is split and continues while at least
  // this many anchors remain.
  int cigar_dp_split_min_anchors = 3;  // minimap2 opt->min_cnt
  // minimap2 min_ksw_len; 0 on RNA presets.
  int cigar_dp_min_ksw_len = kDnaMinKswLen;
  int cigar_local_interval_anchor_occ_cap = 32;
  int cigar_local_interval_anchor_interval_pad = 512;
  // The chains' maximum gap (-g); on a splice preset the fine chain's query
  // gap.
  int cigar_local_interval_anchor_chain_max_gap = 5000;
  // The DNA screening chain's diagonal band.
  int screen_diag_band = 20000;
  bool cigar_local_occ_mode_global = true;
  // Global occurrence cap: seeds whose exact reference occurrence exceeds it
  // are dropped. Only the RNA path reads this field (compose_rna_runtime); the
  // DNA path uses N from dna_chain_occ_thresholds() instead.
  int cigar_local_global_occ = 200;
  // Dense-chain search crossover (kDnaDenseDiagMinRuns); not settable.
  int dna_dense_diag_min_runs = kDnaDenseDiagMinRuns;
  // Tandem release half-window (kDnaTandemWindow); 0 on RNA presets.
  int dna_tandem_window = kDnaTandemWindow;
  // HiFi only (dna/chain_mapq.h): when the rival that owns dp2 has fewer
  // chain anchors than the winner and no other rival is within 5 % of it, MAPQ is
  // BWA-MEM's margin Phred on dp_max - dp_max2. The full guard is in
  // chain_mapq.h.
  bool chain_mapq_hifi_margin = false;
  // HiFi only (dna/inv_local_chain.h): the late inversion probe of a seam or
  // piece runs only where the read's opposite-lane fine seeds in the drop
  // window chain to kDnaInvLocalMinAnchors anchors.
  bool inversion_probe_local_gate = false;
  // Pool occurrence gate (--max-chain-occ); 0 is ungated. DNA presets
  // install kDnaPoolGateVoteCap.
  int dna_pool_gate_occ = kDnaPoolGateOcc;
  // The vote's empty-tile rescue M (kDnaTileRescueOcc); 0 is none. DNA
  // presets install it; not settable.
  int dna_tile_rescue_occ = 0;
  // -N: the alternatives whose whole-query chains enter the ownership
  // selection (dna/placement_chaining.cpp).
  int alternative_realize_max = 1;
};

// The lane bound: the vote's peaks and the catalogue's candidates per strand.
// --max-cands sets it; otherwise kCatalogueLaneBound under ratio admission,
// kCountAdmissionLaneBound under count admission.
inline int dna_chain_max_candidates(const DnaLongOptions& mapping) noexcept {
  if (mapping.max_cands > 0)
    return mapping.max_cands;
  return std::max(mapping.chain_max_candidates_per_window,
                  mapping.vote_admission_ratio > 0.0
                      ? ::fa::cpu::voting::kCatalogueLaneBound
                      : ::fa::cpu::voting::kCountAdmissionLaneBound);
}

// N, the vote's occurrence cap for this run; 0 is no cap. Before an index is
// attached, the Platform rule gives the preset floor.
inline int dna_resolved_vote_occ_cap(const DnaLongOptions& mapping) noexcept {
  LongOccPolicyConfig policy;
  policy.policy = mapping.long_occ_policy;
  policy.primary_occ_cap = mapping.long_primary_occ_cap;
  policy.platform_occ_cap = mapping.long_occ_cap;
  return effective_long_primary_occ_cap(policy);
}

// Chain-side occurrence thresholds (dna/context.h). pool_gate_occ is N, or
// the --max-chain-occ value; 0 is ungated. global_occ is N, or INT_MAX when
// the vote has no cap.
struct DnaChainOccThresholds {
  int pool_gate_occ = 0;
  int global_occ = 0;
};

inline DnaChainOccThresholds dna_chain_occ_thresholds(
    const DnaLongOptions& mapping) noexcept {
  const int vote_cap = dna_resolved_vote_occ_cap(mapping);
  DnaChainOccThresholds thresholds;
  thresholds.pool_gate_occ = mapping.dna_pool_gate_occ == kDnaPoolGateVoteCap
                                 ? vote_cap
                                 : mapping.dna_pool_gate_occ;
  thresholds.global_occ =
      vote_cap > 0 ? vote_cap : std::numeric_limits<int>::max();
  return thresholds;
}

} // namespace lr
} // namespace cpu
} // namespace fa
