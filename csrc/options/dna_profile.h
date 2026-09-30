// DNA mapping options. Seeding (k, s) and the options shared with RNA are in
// options/common_profile.h.
#pragma once

#include "../seeding/context.h"            // LongOccPolicy
#include "../voting/query_partition.h"

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

// Dense-chain pool admission. Postings are harvested under a per-pass budget
// (dna/placement_chaining.h kDnaSkipPoolBudget) shared over the query tiles
// that hold them; a tile over its share keeps only slices at or under the
// global occurrence cap, rarest first, and the chain crosses a tile that keeps
// none as a gap. The budget bounds work rather than evidence: any occurrence
// threshold on runs low enough to shrink a satellite pool also cuts the unique
// flanks a placement rests on.
//
// The pool occurrence gate (--max-chain-occ). At N > 0 the whole-query pass
// drops every seed whose exact reference occurrence
// (KmerPostingIntervalView::global_count) exceeds N, the same test the
// screening pass applies. Like minimap2's cap before chaining, it matters most
// on ONT, whose errors break exact runs so that an ungated pool is mostly
// short, noisy runs. The gate and the run collapse are complementary: the
// high-occurrence postings it removes are the ones that scatter across
// diagonals and collapse poorly. 0 is the ungated pass.
inline constexpr int kDnaPoolGateOcc = 0;

// Every DNA occurrence threshold is one number, N: the vote's seed occurrence
// cap as resolved against the index (seeding/context.h, engine/aligner.h
// resolve_index_occ_cap_into):
//   lr     N = max(200, mm_idx_cal_max_occ(index, 1.81e-4))
//   lr:hq  N = min(500, max(200, mm_idx_cal_max_occ(index, 1.81e-4)))
// where 500 is minimap2's map-hifi max_mid_occ. dna_chain_occ_thresholds()
// hands N to the pool gate (and through it the terminal clip nomination's
// ceiling, min(4095, gate)) and to the global cap read by the screening pass,
// pool admission, residue recovery and terminal-clip recovery.
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
// exactly the seeds the vote rescued (dna/retained_seed_density.h); stage 1,
// the inversion probe and the clip nomination keep their gates. With no cap, or
// N >= M, nothing is rescued. lr's M is minimap2's max_max_occ, a constant
// there too; the RNA presets install 0.
inline constexpr int kDnaTileRescueOcc = 4095;
inline constexpr int kDnaHiFiTileRescueOcc = 1024;

// minimap2's seed-repetitiveness discount (uniq_ratio from rep_len, hit.c) is
// deliberately not applied: pen_s1 is the bare term (dna/chain_mapq.cpp).

struct DnaLongOptions {
  // --vote-ratio R: nominate every vote peak with vote >= R * best, with
  // kCatalogueLaneBound as a cost bound only; 0 admits by count. DNA presets
  // install kDnaProductionVoteAdmissionRatio.
  double vote_admission_ratio = 0.0;
  // Objective of the 128-tile query partition (--tile-score).
  ::fa::cpu::voting::QueryPartitionParameters query_partition;
  float cigar_band_frac = 0.10f;
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
  // dp/params.h bw_eff() and bw_long_eff() apply.
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
  // minimap2's min_chain_score default; not settable.
  int cigar_dp_inversion_min_chain_score = 40;
  // As minimap2, a Z-dropped region is split and continues while at least
  // this many anchors remain.
  int cigar_dp_split_min_anchors = 3;  // minimap2 opt->min_cnt
  // minimap2 min_ksw_len; 0 on RNA presets.
  int cigar_dp_min_ksw_len = kDnaMinKswLen;
  int cigar_local_interval_anchor_occ_cap = 32;
  int cigar_local_interval_anchor_interval_pad = 512;
  int cigar_local_interval_anchor_chain_max_gap = 5000;
  int cigar_local_diag_band = 20000;
  bool cigar_local_occ_mode_global = true;
  // Global occurrence cap: seeds whose exact reference occurrence exceeds it
  // are dropped. Only the RNA path reads this field (compose_rna_runtime); the
  // DNA path uses N from dna_chain_occ_thresholds() instead.
  int cigar_local_global_occ = 200;
  // Dense-chain search crossover (kDnaDenseDiagMinRuns); not settable.
  int dna_dense_diag_min_runs = kDnaDenseDiagMinRuns;
  // Tandem release half-window (kDnaTandemWindow); 0 on RNA presets.
  int dna_tandem_window = kDnaTandemWindow;
  // Residue recovery does not run below this many chain anchors.
  int residue_recovery_anchor_floor = 800;
  // Minimum query length in base pairs of an owned-but-unsupported run before
  // bounded cached-evidence recovery may consider it.
  int residue_min_interval_bp = 200;
  // Minimum anchor density of an admitted chain, in anchors per 100 query bp:
  // chain_anchors * 100 >= D * chain_query_span.
  int residue_min_anchor_density_per_100bp = 9;
  // HiFi only (dna/postdp_scoring.h): the MAPQ's dp1/dp2 come from the
  // log-gap CIGAR sweep rescored by minibwa's b2 formula, with its clip term
  // when the two primary records' query spans compete, instead of the raw DP
  // scores.
  bool postdp_rescoring = false;
  // HiFi only (dna/chain_mapq.h): when the realized rival has fewer chain
  // anchors than the winner and no other rival is within 5 % of it, MAPQ is
  // BWA-MEM's margin Phred on the raw ksw2 margin. The full guard is in
  // chain_mapq.h.
  bool chain_mapq_hifi_margin = false;
  // HiFi only (dna/inv_local_chain.h): the late inversion probe of a seam or
  // piece runs only where the read's opposite-lane fine seeds in the drop
  // window chain to kDnaInvLocalMinAnchors anchors. Bridges keep their probe.
  bool inversion_probe_local_gate = false;
  // Terminal clip nomination: a terminal clip of at least
  // kDnaClipNominateMinIntervalBp unclaimed query gets a second recovery
  // attempt after the production one declines, with an interval-scaled
  // anchor floor, no absolute density bar, a wider per-strand cluster budget,
  // a higher occurrence ceiling when the production cap finds nothing, and
  // cross-contig, cross-strand geometry. It recovers split records whose
  // second locus the whole-read vote never nominates, such as the far
  // breakend of a translocation. A read still gains at most
  // kDnaResidueMaxAdmissionsPerRead records.
  bool dna_clip_nominate = true;
  // Pool occurrence gate (--max-chain-occ); 0 is ungated. DNA presets
  // install kDnaPoolGateVoteCap.
  int dna_pool_gate_occ = kDnaPoolGateOcc;
  // The vote's empty-tile rescue M (kDnaTileRescueOcc); 0 is none. DNA
  // presets install it; not settable.
  int dna_tile_rescue_occ = 0;
};

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
