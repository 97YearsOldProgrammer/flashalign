// Preset names, the DNA and RNA option profiles, and the inputs and output of option resolution.
#pragma once

#include "../seed/seed.h" // LongOccPolicy
#include "../vote/vote.h"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

/**** common_profile ****/
// Options shared by DNA and RNA mapping.

namespace fa {
namespace cpu {
namespace options {

// Seeding parameters. from_index means they were read from a prebuilt index;
// otherwise they come from the preset.
struct IndexIdentity {
  bool from_index = false;
  int k = 15;
  int syncmer_s = 9;
  int syncmer_downsample = 1;
};

// cs:Z form (minimap2 --cs[=short|long]); None emits no cs:Z.
enum class CsMode { None, Short, Long };

struct CommonOptions {
  int num_threads = 1;
  int min_support = 3;
  int max_query_seeds_per_strand = 128;
  bool enable_full_read_cigar = true;
  // minimap2 --cs / --MD; both imply a realized CIGAR.
  CsMode cs = CsMode::None;
  bool emit_md = false;
  // =/X CIGAR operators instead of M (minimap2 --eqx).
  bool emit_eqx = false;
};

} // namespace options
} // namespace cpu
} // namespace fa

/**** dna_profile ****/
// DNA mapping options. Seeding (k, s) and the options shared with RNA are in
// options/options.h.

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

// minimap2's min_ksw_len, in bases: the piece length of its gap-filling loop
// over unverified material, which dna/dna.h
// GeometryPlanner reproduces. A piece ends at the first retained anchor
// where both the query and the reference span since the previous corner reach
// it.
inline constexpr int kDnaMinKswLen = 200;

// The pool occurrence gate (--max-chain-occ). At N > 0 the whole-query pass
// drops every seed whose exact reference occurrence
// (KmerPostingIntervalView::global_count) exceeds N. Like minimap2's cap
// before chaining, it matters most on ONT, whose errors break exact runs so
// that an ungated pool is mostly short, noisy runs. The gate and the run collapse are complementary: the
// high-occurrence postings it removes are the ones that scatter across
// diagonals and collapse poorly. 0 is the ungated pass. Nothing else bounds
// the pool, so at 0 nothing does.
inline constexpr int kDnaPoolGateOcc = 0;

// The within-read repeat filter (--q-occ-frac), minimap2's q_occ_frac: a
// read drops every seed key it holds more than N times and in more than this
// fraction of its seeds (seed/seed.h drop_query_repeats). N is the lane's
// cap: the pool gate's on a DNA preset, the vote's on a splice preset. 0 is
// off.
inline constexpr double kQueryOccFrac = 0.01;

// minimap2's occ_dist (-e) and its seed selection, mm_seed_select: each streak
// of the read's seeds over the lane's N, between two seeds at most N or a read
// end, re-admits its rarest seeds past N, one per occ_dist bases of the streak
// and at most 128, a re-admitted seed over 4,095 dropped (seed/seed.h
// seed_select). Every occurrence gate of the read admits the re-admitted
// seeds. The splice presets install kOccDist; the DNA presets 0, which is off.
inline constexpr int kOccDist = 500;

// Every occurrence threshold of a lane is one number, N: the vote's seed
// occurrence cap as resolved against the index (seed/seed.h, engine/engine.cpp
// resolve_index_occ_cap_into):
//   lr                N = max(200, mm_idx_cal_max_occ(index, 1.81e-4))
//   lr:hq             N = min(500, max(200, mm_idx_cal_max_occ(index, 1.81e-4)))
//   splice, splice:hq N = max(200, mm_idx_cal_max_occ(index, 1.62e-3))
// where 500 is minimap2's map-hifi max_mid_occ. dna_chain_occ_thresholds()
// hands N to the pool gate; on the splice presets N is also the second pass's,
// the fine harvest's and the leftover step's cap (engine/engine.cpp
// compose_rna_runtime).
// --max-vote-occ INT sets N for every gate alike, and 0 removes every cap;
// --max-chain-occ INT replaces N on the pool gate only.
//
// kDnaPoolGateVoteCap is the dna_pool_gate_occ value that means N. It is
// outside --max-chain-occ's domain.
inline constexpr int kDnaPoolGateVoteCap = -1;
// Ceiling of the vote's cap rule (long_occ_ceiling); 0 is none.
inline constexpr int kDnaOccCeiling = 0;
inline constexpr int kDnaHiFiOccCeiling = 500;

// The vote's empty-tile rescue, M. A vote tile (seed/seed.h) whose found
// seeds are all over N votes with the rarest of them at or under M, ranked as a
// tile ranks its centred seed (occurrence, distance to the tile centre,
// read_pos), and that seed is also an exact-refine view. Per read and strand,
// in tile order, a rescue that would take the rescued occurrence past
// kDnaTileRescueBudget is skipped; reads shorter than kDnaTileRescueMinLen bp
// are not rescued (both constants in seed/seed.h). The whole-query pool
// admits exactly the seeds the vote rescued (dna/dna.h); the
// inversion probe keeps its gate. With no cap, or N >= M, nothing is rescued.
// lr's M is minimap2's max_max_occ, a constant there too; the RNA presets
// install 0.
inline constexpr int kDnaTileRescueOcc = 4095;
inline constexpr int kDnaHiFiTileRescueOcc = 1024;

// minimap2's seed-repetitiveness discount (uniq_ratio from rep_len, hit.c) is
// deliberately not applied: pen_s1 is the bare term (dna/mapq.cpp).

struct DnaLongOptions {
  // --vote-ratio R: nominate every vote peak with vote >= R * best, with the
  // lane bound (dna_chain_max_candidates) as a cost bound only; 0 admits by
  // count. DNA presets install kDnaProductionVoteAdmissionRatio.
  double vote_admission_ratio = 0.0;
  // The query tiles per read, the grid of the catalogue's masks and of the
  // read's blocks: the preset's, or -v's NUM2. The splice presets partition
  // into kQueryTileCount.
  int query_tiles = ::fa::cpu::voting::kQueryTileCount;
  // A read that is itself in the index, found by its exact name, is left out
  // of its own vote: a seed's occurrence does not count the key's postings on
  // the read's own contig (seed/seed.h vote_seed_occurrence), and the
  // read's own exact diagonal casts no vote, as minimap2's -D. In the
  // all-chains lane the chain takes no anchor on that diagonal either
  // (dna/dna.h DnaContext::self_contig); everything else after the vote
  // still reads the read's own postings.
  bool skip_self = false;
  // The all-chains lane, minimap2's -P for map-only PAF: the catalogue
  // candidates on each contig and strand share one whole-query chain call and
  // take its chains best first, one each, and each chain at -m or above is
  // printed once, MAPQ 0 and tp:A:S (dna/map.cpp emit_all_chains). No
  // tile mask, block ownership, -p, -N or MAPQ runs, so its lane bound
  // (max_cands) may reach kAllChainsLaneBound. Map-only.
  bool all_chains = false;
  // In the all-chains lane, whether a pair prints from both of its reads.
  // When false (--dual=no) the reads are ordered by length, then name, and a
  // candidate on a contig ordered before the read is not chained, so a pair
  // prints once, from its shorter read (engine/engine.cpp bind_dual_order,
  // dna/dna.h DnaContext::dual_rank).
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
  // (seed/seed.h). --max-vote-occ selects Fixed or Off instead.
  int long_occ_cap = 12;
  // Ceiling of the same rule, 0 = none: the cap is
  // min(ceiling, max(floor, quantile)), and a ceiling below the floor does
  // not apply, as minimap2's mm_mapopt_update clamps mid_occ.
  int long_occ_ceiling = kDnaOccCeiling;
  // The rule's quantile fraction f: kLongOccQuantileF, or kRnaOccQuantileF on
  // the splice presets.
  double long_occ_quantile_f = kLongOccQuantileF;
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
  // dp/dp.h bw_eff() and bw_long_eff() apply. On a splice preset
  // cigar_dp_bw is only the fine chain's band and cigar_dp_bw_long only
  // bounds it.
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
  // minimap2's -p (pri_ratio): the -c lane's secondary records are credible
  // at p times the owner's chain score, or within 2k of it.
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
  // The dense (whole-query) chain passes' diagonal band (-b).
  int chain_band = 20000;
  bool cigar_local_occ_mode_global = true;
  // Dense-chain search crossover (kDnaDenseDiagMinRuns); not settable.
  int dna_dense_diag_min_runs = kDnaDenseDiagMinRuns;
  // HiFi only (dna/dna.h): when the rival that owns dp2 has fewer
  // chain anchors than the winner and no other rival is within 5 % of it, MAPQ is
  // BWA-MEM's margin Phred on dp_max - dp_max2. The full guard is in
  // dna/dna.h.
  bool chain_mapq_hifi_margin = false;
  // HiFi only (dna/dna.h): the late inversion probe of a seam or
  // piece runs only where the read's opposite-lane fine seeds in the drop
  // window chain to kDnaInvLocalMinAnchors anchors.
  bool inversion_probe_local_gate = false;
  // Pool occurrence gate (--max-chain-occ); 0 is ungated. DNA presets
  // install kDnaPoolGateVoteCap.
  int dna_pool_gate_occ = kDnaPoolGateOcc;
  // The within-read repeat filter's fraction (--q-occ-frac); 0 is off.
  // Every preset but the overlap ones installs kQueryOccFrac.
  double query_occ_frac = 0.0;
  // minimap2's occ_dist (-e, kOccDist); 0 is off. The splice presets install
  // kOccDist, the DNA presets 0.
  int occ_dist = 0;
  // The vote's empty-tile rescue M (kDnaTileRescueOcc); 0 is none. DNA
  // presets install it; not settable.
  int dna_tile_rescue_occ = 0;
  // -N: the secondaries both lanes keep, map-only projected and the -c lane
  // realized; one prices the MAPQ.
  int secondary_max = 5;
  // A secondary-role region that splits, on a confirmed Z-drop or a fill over
  // the matrix cap, ends at the cut: no continuation. lr and lr:hq install it;
  // not settable.
  bool rival_split_ends = false;
};

// The lane bound: the vote's peaks and the catalogue's candidates per strand.
// --max-cands sets it; otherwise kCatalogueLaneBound under ratio admission,
// kCountAdmissionLaneBound under count admission.
int dna_chain_max_candidates(const DnaLongOptions& mapping) noexcept;

// Chain-side occurrence thresholds (dna/dna.h). pool_gate_occ is N, or
// the --max-chain-occ value; 0 is ungated.
struct DnaChainOccThresholds {
  int pool_gate_occ = 0;
};

DnaChainOccThresholds dna_chain_occ_thresholds(
    const DnaLongOptions& mapping) noexcept;

} // namespace lr
} // namespace cpu
} // namespace fa

/**** rna_profile ****/
// RNA (spliced) mapping options: the DNA options as `base`, plus intron,
// strand, rival and splice-scoring settings.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct RnaLongOptions {
  // Seeding, occurrence, chaining and DP scoring shared with DNA.
  DnaLongOptions base;
  // --min-intron / -G
  int min_intron = 20;
  int max_intron = 200000;
  int max_locus_chains = 6;
  int max_chain_predecessors = 64;
  // --vote-spacing S: a read of length L votes max(the seeds per strand,
  // ceil(L / S)) seeds per strand; 0 keeps the seeds per strand.
  int vote_seed_spacing = 0;
  // Rival realization (rna/rna.h). A rival is
  // retained when its chain score is within rival_pri_ratio of the primary's,
  // or within rival_min_diff of it absolutely (default 30, minimap2's 2k at
  // k = 15; not rescaled with the index's k).
  double rival_pri_ratio = 0.8;
  int rival_min_diff = 30;
  // At most this many rival realizations per read; -N sets it to the loci
  // kept minus one.
  int rival_realize_max = 5;
  // Splice-DP terms with no DNA analogue, as in minimap2's splice presets.
  int splice_transition = 0;
  int splice_junction_bonus = 9;
  // Read only by the scored-site mode, which no option selects.
  int splice_junction_penalty = 5;
  std::string junction_bed;
  int splice_inversion_zdrop = 100;
  // Transcript strand as rna::StrandMode: -1 unset, 0 auto, 1 forward,
  // 2 reverse, 3 none (no splice motif scored).
  int strand_mode = -1;
  // RNA MAPQ parameters (rna/hit.cpp), set by the preset only.
  //
  // Coverage threshold, in (0,1]: the penalty's coverage term is
  // clip((sel_qcov/read_len) / tau, 0, 1).
  double rna_mapq_qcov_tau = 0.7;
  // In [0,1]: a rival that fails the mask_level test adds
  // min(1, damp * rs_i/rs1) to the ambiguity term; 0 disables this.
  double rna_mapq_disjoint_vote_damp = 0.0;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa

/**** lane_stage ****/
// A stage of the mapping pipeline, which a preset's lane runs or leaves out
// (options.cpp lane_runs_stage). An option that configures a stage the lane
// leaves out is refused.

namespace fa {
namespace cpu {
namespace options {

enum class LaneStage {
  None,            // run by every lane
  BaseLevelOutput, // a requested CIGAR: SAM, or PAF with one
  BaseAlignment,   // base-level alignment: DP scoring, bands, CIGAR operators
  SamOutput,       // the shape of SAM output
  Selection,       // primary and secondary selection
  CandidateChains, // the catalogue candidates' whole-query chains
  OverlapPairs,    // which reads of an overlap pair print it
};

} // namespace options
} // namespace cpu
} // namespace fa

/**** resolved_options ****/
// The resolved mapping configuration. The variant holds the mode's options;
// RNA options embed the DNA ones in RnaLongOptions::base.

namespace fa::cpu::options {

struct ResolvedOptions {
  CommonOptions common;
  IndexIdentity index;
  std::variant<::fa::cpu::lr::DnaLongOptions,
               ::fa::cpu::lr::rna::RnaLongOptions>
      mapping;

  bool is_rna() const noexcept {
    return std::holds_alternative<
        ::fa::cpu::lr::rna::RnaLongOptions>(mapping);
  }

  ::fa::cpu::lr::DnaLongOptions& long_read() noexcept {
    if (auto* rna =
            std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping)) {
      return rna->base;
    }
    return std::get<::fa::cpu::lr::DnaLongOptions>(mapping);
  }

  const ::fa::cpu::lr::DnaLongOptions& long_read() const noexcept {
    if (const auto* rna =
            std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping)) {
      return rna->base;
    }
    return std::get<::fa::cpu::lr::DnaLongOptions>(mapping);
  }

  ::fa::cpu::lr::rna::RnaLongOptions* rna() noexcept {
    return std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping);
  }

  const ::fa::cpu::lr::rna::RnaLongOptions* rna() const noexcept {
    return std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping);
  }
};

}  // namespace fa::cpu::options

/**** types ****/
// Inputs and output of option resolution (options/options.h).

namespace fa {
namespace cpu {
namespace options {

enum class BackendMode { DnaLong, RnaSplice };

// (k, s) of a prebuilt index, which fix what the index contains. Without an
// index they come from the preset and any explicit request.
struct IndexMetadata {
  bool has_index = false;
  int k = 0;
  int syncmer_s = 0;
};

// Explicit overrides; std::nullopt keeps the value from the layer below.
struct UserOverrides {
  // With an index, k and s must match it; otherwise they override the
  // preset. downsample applies to queries only.
  std::optional<int> k;
  std::optional<int> syncmer_s;
  std::optional<int> syncmer_downsample;

  std::optional<int> min_support;
  // -m, minimap2's min_chain_score, >= 1. DNA presets only.
  std::optional<int> min_chain_score;
  // -v NUM1 and the Python Config.max_query_seeds, >= 0. DNA presets only.
  std::optional<int> max_query_seeds_per_strand;
  // -v NUM2, 1..kMaxQueryTiles. DNA presets only.
  std::optional<int> query_tiles;
  std::optional<int>
      num_threads; // the CLI passes the resolved count
  std::optional<bool> enable_full_read_cigar;
  // --cs, --MD and --eqx; apply wherever a CIGAR is realized.
  std::optional<CsMode> cs;
  std::optional<bool> emit_md;
  std::optional<bool> emit_eqx;
  std::optional<int> vote_diag_bin_width;
  // Adaptive vote-width coefficients (--dw-slope-den / --dw-max). Unlike
  // vote_diag_bin_width they keep the model on; den == 0 drops the length
  // term. DNA presets only.
  std::optional<int> vote_diag_slope_den;
  std::optional<int> vote_diag_width_max;
  std::optional<int> primary_occ_cap;
  std::optional<int> long_occ_cap;
  std::optional<std::string> occ_policy;

  // DP scoring (-A, -B, --score-N, -O, -E, -z, --end-bonus).
  std::optional<int> dp_match;
  std::optional<int> dp_mismatch;
  std::optional<int> dp_ambi;
  std::optional<int> dp_gap_open1;
  std::optional<int> dp_gap_open2;
  std::optional<int> dp_gap_extend1;
  std::optional<int> dp_gap_extend2;
  std::optional<int> dp_tail_zdrop;
  std::optional<int> dp_inversion_zdrop; // -z's second value
  std::optional<int> dp_tail_end_bonus;
  std::optional<int> dp_min_dp_max; // -S (minimap2 -s)
  // DP bandwidths (-r INT[,INT]) and maximum gap (-g); -g is also the chains'
  // gap. On a splice preset the first bandwidth is only the fine chain's band
  // and the second only bounds it; both apply after rna_max_intron, which
  // sets both.
  std::optional<int> dp_bw;
  std::optional<int> dp_bw_long;
  std::optional<int> dp_max_gap;
  // -b NUM: the dense chain passes' diagonal band, >= 1. DNA presets only.
  std::optional<int> chain_band;
  std::optional<int> rna_min_intron;
  std::optional<int> rna_max_intron;
  // RNA transcript strand (-u) as rna::StrandMode: 0 auto, 1 forward,
  // 2 reverse, 3 none (no splice motif scored).
  std::optional<int> rna_strand;
  std::optional<std::string> rna_junction_bed;
  std::optional<int> rna_junction_bonus;
  // -p: on a DNA preset the credibility ratio of the -c lane's secondary
  // records, on a splice preset the rival retention ratio.
  std::optional<double> pri_ratio;
  // RNA rival retention band (--rival-min-diff) and catalogue depth (-N + 1),
  // which also sets the realization budget to one fewer.
  std::optional<int> rna_rival_min_diff;
  std::optional<int> rna_max_loci; // >= 1
  // --vote-spacing, >= 0. RNA presets only.
  std::optional<int> rna_vote_seed_spacing;
  // -N on a DNA preset, minimap2's best_n: the cap on secondary alignments,
  // >= 1.
  std::optional<int> dna_secondary_max;
  // --vote-ratio: admit vote peaks by ratio to the read's best vote, in
  // [0,1]; 0 admits by count. DNA presets only.
  std::optional<double> dna_vote_admission_ratio;
  // --max-chain-occ: the dense chain's pool occurrence gate, >= 0; 0 is no
  // gate. Unset, the vote's cap. DNA presets only.
  std::optional<int> dna_pool_gate_occ;
  // --q-occ-frac: the within-read repeat filter's fraction, in [0,1]; 0 is
  // off. Every preset but the overlap ones.
  std::optional<double> query_occ_frac;
  // -e: minimap2's occ_dist (lr::kOccDist), >= 0; 0 is off. Every preset.
  std::optional<int> occ_dist;
  // --max-cands: the lane bound, the vote's peaks and the catalogue's
  // candidates per strand, 1..64, or 1..16383 in the all-chains lane. DNA
  // presets only.
  std::optional<int> max_cands;
  // --dual: whether an overlap preset prints a pair from both of its reads.
  // Overlap presets only.
  std::optional<bool> dual;
};

struct ResolveRequest {
  std::string preset;
  IndexMetadata index;
  UserOverrides user;
};

class ResolvedMapOptions; // forward
ResolvedMapOptions resolve_options(const ResolveRequest&);

// The immutable result of resolve_options().
class ResolvedMapOptions {
public:
  BackendMode mode() const {
    return resolved_.is_rna() ? BackendMode::RnaSplice : BackendMode::DnaLong;
  }
  bool is_rna() const { return resolved_.is_rna(); }
  const ResolvedOptions& resolved() const { return resolved_; }
  int k() const { return resolved_.index.k; }
  int syncmer_s() const { return resolved_.index.syncmer_s; }
  int syncmer_downsample() const { return resolved_.index.syncmer_downsample; }
  int min_support() const { return resolved_.common.min_support; }

private:
  friend ResolvedMapOptions resolve_options(const ResolveRequest&);
  ResolvedMapOptions() = default;
  ResolvedOptions resolved_{};
};

} // namespace options
} // namespace cpu
} // namespace fa

/**** resolve ****/
// The option resolver shared by the CLI and the Python binding. Layers, lowest
// first: defaults, the named preset, the index's (k, s), explicit overrides,
// derived rules, then validation. req.preset is already the final preset name.

namespace fa {
namespace cpu {
namespace options {

// Throws std::invalid_argument on an unknown preset or occurrence policy, an
// explicit k or s that contradicts a prebuilt index, or an invalid value.
ResolvedMapOptions resolve_options(const ResolveRequest &req);

} // namespace options
} // namespace cpu
} // namespace fa

/**** presets ****/
// Preset names and their seeding and mapping defaults, shared by index
// construction, option resolution and the Python binding.

namespace fa {
namespace cpu {
namespace options {

// asm5, asm10 and asm20: experimental.
bool is_assembly_preset(std::string_view preset);

// ava-ont and ava-hifi: experimental, map-only.
bool is_overlap_preset(std::string_view preset);

// Whether the lane `preset` selects runs `stage`. The all-chains lane of an
// overlap row prints map-only PAF and every chain, chains each contig and
// strand once, and alone chooses which reads of a pair print it.
bool lane_runs_stage(LaneStage stage, std::string_view preset);

bool is_rna_preset(std::string_view preset);
bool preset_is_valid(std::string_view preset);

// Preset (k, s) for building an index. A loaded index uses its own values.
struct PresetSeeding { int k; int syncmer_s; };
PresetSeeding resolve_preset_seeding(std::string_view preset);

// The batches in flight in the compute stage, --batch-window's default: 1 on
// the assembly presets, one batch at a time as in minimap2, 3 otherwise.
int preset_batch_window(std::string_view preset);

// Every accepted preset name, for error messages; the experimental ones last.
std::string accepted_preset_names();

} // namespace options
} // namespace cpu
} // namespace fa
