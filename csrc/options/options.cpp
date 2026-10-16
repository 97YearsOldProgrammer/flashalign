// The presets' defaults, and the resolver the CLI and the Python binding share: defaults, preset, index-derived k and s, explicit overrides, derived rules, validation.
#include "options.h"

#include "../index/format.h"
#include "../seed/seed.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

/**** dna_profile ****/

namespace fa {
namespace cpu {
namespace lr {

int dna_chain_max_candidates(const DnaLongOptions& mapping) noexcept {
  if (mapping.max_cands > 0)
    return mapping.max_cands;
  return std::max(mapping.chain_max_candidates_per_window,
                  mapping.vote_admission_ratio > 0.0
                      ? ::fa::cpu::voting::kCatalogueLaneBound
                      : ::fa::cpu::voting::kCountAdmissionLaneBound);
}

// N, the vote's occurrence cap for this run; 0 is no cap. Before an index is
// attached, the Platform rule gives the preset floor.
static int dna_resolved_vote_occ_cap(const DnaLongOptions& mapping) noexcept {
  LongOccPolicyConfig policy;
  policy.policy = mapping.long_occ_policy;
  policy.primary_occ_cap = mapping.long_primary_occ_cap;
  policy.platform_occ_cap = mapping.long_occ_cap;
  return effective_long_primary_occ_cap(policy);
}

DnaChainOccThresholds dna_chain_occ_thresholds(
    const DnaLongOptions& mapping) noexcept {
  const int vote_cap = dna_resolved_vote_occ_cap(mapping);
  DnaChainOccThresholds thresholds;
  thresholds.pool_gate_occ = mapping.dna_pool_gate_occ == kDnaPoolGateVoteCap
                                 ? vote_cap
                                 : mapping.dna_pool_gate_occ;
  return thresholds;
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** presets ****/

namespace fa {
namespace cpu {
namespace options {

// Asm marks the experimental assembly-to-reference rows.
enum class DnaPresetKind { Ont, HiFi, Asm };

// The placement values: the vote's seeds per strand, which --vote-seeds
// replaces, and the query tiles per read, which no option sets.
struct DnaPlacementDefaults {
  int vote_seeds;
  int query_tiles;
};

constexpr DnaPlacementDefaults kReadPlacement{
    128, ::fa::cpu::voting::kQueryTileCount};
// Above 128 vote seeds the nested sampler applies (seed/seed.h).
constexpr DnaPlacementDefaults kAssemblyPlacement{4096, 2048};
constexpr DnaPlacementDefaults kOntOverlapPlacement{
    256, ::fa::cpu::voting::kQueryTileCount};

// An all-vs-all overlap row turns on both of its mechanisms (options.h
// skip_self, all_chains) and replaces --max-cands and --vote-ratio, and -m and
// --min-support where positive. Off on every other row.
struct DnaOverlapDefaults {
  bool on;
  int max_cands;
  double vote_ratio;
  int min_chain_score;
  int min_support;
  // --dual's default (options.h dual).
  bool dual;
};

constexpr DnaOverlapDefaults kNoOverlap{false, 0, 0.0, 0, 0, true};
// The all-chains lane's widest bound and count admission, minimap2's ava -m100
// and --dual=no.
constexpr DnaOverlapDefaults kHiFiOverlap{
    true, ::fa::cpu::voting::kAllChainsLaneBound, 0.0, 100, 0, false};
// The same with --min-support 2.
constexpr DnaOverlapDefaults kOntOverlap{
    true, ::fa::cpu::voting::kAllChainsLaneBound, 0.0, 100, 2, false};

struct DnaPresetProfile {
  std::string_view name;
  DnaPresetKind kind;
  int k;
  int syncmer_s;
  int vote_diag_bin_width;
  // Adaptive vote width: W(L) = clamp(L / vote_diag_slope_den,
  // vote_diag_bin_width, vote_diag_width_max), in integers. A denominator of
  // 0 or less drops the length term.
  int vote_diag_slope_den;
  int vote_diag_width_max;
  int cigar_dp_match;
  int cigar_dp_mismatch;
  int cigar_dp_ambi;
  int cigar_dp_gap_open1;
  int cigar_dp_gap_extend1;
  int cigar_dp_gap_open2;
  int cigar_dp_gap_extend2;
  int cigar_dp_tail_zdrop;
  int cigar_dp_tail_end_bonus;
  int cigar_dp_bw;
  int cigar_dp_bw_long;
  int cigar_dp_max_gap;
  int cigar_dp_min_dp_max;
  // The gap-fill row -A -B -O -E -z --score-N set; the cigar_dp_* row above
  // is the end row, which prices every path. On an Asm row the two are equal
  // and the letters set both (options/options.cpp).
  int fill_dp_match;
  int fill_dp_mismatch;
  int fill_dp_ambi;
  int fill_dp_gap_open1;
  int fill_dp_gap_extend1;
  int fill_dp_gap_open2;
  int fill_dp_gap_extend2;
  int fill_dp_tail_zdrop;
  int fill_dp_inversion_zdrop;
  // The dense chain's maximum gap; the splice presets borrow it as the fine
  // chain's query gap. -g sets it.
  int chain_max_gap;
  // The dense (whole-query) chain passes' diagonal band; -b sets it.
  int chain_band;
  DnaPlacementDefaults placement;
  DnaOverlapDefaults overlap;
};

constexpr std::array<DnaPresetProfile, 7> kDnaPresetProfiles{{
    // End row as minimap2: lr is map-ont (-A2 -B4 -O4,24 -E2,1), lr:hq is
    // map-hifi (-A1 -B4 -O6,26 -E2,1). The fill rows are gentler on gaps.
    {"lr", DnaPresetKind::Ont, 21, 9, 64, 128, 2048,
     2, 4, 1, 4, 2, 24, 1, 400, -1, 500, 20000,
     5000, 80,
     4, 8, 2, 8, 4, 48, 1, 800, 200,
     20000, 20000, kReadPlacement, kNoOverlap},
    {"lr:hq", DnaPresetKind::HiFi, 21, 5, 48, 64, 2048,
     1, 4, 1, 6, 2, 26, 1, 400, -1, 500, 20000,
     10000, 200,
     3, 12, 3, 18, 6, 78, 1, 1200, 600,
     10000, 20000, kReadPlacement, kNoOverlap},
    // minimap2's asm5, asm10 and asm20 scoring with its -r1k,100k -g10k -s200
    // -z200, one row for the fills and the read ends, and its 100k bandwidth
    // as the dense chains' band; lr's seeding, vote width and chain gap.
    {"asm5", DnaPresetKind::Asm, 21, 9, 64, 128, 2048,
     1, 19, 1, 39, 3, 81, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 19, 1, 39, 3, 81, 1, 200, 200,
     20000, 100000, kAssemblyPlacement, kNoOverlap},
    {"asm10", DnaPresetKind::Asm, 21, 9, 64, 128, 2048,
     1, 9, 1, 16, 2, 41, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 9, 1, 16, 2, 41, 1, 200, 200,
     20000, 100000, kAssemblyPlacement, kNoOverlap},
    {"asm20", DnaPresetKind::Asm, 21, 9, 64, 128, 2048,
     1, 4, 1, 6, 2, 26, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 4, 1, 6, 2, 26, 1, 200, 200,
     20000, 100000, kAssemblyPlacement, kNoOverlap},
    // Experimental all-vs-all read overlap, map-only: lr's row at k17 with
    // the vote width L/64, and lr:hq's row, as overlap rows. They follow the
    // row of their kind, which dna_preset_profile returns.
    {"ava-ont", DnaPresetKind::Ont, 17, 9, 64, 64, 2048,
     2, 4, 1, 4, 2, 24, 1, 400, -1, 500, 20000,
     5000, 80,
     4, 8, 2, 8, 4, 48, 1, 800, 200,
     20000, 20000, kOntOverlapPlacement, kOntOverlap},
    {"ava-hifi", DnaPresetKind::HiFi, 21, 5, 48, 64, 2048,
     1, 4, 1, 6, 2, 26, 1, 400, -1, 500, 20000,
     10000, 200,
     3, 12, 3, 18, 6, 78, 1, 1200, 600,
     10000, 20000, kReadPlacement, kHiFiOverlap},
}};

static const DnaPresetProfile* find_dna_preset_profile(
    std::string_view preset) {
  for (const auto& profile : kDnaPresetProfiles) {
    if (profile.name == preset) return &profile;
  }
  return nullptr;
}

static const DnaPresetProfile& dna_preset_profile(DnaPresetKind kind) {
  for (const auto& profile : kDnaPresetProfiles) {
    if (profile.kind == kind) return profile;
  }
  throw std::logic_error("missing DNA preset profile");
}

// Installs a DNA preset's mapping defaults. (k, s) are resolved separately.
static void set_dna_long_platform_fields(
    lr::DnaLongOptions& mapping, CommonOptions& common,
    const DnaPresetProfile& profile) {
    mapping.vote_diag_bin_width = profile.vote_diag_bin_width;
    mapping.vote_diag_bin_width_adaptive = true;
    mapping.vote_diag_slope_den = profile.vote_diag_slope_den;
    mapping.vote_diag_width_max = profile.vote_diag_width_max;
    // The floor of the vote occurrence cap. On a DNA run the engine raises
    // it to the index's mid-occurrence quantile when larger
    // (seed/seed.h), up to long_occ_ceiling.
    mapping.long_occ_cap = 200;
    mapping.query_tiles = profile.placement.query_tiles;
    mapping.vote_admission_ratio = lr::kDnaProductionVoteAdmissionRatio;
    common.min_support = 3;
    if (profile.overlap.on) {
      mapping.skip_self = true;
      mapping.all_chains = true;
      mapping.max_cands = profile.overlap.max_cands;
      mapping.vote_admission_ratio = profile.overlap.vote_ratio;
      if (profile.overlap.min_chain_score > 0)
        mapping.min_chain_score = profile.overlap.min_chain_score;
      if (profile.overlap.min_support > 0)
        common.min_support = profile.overlap.min_support;
      mapping.dual = profile.overlap.dual;
    }
    // minimap2's min_ksw_len, the piece length of its gap-filling loop.
    mapping.cigar_dp_min_ksw_len = lr::kDnaMinKswLen;
    const bool hifi_row = profile.kind == DnaPresetKind::HiFi;
    // The vote's cap N is also the dense chain's pool gate
    // (engine/engine.cpp make_dna_context). lr:hq caps N at minimap2's
    // map-hifi 500. --max-chain-occ replaces N on the gate only.
    mapping.long_occ_ceiling =
        hifi_row ? lr::kDnaHiFiOccCeiling : lr::kDnaOccCeiling;
    mapping.dna_pool_gate_occ = lr::kDnaPoolGateVoteCap;
    // The within-read repeat filter; the all-chains lane keeps every seed.
    mapping.query_occ_frac =
        profile.overlap.on ? 0.0 : lr::kQueryOccFrac;
    // The vote's empty-tile rescue past N, up to M (seed/seed.h).
    mapping.dna_tile_rescue_occ =
        hifi_row ? lr::kDnaHiFiTileRescueOcc : lr::kDnaTileRescueOcc;
    // HiFi only: the chain MAPQ margin rule (dna/dna.h).
    mapping.chain_mapq_hifi_margin = profile.kind == DnaPresetKind::HiFi;
    // HiFi only: the late inversion probe's local-chain gate
    // (dna/dna.h).
    mapping.inversion_probe_local_gate = profile.kind == DnaPresetKind::HiFi;
    common.max_query_seeds_per_strand = profile.placement.vote_seeds;
    mapping.cigar_dp_match = profile.cigar_dp_match;
    mapping.cigar_dp_mismatch = profile.cigar_dp_mismatch;
    mapping.cigar_dp_ambi = profile.cigar_dp_ambi;
    mapping.cigar_dp_gap_open1 = profile.cigar_dp_gap_open1;
    mapping.cigar_dp_gap_extend1 = profile.cigar_dp_gap_extend1;
    mapping.cigar_dp_gap_open2 = profile.cigar_dp_gap_open2;
    mapping.cigar_dp_gap_extend2 = profile.cigar_dp_gap_extend2;
    mapping.cigar_dp_tail_zdrop = profile.cigar_dp_tail_zdrop;
    mapping.cigar_dp_tail_end_bonus = profile.cigar_dp_tail_end_bonus;
    mapping.cigar_dp_bw = profile.cigar_dp_bw;
    mapping.cigar_dp_bw_long = profile.cigar_dp_bw_long;
    mapping.chain_band = profile.chain_band;
    mapping.cigar_dp_max_gap = profile.cigar_dp_max_gap;
    mapping.cigar_dp_min_dp_max = profile.cigar_dp_min_dp_max;
    mapping.fill_dp_match = profile.fill_dp_match;
    mapping.fill_dp_mismatch = profile.fill_dp_mismatch;
    mapping.fill_dp_ambi = profile.fill_dp_ambi;
    mapping.fill_dp_gap_open1 = profile.fill_dp_gap_open1;
    mapping.fill_dp_gap_extend1 = profile.fill_dp_gap_extend1;
    mapping.fill_dp_gap_open2 = profile.fill_dp_gap_open2;
    mapping.fill_dp_gap_extend2 = profile.fill_dp_gap_extend2;
    mapping.fill_dp_tail_zdrop = profile.fill_dp_tail_zdrop;
    mapping.fill_dp_inversion_zdrop = profile.fill_dp_inversion_zdrop;
    // 20000 under lr and the asm rows (minimap2's bw_long), 10000 under lr:hq
    // (map-hifi's max_gap).
    mapping.cigar_local_interval_anchor_chain_max_gap = profile.chain_max_gap;
}

static bool is_hifi_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->kind == DnaPresetKind::HiFi;
}

bool is_assembly_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->kind == DnaPresetKind::Asm;
}

bool is_overlap_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->overlap.on;
}

bool lane_runs_stage(LaneStage stage, std::string_view preset) {
    const bool overlap = is_overlap_preset(preset);
    switch (stage) {
        case LaneStage::None:
            return true;
        case LaneStage::BaseLevelOutput:
        case LaneStage::BaseAlignment:
        case LaneStage::SamOutput:
        case LaneStage::Selection:
        case LaneStage::CandidateChains:
            return !overlap;
        case LaneStage::OverlapPairs:
            return overlap;
    }
    return true;
}

bool is_rna_preset(std::string_view preset) {
    return preset == "splice" || preset == "splice:hq";
}
static bool is_rna_hifi_preset(std::string_view preset) {
    return preset == "splice:hq";
}

static bool is_dna_long_preset(std::string_view preset) {
    return find_dna_preset_profile(preset) != nullptr;
}
bool preset_is_valid(std::string_view preset) {
    return is_dna_long_preset(preset) || is_rna_preset(preset);
}

PresetSeeding resolve_preset_seeding(std::string_view preset) {
    if (const auto* profile = find_dna_preset_profile(preset)) {
      return PresetSeeding{profile->k, profile->syncmer_s};
    }
    // Matches the seed density of minimap2's k=15, w=5 minimizers: closed
    // syncmers with k=15, s=10 have density 2/(k-s+1) = 1/3.
    return PresetSeeding{15, 10};
}

std::string accepted_preset_names() {
  std::string names;
  const auto append = [&names](std::string_view name) {
    if (!names.empty()) names += ", ";
    names += name;
  };
  const auto experimental = [](const DnaPresetProfile& profile) {
    return profile.kind == DnaPresetKind::Asm || profile.overlap.on;
  };
  for (const auto& profile : kDnaPresetProfiles)
    if (!experimental(profile)) append(profile.name);
  append("splice");
  append("splice:hq");
  for (const auto& profile : kDnaPresetProfiles)
    if (experimental(profile)) append(profile.name);
  return names;
}

// Splice presets start from the matching DNA preset's placement defaults,
// then install minimap2 2.30's splice scoring. Explicit overrides come after.
static void set_splice_fields(
    lr::rna::RnaLongOptions& mapping, CommonOptions& common, bool hq) {
    set_dna_long_platform_fields(
        mapping.base, common,
        dna_preset_profile(hq ? DnaPresetKind::HiFi : DnaPresetKind::Ont));
    // RNA diagonals move with splicing, not indel drift, so RNA uses the
    // fixed vote width.
    mapping.base.vote_diag_bin_width_adaptive = false;
    mapping.base.vote_diag_slope_den = 0;
    mapping.base.vote_diag_width_max = mapping.base.vote_diag_bin_width;
    // Clear DNA-only settings inherited from the DNA preset. RNA never reads
    // them, but the resolved configuration should not show them.
    mapping.base.vote_admission_ratio = 0.0;
    mapping.base.cigar_dp_min_ksw_len = 0;
    mapping.base.chain_mapq_hifi_margin = false;
    mapping.base.inversion_probe_local_gate = false;
    mapping.base.dna_pool_gate_occ = 0;
    mapping.base.dna_tile_rescue_occ = 0;
    mapping.base.long_occ_ceiling = 0;
    // The within-read repeat filter and the seed selection, at the vote's cap.
    mapping.base.query_occ_frac = lr::kQueryOccFrac;
    mapping.base.occ_dist = lr::kOccDist;
    // RNA's own support floor.
    common.min_support = 2;

    // minimap2's splice: a=1,b=2,q=2,e=1,q2=32,e2=0; splice:hq: b=4,q=6,q2=24.
    mapping.base.cigar_dp_match       = 1;
    mapping.base.cigar_dp_mismatch    = hq ? 4 : 2;
    mapping.base.cigar_dp_ambi        = 1;
    mapping.base.cigar_dp_gap_open1   = hq ? 6 : 2;
    mapping.base.cigar_dp_gap_extend1 = 1;
    mapping.base.cigar_dp_gap_open2   = hq ? 24 : 32;
    mapping.base.cigar_dp_gap_extend2 = 0;
    mapping.base.cigar_dp_tail_zdrop = 200;
    mapping.base.cigar_dp_tail_end_bonus = -1;
    mapping.base.cigar_dp_min_dp_max = 80;
    // -r: the fine chain's band, and a bound on it; -G sets both.
    mapping.base.cigar_dp_bw = 200000;
    mapping.base.cigar_dp_bw_long = 200000;
    mapping.base.cigar_dp_max_gap = 2000;
    // minimap2's max_gap under splice is also the fine chain's query gap.
    mapping.base.cigar_local_interval_anchor_chain_max_gap = 2000;
    mapping.splice_transition = 0;
    // The penalty is used only by ksw2's scored-splice-site mode, which is
    // not exposed.
    mapping.splice_junction_bonus = 9;
    mapping.splice_junction_penalty = 5;
    mapping.splice_inversion_zdrop = 100;
    // Both strands, the GT-AG motif deciding: minimap2's -ub.
    mapping.strand_mode = 0;
    mapping.min_intron = 20;
    mapping.max_intron = 200000;
    // One vote seed per this many query bases, above the 128 per strand.
    mapping.vote_seed_spacing = hq ? 16 : 8;

    // RNA MAPQ parameters (rna/hit.cpp). ONT cDNA alignments often
    // cover less of the read, so `splice` uses a lower coverage threshold.
    mapping.rna_mapq_qcov_tau = hq ? 0.70 : 0.55;
    // Rivals that mask_level drops (span-disjoint loci the vote found) still
    // count toward ambiguity, discounted because a non-overlapping rival is
    // weaker evidence.
    mapping.rna_mapq_disjoint_vote_damp = 0.90;
}

} // namespace options
} // namespace cpu
} // namespace fa

/**** resolve ****/

namespace fa::cpu::options {

using ::fa::cpu::lr::LongOccPolicy;

namespace {

// -p and -N, and the splice presets' --junc-bed, --junc-bonus and
// --rival-min-diff; each refused on the lane it does not apply to.
void apply_secondary_and_junction_overrides(bool rna_mode,
    ::fa::cpu::lr::DnaLongOptions& mapping,
    ::fa::cpu::lr::rna::RnaLongOptions* rna, const UserOverrides& user) {
  if (!rna_mode && user.rna_junction_bed && !user.rna_junction_bed->empty())
    throw std::invalid_argument(
        "--junc-bed is valid only with an RNA preset");
  if (!rna_mode && user.rna_junction_bonus)
    throw std::invalid_argument(
        "--junc-bonus is valid only with an RNA preset");
  if (!rna_mode && user.rna_rival_min_diff)
    throw std::invalid_argument(
        "--rival-min-diff is valid only with an RNA preset");
  // -p, minimap2's pri_ratio. Validated here as well as in the CLI, for the
  // Python binding.
  if (user.pri_ratio) {
    if (!(*user.pri_ratio >= 0.0) || *user.pri_ratio > 1.0)
      throw std::invalid_argument("-p must be within [0,1]");
    if (rna_mode)
      rna->rival_pri_ratio = *user.pri_ratio;
    else
      mapping.pri_ratio = *user.pri_ratio;
  }
  if (!rna_mode && user.rna_max_loci)
    throw std::invalid_argument(
        "rna_max_loci is valid only with an RNA preset");
  if (rna_mode && user.dna_secondary_max)
    throw std::invalid_argument(
        "dna_secondary_max is valid only with a DNA preset");
  if (user.dna_secondary_max) {
    if (*user.dna_secondary_max < 1)
      throw std::invalid_argument("-N must be >= 1");
    mapping.secondary_max = *user.dna_secondary_max;
  }
  if (rna_mode) {
    if (user.rna_junction_bed)
      rna->junction_bed = *user.rna_junction_bed;
    if (user.rna_junction_bonus) {
      if (*user.rna_junction_bonus < 0 || *user.rna_junction_bonus > 127)
        throw std::invalid_argument(
            "--junc-bonus must be within [0,127]");
      rna->splice_junction_bonus = *user.rna_junction_bonus;
    }
    if (user.rna_rival_min_diff) {
      if (*user.rna_rival_min_diff < 0)
        throw std::invalid_argument(
            "--rival-min-diff must be >= 0");
      rna->rival_min_diff = *user.rna_rival_min_diff;
    }
    // -N n keeps n+1 loci and realizes up to n of them.
    if (user.rna_max_loci) {
      if (*user.rna_max_loci < 1)
        throw std::invalid_argument("-N must be >= 1");
      rna->max_locus_chains = *user.rna_max_loci;
      rna->rival_realize_max = *user.rna_max_loci - 1;
    }
  }
}

// The DNA presets' placement: --max-cands.
void apply_placement_overrides(bool rna_mode,
    ::fa::cpu::lr::DnaLongOptions& mapping, const UserOverrides& user) {
  if (user.max_cands) {
    if (rna_mode)
      throw std::invalid_argument(
          "--max-cands is valid only with a DNA preset");
    // The all-chains lane chains each contig and strand once and takes a
    // wider bound.
    const int ceiling = mapping.all_chains
                            ? ::fa::cpu::voting::kAllChainsLaneBound
                            : ::fa::cpu::voting::kMaxCatalogueLaneBound;
    if (*user.max_cands < 1 || *user.max_cands > ceiling)
      throw std::invalid_argument("--max-cands must be within [1," +
                                  std::to_string(ceiling) + "]");
    mapping.max_cands = *user.max_cands;
  }
}

// k and s: a prebuilt index's own, which an explicit value may not
// contradict; otherwise the explicit value, then the preset's.
void resolve_index_identity(const ResolveRequest& request, IndexIdentity& index,
    const PresetSeeding& seeding, const UserOverrides& user) {
  // 1 <= s <= k; s == k is legal, like minimap2's -w 1.
  if (user.k && *user.k < 1)
    throw std::invalid_argument("-k must be at least 1");
  const int resolved_k = request.index.has_index
                             ? request.index.k
                             : (user.k ? *user.k : seeding.k);
  if (user.syncmer_s && (*user.syncmer_s < 1 || *user.syncmer_s > resolved_k))
    throw std::invalid_argument(
        "-s must be within [1,k] (k=" +
        std::to_string(resolved_k) + ")");
  // An in-memory index has the same k limit as a .faix. A prebuilt index
  // checked its k at load time.
  if (!request.index.has_index && resolved_k > ::fa::cpu::kFaixMaxK)
    throw std::invalid_argument(
        "k=" + std::to_string(resolved_k) +
        " is too large; the maximum is " + std::to_string(::fa::cpu::kFaixMaxK));
  if (request.index.has_index) {
    if (user.k && *user.k != request.index.k) {
      throw std::invalid_argument(
          "explicit -k=" + std::to_string(*user.k) +
          " contradicts the prebuilt index (k=" +
          std::to_string(request.index.k) +
          "); seeding is fixed at build time -- rebuild the index to change k");
    }
    if (user.syncmer_s && *user.syncmer_s != request.index.syncmer_s) {
      throw std::invalid_argument(
          "explicit -s=" + std::to_string(*user.syncmer_s) +
          " contradicts the prebuilt index (s=" +
          std::to_string(request.index.syncmer_s) +
          "); seeding is fixed at build time -- rebuild the index to change s");
    }
    index.k = request.index.k;
    index.syncmer_s = request.index.syncmer_s;
  } else {
    index.k = resolved_k;
    index.syncmer_s = user.syncmer_s ? *user.syncmer_s : seeding.syncmer_s;
  }
  index.syncmer_downsample =
      user.syncmer_downsample && *user.syncmer_downsample > 0
          ? *user.syncmer_downsample
          : 1;
}

// Threads and the requested CIGAR with --cs, --MD and --eqx.
void apply_output_overrides(const ResolveRequest& request, bool overlap_mode,
    CommonOptions& common, const UserOverrides& user) {
  if (user.num_threads)
    common.num_threads = *user.num_threads;
  if (user.enable_full_read_cigar) {
    common.enable_full_read_cigar = *user.enable_full_read_cigar;
  }
  // --cs, --MD and --eqx come from core/core.cpp's replay, which both the
  // DNA and the splice realizers use, so every preset accepts them.
  if (user.cs)
    common.cs = *user.cs;
  if (user.emit_md)
    common.emit_md = *user.emit_md;
  if (user.emit_eqx)
    common.emit_eqx = *user.emit_eqx;
  // --cs and --MD imply a CIGAR, as minimap2's --cs does; --eqx does not.
  if (common.cs != CsMode::None || common.emit_md) {
    common.enable_full_read_cigar = true;
  }
  // The CLI refuses this first, as a usage error.
  if (overlap_mode && common.enable_full_read_cigar)
    throw std::invalid_argument(
        "-x " + request.preset +
        " prints placements only; it is refused with a requested CIGAR");
}

// The vote and its admission: --max-chain-occ, --q-occ-frac, --vote-ratio,
// --dual, the vote seeds, --min-support, -m, the occurrence policy and caps.
void apply_vote_overrides(bool rna_mode, bool overlap_mode,
    CommonOptions& common, ::fa::cpu::lr::DnaLongOptions& mapping,
    ::fa::cpu::lr::rna::RnaLongOptions* rna, const UserOverrides& user) {
  // --max-chain-occ: DNA only. The CLI relies on this check.
  if (rna_mode && user.dna_pool_gate_occ)
    throw std::invalid_argument(
        "--max-chain-occ is valid only with a DNA preset");
  // Replaces the vote's cap on the pool gate only; 0 is ungated. A negative
  // value would alias kDnaPoolGateVoteCap.
  if (!rna_mode && user.dna_pool_gate_occ) {
    if (*user.dna_pool_gate_occ < 0)
      throw std::invalid_argument(
          "--max-chain-occ must be >= 0");
    mapping.dna_pool_gate_occ = *user.dna_pool_gate_occ;
  }
  // --q-occ-frac: read placement, DNA and RNA; the all-chains lane keeps
  // every seed. The CLI relies on this check.
  if (overlap_mode && user.query_occ_frac)
    throw std::invalid_argument(
        "--q-occ-frac is valid only with lr, lr:hq, asm5, asm10, asm20, "
        "splice or splice:hq");
  if (user.query_occ_frac) {
    if (!(*user.query_occ_frac >= 0.0) || *user.query_occ_frac > 1.0)
      throw std::invalid_argument("--q-occ-frac must be within [0,1]");
    mapping.query_occ_frac = *user.query_occ_frac;
  }
  // -e: the seed selection's streak spacing, on every preset; 0 is off.
  if (user.occ_dist) {
    if (*user.occ_dist < 0)
      throw std::invalid_argument("-e must be >= 0");
    mapping.occ_dist = *user.occ_dist;
  }
  // --vote-ratio overrides the preset's ratio; an explicit 0 admits by count.
  if (rna_mode && user.dna_vote_admission_ratio)
    throw std::invalid_argument(
        "--vote-ratio is valid only with a DNA preset");
  if (!rna_mode && user.dna_vote_admission_ratio) {
    if (!(*user.dna_vote_admission_ratio >= 0.0) ||
        *user.dna_vote_admission_ratio > 1.0)
      throw std::invalid_argument(
          "--vote-ratio must be within [0,1]");
    mapping.vote_admission_ratio = *user.dna_vote_admission_ratio;
  }
  // --dual: the overlap presets only.
  if (user.dual) {
    if (!overlap_mode)
      throw std::invalid_argument(
          "--dual is valid only with ava-ont or ava-hifi");
    mapping.dual = *user.dual;
  }
  // The splice presets' vote seeds follow --vote-spacing.
  if (rna_mode && user.max_query_seeds_per_strand)
    throw std::invalid_argument(
        "--vote-seeds is valid only with a DNA preset; the splice presets "
        "take --vote-spacing");
  if (user.max_query_seeds_per_strand) {
    common.max_query_seeds_per_strand =
        std::max(0, *user.max_query_seeds_per_strand);
  }
  if (!rna_mode && user.rna_vote_seed_spacing)
    throw std::invalid_argument(
        "--vote-spacing is valid only with an RNA preset");
  if (user.rna_vote_seed_spacing) {
    if (*user.rna_vote_seed_spacing < 0)
      throw std::invalid_argument("--vote-spacing must be >= 0");
    rna->vote_seed_spacing = *user.rna_vote_seed_spacing;
  }
  if (user.min_support) {
    common.min_support = std::max(0, *user.min_support);
  }
  if (rna_mode && user.min_chain_score)
    throw std::invalid_argument("-m is valid only with a DNA preset");
  if (user.min_chain_score) {
    if (*user.min_chain_score < 1)
      throw std::invalid_argument("-m must be at least 1");
    mapping.min_chain_score = *user.min_chain_score;
  }
  if (user.occ_policy) {
    if (!::fa::cpu::lr::long_occ_policy_string_valid(*user.occ_policy)) {
      throw std::invalid_argument(
          "unknown DNA-long occurrence policy: " +
          *user.occ_policy);
    }
    mapping.long_occ_policy = ::fa::cpu::lr::long_occ_policy_from_string(
        *user.occ_policy, mapping.long_occ_policy);
  }
  if (user.primary_occ_cap) {
    mapping.long_primary_occ_cap = std::max(0, *user.primary_occ_cap);
    if (mapping.long_primary_occ_cap > 0 &&
        mapping.long_occ_policy == LongOccPolicy::Off) {
      mapping.long_occ_policy = LongOccPolicy::Fixed;
    }
  }
  if (user.long_occ_cap) {
    if (*user.long_occ_cap < 0) {
      throw std::invalid_argument(
          "long occurrence cap must be >= 0 (got " +
          std::to_string(*user.long_occ_cap) + ")");
    }
    // Under the default policy this is the floor, which the index quantile
    // can still raise; the fixed policy pins it.
    mapping.long_occ_cap = *user.long_occ_cap;
  }
}

// The vote's diagonal width: a typed width turns the adaptive model off.
void apply_vote_width_overrides(bool rna_mode,
    ::fa::cpu::lr::DnaLongOptions& mapping, const UserOverrides& user) {
  if (user.vote_diag_bin_width) {
    if (*user.vote_diag_bin_width <= 0) {
      throw std::invalid_argument(
          "vote diagonal-bin width must be > 0 (got " +
          std::to_string(*user.vote_diag_bin_width) + ")");
    }
    mapping.vote_diag_bin_width = *user.vote_diag_bin_width;
    mapping.vote_diag_bin_width_adaptive = false;
  }
  // The adaptive model's coefficients leave the model on. RNA presets have
  // no adaptive model.
  if (rna_mode && (user.vote_diag_slope_den || user.vote_diag_width_max))
    throw std::invalid_argument(
        "--dw-slope-den and --dw-max are valid only with "
        "a DNA preset");
  if (user.vote_diag_slope_den) {
    if (*user.vote_diag_slope_den < 0) {
      throw std::invalid_argument(
          "vote diagonal-width slope denominator must be "
          ">= 0 (got " +
          std::to_string(*user.vote_diag_slope_den) + ")");
    }
    mapping.vote_diag_slope_den = *user.vote_diag_slope_den;
  }
  if (user.vote_diag_width_max) {
    if (*user.vote_diag_width_max <= 0) {
      throw std::invalid_argument(
          "vote diagonal-width cap must be > 0 (got " +
          std::to_string(*user.vote_diag_width_max) + ")");
    }
    mapping.vote_diag_width_max = *user.vote_diag_width_max;
  }
}

// -r and -G's bands, -g and -b.
void apply_band_and_gap_overrides(bool rna_mode,
    ::fa::cpu::lr::DnaLongOptions& mapping, const UserOverrides& user) {
  // On a splice preset -G sets both -r values, as minimap2's
  // mm_mapopt_max_intron_len, and -r then applies on top: the parser drops a
  // -r typed before -G, and the API's dp_bw acts as typed after it.
  if (rna_mode && user.rna_max_intron && *user.rna_max_intron > 0) {
    mapping.cigar_dp_bw = *user.rna_max_intron;
    mapping.cigar_dp_bw_long = *user.rna_max_intron;
  }
  if (user.dp_bw)
    mapping.cigar_dp_bw = std::max(1, *user.dp_bw);
  if (user.dp_bw_long)
    mapping.cigar_dp_bw_long = std::max(1, *user.dp_bw_long);
  // minimap2's mm_check_opt.
  if (mapping.cigar_dp_bw > mapping.cigar_dp_bw_long)
    throw std::invalid_argument(
        "with '-rNUM1,NUM2', NUM1 (" + std::to_string(mapping.cigar_dp_bw) +
        ") can't be larger than NUM2 (" +
        std::to_string(mapping.cigar_dp_bw_long) + ")");
  // -g, minimap2's max_gap, sets the DP's gap and the chains' gap, which on a
  // splice preset is the fine chain's query gap. Untyped, each keeps its own
  // preset value.
  if (user.dp_max_gap) {
    mapping.cigar_dp_max_gap = std::max(1, *user.dp_max_gap);
    mapping.cigar_local_interval_anchor_chain_max_gap =
        mapping.cigar_dp_max_gap;
  }
  // -b: the dense chain passes' band.
  if (user.chain_band) {
    if (rna_mode)
      throw std::invalid_argument("-b is valid only with a DNA preset");
    mapping.chain_band = std::max(1, *user.chain_band);
  }
}

// The splice presets' introns and strand, then minimap2's checks on their
// one scoring row.
void apply_splice_overrides(::fa::cpu::lr::DnaLongOptions& mapping,
    ::fa::cpu::lr::rna::RnaLongOptions* rna, const UserOverrides& user) {
  if (user.rna_min_intron && *user.rna_min_intron <= 0) {
    throw std::invalid_argument(
        "--min-intron must be > 0 (got " +
        std::to_string(*user.rna_min_intron) + ")");
  }
  if (user.rna_max_intron && *user.rna_max_intron <= 0) {
    throw std::invalid_argument(
        "-G must be > 0 (got " +
        std::to_string(*user.rna_max_intron) + ")");
  }
  if (user.rna_min_intron && *user.rna_min_intron > 0)
    rna->min_intron = *user.rna_min_intron;
  if (user.rna_max_intron && *user.rna_max_intron > 0)
    rna->max_intron = *user.rna_max_intron;
  // -1 unset, 0 auto (both strands), 1 forward, 2 reverse, 3 none (no
  // splice motif, no ts:A). The API passes the integer straight in.
  if (user.rna_strand) {
    if (*user.rna_strand < -1 || *user.rna_strand > 3)
      throw std::invalid_argument(
          "rna_strand_mode must be -1 (unset), 0 (auto), "
          "1 (forward), 2 (reverse) or 3 (none) (got " +
          std::to_string(*user.rna_strand) + ")");
    rna->strand_mode = *user.rna_strand;
  }
  if (rna->max_intron < rna->min_intron) {
    throw std::invalid_argument("-G (" +
                                std::to_string(rna->max_intron) +
                                ") must be >= --min-intron (" +
                                std::to_string(rna->min_intron) + ")");
  }
  if (mapping.cigar_dp_match <= 0 || mapping.cigar_dp_match > 127) {
    throw std::invalid_argument(
        "RNA match score (-A) must be within [1,127]");
  }
  if (mapping.cigar_dp_mismatch <= 0 || mapping.cigar_dp_mismatch > 127) {
    throw std::invalid_argument(
        "RNA mismatch penalty (-B) must be within "
        "[1,127]");
  }
  if (mapping.cigar_dp_ambi < 0 ||
      mapping.cigar_dp_ambi >= mapping.cigar_dp_mismatch) {
    throw std::invalid_argument(
        "RNA --score-N must be within [0,-B)");
  }
  if (mapping.cigar_dp_gap_extend2 != 0) {
    throw std::invalid_argument(
        "the splice presets require the second "
        "gap-extension penalty E2 to be 0");
  }
  const std::int64_t q1e1 =
      static_cast<std::int64_t>(mapping.cigar_dp_gap_open1) +
      mapping.cigar_dp_gap_extend1;
  const std::int64_t q2e2 =
      static_cast<std::int64_t>(mapping.cigar_dp_gap_open2) +
      mapping.cigar_dp_gap_extend2;
  const bool one_gap_model =
      mapping.cigar_dp_gap_open1 == mapping.cigar_dp_gap_open2 &&
      mapping.cigar_dp_gap_extend1 == mapping.cigar_dp_gap_extend2;
  const bool valid_dual_gap =
      mapping.cigar_dp_gap_extend1 > mapping.cigar_dp_gap_extend2 &&
      q1e1 < q2e2;
  if (!one_gap_model && !valid_dual_gap) {
    throw std::invalid_argument(
        "RNA dual-gap penalties must satisfy "
        "E1>E2 and O1+E1<O2+E2");
  }
  if (q1e1 + q2e2 > 127) {
    throw std::invalid_argument("RNA scoring requires "
                                "(O1+E1)+(O2+E2) <= 127");
  }
  if (mapping.cigar_dp_tail_zdrop < 0 || rna->splice_inversion_zdrop < 0 ||
      mapping.cigar_dp_tail_zdrop < rna->splice_inversion_zdrop) {
    throw std::invalid_argument(
        "RNA Z-drop must be nonnegative and should not be less than "
        "inversion-Z-drop");
  }
  if (mapping.cigar_dp_tail_end_bonus < -1) {
    throw std::invalid_argument(
        "RNA --end-bonus must be >= -1");
  }
}

// While the adaptive vote-width model is on, a cap below the base width would
// silently reduce it to the fixed base width, so it is refused.
void check_vote_diag_model(const ::fa::cpu::lr::DnaLongOptions& mapping) {
  if (mapping.vote_diag_bin_width_adaptive &&
      mapping.vote_diag_width_max < std::max(1, mapping.vote_diag_bin_width)) {
    throw std::invalid_argument(
        "vote diagonal-width cap must be >= the base bin "
        "width (got cap " +
        std::to_string(mapping.vote_diag_width_max) + " < width " +
        std::to_string(mapping.vote_diag_bin_width) + ")");
  }
}

} // namespace

ResolvedMapOptions resolve_options(const ResolveRequest& request) {
  if (!preset_is_valid(request.preset)) {
    throw std::invalid_argument("unknown preset '" +
                                request.preset + "' (expected " +
                                accepted_preset_names() + ")");
  }
  const bool rna_mode = is_rna_preset(request.preset);
  const bool asm_mode = is_assembly_preset(request.preset);
  const bool overlap_mode = is_overlap_preset(request.preset);
  const bool hifi_class = rna_mode ? is_rna_hifi_preset(request.preset)
                                   : is_hifi_preset(request.preset);

  ResolvedMapOptions output;
  output.resolved_.index.from_index = request.index.has_index;
  if (rna_mode) {
    output.resolved_.mapping = ::fa::cpu::lr::rna::RnaLongOptions{};
    set_splice_fields(*output.resolved_.rna(), output.resolved_.common,
                      hifi_class);
  } else {
    output.resolved_.mapping = ::fa::cpu::lr::DnaLongOptions{};
    set_dna_long_platform_fields(output.resolved_.long_read(),
                                 output.resolved_.common,
                                 *find_dna_preset_profile(request.preset));
  }

  auto& common = output.resolved_.common;
  auto& index = output.resolved_.index;
  auto& mapping = output.resolved_.long_read();
  auto* rna = output.resolved_.rna();
  const PresetSeeding seeding = resolve_preset_seeding(request.preset);
  const UserOverrides& user = request.user;

  apply_secondary_and_junction_overrides(rna_mode, mapping, rna, user);

  apply_placement_overrides(rna_mode, mapping, user);
  resolve_index_identity(request, index, seeding, user);

  apply_output_overrides(request, overlap_mode, common, user);
  apply_vote_overrides(rna_mode, overlap_mode, common, mapping, rna, user);
  apply_vote_width_overrides(rna_mode, mapping, user);
  // On lr and lr:hq -A -B -O -E -z --score-N set the gap-fill row and the
  // end row stays the preset's; a splice or assembly preset has one row.
  const bool fill_row = !rna_mode && !asm_mode;
  int& dp_match = fill_row ? mapping.fill_dp_match : mapping.cigar_dp_match;
  int& dp_mismatch =
      fill_row ? mapping.fill_dp_mismatch : mapping.cigar_dp_mismatch;
  int& dp_ambi = fill_row ? mapping.fill_dp_ambi : mapping.cigar_dp_ambi;
  int& dp_gap_open1 =
      fill_row ? mapping.fill_dp_gap_open1 : mapping.cigar_dp_gap_open1;
  int& dp_gap_open2 =
      fill_row ? mapping.fill_dp_gap_open2 : mapping.cigar_dp_gap_open2;
  int& dp_gap_extend1 =
      fill_row ? mapping.fill_dp_gap_extend1 : mapping.cigar_dp_gap_extend1;
  int& dp_gap_extend2 =
      fill_row ? mapping.fill_dp_gap_extend2 : mapping.cigar_dp_gap_extend2;
  int& dp_tail_zdrop =
      fill_row ? mapping.fill_dp_tail_zdrop : mapping.cigar_dp_tail_zdrop;
  int& dp_inversion_zdrop = fill_row ? mapping.fill_dp_inversion_zdrop
                                     : mapping.cigar_dp_inversion_zdrop;
  // A splice preset clamps a typed value; a DNA preset takes it as typed and
  // checks the row below.
  const auto typed = [rna_mode](int value, int floor) {
    return rna_mode ? std::max(floor, value) : value;
  };
  if (user.dp_match)
    dp_match = typed(*user.dp_match, 1);
  if (user.dp_mismatch)
    dp_mismatch = typed(*user.dp_mismatch, 0);
  if (user.dp_ambi)
    dp_ambi = typed(*user.dp_ambi, 0);
  if (user.dp_gap_open1)
    dp_gap_open1 = typed(*user.dp_gap_open1, 1);
  if (user.dp_gap_open2)
    dp_gap_open2 = typed(*user.dp_gap_open2, 0);
  if (user.dp_gap_extend1)
    dp_gap_extend1 = typed(*user.dp_gap_extend1, 1);
  if (user.dp_gap_extend2)
    dp_gap_extend2 = typed(*user.dp_gap_extend2, 0);
  if (user.dp_tail_zdrop) {
    dp_tail_zdrop = *user.dp_tail_zdrop;
    // As in minimap2, a scalar -z also sets the inversion Z-drop.
    const int inversion_zdrop =
        user.dp_inversion_zdrop.value_or(*user.dp_tail_zdrop);
    dp_inversion_zdrop = inversion_zdrop;
    if (rna)
      rna->splice_inversion_zdrop = inversion_zdrop;
  }
  if (user.dp_tail_end_bonus)
    mapping.cigar_dp_tail_end_bonus = *user.dp_tail_end_bonus;
  // An assembly preset fills its gaps under that one row.
  if (asm_mode) {
    mapping.fill_dp_match = mapping.cigar_dp_match;
    mapping.fill_dp_mismatch = mapping.cigar_dp_mismatch;
    mapping.fill_dp_ambi = mapping.cigar_dp_ambi;
    mapping.fill_dp_gap_open1 = mapping.cigar_dp_gap_open1;
    mapping.fill_dp_gap_extend1 = mapping.cigar_dp_gap_extend1;
    mapping.fill_dp_gap_open2 = mapping.cigar_dp_gap_open2;
    mapping.fill_dp_gap_extend2 = mapping.cigar_dp_gap_extend2;
    mapping.fill_dp_tail_zdrop = mapping.cigar_dp_tail_zdrop;
    mapping.fill_dp_inversion_zdrop = mapping.cigar_dp_inversion_zdrop;
  }
  if (user.dp_min_dp_max) {
    mapping.cigar_dp_min_dp_max = std::max(0, *user.dp_min_dp_max);
  }
  apply_band_and_gap_overrides(rna_mode, mapping, user);

  if (rna_mode) {
    apply_splice_overrides(mapping, rna, user);
  }
  // minimap2's mm_check_opt on the row -A -B -O -E --score-N set; the splice
  // presets check theirs above.
  if (!rna_mode) {
    if (dp_match < 1 || dp_match > 127)
      throw std::invalid_argument("-A must be within [1,127]");
    if (dp_mismatch < 1 || dp_mismatch > 127)
      throw std::invalid_argument("-B must be within [1,127]");
    if (dp_gap_open1 <= 0 || dp_gap_extend1 <= 0 || dp_gap_open2 < 0 ||
        dp_gap_extend2 < 0)
      throw std::invalid_argument("-O and -E must be positive");
    const std::int64_t q1e1 =
        static_cast<std::int64_t>(dp_gap_open1) + dp_gap_extend1;
    const std::int64_t q2e2 =
        static_cast<std::int64_t>(dp_gap_open2) + dp_gap_extend2;
    if ((dp_gap_open1 != dp_gap_open2 || dp_gap_extend1 != dp_gap_extend2) &&
        !(dp_gap_extend1 > dp_gap_extend2 && q1e1 < q2e2))
      throw std::invalid_argument(
          "dual gap penalties violating E1>E2 and O1+E1<O2+E2");
    if (q1e1 + q2e2 > 127)
      throw std::invalid_argument(
          "scoring system violating ({-O}+{-E})+({-O2}+{-E2}) <= 127");
    if (dp_ambi < 0 || dp_ambi >= dp_mismatch)
      throw std::invalid_argument("--score-N should be within [0,{-B})");
    if (mapping.cigar_dp_tail_end_bonus < -1)
      throw std::invalid_argument("--end-bonus must be >= -1");
  }
  if (!rna_mode &&
      mapping.fill_dp_tail_zdrop < mapping.fill_dp_inversion_zdrop) {
    throw std::invalid_argument(
        "Z-drop should not be less than inversion-Z-drop");
  }
  if (!rna_mode && mapping.fill_dp_inversion_zdrop < 0) {
    throw std::invalid_argument("inversion Z-drop should not be negative");
  }
  if (mapping.vote_diag_bin_width <= 0) {
    throw std::invalid_argument(
        "resolved vote diagonal-bin width must be > 0 "
        "(got " +
        std::to_string(mapping.vote_diag_bin_width) + ")");
  }
  // -S is in the end row's units; the inversion gates of a fill compare
  // fill-row scores. Scaled after the checks above, which keep a typed end
  // row's match score positive.
  mapping.fill_dp_min_dp_max = static_cast<int>(std::min<std::int64_t>(
      std::numeric_limits<int>::max(),
      static_cast<std::int64_t>(mapping.cigar_dp_min_dp_max) *
          mapping.fill_dp_match / mapping.cigar_dp_match));
  check_vote_diag_model(mapping);
  return output;
}

} // namespace fa::cpu::options
