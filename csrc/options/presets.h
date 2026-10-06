// Preset names and their seeding and mapping defaults, shared by index
// construction, option resolution and the Python binding.
#pragma once

#include "common_profile.h"
#include "dna_profile.h"
#include "lane_stage.h"
#include "rna_profile.h"

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>

namespace fa {
namespace cpu {
namespace options {

// Asm marks the experimental assembly-to-reference rows.
enum class DnaPresetKind { Ont, HiFi, Asm };

// The placement values --vote-seeds, --tiles and --tile-owner replace.
struct DnaPlacementDefaults {
  int vote_seeds;
  int query_tiles;
  bool tile_owner_anchors;
};

inline constexpr DnaPlacementDefaults kReadPlacement{
    128, ::fa::cpu::voting::kQueryTileCount, false};
// Above 128 vote seeds the nested sampler applies (seeding/syncmer.h).
inline constexpr DnaPlacementDefaults kAssemblyPlacement{4096, 2048, true};
inline constexpr DnaPlacementDefaults kOntOverlapPlacement{
    256, ::fa::cpu::voting::kQueryTileCount, false};

// An all-vs-all overlap row turns on both of its mechanisms (dna_profile.h
// skip_self, all_chains) and replaces --max-cands and --vote-ratio, and -m
// where min_chain_score is positive. Off on every other row.
struct DnaOverlapDefaults {
  bool on;
  int max_cands;
  double vote_ratio;
  int min_chain_score;
  // --dual's default (dna_profile.h dual).
  bool dual;
};

inline constexpr DnaOverlapDefaults kNoOverlap{false, 0, 0.0, 0, true};
// The all-chains lane's widest bound and count admission, and minimap2's
// --dual=no.
inline constexpr DnaOverlapDefaults kHiFiOverlap{
    true, ::fa::cpu::voting::kAllChainsLaneBound, 0.0, 0, false};
// The same with minimap2's ava -m100.
inline constexpr DnaOverlapDefaults kOntOverlap{
    true, ::fa::cpu::voting::kAllChainsLaneBound, 0.0, 100, false};

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
  ::fa::cpu::voting::QueryPartitionParameters query_partition;
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
  // and the letters set both (options/resolve.cpp).
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
  int residue_min_interval_bp;
  int residue_min_anchor_density_per_100bp;
  DnaPlacementDefaults placement;
  DnaOverlapDefaults overlap;
};

inline constexpr std::array<DnaPresetProfile, 7> kDnaPresetProfiles{{
    // End row as minimap2: lr is map-ont (-A2 -B4 -O4,24 -E2,1), lr:hq is
    // map-hifi (-A1 -B4 -O6,26 -E2,1). The fill rows are gentler on gaps.
    {"lr", DnaPresetKind::Ont, 21, 9, 64, 128, 2048, {4, 12, 0, 1, 2},
     2, 4, 1, 4, 2, 24, 1, 400, -1, 500, 20000,
     5000, 80,
     4, 8, 2, 8, 4, 48, 1, 800, 200,
     20000, 200, 9, kReadPlacement, kNoOverlap},
    {"lr:hq", DnaPresetKind::HiFi, 21, 5, 48, 64, 2048, {4, 12, 0, 1, 2},
     1, 4, 1, 6, 2, 26, 1, 400, -1, 500, 20000,
     10000, 200,
     3, 12, 3, 18, 6, 78, 1, 1200, 600,
     10000, 100, 9, kReadPlacement, kNoOverlap},
    // minimap2's asm5, asm10 and asm20 scoring with its -r1k,100k -g10k -s200
    // -z200, one row for the fills and the read ends; lr's seeding, vote
    // width, chain gap and terminal-clip bounds.
    {"asm5", DnaPresetKind::Asm, 21, 9, 64, 128, 2048, {4, 12, 0, 1, 2},
     1, 19, 1, 39, 3, 81, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 19, 1, 39, 3, 81, 1, 200, 200,
     20000, 200, 9, kAssemblyPlacement, kNoOverlap},
    {"asm10", DnaPresetKind::Asm, 21, 9, 64, 128, 2048, {4, 12, 0, 1, 2},
     1, 9, 1, 16, 2, 41, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 9, 1, 16, 2, 41, 1, 200, 200,
     20000, 200, 9, kAssemblyPlacement, kNoOverlap},
    {"asm20", DnaPresetKind::Asm, 21, 9, 64, 128, 2048, {4, 12, 0, 1, 2},
     1, 4, 1, 6, 2, 26, 1, 200, -1, 1000, 100000,
     10000, 200,
     1, 4, 1, 6, 2, 26, 1, 200, 200,
     20000, 200, 9, kAssemblyPlacement, kNoOverlap},
    // Experimental all-vs-all read overlap, map-only: lr's row at k17 with
    // the vote width L/64, and lr:hq's row, as overlap rows. They follow the
    // row of their kind, which dna_preset_profile returns.
    {"ava-ont", DnaPresetKind::Ont, 17, 9, 64, 64, 2048, {4, 12, 0, 1, 2},
     2, 4, 1, 4, 2, 24, 1, 400, -1, 500, 20000,
     5000, 80,
     4, 8, 2, 8, 4, 48, 1, 800, 200,
     20000, 200, 9, kOntOverlapPlacement, kOntOverlap},
    {"ava-hifi", DnaPresetKind::HiFi, 21, 5, 48, 64, 2048, {4, 12, 0, 1, 2},
     1, 4, 1, 6, 2, 26, 1, 400, -1, 500, 20000,
     10000, 200,
     3, 12, 3, 18, 6, 78, 1, 1200, 600,
     10000, 100, 9, kReadPlacement, kHiFiOverlap},
}};

inline const DnaPresetProfile* find_dna_preset_profile(
    std::string_view preset) {
  for (const auto& profile : kDnaPresetProfiles) {
    if (profile.name == preset) return &profile;
  }
  return nullptr;
}

inline const DnaPresetProfile& dna_preset_profile(DnaPresetKind kind) {
  for (const auto& profile : kDnaPresetProfiles) {
    if (profile.kind == kind) return profile;
  }
  throw std::logic_error("missing DNA preset profile");
}

// Installs a DNA preset's mapping defaults. (k, s) are resolved separately.
inline void set_dna_long_platform_fields(
    lr::DnaLongOptions& mapping, CommonOptions& common,
    const DnaPresetProfile& profile) {
    mapping.vote_diag_bin_width = profile.vote_diag_bin_width;
    mapping.vote_diag_bin_width_adaptive = true;
    mapping.vote_diag_slope_den = profile.vote_diag_slope_den;
    mapping.vote_diag_width_max = profile.vote_diag_width_max;
    // The floor of the vote occurrence cap. On a DNA run the engine raises
    // it to the index's mid-occurrence quantile when larger
    // (seeding/context.h), up to long_occ_ceiling.
    mapping.long_occ_cap = 200;
    mapping.query_partition = profile.query_partition;
    mapping.query_tiles = profile.placement.query_tiles;
    mapping.tile_owner_anchors = profile.placement.tile_owner_anchors;
    mapping.vote_admission_ratio = lr::kDnaProductionVoteAdmissionRatio;
    if (profile.overlap.on) {
      mapping.skip_self = true;
      mapping.all_chains = true;
      mapping.max_cands = profile.overlap.max_cands;
      mapping.vote_admission_ratio = profile.overlap.vote_ratio;
      if (profile.overlap.min_chain_score > 0)
        mapping.min_chain_score = profile.overlap.min_chain_score;
      mapping.dual = profile.overlap.dual;
    }
    mapping.dna_tandem_window = lr::kDnaTandemWindow;
    // minimap2's min_ksw_len, the piece length of its gap-filling loop.
    mapping.cigar_dp_min_ksw_len = lr::kDnaMinKswLen;
    const bool hifi_row = profile.kind == DnaPresetKind::HiFi;
    // The vote's cap N is also the dense chain's pool gate and
    // cigar_local_global_occ (engine/aligner.h make_dna_context). lr:hq caps
    // N at minimap2's map-hifi 500. --max-chain-occ replaces N on the gate
    // only.
    mapping.long_occ_ceiling =
        hifi_row ? lr::kDnaHiFiOccCeiling : lr::kDnaOccCeiling;
    mapping.dna_pool_gate_occ = lr::kDnaPoolGateVoteCap;
    // The vote's empty-tile rescue past N, up to M (seeding/syncmer.h).
    mapping.dna_tile_rescue_occ =
        hifi_row ? lr::kDnaHiFiTileRescueOcc : lr::kDnaTileRescueOcc;
    // HiFi only: the MAPQ's dp1/dp2 come from minibwa-style rescoring of the
    // CIGAR (dna/postdp_scoring.h) instead of the raw DP scores.
    mapping.postdp_rescoring = profile.kind == DnaPresetKind::HiFi;
    // HiFi only: the chain MAPQ margin rule (dna/chain_mapq.h).
    mapping.chain_mapq_hifi_margin = profile.kind == DnaPresetKind::HiFi;
    // HiFi only: the late inversion probe's local-chain gate
    // (dna/inv_local_chain.h).
    mapping.inversion_probe_local_gate = profile.kind == DnaPresetKind::HiFi;
    common.min_support = 3;
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
    mapping.residue_min_interval_bp = profile.residue_min_interval_bp;
    mapping.residue_min_anchor_density_per_100bp =
        profile.residue_min_anchor_density_per_100bp;
}

inline bool is_hifi_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->kind == DnaPresetKind::HiFi;
}

// asm5, asm10 and asm20: experimental, not qualified on intact chromosomes.
inline bool is_assembly_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->kind == DnaPresetKind::Asm;
}

// ava-ont and ava-hifi: experimental, map-only.
inline bool is_overlap_preset(std::string_view preset) {
    const auto* profile = find_dna_preset_profile(preset);
    return profile && profile->overlap.on;
}

// Whether the lane `preset` selects runs `stage`. The all-chains lane of an
// overlap row prints map-only PAF and every chain, solves no partition, and
// alone chooses which reads of a pair print it.
inline bool lane_runs_stage(LaneStage stage, std::string_view preset) {
    const bool overlap = is_overlap_preset(preset);
    switch (stage) {
        case LaneStage::None:
            return true;
        case LaneStage::BaseLevelOutput:
        case LaneStage::BaseAlignment:
        case LaneStage::SamOutput:
        case LaneStage::Selection:
        case LaneStage::Partition:
            return !overlap;
        case LaneStage::OverlapPairs:
            return overlap;
    }
    return true;
}

inline bool is_rna_preset(std::string_view preset) {
    return preset == "splice" || preset == "splice:hq";
}
inline bool is_rna_hifi_preset(std::string_view preset) {
    return preset == "splice:hq";
}

inline bool is_dna_long_preset(std::string_view preset) {
    return find_dna_preset_profile(preset) != nullptr;
}
inline bool preset_is_valid(std::string_view preset) {
    return is_dna_long_preset(preset) || is_rna_preset(preset);
}

// Preset (k, s) for building an index. A loaded index uses its own values.
struct PresetSeeding { int k; int syncmer_s; };
inline PresetSeeding resolve_preset_seeding(std::string_view preset) {
    if (const auto* profile = find_dna_preset_profile(preset)) {
      return PresetSeeding{profile->k, profile->syncmer_s};
    }
    // Matches the seed density of minimap2's k=15, w=5 minimizers: closed
    // syncmers with k=15, s=10 have density 2/(k-s+1) = 1/3.
    return PresetSeeding{15, 10};
}

// Every accepted preset name, for error messages; the experimental ones last.
inline std::string accepted_preset_names() {
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
inline void set_splice_fields(
    lr::rna::RnaLongOptions& mapping, CommonOptions& common, bool hq) {
    set_dna_long_platform_fields(
        mapping.base, common,
        dna_preset_profile(hq ? DnaPresetKind::HiFi : DnaPresetKind::Ont));
    // The query-tile objective is DNA-only.
    mapping.base.query_partition =
        ::fa::cpu::voting::QueryPartitionParameters{};
    // RNA diagonals move with splicing, not indel drift, so RNA uses the
    // fixed vote width.
    mapping.base.vote_diag_bin_width_adaptive = false;
    mapping.base.vote_diag_slope_den = 0;
    mapping.base.vote_diag_width_max = mapping.base.vote_diag_bin_width;
    // Clear DNA-only settings inherited from the DNA preset. RNA never reads
    // them, but the resolved configuration should not show them.
    mapping.base.vote_admission_ratio = 0.0;
    mapping.base.dna_tandem_window = 0;
    mapping.base.cigar_dp_min_ksw_len = 0;
    mapping.base.postdp_rescoring = false;
    mapping.base.chain_mapq_hifi_margin = false;
    mapping.base.inversion_probe_local_gate = false;
    mapping.base.dna_pool_gate_occ = 0;
    mapping.base.dna_tile_rescue_occ = 0;
    mapping.base.long_occ_ceiling = 0;
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
    mapping.splice_transition = 0;
    // The penalty is used only by ksw2's scored-splice-site mode, which is
    // not exposed.
    mapping.splice_junction_bonus = 9;
    mapping.splice_junction_penalty = 5;
    mapping.splice_inversion_zdrop = 100;
    mapping.min_intron = 20;
    mapping.max_intron = 200000;

    // RNA MAPQ parameters (rna/chain_mapq.h). ONT cDNA alignments often
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
