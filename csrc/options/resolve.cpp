#include "resolve.h"

#include "presets.h"
#include "../index/format.h"
#include "../seeding/context.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace fa::cpu::options {

using ::fa::cpu::lr::LongOccPolicy;

namespace {

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

  if (!rna_mode && user.rna_junction_bed && !user.rna_junction_bed->empty())
    throw std::invalid_argument(
        "--junc-bed is valid only with an RNA preset");
  if (!rna_mode && user.rna_junction_bonus)
    throw std::invalid_argument(
        "--junc-bonus is valid only with an RNA preset");
  if (!rna_mode &&
      (user.rna_pri_ratio || user.rna_rival_min_diff || user.rna_realize_max ||
       user.rna_max_loci))
    throw std::invalid_argument(
        "-p, -N, --realize-max and --rival-min-diff "
        "are valid only with an RNA preset");
  if (rna_mode) {
    if (user.rna_junction_bed)
      rna->junction_bed = *user.rna_junction_bed;
    if (user.rna_junction_bonus) {
      if (*user.rna_junction_bonus < 0 || *user.rna_junction_bonus > 127)
        throw std::invalid_argument(
            "--junc-bonus must be within [0,127]");
      rna->splice_junction_bonus = *user.rna_junction_bonus;
    }
    // Validated here as well as in the CLI, for the Python binding.
    if (user.rna_pri_ratio) {
      if (!(*user.rna_pri_ratio >= 0.0) || *user.rna_pri_ratio > 1.0)
        throw std::invalid_argument("-p must be within [0,1]");
      rna->rival_pri_ratio = *user.rna_pri_ratio;
    }
    if (user.rna_rival_min_diff) {
      if (*user.rna_rival_min_diff < 0)
        throw std::invalid_argument(
            "--rival-min-diff must be >= 0");
      rna->rival_min_diff = *user.rna_rival_min_diff;
    }
    if (user.rna_realize_max) {
      if (*user.rna_realize_max < 0)
        throw std::invalid_argument(
            "--realize-max must be >= 0");
      rna->rival_realize_max = *user.rna_realize_max;
    }
    if (user.rna_max_loci) {
      if (*user.rna_max_loci < 1)
        throw std::invalid_argument("-N must be >= 1");
      rna->max_locus_chains = *user.rna_max_loci;
    }
  }

  const bool has_tile_override =
      user.query_tile_supported_reward || user.query_tile_block_open_cost ||
      user.query_tile_null_cost || user.query_tile_unsupported_cost;
  if (rna_mode && has_tile_override) {
    throw std::invalid_argument(
        "--tile-score is incompatible with RNA preset '" +
        request.preset + "'");
  }
  const auto validate_tile_cost = [](int value, const char* key) {
    if (value < 0) {
      throw std::invalid_argument(
          std::string("tile score '") + key +
          "' must be nonnegative");
    }
  };
  if (user.query_tile_supported_reward) {
    if (*user.query_tile_supported_reward <= 0) {
      throw std::invalid_argument(
          "tile score 'hit' must be positive");
    }
    mapping.query_partition.supported_tile_reward =
        *user.query_tile_supported_reward;
  }
  if (user.query_tile_block_open_cost) {
    validate_tile_cost(*user.query_tile_block_open_cost, "block");
    mapping.query_partition.block_open_cost = *user.query_tile_block_open_cost;
  }
  if (user.query_tile_null_cost) {
    validate_tile_cost(*user.query_tile_null_cost, "null");
    mapping.query_partition.null_tile_cost = *user.query_tile_null_cost;
  }
  if (user.query_tile_unsupported_cost) {
    validate_tile_cost(*user.query_tile_unsupported_cost, "miss");
    mapping.query_partition.unsupported_ownership_cost =
        *user.query_tile_unsupported_cost;
  }
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

  if (rna_mode && !native_splice_k_supported(index.k)) {
    throw std::invalid_argument(
        "the splice presets require k=" +
        std::to_string(kNativeSpliceK) + " (got k=" + std::to_string(index.k) +
        "); rebuild the splice index with -k " +
        std::to_string(kNativeSpliceK) + " or choose a DNA preset");
  }

  // The splice kernel has no bandwidth parameter, so -r is refused rather
  // than ignored.
  if (rna_mode && (user.dp_bw || user.dp_bw_long)) {
    throw std::invalid_argument(
        "-r is not supported under the splice presets");
  }

  if (user.chain_max_gap) {
    mapping.cigar_local_interval_anchor_chain_max_gap =
        std::max(1, *user.chain_max_gap);
  }
  if (user.interval_pad) {
    mapping.cigar_local_interval_anchor_interval_pad =
        std::max(0, *user.interval_pad);
  }
  if (user.num_threads)
    common.num_threads = *user.num_threads;
  if (user.enable_full_read_cigar) {
    common.enable_full_read_cigar = *user.enable_full_read_cigar;
  }
  // --cs, --MD and --eqx come from core/cigar.cpp's replay, which both the
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
  if (user.max_query_seeds_per_strand) {
    common.max_query_seeds_per_strand =
        std::max(0, *user.max_query_seeds_per_strand);
  }
  if (user.min_support) {
    common.min_support = std::max(0, *user.min_support);
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
  if (user.cigar_band_frac) {
    if (!std::isfinite(*user.cigar_band_frac) ||
        *user.cigar_band_frac <= 0.0f) {
      throw std::invalid_argument(
          "CIGAR band fraction must be finite and > 0");
    }
    mapping.cigar_band_frac = *user.cigar_band_frac;
  }
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
  if (user.dp_match)
    mapping.cigar_dp_match = std::max(1, *user.dp_match);
  if (user.dp_mismatch)
    mapping.cigar_dp_mismatch = std::max(0, *user.dp_mismatch);
  if (user.dp_ambi)
    mapping.cigar_dp_ambi = std::max(0, *user.dp_ambi);
  if (user.dp_gap_open1)
    mapping.cigar_dp_gap_open1 = std::max(1, *user.dp_gap_open1);
  if (user.dp_gap_open2)
    mapping.cigar_dp_gap_open2 = std::max(0, *user.dp_gap_open2);
  if (user.dp_gap_extend1)
    mapping.cigar_dp_gap_extend1 = std::max(1, *user.dp_gap_extend1);
  if (user.dp_gap_extend2)
    mapping.cigar_dp_gap_extend2 = std::max(0, *user.dp_gap_extend2);
  if (user.dp_tail_zdrop) {
    mapping.cigar_dp_tail_zdrop = *user.dp_tail_zdrop;
    // As in minimap2, a scalar -z also sets the inversion Z-drop.
    if (rna)
      rna->splice_inversion_zdrop = *user.dp_tail_zdrop;
  }
  if (user.dp_tail_end_bonus)
    mapping.cigar_dp_tail_end_bonus = *user.dp_tail_end_bonus;
  if (user.dp_min_dp_max) {
    mapping.cigar_dp_min_dp_max = std::max(0, *user.dp_min_dp_max);
  }
  if (user.dp_bw)
    mapping.cigar_dp_bw = std::max(1, *user.dp_bw);
  if (user.dp_bw_long)
    mapping.cigar_dp_bw_long = std::max(1, *user.dp_bw_long);
  if (user.dp_max_gap)
    mapping.cigar_dp_max_gap = std::max(1, *user.dp_max_gap);

  if (rna_mode) {
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
          "RNA Z-drop must be nonnegative and not less "
          "than inversion Z-drop");
    }
    if (mapping.cigar_dp_tail_end_bonus < -1) {
      throw std::invalid_argument(
          "RNA --end-bonus must be >= -1");
    }
  }
  if (mapping.vote_diag_bin_width <= 0) {
    throw std::invalid_argument(
        "resolved vote diagonal-bin width must be > 0 "
        "(got " +
        std::to_string(mapping.vote_diag_bin_width) + ")");
  }
  check_vote_diag_model(mapping);
  return output;
}

} // namespace fa::cpu::options
