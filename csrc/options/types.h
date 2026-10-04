// Inputs and output of option resolution (options/resolve.h).
#pragma once

#include "resolved_options.h"

#include <optional>
#include <string>

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
  std::optional<int> max_query_seeds_per_strand;
  std::optional<int> query_tile_supported_reward;
  std::optional<int> query_tile_block_open_cost;
  std::optional<int> query_tile_null_cost;
  std::optional<int> query_tile_unsupported_cost;
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
  // DP bandwidths (-r INT[,INT]) and maximum gap (-g), which on a DNA preset
  // is also the dense chain's.
  std::optional<int> dp_bw;
  std::optional<int> dp_bw_long;
  std::optional<int> dp_max_gap;
  std::optional<int> rna_min_intron;
  std::optional<int> rna_max_intron;
  // RNA transcript strand (-u) as rna::StrandMode: 0 auto, 1 forward,
  // 2 reverse, 3 none (no splice motif scored).
  std::optional<int> rna_strand;
  std::optional<std::string> rna_junction_bed;
  std::optional<int> rna_junction_bonus;
  // -p: on a DNA preset the credibility ratio of the alternative and of a
  // block's rival, on a splice preset the rival retention ratio.
  std::optional<double> pri_ratio;
  // RNA rival retention band (--rival-min-diff) and catalogue depth (-N + 1),
  // which also sets the realization budget to one fewer.
  std::optional<int> rna_rival_min_diff;
  std::optional<int> rna_max_loci; // >= 1
  // -N on a DNA preset: alternatives realized per read, >= 1.
  std::optional<int> dna_alternative_realize_max;
  // --vote-ratio: admit vote peaks by ratio to the read's best vote, in
  // [0,1]; 0 admits by count. DNA presets only.
  std::optional<double> dna_vote_admission_ratio;
  // --max-chain-occ: the dense chain's pool occurrence gate, >= 0; 0 is no
  // gate. Unset, the vote's cap. DNA presets only.
  std::optional<int> dna_pool_gate_occ;
  // --max-cands: the lane bound, the vote's peaks and the catalogue's
  // candidates per strand, 1..64. DNA presets only.
  std::optional<int> max_cands;
  // --tiles: the query partition's tiles per read, 2..4096. DNA presets only.
  std::optional<int> query_tiles;
  // --tile-owner: true for anchors, false for span. DNA presets only.
  std::optional<bool> tile_owner_anchors;
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
