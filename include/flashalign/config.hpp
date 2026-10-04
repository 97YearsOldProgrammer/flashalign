#pragma once

#include <optional>
#include <string>
#include <utility>

namespace flashalign {

/** Form of the cs:Z difference string (minimap2 --cs[=short|long]); None emits none. */
enum class CsMode { None, Short, Long };

/**
 * Mapping options, named as `flashalign align --show-config` prints them (a dot there is an
 * underscore here). -1, an empty optional or an empty string means "take the preset's
 * value"; anything else is an explicit setting, which the option resolver may reject with
 * std::invalid_argument. The bracketed names are the matching CLI options.
 */
struct Config {
  std::string preset = "lr";    ///< 'lr', 'lr:hq', 'splice' or 'splice:hq' [-x]
  int k = -1;                   ///< seed k-mer length [-k]; fixed by a loaded index
  int min_support = -1;         ///< minimum anchor support [--min-support]
  int max_query_seeds = -1;     ///< query seeds kept per strand
  bool full_read_cigar = true;  ///< realize base-level CIGARs; false maps only (plain PAF)
  int syncmer_s = -1;           ///< closed-syncmer s [-s]; fixed by a loaded index
  int syncmer_downsample = 1;   ///< query seed downsampling stride; 1 = none
  int vote_diag_bin_width = -1; ///< vote diagonal bin width [--dw]
  int long_occ_cap = -1;        ///< seed occurrence cap [--max-vote-occ]
  int long_primary_occ_cap = -1;  ///< primary-pass occurrence cap
  int threads = 0;              ///< worker threads [-t]; 0 = one per available core
  int dp_match = -1;            ///< match score [-A]
  int dp_mismatch = -1;         ///< mismatch penalty [-B]
  int dp_score_n = -1;          ///< score against an N [--score-N]
  int dp_gap_open1 = -1;        ///< first gap-open penalty [-O]
  int dp_gap_extend1 = -1;      ///< first gap-extension penalty [-E]
  int dp_gap_open2 = -1;        ///< second gap-open penalty [-O]
  int dp_gap_extend2 = -1;      ///< second gap-extension penalty [-E]
  int dp_zdrop = -1;            ///< Z-drop [-z]
  /// Inversion Z-drop, -z's second value (the dp_inversion_zdrop row; splice_inv_zdrop under a
  /// splice preset). -1 follows dp_zdrop, as a lone -z does. config() reports both, so after
  /// editing dp_zdrop on a copy of it set this too, or -1.
  int dp_zdrop_inv = -1;
  int dp_end_bonus = -2;        ///< end bonus [--end-bonus]; -2 = preset's value
  int rna_min_intron = 20;      ///< splice: minimum intron length [--min-intron]
  int rna_max_intron = 200000;  ///< splice: maximum intron length [-G]
  /// splice: transcript strand [-u]: -1 unset, 0 auto, 1 forward, 2 reverse, 3 none (no
  /// splice-motif scoring).
  int rna_strand_mode = -1;
  /// Query-tile partition scores [--tile-score]: supported-tile reward, block-open cost,
  /// null-tile cost and unsupported-tile cost.
  std::optional<int> tile_score_hit;
  std::optional<int> tile_score_block;
  std::optional<int> tile_score_null;
  std::optional<int> tile_score_miss;
  std::string rna_junction_bed;  ///< splice: known junctions, BED6/BED12 [--junc-bed]
  int rna_junction_bonus = -1;   ///< splice: known-junction bonus [--junc-bonus]
  double rna_rival_pri_ratio = -1.0;  ///< splice: secondary-to-primary score ratio [-p], [0,1]
  int rna_max_loci = -1;              ///< splice: candidate loci kept per read [-N], >= 1
  int dp_min_score = -1;        ///< minimum DP alignment score [-S] (minimap2 -s)
  /// Bandwidths [-r INT[,INT]]. A splice preset rejects them, and config() reports -1 there.
  int dp_bw = -1;
  int dp_bw_long = -1;
  /// cs and MD output [--cs, --MD]. Either one turns on full_read_cigar.
  CsMode cs = CsMode::None;
  bool emit_md = false;
  /// =/X instead of M in every realized CIGAR [--eqx]. NM, MD and cs are unchanged.
  bool emit_eqx = false;
};

/**
 * The (k, syncmer_s) pair a preset builds its index with. Mapping always uses the seeding
 * stored in the index, so pass this pair to Index::build or Index::build_from_fasta when the
 * index is meant for a preset other than lr.
 * @throws std::invalid_argument for an unknown preset
 */
std::pair<int, int> preset_seeding(const std::string& preset);

}  // namespace flashalign
