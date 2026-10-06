#pragma once

#include "cli/option_registry.h"

#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace fa::cpu::cli {

// Parsed align command line. An unset optional means the option was not given
// and the preset value stands. Seeding (k, s) is not here: align takes it from
// the target .faix or the preset.
struct AlignOptions {
    // The @PG CL: value: "flashalign" (not argv[0], which carries the install
    // path) and argv[1..] joined by spaces, as minimap2 writes it.
    std::string command_line;
    // <ref.fa|ref.faix>. The parser does no I/O; run_align decides which of
    // the two it is from the file's first bytes.
    std::string target_path;
    std::vector<std::string> reads_paths;
    std::string output_path = "-";
    std::string format = "paf";
    bool no_header = false;
    // --sam-hit-only: omit unmapped reads from SAM. PAF has them only
    // with --paf-no-hit.
    bool sam_hit_only = false;
    bool paf_no_hit = false;
    bool paf_cigar = false;
    // --cs[=short|long] / --MD: minimap2 difference strings, DNA presets only.
    // Both imply CIGAR realization, so PAF gets cg:Z too. cs is "" when not
    // requested, otherwise "short" or "long".
    std::string cs;
    bool emit_md = false;
    // --eqx: =/X CIGAR operators instead of M, as minimap2 --eqx.
    bool emit_eqx = false;
    // -Y/--soft-clip-supp: soft-clip supplementary records, keeping the full
    // SEQ. The default is hard clipping, as in minimap2.
    bool soft_clip_supp = false;
    // -y: copy the FASTA/Q comment onto the output records (minimap2 -y),
    // verbatim. uBAM input has no comment.
    bool copy_comment = false;
    // --secondary yes|no: emit one record per mapped alternative (FLAG 0x100,
    // tp:A:S) in every format. The primary's XA:Z lists them either way,
    // output-only DNA secondaries (-N >= 2) aside.
    bool output_secondary = false;
    // -R/--rg: the unescaped @RG header line and its ID; empty when not given.
    std::string read_group_line;
    std::string read_group_id;
    bool progress = false;
    bool quiet = false;
    bool stats = false;        // --stats: print the human run-summary block to stderr
    bool show_config = false;  // --show-config: print the resolved config and exit 0
    // The resolved preset: the builtin "lr", overridden by the preset the
    // target .faix records, overridden by an explicit -x. preset_source says
    // which, for --show-config.
    std::string preset = "lr";
    std::string preset_source = "builtin";
    std::optional<int> threads;
    // The reader closes a batch after the cumulative input length reaches this
    // budget. A single over-budget read forms a singleton batch.
    int64_t batch_bp = 500'000'000;   // -K: bp per alignment batch [default 500M]
    // --batch-window: batches in flight in the compute stage. With two or
    // more, a worker that finishes batch k starts on k+1 instead of waiting.
    // Output is identical for every value.
    int batch_window = 3;
    // --io-staging: MiB of decompressed input the parallel gzip decoder may
    // hold ahead of the reader, as a ceiling; 0 scales it with -t.
    int io_staging_mib = 0;
    std::optional<int> min_support;
    // --vote-seeds INT: the vote's seeds per strand, >= 0; 0 is every seed.
    std::optional<int> vote_seeds;
    // --max-cands INT: the vote's peaks and the catalogue's candidates per
    // strand, 1..64, or 1..16383 in the all-chains lane
    // (check_max_cands_range). DNA presets only.
    std::optional<int> max_cands;
    // --tiles INT: the query partition's tiles per read, 2..4096. DNA
    // presets only.
    std::optional<int> tiles;
    // --tile-owner span|anchors: placement's first tile-ownership rule,
    // stored as "span" or "anchors". DNA presets only.
    std::optional<std::string> tile_owner;
    // -m: minimap2's minimal chain score, >= 1. DNA presets only.
    std::optional<int> min_chain_score;
    // --dual yes|no: whether an overlap preset prints a pair from both of its
    // reads or, with no, only from the read whose name sorts first. Unset, the
    // preset's value. Overlap presets only.
    std::optional<bool> dual;
    std::optional<int> tile_supported_reward;
    std::optional<int> tile_block_open_cost;
    std::optional<int> tile_null_cost;
    std::optional<int> tile_unsupported_cost;
    std::optional<int> vote_diag_bin_width;
    // --dw-slope-den / --dw-max: coefficients of the adaptive vote-width model.
    // Unlike --dw they do not switch it off. DNA presets only.
    std::optional<int> vote_diag_slope_den;
    std::optional<int> vote_diag_width_max;
    // RNA intron bounds; -G takes a k/m/g suffix.
    std::optional<int> min_intron;
    std::optional<int> max_intron;
    // -u {f,b,r,n} (RNA only), stored as forward, auto, reverse or none.
    // auto is the preset default; none scores no splice motif and writes no
    // ts:A.
    std::optional<std::string> splice_strand;
    std::optional<std::string> rna_junction_bed;
    std::optional<int> rna_junction_bonus;
    // -p: the DNA alternative's credibility ratio, the RNA rival retention
    // ratio.
    std::optional<double> pri_ratio;
    // RNA: the secondary retention band (--rival-min-diff) and the catalogue
    // depth (-N + 1), which also sets the realization budget.
    std::optional<int> rna_rival_min_diff;
    std::optional<int> rna_max_loci;
    // DNA: -N's count of alternatives realized per read; unset without
    // secondary output.
    std::optional<int> dna_alternative_realize_max;
    // --vote-ratio: admit vote peaks by their ratio to the read's best vote,
    // in [0,1]; 0 admits by count. DNA presets only.
    std::optional<double> dna_vote_admission_ratio;
    // --max-chain-occ N: global occurrence gate on the dense chain's anchor
    // pool; 0 is no gate. Unset, the vote's cap applies. DNA presets only.
    std::optional<int> max_chain_occ;
    // --max-vote-occ INT: the vote's seed occurrence cap. Unset, the cap is
    // max(the preset's value, the index's mid-occurrence quantile); 0 turns
    // filtering off and >0 fixes the cap. On a DNA preset the chain's
    // occurrence thresholds follow it.
    std::optional<int> max_vote_occ;
    // DP scoring. A single -O or -E value sets both affine gap costs.
    std::optional<int> dp_match;       // -A          -> cigar_dp_match_
    std::optional<int> dp_mismatch;    // -B          -> cigar_dp_mismatch_
    std::optional<int> dp_score_n;     // --score-N   -> cigar_dp_ambi_
    std::optional<int> dp_gap_open1;   // -O (1st)    -> cigar_dp_gap_open1_
    std::optional<int> dp_gap_open2;   // -O (2nd)    -> cigar_dp_gap_open2_
    std::optional<int> dp_gap_extend1; // -E (1st)    -> cigar_dp_gap_extend1_
    std::optional<int> dp_gap_extend2; // -E (2nd)    -> cigar_dp_gap_extend2_
    std::optional<int> dp_zdrop;       // -z (1st)    -> cigar_dp_tail_zdrop_
    std::optional<int> dp_zdrop_inv;   // -z (2nd)    -> cigar_dp_inversion_zdrop_
    std::optional<int> dp_end_bonus;   // --end-bonus -> cigar_dp_tail_end_bonus_
    std::optional<int> dp_min_score;   // -S (mm2 -s) -> cigar_dp_min_dp_max_
    // DP bandwidths. A single -r value leaves the long-join bandwidth unset.
    std::optional<int> dp_bw;          // -r (1st)    -> cigar_dp_bw_
    std::optional<int> dp_bw_long;     // -r (2nd)    -> cigar_dp_bw_long_
    // -g: the DP's maximum gap, and on a DNA preset the chain's too.
    std::optional<int> dp_max_gap;
    // Every option the command line gave, filled where the parser dispatches
    // a spelling.
    std::set<OptionId> given_options;
};

struct IndexOptions {
    std::string ref_path;
    std::string out_path;
    std::string preset = "lr";
    std::optional<int> k;          // unset => k comes from the preset
    std::optional<int> threads;    // unset => build with all available
    std::optional<int> syncmer_s;  // >0 overrides the preset-derived closed-syncmer s
    bool no_seq = false;           // --idx-no-seq: write no reference payload
    // -I: at most this many reference bases per index part; unset is one part.
    std::optional<int64_t> batch_bp;
    bool quiet = false;
};

AlignOptions parse_align_args(int argc, char** argv, int start);
// Throws UsageError unless --max-cands is within the lane bound's range for
// opt.preset: 1..kMaxCatalogueLaneBound, or 1..kAllChainsLaneBound where
// an overlap preset selects the all-chains lane.
void check_max_cands_range(const AlignOptions& opt);
// The refusal of an option in `given` whose stage (OptionSpec::stage) the lane
// of `preset` does not run, or none. Every option of that stage is named, in
// registry order.
std::optional<std::string> option_stage_refusal(
    const std::set<OptionId>& given, const std::string& preset);
IndexOptions parse_index_args(int argc, char** argv, int start);

}  // namespace fa::cpu::cli
