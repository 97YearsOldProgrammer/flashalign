// The CLI's parsed command lines, option registry, errors, output writer and report declarations.
#pragma once

#include "../core/core.h"
#include "../io/io.h"
#include "api/api.h" // LaneStage
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

/**** errors ****/
// UsageError means the command line is wrong (main() exits 2 with a hint);
// any other exception means the run failed (main() exits 1).

namespace fa::cpu::cli {

class UsageError : public std::invalid_argument {
public:
  explicit UsageError(const std::string& what) : std::invalid_argument(what) {}
  explicit UsageError(const char* what) : std::invalid_argument(what) {}
};

} // namespace fa::cpu::cli

/**** option_registry ****/
// The CLI option registry: one table drives both the parser (parse.cpp) and
// the help screens (parse.cpp).
//
// Single letters keep minimap2's meanings, with one exception: -s is the
// closed-syncmer size, so minimap2's minimum DP score is -S here.
// Flash-specific tuning is preset-owned and reported by --show-config.

namespace fa::cpu::cli {

// The field a spelling assigns into. Several spellings may share one id (e.g.
// -x and --preset).
enum class OptionId {
  // clang-format off: the ids are grouped by CLI section.
    Help,
    // align: input/output
    Output, PafCigar, NoHeader, SamHitOnly, PafNoHit, SoftClipSupp,
    CopyComment,
    Secondary, ReadGroup,
    // align: the overlap presets' --dual, minimap2's spelling
    Dual,
    // -a, minimap2's "output SAM".
    OutputSam,
    Cs, Md, Eqx,
    // seeding (align and index)
    Preset, K, SyncmerS,
    // index: build the .faix without reference bases. Names and lengths are
    // still written, so the index supports plain PAF only.
    IdxNoSeq,
    // index: -I, at most NUM reference bases (whole contigs) per index part;
    // align maps one part at a time.
    IndexBatchBp,
    // align: runtime
    Threads, BatchBp, BatchWindow, IoStaging, Progress, Quiet,
    // align: seed occurrence caps for the vote and the dense chain, and the
    // read's own repeats
    MaxVoteOcc, MaxChainOcc, QOccFrac, OccDist,
    // align: minimal chain score (minimap2 -m)
    MinChainScore,
    // align: vote peak admission
    VoteRatio,
    VoteDiagBinWidth,
    VoteDiagSlopeDen, VoteDiagWidthMax,
    MinSupport,
    // align: the vote's seeds per strand and its candidates per strand
    VoteSeeds, MaxCands,
    // align: DP scoring (minimap2 -A/-B/-O/-E/-z, --score-N, --end-bonus).
    // DpMinScore is minimap2's -s, spelled -S here.
    DpMatch, DpMismatch, DpGapOpen, DpGapExtend, DpZdrop, DpScoreN, DpEndBonus,
    DpMinScore,
    // align: band geometry (minimap2 -r) and maximum gap (minimap2 -g).
    DpBw, DpMaxGap,
    // the DNA dense chain passes' band (-b)
    ChainBand,
    // RNA intron bounds, strand and junction annotation
    MinIntron, MaxIntron, SpliceStrand, RnaJunctionBed, RnaJunctionBonus,
    // Secondary retention ratio (-p), the RNA retention band
    // (--rival-min-diff) and the catalogue depth and realization budget (-N)
    PriRatio, SpliceRivalMinDiff, SpliceMaxLoci,
    // RNA: one vote seed per this many query bases
    VoteSpacing,
    // align: reporting
    Stats, ShowConfig,
  // clang-format on
};

// Whether a spelling consumes an argument. The typed parse is done by
// parse.cpp's dispatch switch.
enum class ValueKind {
  None,    // boolean flag, no value
  Int,     // single integer
  IntPair, // INT or INT,INT
  Size,    // base-pair size with optional k/m/g suffix
  Str,     // free string
};

// --help shows Stable rows. Dev rows parse but are never shown.
enum class HelpTier {
  Stable,
  Dev,
};

// Which subcommands accept a spelling.
enum ModeMask : unsigned {
  ModeNone = 0,
  ModeAlign = 1u << 0,
  ModeIndex = 1u << 1,
};

using ::fa::cpu::api::LaneStage;

// One row of the registry. Either spelling may be absent.
struct OptionSpec {
  OptionId id;
  char short_name;            // '\0' if none
  std::string_view long_name; // includes "--"; "" if none
  ValueKind kind;
  HelpTier tier;
  unsigned modes; // bitmask of ModeMask
  std::string_view metavar; // value placeholder, "" if none
  std::string_view section; // help section header, "" for Dev/hidden rows
  std::string_view help;    // description; '\n' forces a wrapped help line
  // Dev rows only: why the row is off the help screen.
  std::string_view study = "";
  // The stage the option configures; refused where the preset's lane does
  // not run it (cli/parse.cpp option_stage_refusal).
  LaneStage stage = LaneStage::None;
};

} // namespace fa::cpu::cli

/**** parse ****/

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
    // tp:A:S) in every format. The primary's XA:Z lists them either way.
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
    // --batch-window: batches in flight in the compute stage; unset, the
    // preset's (options/options.h preset_batch_window). With two or more, a
    // worker that finishes batch k starts on k+1 instead of waiting. Output is
    // identical for every value.
    std::optional<int> batch_window;
    // --io-staging: MiB of decompressed input the parallel gzip decoder may
    // hold ahead of the reader, as a ceiling; 0 scales it with -t.
    int io_staging_mib = 0;
    std::optional<int> min_support;
    // --vote-seeds INT: the vote's seeds per strand, >= 0; 0 is every seed.
    // DNA presets only.
    std::optional<int> vote_seeds;
    // --vote-spacing INT: one vote seed per INT query bases above the seeds
    // per strand, >= 0; 0 keeps the seeds per strand. RNA presets only.
    std::optional<int> vote_spacing;
    // --max-cands INT: the vote's peaks and the catalogue's candidates per
    // strand, 1..64, or 1..16383 in the all-chains lane
    // (check_max_cands_range). DNA presets only.
    std::optional<int> max_cands;
    // -m: minimap2's minimal chain score, >= 1. DNA presets only.
    std::optional<int> min_chain_score;
    // --dual yes|no: whether an overlap preset prints a pair from both of its
    // reads or, with no, only from its shorter read, ties broken by name.
    // Unset, the preset's value. Overlap presets only.
    std::optional<bool> dual;
    std::optional<int> vote_diag_bin_width;
    // --dw-slope-den / --dw-max: coefficients of the adaptive vote-width model.
    // Unlike --dw they do not switch it off. DNA presets only.
    std::optional<int> vote_diag_slope_den;
    std::optional<int> vote_diag_width_max;
    // RNA intron bounds; -G takes a k/m/g suffix and drops a -r typed before
    // it.
    std::optional<int> min_intron;
    std::optional<int> max_intron;
    // -u {f,b,r,n} (RNA only), stored as forward, auto, reverse or none.
    // auto is the preset default; none scores no splice motif and writes no
    // ts:A.
    std::optional<std::string> splice_strand;
    std::optional<std::string> rna_junction_bed;
    std::optional<int> rna_junction_bonus;
    // -p: the DNA -c lane's secondary credibility ratio, the RNA rival
    // retention ratio.
    std::optional<double> pri_ratio;
    // RNA: the secondary retention band (--rival-min-diff) and the catalogue
    // depth (-N + 1), which also sets the realization budget.
    std::optional<int> rna_rival_min_diff;
    std::optional<int> rna_max_loci;
    // DNA: -N, minimap2's best_n, the cap on secondary alignments.
    std::optional<int> dna_secondary_max;
    // --vote-ratio: admit vote peaks by their ratio to the read's best vote,
    // in [0,1]; 0 admits by count. DNA presets only.
    std::optional<double> dna_vote_admission_ratio;
    // --max-chain-occ N: global occurrence gate on the dense chain's anchor
    // pool; 0 is no gate. Unset, the vote's cap applies. DNA presets only.
    std::optional<int> max_chain_occ;
    // --q-occ-frac F: minimap2's q_occ_frac, the within-read repeat
    // filter's fraction of the read's seeds, in [0,1]; 0 turns the filter
    // off. Every preset but the overlap ones.
    std::optional<double> q_occ_frac;
    // -e INT: minimap2's occ_dist, the seed selection's streak spacing, >= 0;
    // 0 turns the selection off. Every preset.
    std::optional<int> occ_dist;
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
    // -g: the DP's maximum gap and the chains' gap (on a splice preset the
    // fine chain's query gap).
    std::optional<int> dp_max_gap;
    // -b NUM: the dense chain passes' band, >= 1. DNA presets only.
    std::optional<int> chain_band;
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

/**** help ****/

namespace fa::cpu::cli {

void print_top_help(std::ostream& out);
void print_align_help(std::ostream& out);
void print_index_help(std::ostream& out);

}  // namespace fa::cpu::cli

/**** read_group ****/
// -R/--rg validation and unescaping, as minimap2's sam_write_rg_line and
// mm_escape, except that an unrecognized escape is an error instead of being
// silently dropped.

namespace fa::cpu::cli {

struct ReadGroup {
  std::string line; // unescaped @RG header line (real tab characters)
  std::string id;   // the ID: field value, for RG:Z and duplicate detection
};

// The ID declared by an @RG header line (e.g. a uBAM passthrough line), or an
// empty string when the line declares none.
std::string read_group_line_id(const std::string& line);

} // namespace fa::cpu::cli

/**** output_writer ****/

namespace fa::cpu::cli {

// The tag keys FlashAlign writes on SAM records. An input tag with one of
// them is dropped from the passthrough in favour of ours, with one warning
// per key.
class OwnTagKeys {
 public:
    explicit OwnTagKeys(std::vector<std::string> keys);
    // True when the two characters at `key` are one of ours.
    bool drop(const char* key) const;

 private:
    std::vector<std::string> keys_;
    std::unique_ptr<std::atomic<bool>[]> warned_;
};

struct AlignmentOutputContext {
    std::ostream& text;
    const std::unordered_map<std::string, int64_t>* reference_lengths = nullptr;
    std::string format;
    bool include_unmapped = false;
    bool paf_cigar = false;
    // -y: copy the read's FASTA/Q comment onto the record (minimap2 -y).
    bool copy_comment = false;
    output::SamEmitOptions sam;
    const OwnTagKeys* own_tag_keys = nullptr;
};

// The header lines a run carries over from its BAM inputs, in input order:
// every @PG line (the first of each ID), the union of the @RG lines by ID,
// and every @CO line.
struct InputHeader {
    std::vector<std::string> program_lines;
    // @RG and @CO lines, in the order read.
    std::vector<std::string> passthrough;
    // The ID of the last @PG line of the first input that has one: our @PG's PP.
    std::string previous_program;
};

class AlignmentOutputWriter {
 public:
    AlignmentOutputWriter(
        const std::string& path,
        const std::string& format,
        const io::GenomeLoad& genome,
        // @PG CL: value; empty omits the field.
        const std::string& command_line = {},
        // @RG line to declare, already validated and unescaped; empty means
        // no read group.
        const std::string& read_group_line = {});
    ~AlignmentOutputWriter();

    AlignmentOutputWriter(const AlignmentOutputWriter&) = delete;
    AlignmentOutputWriter& operator=(const AlignmentOutputWriter&) = delete;

    void write_header(
        const InputHeader& input_header,
        bool suppress_sam_header);
    // `own_tag_keys`: the SAM tag keys this run writes (see OwnTagKeys).
    AlignmentOutputContext context(
        bool include_unmapped,
        bool paf_cigar,
        bool copy_comment,
        output::SamEmitOptions sam,
        std::vector<std::string> own_tag_keys);
    void close();
    void flush();

 private:
    std::string format_;
    std::string command_line_;
    std::string read_group_line_;
    std::string read_group_id_;
    // Declared first, so the stream below is destroyed (and flushed) before it.
    std::unique_ptr<std::streambuf> stdout_sink_;
    std::unique_ptr<std::ostream> owner_;
    std::ostream* text_ = nullptr;
    std::vector<output::SamReference> references_;
    std::unordered_map<std::string, int64_t> reference_lengths_;
    std::unique_ptr<OwnTagKeys> own_tag_keys_;
};

}  // namespace fa::cpu::cli

/**** report ****/

namespace fa::cpu::cli {

// What --stats prints. Over a multi-part index the read counts sum over the
// passes; bases are counted once.
struct RunSummary {
  std::string preset;
  std::string output_format;
  bool rna = false;
  int64_t reads = 0;
  int64_t bases = 0;
  int64_t mapped_reads = 0;
  int64_t records = 0; // every record written, unmapped SAM rows included
  int64_t supplementary_records = 0;
  int64_t secondary_records = 0; // written only with --secondary yes
  // RNA: reads whose primary alignment has an intron, and the introns (N
  // operators) of all primary alignments.
  int64_t spliced_reads = 0;
  int64_t junctions = 0;
  double wall_seconds = 0.0; // from before the index is loaded
  double map_seconds = 0.0;  // the mapping passes alone, for throughput
  int64_t peak_rss_bytes = 0;
};

// The closing stderr lines of every align and index run: version, the command
// line as typed, and real time, CPU time and peak RSS.
void write_run_footer(int argc, char** argv, double real_seconds);

}  // namespace fa::cpu::cli

/**** align_command ****/

namespace fa::cpu::cli {

int run_align(const AlignOptions& opt);

} // namespace fa::cpu::cli
