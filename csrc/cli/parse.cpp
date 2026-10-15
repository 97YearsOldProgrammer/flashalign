// Option parsing: the registry table that drives the parser and the help screens, and -R/--rg validation.
#include "cli/cli.h"

#include "api/api.h" // preset_is_valid

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**** option_registry ****/

namespace fa::cpu::cli {

namespace {

// Row fields, as in OptionSpec:
//   { id, short, long, kind, tier, modes, metavar, section, help, study,
//     stage }
constexpr ValueKind kNone = ValueKind::None;
constexpr ValueKind kInt = ValueKind::Int;
constexpr ValueKind kPair = ValueKind::IntPair;
constexpr ValueKind kSize = ValueKind::Size;
constexpr ValueKind kStr = ValueKind::Str;

constexpr HelpTier kStable = HelpTier::Stable;
constexpr HelpTier kDev = HelpTier::Dev;

constexpr unsigned A = ModeAlign;
constexpr unsigned I = ModeIndex;

constexpr LaneStage kBaseOutput = LaneStage::BaseLevelOutput;
constexpr LaneStage kBaseAlign = LaneStage::BaseAlignment;
constexpr LaneStage kSamOutput = LaneStage::SamOutput;
constexpr LaneStage kSelection = LaneStage::Selection;
constexpr LaneStage kCandidateChains = LaneStage::CandidateChains;

// Section headers of the align screen, in screen order. The index screen is
// one list under "Options:".
constexpr std::string_view SEC_PLACEMENT = "Placement:";
constexpr std::string_view SEC_ALIGN = "Alignment:";
constexpr std::string_view SEC_SPLICE = "Splice:";
constexpr std::string_view SEC_IO = "Input/Output:";
constexpr std::string_view SEC_PRESET = "Preset:";
constexpr std::string_view SEC_IOPT = "options:"; // index subcommand
constexpr std::string_view HIDDEN = "";           // Dev rows render nowhere

// Why a Dev row is off the help screen (OptionSpec::study).
constexpr std::string_view STUDY_CLI_MANUAL =
    "documented in flashalign.1, not on the help screen";

const std::vector<OptionSpec>& specs_table() {
  // clang-format off: preserve the hand-aligned option table.
    static const std::vector<OptionSpec> kSpecs = {
        // align, shown by --help
        // minimap2's -g. On a DNA preset it also bounds the chains' reach,
        // which never falls below a chain's band (-b); on a splice preset it
        // also sets the fine chain's query gap, which the band does not
        // raise.
        {OptionId::DpMaxGap, 'g', "", kSize, kStable, A,
         "NUM", SEC_PLACEMENT,
         "stop alignment elongation if there are no seeds in NUM-bp [5000]"},
        // minimap2's -r. On a DNA preset NUM1 bands the read-end extensions
        // and NUM2 the gap fills; neither bands a chain. On a splice preset
        // NUM1 bands the fine chain and NUM2 only bounds NUM1.
        {OptionId::DpBw, 'r', "", kPair, kStable, A,
         "NUM[,NUM]", SEC_PLACEMENT,
         "alignment bandwidth and gap-fill bandwidth [500,20000]", "", kBaseAlign},
        // NUM bands each catalogue candidate's dense (whole-query) chain
        // pass. A size, as -r's. DNA presets but the overlap ones; the asm
        // presets' NUM is 100000.
        {OptionId::ChainBand, 'b', "--chain-band", kSize, kStable, A,
         "NUM", SEC_PLACEMENT,
         "chaining bandwidth of the dense chains [20000]", "",
         kCandidateChains},
        // minimap2's -m. DNA presets only; the splice presets keep their own
        // chain floors.
        {OptionId::MinChainScore, 'm', "", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "minimal chaining score (matching bases minus log gap penalty) [40]"},
        {OptionId::PriRatio, 'p', "", kStr, kStable, A,
         "FLOAT", SEC_PLACEMENT,
         "min secondary-to-primary score ratio [0.8]", "", kSelection},
        // As in minimap2: -N 5 keeps the best locus and at most five rivals,
        // and -N 0 is --secondary no. On a splice preset up to INT rivals are
        // realized; on a DNA preset the -c lane retains at most INT
        // secondary alignments.
        {OptionId::SpliceMaxLoci, 'N', "", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "retain at most INT secondary alignments [1; splice 5]", "", kSelection},
        {OptionId::MinSupport, '\0', "--min-support", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "minimal number of seeds on a vote peak [3]"},
        // Unset, the cap is max(the preset's cap, the index's mid-occurrence
        // quantile); 0 disables it and a positive value fixes it. On a DNA
        // preset the chain's occurrence thresholds follow it.
        {OptionId::MaxVoteOcc, '\0', "--max-vote-occ", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "max seed occurrences in voting [200]"},
        // Keeps every vote peak within ratio R of the read's best vote; 0
        // admits by count instead. A float, so parsed from a string like -p.
        // DNA presets only.
        {OptionId::VoteRatio, '\0', "--vote-ratio", kStr, kStable, A,
         "FLOAT", SEC_PLACEMENT,
         "min vote-to-best-vote ratio to keep a peak [0.25]"},
        // Drops seeds with more than N reference occurrences from the dense
        // chain's anchor pool, except the seeds the vote rescued. Unset, N is
        // the vote's cap; 0 is no gate. DNA presets only. Only this gate
        // bounds the pool; at 0 nothing does.
        {OptionId::MaxChainOcc, '\0', "--max-chain-occ", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "max seed occurrences in chaining; 0 to disable, leaving repeats unbounded [200]"},
        // The vote bin width is W = clamp(L / den, dw, dw-max) for a read of
        // length L. --dw fixes W; --dw-slope-den and --dw-max are DNA only.
        {OptionId::VoteDiagBinWidth, '\0', "--dw", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "min diagonal bin width in voting; set to fix the width [64]"},
        {OptionId::VoteDiagSlopeDen, '\0', "--dw-slope-den", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "bin width grows as read length / INT; 0 to disable [128]"},
        {OptionId::VoteDiagWidthMax, '\0', "--dw-max", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "max diagonal bin width in voting [2048]"},
        // The Python Config.max_query_seeds. DNA presets only; the splice
        // presets take --vote-spacing. Above 128 a DNA preset keeps the
        // 128-seed selection and fills it from the larger one
        // (seed/seed.h).
        {OptionId::VoteSeeds, '\0', "--vote-seeds", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "seeds per strand in voting; 0 for all [128]"},
        // The lane bound, DNA presets only, 1..64: the vote's peaks and the
        // catalogue's candidates per strand. Unset, 16 under ratio admission
        // and 4 under count admission (--vote-ratio 0).
        {OptionId::MaxCands, '\0', "--max-cands", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "max placement candidates per strand [16]"},

        // On lr and lr:hq -A -B -O -E -z --score-N set the DP row of the gap
        // fills between anchors; the read ends keep the preset's own row,
        // which also prices every path. On a splice or assembly preset they
        // set its one row. Defaults are lr's.
        {OptionId::DpMatch, 'A', "", kInt, kStable, A,
         "INT", SEC_ALIGN, "matching score [4]", "", kBaseAlign},
        {OptionId::DpMismatch, 'B', "", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "mismatch penalty (larger value for lower divergence) [8]", "", kBaseAlign},
        {OptionId::DpGapOpen, 'O', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN, "gap open penalty [8,48]", "", kBaseAlign},
        {OptionId::DpGapExtend, 'E', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN,
         "gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [4,1]", "", kBaseAlign},
        {OptionId::DpZdrop, 'z', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN,
         "Z-drop score and inversion Z-drop score [800,200]", "", kBaseAlign},
        // minimap2's -s. RNA prints a supplementary piece only when its DP
        // maximum reaches it; a primary below it keeps its placement. DNA uses it as
        // the per-record emission floor when a CIGAR is realized, as
        // minimap2's mm_filter_regs; map-only ignores it.
        {OptionId::DpMinScore, 'S', "", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "minimal peak DP alignment score [80]", "", kBaseAlign},
        {OptionId::DpScoreN, '\0', "--score-N", kInt, kStable, A,
         "INT", SEC_ALIGN, "penalty of a mismatch involving ambiguous bases [2]", "", kBaseAlign},
        {OptionId::DpEndBonus, '\0', "--end-bonus", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "score bonus when alignment extends to the end of the query sequence [-1]", "", kBaseAlign},

        // RNA presets only.
        {OptionId::MinIntron, '\0', "--min-intron", kInt, kStable, A,
         "INT", SEC_SPLICE,
         "min intron length [20]"},
        // As minimap2's -G, it also sets both -r values; between -G and -r
        // the later typed wins.
        {OptionId::MaxIntron, 'G', "", kSize, kStable, A,
         "NUM", SEC_SPLICE,
         "max intron length (effective with -xsplice; changing -r) [200k]"},
        {OptionId::SpliceStrand, 'u', "", kStr, kStable, A,
         "CHAR", SEC_SPLICE,
         "how to find GT-AG. f:transcript strand, b:both strands, r:reverse strand, n:don't match GT-AG [b]"},
        {OptionId::RnaJunctionBed, '\0', "--junc-bed", kStr,
         kStable, A, "FILE", SEC_SPLICE,
         "junctions to prefer during base alignment []"},
        {OptionId::RnaJunctionBonus, '\0', "--junc-bonus", kInt,
         kStable, A, "INT", SEC_SPLICE,
         "score bonus for a splice donor or acceptor found in annotation [9]"},
        // A read of length L votes max(128, ceil(L / INT)) seeds per strand.
        // The Python Config.vote_seed_spacing.
        {OptionId::VoteSpacing, '\0', "--vote-spacing", kInt, kStable, A,
         "INT", SEC_SPLICE,
         "one vote seed per INT query bases, at least 128 per strand; 0 for a fixed 128 [8; splice:hq 16]"},

        {OptionId::Output,   'o', "--output", kStr, kStable, A,
         "FILE", SEC_IO, "output alignments to FILE [stdout]"},
        {OptionId::OutputSam, 'a', "", kNone, kStable, A,
         "", SEC_IO, "output in the SAM format (PAF by default)", "",
         kBaseOutput},
        {OptionId::PafCigar, 'c', "", kNone, kStable, A,
         "", SEC_IO, "output CIGAR in PAF", "", kBaseOutput},
        // As in minimap2 the value can only be attached (--cs=long); a bare
        // --cs means short.
        {OptionId::Cs, '\0', "--cs", kNone, kStable, A,
         "[=STR]", SEC_IO,
         "output the cs tag; STR is 'short' (if absent) or 'long' [none]", "",
         kBaseOutput},
        {OptionId::Md, '\0', "--MD", kNone, kStable, A,
         "", SEC_IO,
         "output the MD tag", "", kBaseOutput},
        // Applies to every realized CIGAR (SAM, SA:Z, cg:Z); plain PAF
        // has none.
        {OptionId::Eqx, '\0', "--eqx", kNone, kStable, A,
         "", SEC_IO,
         "write =/X CIGAR operators", "", kBaseAlign},
        {OptionId::SoftClipSupp, 'Y', "--soft-clip-supp", kNone, kStable, A,
         "", SEC_IO, "use soft clipping for supplementary alignments", "", kSamOutput},
        {OptionId::ReadGroup, 'R', "--rg", kStr, kStable, A,
         "STR", SEC_IO,
         "SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []", "", kSamOutput},
        {OptionId::SamHitOnly, '\0', "--sam-hit-only", kNone, kStable, A,
         "", SEC_IO,
         "in SAM, don't output unmapped reads", "", kSamOutput},
        {OptionId::PafNoHit, '\0', "--paf-no-hit", kNone, kStable, A,
         "", SEC_IO,
         "in PAF, output unmapped queries; the strand and the reference name fields are set to '*'"},
        // Also --secondary yes, the value as a separate word. yes emits one
        // record per mapped alternative (FLAG 0x100, tp:A:S, MAPQ 0, no SEQ); either
        // way the alternatives appear on the primary's XA:Z, md:i and s2:i.
        {OptionId::Secondary, '\0', "--secondary", kStr, kStable, A,
         "=yes|no", SEC_IO,
         "whether to output secondary alignments [no]", "", kSelection},
        {OptionId::NoHeader, '\0', "--no-header", kNone, kStable, A,
         "", SEC_IO,
         "don't output the SAM header", "", kSamOutput},
        // PAF and SAM copy the comment verbatim, as minimap2 does.
        {OptionId::CopyComment, 'y', "", kNone, kStable, A,
         "", SEC_IO,
         "copy FASTA/Q comments to output SAM"},
        {OptionId::Threads, 't', "", kInt, kStable, A,
         // Unset, two cores are left for the reader and writer. Parallel
         // decoding shares the -t budget.
         "INT", SEC_IO,
         "number of threads [all-2]"},
        {OptionId::BatchBp, 'K', "--batch-bp", kSize, kStable, A,
         "NUM", SEC_IO,
         "minibatch size for mapping [500M]"},
        {OptionId::ShowConfig, '\0', "--show-config", kNone, kStable, A,
         "", SEC_IO,
         "print the resolved configuration and exit"},
        {OptionId::Stats, '\0', "--stats", kNone, kStable, A,
         "", SEC_IO, "print a run summary to stderr"},
        {OptionId::Progress, '\0', "--progress", kNone, kStable, A,
         "", SEC_IO, "print progress to stderr"},
        {OptionId::Quiet, '\0', "--quiet", kNone, kStable, A,
         "", SEC_IO, "suppress non-error logging"},
        {OptionId::Help, '\0', "--help", kNone, kStable, A,
         "", SEC_IO, "print this help message"},

        {OptionId::Preset, 'x', "--preset", kStr, kStable, A,
         "STR", SEC_PRESET,
         "preset: lr, lr:hq, splice, splice:hq; experimental: asm5, asm10, asm20, ava-ont, ava-hifi (see 'flashalign index') [lr]"},

        // align, not on the help screen; the manual page lists these.
        // --batch-window INT [3]: batches in flight in the compute stage.
        // Output is identical for every value.
        {OptionId::BatchWindow, '\0', "--batch-window", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // --io-staging MIB [0]: ceiling on the parallel gzip decoder's
        // read-ahead; 0 scales it with -t.
        {OptionId::IoStaging, '\0', "--io-staging", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // minimap2's long forms of -G, -S (its -s) and -m.
        {OptionId::MaxIntron, '\0', "--max-intron-len", kSize, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        {OptionId::DpMinScore, '\0', "--min-dp-score", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        {OptionId::MinChainScore, '\0', "--min-chain-score", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // RNA presets only: --rival-min-diff INT [30], the absolute retention
        // band beside -p.
        {OptionId::SpliceRivalMinDiff, '\0', "--rival-min-diff", kInt, kDev,
         A, "", HIDDEN, "", STUDY_CLI_MANUAL},

        // --dual yes|no, also --dual=yes|no: the overlap presets only. no
        // prints a pair once, from its shorter read, ties broken by name
        // (options/options.h dual).
        {OptionId::Dual, '\0', "--dual", kStr, kDev, A, "", HIDDEN, "",
         STUDY_CLI_MANUAL, LaneStage::OverlapPairs},

        // index, shown by --help
        {OptionId::Preset, 'x', "--preset", kStr, kStable, I,
         "STR", SEC_IOPT,
         // (k, s) per preset: resolve_preset_seeding() in options/options.cpp.
         "preset [lr]\n"
         "  lr         noisy long reads (Nanopore) vs reference mapping                 [-k21 -s9]\n"
         "  lr:hq      accurate long reads (HiFi, error rate <1%) vs reference mapping  [-k21 -s5]\n"
         "  splice     spliced alignment for long RNA reads                             [-k15 -s10]\n"
         "  splice:hq  spliced alignment for accurate long RNA reads                    [-k15 -s10]\n"
         "  asm5       experimental: assembly vs reference, ~0.1% divergence            [-k21 -s9]\n"
         "  asm10      experimental: assembly vs reference, ~1% divergence              [-k21 -s9]\n"
         "  asm20      experimental: assembly vs reference, several % divergence        [-k21 -s9]\n"
         "  ava-ont    experimental: all-vs-all overlap of noisy long reads (Nanopore)  [-k17 -s9]\n"
         "  ava-hifi   experimental: all-vs-all overlap of accurate long reads (HiFi)   [-k21 -s5]"},
        {OptionId::K, 'k', "", kInt, kStable, I,
         "INT", SEC_IOPT, "k-mer size (no larger than 23) [21]"},
        {OptionId::SyncmerS, 's', "", kInt, kStable, I,
         "INT", SEC_IOPT, "closed-syncmer s-mer size [9]"},
        {OptionId::IndexBatchBp, 'I', "", kSize, kStable, I,
         "NUM", SEC_IOPT,
         "split index for every ~NUM input bases [no split]"},
        {OptionId::Threads, 't', "", kInt, kStable, I,
         "INT", SEC_IOPT, "number of threads [all]"},
        {OptionId::Quiet, '\0', "--quiet", kNone, kStable, I,
         "", SEC_IOPT, "suppress non-error logging"},
        {OptionId::Help, '\0', "--help", kNone, kStable, I,
         "", SEC_IOPT, "print this help message"},

        // index, not on the help screen (as in minimap2)
        {OptionId::IdxNoSeq, '\0', "--idx-no-seq", kNone, kDev, I, "", HIDDEN,
         "", STUDY_CLI_MANUAL},
    };
  // clang-format on
  return kSpecs;
}

// The registry, ordered align-Stable (grouped by section), align-Dev,
// index-Stable, index-Dev, so help can print each section header once.
const std::vector<OptionSpec>& option_specs() { return specs_table(); }

// The spec a token selects in `mode`, or nullptr if none matches.
const OptionSpec* find_option_spec(std::string_view token, unsigned mode) {
  for (const OptionSpec& s : specs_table()) {
    if (!(s.modes & mode))
      continue;
    if (s.short_name != '\0' && token.size() == 2 && token[0] == '-' &&
        token[1] == s.short_name) {
      return &s;
    }
    if (!s.long_name.empty() && token == s.long_name)
      return &s;
  }
  return nullptr;
}

} // namespace

} // namespace fa::cpu::cli

/**** help ****/

namespace fa::cpu::cli {

namespace {

// Description column of each screen. A left block that reaches it is followed
// by a single space instead.
constexpr std::size_t kAlignDescCol = 28;
constexpr std::size_t kIndexDescCol = 20;

// Row indent under section headers, and for the index screen's flat list.
constexpr std::size_t kSectionIndent = 4;
constexpr std::size_t kFlatIndent = 2;

// The left-hand block "-s, --long METAVAR" of one option. long_name is empty
// for a flag shown short-only.
std::string left_column(const OptionSpec& s, std::size_t indent) {
    std::string left(indent, ' ');
    bool have_tok = false;
    if (s.short_name != '\0') {
        left += '-';
        left += s.short_name;
        have_tok = true;
    }
    if (!s.long_name.empty()) {
        if (have_tok) left += ", ";
        left += std::string(s.long_name);
        have_tok = true;
    }
    // A metavar starting with '=' or '[' is an attached value (--cs[=STR]).
    if (!s.metavar.empty()) {
        if (s.metavar.front() != '=' && s.metavar.front() != '[') left += ' ';
        left += std::string(s.metavar);
    }
    return left;
}

// One option row. A '\n' in the help text starts a continuation line indented
// to the description column.
void render_row(std::ostream& out, const OptionSpec& s, std::size_t desc_col,
                std::size_t indent) {
    const std::string left = left_column(s, indent);
    const std::size_t gap =
        (left.size() + 1 > desc_col) ? 1 : desc_col - left.size();
    out << left << std::string(gap, ' ');
    const std::string_view help = s.help;
    std::size_t pos = 0;
    bool first = true;
    while (pos <= help.size()) {
        const std::size_t nl = help.find('\n', pos);
        const std::string_view line = help.substr(
            pos, nl == std::string_view::npos ? std::string_view::npos : nl - pos);
        if (!first) out << std::string(desc_col, ' ');
        out << line << '\n';
        first = false;
        if (nl == std::string_view::npos) break;
        pos = nl + 1;
    }
}

// The registry rows for `mode` and `tier`, in order, under their section
// headers. Rows with an empty section are parseable but never shown.
void render_sections(std::ostream& out, unsigned mode, HelpTier tier,
                     bool headers, std::size_t desc_col, std::size_t indent) {
    std::string_view cur_section;
    for (const OptionSpec& s : option_specs()) {
        if (!(s.modes & mode)) continue;
        if (s.tier != tier || s.section.empty()) continue;
        if (headers && s.section != cur_section) {
            out << std::string(indent - 2, ' ') << std::string(s.section)
                << '\n';
            cur_section = s.section;
        }
        render_row(out, s, desc_col, indent);
    }
}

}  // namespace

void print_top_help(std::ostream& out) {
    out
        << "Usage: flashalign <command> <arguments>\n"
        << "Commands:\n"
        << "  index      index reference FASTA\n"
        << "  align      read alignment\n"
        << "  version    print the version number\n";
}

void print_align_help(std::ostream& out) {
    out
        << "Usage: flashalign align [options] <ref.fa|ref.faix> <reads.fq> [...]\n"
        << "Options:\n";
    render_sections(out, ModeAlign, HelpTier::Stable, /*headers=*/true,
                    kAlignDescCol, kSectionIndent);
}

void print_index_help(std::ostream& out) {
    out
        << "Usage: flashalign index [options] <ref.fa> [out.faix]\n"
        << "Options:\n";
    // One flat list under "Options:", as minibwa index prints it.
    render_sections(out, ModeIndex, HelpTier::Stable, /*headers=*/false,
                    kIndexDescCol, kFlatIndent);
}

}  // namespace fa::cpu::cli

/**** read_group ****/

namespace fa::cpu::cli {

// Throws UsageError with a diagnostic naming the defect. An empty
// spec is rejected like any other non-@RG line.
static ReadGroup parse_read_group(const std::string& spec) {
  if (spec.rfind("@RG", 0) != 0) {
    throw UsageError(
        "-R/--rg: the read group line is not started with @RG");
  }
  if (spec.find('\t') != std::string::npos) {
    throw UsageError(
        "-R/--rg: the read group line contained literal <tab> characters -- "
        "replace literal tabs with escaped tabs: \\t");
  }
  ReadGroup group;
  for (std::size_t i = 0; i < spec.size(); ++i) {
    if (spec[i] != '\\') {
      group.line.push_back(spec[i]);
      continue;
    }
    if (i + 1 == spec.size()) {
      throw UsageError(
          "-R/--rg: the read group line ends with a dangling backslash");
    }
    const char escaped = spec[i + 1];
    if (escaped == 't')
      group.line.push_back('\t');
    else if (escaped == '\\')
      group.line.push_back('\\');
    else
      throw UsageError(
          std::string("-R/--rg: unsupported escape '\\") + escaped +
          "' in the read group line; only \\t and \\\\ are recognized");
    ++i;
  }
  const std::size_t id_field = group.line.find("\tID:");
  if (id_field == std::string::npos) {
    throw UsageError("-R/--rg: no ID within the read group line");
  }
  const std::size_t begin = id_field + 4;
  std::size_t end = begin;
  while (end < group.line.size() && group.line[end] != '\t' &&
         group.line[end] != '\n') {
    ++end;
  }
  group.id = group.line.substr(begin, end - begin);
  if (group.id.size() > 255) {
    throw UsageError(
        "-R/--rg: @RG:ID is longer than 255 characters");
  }
  if (group.id.empty()) {
    throw UsageError("-R/--rg: @RG:ID is empty");
  }
  return group;
}

std::string read_group_line_id(const std::string& line) {
  const std::size_t id_field = line.find("\tID:");
  if (id_field == std::string::npos)
    return std::string();
  const std::size_t begin = id_field + 4;
  std::size_t end = begin;
  while (end < line.size() && line[end] != '\t' && line[end] != '\n')
    ++end;
  return line.substr(begin, end - begin);
}

} // namespace fa::cpu::cli

/**** option_values ****/
// Option values: integers, k/m/g sizes, ratios and INT[,INT] pairs; a bad value
// throws UsageError naming the option.

namespace fa::cpu::cli {

namespace {

int parse_int(const std::string& value, const std::string& name) {
    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0') {
        throw UsageError("invalid integer for " + name + ": " + value);
    }
    // Reject overflow of long and of int, so the narrowing cast cannot wrap.
    if (errno == ERANGE || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        throw UsageError(
            "integer out of range for " + name + ": " + value);
    }
    return static_cast<int>(parsed);
}

// Parse INT or INT,INT. Without a comma, return the value in both elements.
std::pair<int, int> parse_int_pair(const std::string& value, const std::string& name) {
    const auto comma = value.find(',');
    if (comma == std::string::npos) {
        const int v = parse_int(value, name);
        return {v, v};
    }
    const int a = parse_int(value.substr(0, comma), name);
    const int b = parse_int(value.substr(comma + 1), name);
    return {a, b};
}

// A base-pair size with an optional decimal k/m/g suffix, as minimap2 -K
// ("100m", "1.5g").
int64_t parse_size(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    int64_t mult = 1;
    if (end != value.c_str() && *end != '\0') {
        switch (*end) {
            case 'k': case 'K': mult = 1000LL; ++end; break;
            case 'm': case 'M': mult = 1000000LL; ++end; break;
            case 'g': case 'G': mult = 1000000000LL; ++end; break;
            default: break;
        }
    }
    if (end == value.c_str() || *end != '\0' || parsed < 0.0) {
        throw UsageError("invalid size for " + name + ": " + value);
    }
    // Range-check before the cast, which is undefined out of range. This also
    // catches strtod's HUGE_VAL.
    const double scaled = parsed * static_cast<double>(mult);
    if (!(scaled <= static_cast<double>(std::numeric_limits<int64_t>::max()))) {
        throw UsageError(
            "integer out of range for " + name + ": " + value);
    }
    return static_cast<int64_t>(scaled);
}

// A size that must fit an int (-r, -b, -g, -G).
int parse_int_size(const std::string& value, const std::string& name) {
    const int64_t size = parse_size(value, name);
    if (size > std::numeric_limits<int>::max()) {
        throw UsageError(
            "integer out of range for " + name + ": " + value);
    }
    return static_cast<int>(size);
}

// A real ratio in [0,1], for -p and --vote-ratio.
double parse_ratio(const std::string& value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0' || !(parsed >= 0.0) || parsed > 1.0)
        throw UsageError("invalid ratio for " + name + ": " + value);
    return parsed;
}

}  // namespace

}  // namespace fa::cpu::cli

/**** minimap2_hints ****/
// minimap2 spellings and preset names this tool does not take, and the hint each
// refusal carries.

namespace fa::cpu::cli {

namespace {

// minimap2 spellings this tool does not take, and what to use instead. The
// refusal stands; the hint only names the way here.
struct Minimap2Hint {
    std::string_view spelling;
    std::string_view hint;
};

constexpr std::string_view kSyncmerHint =
    "seeds are closed syncmers, set when the index is built: flashalign "
    "index -k -s";
constexpr std::string_view kSpliceModelHint =
    "FlashAlign always runs minimap2's default splice model (-J1), where it "
    "has no effect";
constexpr std::string_view kSpliceScoreHint =
    "--junc-bed gives annotated junctions a bonus";
constexpr std::string_view kEndFilterHint =
    "the end filters are always on, at minimap2's defaults: the bad-end trim "
    "on DNA, the terminal-exon filter on RNA";
constexpr std::string_view kAltHint = "map to a reference without ALT contigs";
constexpr std::string_view kOverlapHint =
    "-x ava-ont and -x ava-hifi already leave out each read's own diagonal "
    "(-D), set no primary chain (-P) and chain each candidate once, with no "
    "long-join re-chain (--no-long-join)";

constexpr Minimap2Hint kAlignHints[] = {
    {"-s", "minimap2's -s, the minimal peak DP score, is -S here"},
    {"-f", "minimap2's -f FLOAT is --max-vote-occ INT here, a seed "
           "occurrence cap"},
    {"-k", "k is set when the index is built: flashalign index -k"},
    {"-w", kSyncmerHint},
    {"-U", "the occurrence cap's floor (200) and lr:hq's ceiling (500) are "
           "built in; --max-vote-occ fixes the cap"},
    {"-I", "-I is an index option: flashalign index -I"},
    {"-d", "build the index with flashalign index ref.fa out.faix"},
    {"--split-prefix",
     "map against a one-part index (flashalign index without -I)"},
    {"-L", "samtools moves a CIGAR of over 65535 operations to CG:B,I "
           "when it writes BAM"},
    {"-2", "drop it: input and output already run on their own threads"},
    {"--cap-kalloc", "drop it: there is no such memory cap"},
    {"--cap-sw-mem", "the DP matrix cap is fixed at minimap2's default, 100M "
                     "cells"},
    {"--alt", kAltHint},
    {"--alt-drop", kAltHint},
    {"-J", "FlashAlign always runs minimap2's default splice model, -J1"},
    {"-C", kSpliceModelHint},
    {"--splice-flank", kSpliceModelHint},
    {"--end-seed-pen", kEndFilterHint},
    {"--no-end-flt", kEndFilterHint},
    {"--spsc", kSpliceScoreHint},
    {"--spsc0", kSpliceScoreHint},
    {"--junc-pen", kSpliceScoreHint},
    {"--spsc-scale", kSpliceScoreHint},
    {"--write-junc",
     "paftools.js splice2bed writes the junctions of SAM or PAF -c output"},
    {"-X", "-x ava-ont and -x ava-hifi with --dual=no do what -X does"},
    {"-D", kOverlapHint},
    {"-P", kOverlapHint},
    {"--no-long-join", kOverlapHint},
    {"-H", "FlashAlign's seeds are not homopolymer-compressed"},
    {"-e", "FlashAlign has no such option"},
    {"-n", "FlashAlign has no such option"},
};

// minimap2 preset names (-x) this tool does not take.
constexpr Minimap2Hint kPresetHints[] = {
    {"map-ont", "minimap2's map-ont is -x lr here"},
    {"map-hifi", "minimap2's map-hifi is -x lr:hq here"},
    {"map-ccs", "minimap2's map-ccs is -x lr:hq here"},
    {"cdna", "minimap2's cdna is -x splice here"},
    {"ava-pb", "minimap2's ava-pb overlaps PacBio CLR reads, which have no "
               "preset here; -x ava-hifi overlaps HiFi reads"},
};

constexpr Minimap2Hint kIndexHints[] = {
    {"-w", "seeds are closed syncmers: -k and -s set their density"},
    {"-d", "the output is the second operand: flashalign index ref.fa out.faix"},
};

// "; <hint>" when minimap2 has the spelling, else "". A long spelling may
// carry its value (--splice-flank=no).
std::string minimap2_hint(std::string_view token, unsigned mode) {
    const std::string_view spelling = token.substr(0, token.find('='));
    if (mode == ModeAlign) {
        for (const Minimap2Hint& h : kAlignHints)
            if (h.spelling == spelling) return "; " + std::string(h.hint);
    } else {
        for (const Minimap2Hint& h : kIndexHints)
            if (h.spelling == spelling) return "; " + std::string(h.hint);
    }
    return {};
}

// The refusal of an unknown -x, with its hint when minimap2 has the name.
std::string preset_refusal(const std::string& preset) {
    std::string message =
        "-x preset must be one of: " + fa::cpu::api::accepted_preset_names();
    for (const Minimap2Hint& h : kPresetHints)
        if (h.spelling == preset) message += "; " + std::string(h.hint);
    return message;
}

}  // namespace

}  // namespace fa::cpu::cli

/**** tokens ****/
// argv as tokens: getopt-style short-option bundles expanded, an option's value
// consumed, SPELLING=VALUE.

namespace fa::cpu::cli {

namespace {

// argv with short-option bundles expanded. `attached` marks a value that was
// glued to its letter ("8" from "-t8").
struct TokenStream {
    std::vector<std::string> tokens;
    std::vector<char> attached;

    void push(std::string token, bool was_attached) {
        tokens.push_back(std::move(token));
        attached.push_back(was_attached ? 1 : 0);
    }
    std::size_t size() const { return tokens.size(); }
};

// A negative number such as "-12" is never expanded as a bundle, so a stray
// one is reported whole.
bool looks_like_negative_number(const std::string& token) {
    if (token.size() < 2 || token[0] != '-') return false;
    for (std::size_t i = 1; i < token.size(); ++i) {
        const unsigned char ch = static_cast<unsigned char>(token[i]);
        if (!std::isdigit(ch) && ch != '.') return false;
    }
    return true;
}

// getopt-style short-option bundling (`-ax lr:hq`, `-ct8`). A bundle is a
// single-dash token longer than two characters that is not a registered
// spelling or a negative number. The first value-taking letter ends it; the
// rest of the token, if any, is its value.
TokenStream expand_short_bundles(
    int argc, char** argv, int start, unsigned mode, const std::string& verb) {
    TokenStream out;
    bool end_of_options = false;
    for (int i = start; i < argc; ++i) {
        const std::string arg = argv[i];
        const bool bundle = !end_of_options && arg.size() > 2 &&
                            arg[0] == '-' && arg[1] != '-' &&
                            !looks_like_negative_number(arg) &&
                            find_option_spec(arg, mode) == nullptr;
        if (arg == "--") end_of_options = true;
        if (!bundle) {
            out.push(arg, false);
            continue;
        }
        for (std::size_t j = 1; j < arg.size(); ++j) {
            const std::string letter = std::string("-") + arg[j];
            const OptionSpec* spec = find_option_spec(letter, mode);
            if (spec == nullptr) {
                throw UsageError(
                    "unknown " + verb + " option: " + letter + " (in " + arg +
                    ")" + minimap2_hint(letter, mode));
            }
            out.push(letter, false);
            if (spec->kind == ValueKind::None) continue;
            const std::string rest = arg.substr(j + 1);
            if (!rest.empty()) out.push(rest, true);
            break;
        }
    }
    return out;
}

// Consumes the value after `opt`. A registered spelling or "--" in that
// position means the value is missing (`-t --quiet`); anything else, such as
// a negative number or "-", is a value. An attached value is always accepted.
std::string require_value(
    std::size_t& i, const TokenStream& ts, const std::string& opt,
    unsigned mode) {
    if (i + 1 >= ts.size()) throw UsageError("missing value for " + opt);
    const std::string& next = ts.tokens[i + 1];
    if (!ts.attached[i + 1] &&
        (next == "--" || find_option_spec(next, mode) != nullptr)) {
        throw UsageError(opt + " expects a value");
    }
    return ts.tokens[++i];
}

// The align spellings that take a value attached with '=', as in minimap2:
// --secondary=yes|no and --dual=yes|no, which also take it separate, and
// --cs=short|long, which takes it only attached (a bare --cs is short).
constexpr std::string_view kAttachedValueSpellings[] = {"--secondary", "--cs",
                                                        "--dual"};

}  // namespace

}  // namespace fa::cpu::cli

/**** preset_checks ****/
// The checks that read the preset; cli/align.cpp repeats them once an index has
// settled a preset not given with -x.

namespace fa::cpu::cli {

void check_max_cands_range(const AlignOptions& opt) {
    if (!opt.max_cands) return;
    const int ceiling = fa::cpu::api::is_overlap_preset(opt.preset)
                            ? ::fa::cpu::voting::kAllChainsLaneBound
                            : ::fa::cpu::voting::kMaxCatalogueLaneBound;
    if (*opt.max_cands < 1 || *opt.max_cands > ceiling)
        throw UsageError("--max-cands must be within [1," +
                         std::to_string(ceiling) + "]");
}

namespace {

// A stage's refusal, checked in this order. {preset} is the preset, {options}
// every spelling of the stage and {is} "is" or "are" to agree with them.
struct StageRefusal {
    LaneStage stage;
    std::string_view text;
};

constexpr StageRefusal kStageRefusals[] = {
    {LaneStage::BaseLevelOutput,
     "-x {preset} prints map-only PAF; it is refused with {options}"},
    {LaneStage::BaseAlignment,
     "-x {preset} prints map-only PAF; {options} set base-level alignment and "
     "are refused"},
    {LaneStage::SamOutput,
     "-x {preset} prints map-only PAF; {options} shape SAM output and are "
     "refused"},
    {LaneStage::Selection,
     "-x {preset} prints every chain; {options} are refused"},
    {LaneStage::CandidateChains,
     "-x {preset} chains each contig and strand once; {options} {is} "
     "refused"},
    {LaneStage::OverlapPairs,
     "{options} {is} valid only with ava-ont or ava-hifi"},
};

void replace_all(std::string& text, std::string_view key,
                 const std::string& value) {
    for (std::size_t at = text.find(key); at != std::string::npos;
         at = text.find(key, at + value.size()))
        text.replace(at, key.size(), value);
}

}  // namespace

std::optional<std::string> option_stage_refusal(
    const std::set<OptionId>& given, const std::string& preset) {
    for (const StageRefusal& refusal : kStageRefusals) {
        if (fa::cpu::api::lane_runs_stage(refusal.stage, preset)) continue;
        // One spelling per option, the short one where there is one.
        std::set<OptionId> listed;
        std::vector<std::string> spellings;
        bool was_given = false;
        for (const OptionSpec& spec : option_specs()) {
            if (spec.stage != refusal.stage || !(spec.modes & ModeAlign) ||
                !listed.insert(spec.id).second)
                continue;
            was_given = was_given || given.count(spec.id) != 0;
            spellings.push_back(spec.short_name != '\0'
                                    ? std::string{'-', spec.short_name}
                                    : std::string(spec.long_name));
        }
        if (!was_given) continue;
        std::string options;
        for (std::size_t i = 0; i < spellings.size(); ++i) {
            if (i > 0) options += i + 1 == spellings.size() ? " and " : ", ";
            options += spellings[i];
        }
        std::string text(refusal.text);
        replace_all(text, "{preset}", preset);
        replace_all(text, "{options}", options);
        replace_all(text, "{is}", spellings.size() == 1 ? "is" : "are");
        return text;
    }
    return std::nullopt;
}

}  // namespace fa::cpu::cli

/**** parse ****/

namespace fa::cpu::cli {

namespace {

// One parsed option's value into `opt`: the typed parse, its range, and the
// field it sets. --secondary and --dual are taken before the dispatch;
// --help prints and exits.
void dispatch_align_option(AlignOptions& opt, bool& secondary_off_by_n,
                           const std::string& arg,
                           const std::optional<std::string>& attached,
                           const OptionSpec* spec, const std::string& val) {
    switch (spec->id) {
        case OptionId::Help:
            print_align_help(std::cout);
            std::exit(0);
        case OptionId::Output:     opt.output_path = val; break;
        case OptionId::PafCigar:   opt.paf_cigar = true; break;
        case OptionId::OutputSam:  opt.format = "sam"; break;
        case OptionId::Cs: {
          const std::string form = attached.value_or("short");
          if (form != "short" && form != "long") {
            throw UsageError(
                "--cs expects 'short' or 'long' (got '" + form + "')");
          }
          opt.cs = form;
          break;
        }
        case OptionId::Md:
          opt.emit_md = true;
          break;
        case OptionId::Eqx:
          opt.emit_eqx = true;
          break;
        case OptionId::NoHeader:   opt.no_header = true; break;
        case OptionId::SamHitOnly: opt.sam_hit_only = true; break;
        case OptionId::PafNoHit:   opt.paf_no_hit = true; break;
        case OptionId::SoftClipSupp: opt.soft_clip_supp = true; break;
        case OptionId::CopyComment: opt.copy_comment = true; break;
        case OptionId::ReadGroup: {
            const ReadGroup group = parse_read_group(val);
            opt.read_group_line = group.line;
            opt.read_group_id = group.id;
            break;
        }
        // An explicit -x outranks the preset an index records; run_align
        // only replaces a "builtin" preset.
        case OptionId::Preset:
            opt.preset = val;
            opt.preset_source = "explicit";
            break;
        case OptionId::Stats:      opt.stats = true; break;
        case OptionId::ShowConfig: opt.show_config = true; break;
        case OptionId::Threads:    opt.threads = parse_int(val, arg); break;
        case OptionId::BatchBp:    opt.batch_bp = parse_size(val, arg); break;
        case OptionId::BatchWindow: opt.batch_window = parse_int(val, arg); break;
        case OptionId::IoStaging:  opt.io_staging_mib = parse_int(val, arg); break;
        case OptionId::Progress:   opt.progress = true; break;
        case OptionId::Quiet:      opt.quiet = true; break;
        case OptionId::VoteRatio:
            opt.dna_vote_admission_ratio = parse_ratio(val, arg);
            break;
        case OptionId::MaxChainOcc:
            opt.max_chain_occ = parse_int(val, arg);
            break;
        case OptionId::MinSupport:    opt.min_support = parse_int(val, arg); break;
        case OptionId::VoteSeeds:
            opt.vote_seeds = parse_int(val, arg);
            if (*opt.vote_seeds < 0)
                throw UsageError("--vote-seeds must be >= 0");
            break;
        case OptionId::VoteSpacing:
            opt.vote_spacing = parse_int(val, arg);
            if (*opt.vote_spacing < 0)
                throw UsageError("--vote-spacing must be >= 0");
            break;
        // Its range depends on the lane (check_max_cands_range).
        case OptionId::MaxCands:
            opt.max_cands = parse_int(val, arg);
            break;
        case OptionId::MinChainScore:
            opt.min_chain_score = parse_int(val, arg);
            if (*opt.min_chain_score < 1)
                throw UsageError("-m must be at least 1");
            break;
        case OptionId::DpMatch:    opt.dp_match = parse_int(val, arg); break;
        case OptionId::DpMismatch: opt.dp_mismatch = parse_int(val, arg); break;
        case OptionId::DpGapOpen: {
            const auto p = parse_int_pair(val, arg);
            opt.dp_gap_open1 = p.first;
            opt.dp_gap_open2 = p.second;
            break;
        }
        case OptionId::DpGapExtend: {
            const auto p = parse_int_pair(val, arg);
            opt.dp_gap_extend1 = p.first;
            opt.dp_gap_extend2 = p.second;
            break;
        }
        case OptionId::DpZdrop: {
            // The second value is the inversion Z-drop; a lone value sets
            // it too, as in minimap2 (options/options.cpp).
            const auto comma = val.find(',');
            opt.dp_zdrop = parse_int(
                comma == std::string::npos ? val : val.substr(0, comma), arg);
            if (comma != std::string::npos)
                opt.dp_zdrop_inv = parse_int(val.substr(comma + 1), arg);
            else
                opt.dp_zdrop_inv.reset();
            break;
        }
        case OptionId::DpScoreN:   opt.dp_score_n = parse_int(val, arg); break;
        case OptionId::DpEndBonus: opt.dp_end_bonus = parse_int(val, arg); break;
        case OptionId::DpMinScore: opt.dp_min_score = parse_int(val, arg); break;
        case OptionId::DpBw: {
            // Sizes ("20k"), and a lone value leaves the long-join
            // bandwidth at its default, as in minimap2.
            const auto comma = val.find(',');
            opt.dp_bw = parse_int_size(
                comma == std::string::npos ? val : val.substr(0, comma), arg);
            if (comma != std::string::npos) {
                opt.dp_bw_long = parse_int_size(val.substr(comma + 1), arg);
            }
            break;
        }
        // -g is a size ("10k"), as in minimap2.
        case OptionId::DpMaxGap: {
            const int gap = parse_int_size(val, arg);
            if (gap <= 0) throw UsageError("-g must be > 0");
            opt.dp_max_gap = gap;
            break;
        }
        // A size, as -r.
        case OptionId::ChainBand: {
            const int band = parse_int_size(val, arg);
            if (band < 1)
                throw UsageError(arg + " must be at least 1");
            opt.chain_band = band;
            break;
        }
        case OptionId::Secondary:
        case OptionId::Dual:
            // Taken by parse_align_args before the dispatch.
            break;
        case OptionId::MaxVoteOcc:    opt.max_vote_occ = parse_int(val, arg); break;
        case OptionId::VoteDiagBinWidth: {
            const int width = parse_int(val, arg);
            if (width <= 0) {
                throw UsageError(
                    arg + " must be > 0 (got " +
                    std::to_string(width) + ")");
            }
            opt.vote_diag_bin_width = width;
            break;
        }
        case OptionId::VoteDiagSlopeDen: {
            const int den = parse_int(val, arg);
            if (den < 0) {
                throw UsageError(
                    arg + " must be >= 0 (got " +
                    std::to_string(den) + ")");
            }
            opt.vote_diag_slope_den = den;
            break;
        }
        case OptionId::VoteDiagWidthMax: {
            const int width_max = parse_int(val, arg);
            if (width_max <= 0) {
                throw UsageError(
                    arg + " must be > 0 (got " +
                    std::to_string(width_max) + ")");
            }
            opt.vote_diag_width_max = width_max;
            break;
        }
        case OptionId::MinIntron:  opt.min_intron = parse_int(val, arg); break;
        // -G is a size ("200k"), as in minimap2. It also sets both -r
        // values (options/options.cpp), so a -r typed before it is
        // dropped and only a later one applies.
        case OptionId::MaxIntron:
            opt.max_intron = parse_int_size(val, arg);
            opt.dp_bw.reset();
            opt.dp_bw_long.reset();
            break;
        case OptionId::RnaJunctionBed:
            opt.rna_junction_bed = val;
            break;
        case OptionId::RnaJunctionBonus:
            opt.rna_junction_bonus = parse_int(val, arg);
            break;
        case OptionId::PriRatio:
            opt.pri_ratio = parse_ratio(val, arg);
            break;
        case OptionId::SpliceRivalMinDiff:
            opt.rna_rival_min_diff = parse_int(val, arg);
            break;
        // minimap2's -N counts secondaries, one fewer than the loci kept.
        // -N 0 keeps the loci and turns secondary output off, as minimap2
        // rewrites it to --secondary=no.
        case OptionId::SpliceMaxLoci: {
            const int n = parse_int(val, arg);
            if (n < 0)
                throw UsageError("-N must be >= 0");
            if (n == std::numeric_limits<int>::max())
                throw UsageError(
                    "integer out of range for " + arg + ": " + val);
            secondary_off_by_n = n == 0;
            if (n > 0) opt.rna_max_loci = n + 1;
            break;
        }
        // minimap2's letters, stored as the words rna::parse_strand_mode
        // reads.
        case OptionId::SpliceStrand:
            if (val == "f")      opt.splice_strand = "forward";
            else if (val == "b") opt.splice_strand = "auto";
            else if (val == "r") opt.splice_strand = "reverse";
            else if (val == "n") opt.splice_strand = "none";
            else throw UsageError(
                "-u must be f (forward), b (both/auto), r (reverse) "
                "or n (no splice motif)");
            break;
    }
}

// After the full parse, so argument order does not matter: the operands, the
// ranges, the options valid with one kind of preset only, and the rules that
// combine options (-N 0, --cs and --MD imply -c, -N on a DNA preset).
void settle_align_options(AlignOptions& opt, bool secondary_off_by_n,
                          const std::vector<std::string>& positional) {
    // Without -x the preset may come from the index: cli/align.cpp
    // checks the range once it is settled.
    if (opt.preset_source != "builtin")
        check_max_cands_range(opt);
    // --show-config reads no sequence, so its operands are optional.
    if (positional.size() < 2 && !opt.show_config) {
        throw UsageError(
            "align expects <ref.fa|ref.faix> and one or more <reads.fq[.gz]> "
            "('-' for stdin)");
    }
    if (opt.batch_bp <= 0) {
        throw UsageError("--batch-bp must be positive");
    }
    if (opt.batch_window < 1) {
        throw UsageError("--batch-window must be at least 1");
    }
    if (opt.threads && *opt.threads < 1) {
        throw UsageError("-t must be at least 1");
    }
    // 0 is the auto scale; a smaller ceiling would make the decoder thrash.
    if (opt.io_staging_mib != 0 && opt.io_staging_mib < 64) {
        throw UsageError(
            "--io-staging must be 0 (auto) or at least 64 (MiB)");
    }
    // <ref.fa|ref.faix> <reads>...; run_align tells FASTA from index by
    // content.
    if (!positional.empty()) {
        opt.target_path = positional[0];
        opt.reads_paths.assign(positional.begin() + 1, positional.end());
    }
    if (!fa::cpu::api::preset_is_valid(opt.preset)) {
        throw UsageError(preset_refusal(opt.preset));
    }
    // Options valid for only one mode, checked after the full parse so
    // argument order does not matter. options/options.cpp repeats most of
    // them for the API.
    const bool rna = opt.preset == "splice" || opt.preset == "splice:hq";
    // RNA presets have no adaptive vote-width model; --dw still sets their
    // fixed width.
    if (rna && (opt.vote_diag_slope_den || opt.vote_diag_width_max)) {
        throw UsageError(
            "--dw-slope-den and --dw-max are valid only with a DNA preset");
    }
    if (secondary_off_by_n)
        opt.output_secondary = false;
    // --cs and --MD imply a CIGAR, as minimap2's --cs does.
    if (!opt.cs.empty() || opt.emit_md)
      opt.paf_cigar = true;
    // The RNA-preset check for --max-chain-occ is in options/options.cpp.
    if (opt.max_chain_occ && *opt.max_chain_occ < 0)
        throw UsageError("--max-chain-occ must be >= 0");
    if (rna && opt.dna_vote_admission_ratio) {
        throw UsageError("--vote-ratio is valid only with a DNA preset");
    }
    if (rna && opt.min_chain_score)
        throw UsageError("-m is valid only with a DNA preset");
    if (rna && opt.max_cands)
        throw UsageError("--max-cands is valid only with a DNA preset");
    if (rna && opt.chain_band)
        throw UsageError("-b is valid only with a DNA preset");
    if (rna && opt.vote_seeds)
        throw UsageError("--vote-seeds is valid only with a DNA preset; the "
                         "splice presets take --vote-spacing");
    // Without -x an index may record a splice preset; options/options.cpp
    // refuses then.
    if (!rna && opt.vote_spacing && opt.preset_source != "builtin")
        throw UsageError("--vote-spacing is valid only with an RNA preset");
    // Without -x the preset may come from the index: cli/align.cpp
    // asks again once it is settled.
    if (opt.preset_source != "builtin") {
        if (const auto refusal =
                option_stage_refusal(opt.given_options, opt.preset))
            throw UsageError(*refusal);
    }
    // The 128-tile and 128-seed trace refusals are in cli/align.cpp.
    if (!rna) {
        if (opt.min_intron) {
            throw UsageError(
                "--min-intron is valid only with splice or splice:hq");
        }
        if (opt.max_intron) {
            throw UsageError(
                "-G is valid only with splice or splice:hq");
        }
        if (opt.splice_strand) {
            throw UsageError(
                "-u is valid only with splice or splice:hq");
        }
        if (opt.rna_junction_bed)
            throw UsageError(
                "--junc-bed is valid only with splice or splice:hq");
        if (opt.rna_junction_bonus)
            throw UsageError(
                "--junc-bonus is valid only with splice or splice:hq");
        if (opt.rna_rival_min_diff)
          throw UsageError(
              "--rival-min-diff is valid only with splice or splice:hq");
        // -N n is minimap2's best_n, whatever --secondary says.
        if (opt.rna_max_loci)
          opt.dna_secondary_max = *opt.rna_max_loci - 1;
        opt.rna_max_loci.reset();
    }
    // Intron bounds must satisfy 0 < min <= max.
    if (opt.min_intron && *opt.min_intron <= 0) {
        throw UsageError("--min-intron must be > 0");
    }
    if (opt.max_intron && *opt.max_intron <= 0) {
        throw UsageError("-G must be > 0");
    }
    if (opt.rna_junction_bonus &&
        (*opt.rna_junction_bonus < 0 || *opt.rna_junction_bonus > 127))
        throw UsageError("--junc-bonus must be within [0,127]");
    if (opt.rna_rival_min_diff && *opt.rna_rival_min_diff < 0)
        throw UsageError("--rival-min-diff must be >= 0");
    if (opt.rna_max_loci && *opt.rna_max_loci < 1)
        throw UsageError("-N must be >= 1");
    if (opt.min_intron && *opt.min_intron > 0 && opt.max_intron && *opt.max_intron > 0 &&
        *opt.max_intron < *opt.min_intron) {
        throw UsageError("-G must be >= --min-intron");
    }
}

}  // namespace

AlignOptions parse_align_args(int argc, char** argv, int start) {
    AlignOptions opt;
    // @PG CL: as minimap2 writes it, but with the program name instead of
    // argv[0] so the header does not carry the install path.
    opt.command_line = "flashalign";
    for (int i = 1; i < argc; ++i) {
        opt.command_line.push_back(' ');
        opt.command_line += argv[i];
    }
    bool secondary_off_by_n = false;
    std::vector<std::string> positional;
    // Everything after "--" is an operand.
    bool end_of_options = false;
    const TokenStream ts =
        expand_short_bundles(argc, argv, start, ModeAlign, "align");
    for (std::size_t i = 0; i < ts.size(); ++i) {
        const std::string arg = ts.tokens[i];
        if (end_of_options) {
            positional.push_back(arg);
            continue;
        }
        if (arg == "--") {
            end_of_options = true;
            continue;
        }
        // SPELLING=VALUE, for the spellings of kAttachedValueSpellings.
        std::string_view spelling = arg;
        std::optional<std::string> attached;
        for (const std::string_view name : kAttachedValueSpellings) {
            if (arg.size() > name.size() && arg[name.size()] == '=' &&
                arg.compare(0, name.size(), name) == 0) {
                spelling = name;
                attached = arg.substr(name.size() + 1);
                break;
            }
        }
        const OptionSpec* spec = find_option_spec(spelling, ModeAlign);
        if (spec == nullptr) {
            if (arg.size() > 1 && arg[0] == '-') {
                throw UsageError("unknown align option: " + arg +
                                 minimap2_hint(arg, ModeAlign));
            }
            // A bare "-" is the stdin marker, not an option -> positional.
            positional.push_back(arg);
            continue;
        }
        opt.given_options.insert(spec->id);
        // --secondary accepts its value attached or separate.
        if (spec->id == OptionId::Secondary) {
            const std::string answer =
                attached ? *attached : require_value(i, ts, arg, ModeAlign);
            if (answer != "yes" && answer != "no")
                throw UsageError("--secondary must be yes or no");
            opt.output_secondary = answer == "yes";
            continue;
        }
        if (spec->id == OptionId::Dual) {
            const std::string answer =
                attached ? *attached : require_value(i, ts, arg, ModeAlign);
            if (answer != "yes" && answer != "no")
                throw UsageError("--dual must be yes or no");
            opt.dual = answer == "yes";
            continue;
        }
        std::string val;
        if (spec->kind != ValueKind::None) {
            val = require_value(i, ts, arg, ModeAlign);
        }
        dispatch_align_option(opt, secondary_off_by_n, arg, attached, spec, val);
    }
    settle_align_options(opt, secondary_off_by_n, positional);
    return opt;
}

IndexOptions parse_index_args(int argc, char** argv, int start) {
    IndexOptions opt;
    std::vector<std::string> positional;
    bool end_of_options = false;
    const TokenStream ts =
        expand_short_bundles(argc, argv, start, ModeIndex, "index");
    for (std::size_t i = 0; i < ts.size(); ++i) {
        const std::string arg = ts.tokens[i];
        if (end_of_options) {
            positional.push_back(arg);
            continue;
        }
        if (arg == "--") {
            end_of_options = true;
            continue;
        }
        const OptionSpec* spec = find_option_spec(arg, ModeIndex);
        if (spec == nullptr) {
            // A bare "-" is stdin.
            if (arg.size() > 1 && arg[0] == '-') {
                throw UsageError("unknown index option: " + arg +
                                 minimap2_hint(arg, ModeIndex));
            }
            positional.push_back(arg);
            continue;
        }
        std::string val;
        if (spec->kind != ValueKind::None)
            val = require_value(i, ts, arg, ModeIndex);
        switch (spec->id) {
            case OptionId::Help:
                print_index_help(std::cout);
                std::exit(0);
            case OptionId::K:        opt.k = parse_int(val, arg); break;
            case OptionId::Preset:   opt.preset = val; break;
            // Unlike minimap2, where -s is the minimum DP score.
            case OptionId::SyncmerS: opt.syncmer_s = parse_int(val, arg); break;
            case OptionId::Threads:  opt.threads = parse_int(val, arg); break;
            case OptionId::IdxNoSeq: opt.no_seq = true; break;
            case OptionId::IndexBatchBp:
                opt.batch_bp = parse_size(val, arg);
                break;
            case OptionId::Quiet:    opt.quiet = true; break;
            default:
                // Unreachable: the index mode admits only the ids above.
                throw UsageError("unknown index option: " + arg);
        }
    }
    // The output defaults to the input path plus ".faix".
    if (positional.empty() || positional.size() > 2) {
        throw UsageError("index expects <ref.fa> [out.faix]");
    }
    opt.ref_path = positional[0];
    opt.out_path =
        positional.size() == 2 ? positional[1] : positional[0] + ".faix";
    if (opt.threads && *opt.threads < 1) {
        throw UsageError("-t must be at least 1");
    }
    if (opt.batch_bp && *opt.batch_bp < 1) {
        throw UsageError("-I must be at least 1 base");
    }
    // A multi-part build reads the FASTA more than once.
    if (opt.batch_bp && opt.ref_path == "-") {
        throw UsageError(
            "-I needs the reference as a file: a multi-part index reads the "
            "FASTA once per part, which stdin cannot repeat");
    }
    if (!fa::cpu::api::preset_is_valid(opt.preset)) {
        throw UsageError(preset_refusal(opt.preset));
    }
    return opt;
}

}  // namespace fa::cpu::cli
