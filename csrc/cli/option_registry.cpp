#include "cli/option_registry.h"

namespace fa::cpu::cli {

namespace {

// Row fields, as in OptionSpec:
//   { id, short, long, kind, tier, modes, metavar, section, help, study }
constexpr ValueKind kNone = ValueKind::None;
constexpr ValueKind kInt = ValueKind::Int;
constexpr ValueKind kPair = ValueKind::IntPair;
constexpr ValueKind kSize = ValueKind::Size;
constexpr ValueKind kStr = ValueKind::Str;

constexpr HelpTier kStable = HelpTier::Stable;
constexpr HelpTier kDev = HelpTier::Dev;

constexpr unsigned A = ModeAlign;
constexpr unsigned I = ModeIndex;

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
constexpr std::string_view STUDY_FLASH_NATIVE =
    "FlashAlign-specific tuning, not on the help screen";
constexpr std::string_view STUDY_CLI_MANUAL =
    "documented in flashalign.1, not on the help screen";

const std::vector<OptionSpec>& specs_table() {
  // clang-format off: preserve the hand-aligned option table.
    static const std::vector<OptionSpec> kSpecs = {
        // align, shown by --help
        {OptionId::DpBw, 'r', "", kPair, kStable, A,
         "INT[,INT]", SEC_PLACEMENT,
         "chaining/alignment bandwidth and long-join bandwidth [500,20000]"},
        {OptionId::SplicePriRatio, 'p', "", kStr, kStable, A,
         "FLOAT", SEC_PLACEMENT,
         "min secondary-to-primary score ratio [0.8]"},
        // The count includes the best locus: -N 6 keeps it plus at most five
        // rivals, where minimap2's -N 6 means six secondaries.
        {OptionId::SpliceMaxLoci, 'N', "", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "retain at most INT candidate loci per read [6]"},
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

        // On a DNA preset -A -B -O -E -z --score-N set the DP row of the gap
        // fills between anchors; the read ends keep the preset's own row,
        // which also prices every path. Defaults are lr's.
        {OptionId::DpMatch, 'A', "", kInt, kStable, A,
         "INT", SEC_ALIGN, "matching score [4]"},
        {OptionId::DpMismatch, 'B', "", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "mismatch penalty (larger value for lower divergence) [8]"},
        {OptionId::DpGapOpen, 'O', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN, "gap open penalty [8,48]"},
        {OptionId::DpGapExtend, 'E', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN,
         "gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [4,1]"},
        {OptionId::DpZdrop, 'z', "", kPair, kStable, A,
         "INT[,INT]", SEC_ALIGN,
         "Z-drop score and inversion Z-drop score [800,200]"},
        // RNA refuses a realization whose DP maximum is below it (minimap2
        // -s). DNA uses it as the per-record emission floor when a CIGAR is
        // realized, as minimap2's mm_filter_regs; map-only ignores it.
        {OptionId::DpMinScore, 'S', "", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "minimal peak DP alignment score [80]"},
        {OptionId::DpScoreN, '\0', "--score-N", kInt, kStable, A,
         "INT", SEC_ALIGN, "score of a mismatch involving ambiguous bases [2]"},
        {OptionId::DpEndBonus, '\0', "--end-bonus", kInt, kStable, A,
         "INT", SEC_ALIGN,
         "score bonus when alignment extends to the end of the query [-1]"},

        // RNA presets only.
        {OptionId::MinIntron, '\0', "--min-intron", kInt, kStable, A,
         "INT", SEC_SPLICE,
         "min intron length [20]"},
        // Unlike minimap2's -G, this does not also change -r.
        {OptionId::MaxIntron, 'G', "", kSize, kStable, A,
         "NUM", SEC_SPLICE,
         "max intron length [200k]"},
        {OptionId::SpliceStrand, 'u', "", kStr, kStable, A,
         "CHAR", SEC_SPLICE,
         "how to find GT-AG. f:transcript strand, b:both strands, r:reverse strand, n:don't match GT-AG [b]"},
        {OptionId::RnaJunctionBed, '\0', "--junc-bed", kStr,
         kStable, A, "FILE", SEC_SPLICE,
         "junctions in BED6/BED12 to guide spliced alignment []"},
        {OptionId::RnaJunctionBonus, '\0', "--junc-bonus", kInt,
         kStable, A, "INT", SEC_SPLICE,
         "score bonus for a known junction [9]"},

        {OptionId::Output,   'o', "--output", kStr, kStable, A,
         "FILE", SEC_IO, "output alignments to FILE [stdout]"},
        {OptionId::Format,   'f', "--format", kStr, kStable, A,
         "STR", SEC_IO, "output format: sam, bam or paf [sam]"},
        {OptionId::PafCigar, 'c', "", kNone, kStable, A,
         "", SEC_IO, "output CIGAR in PAF"},
        // As in minimap2 the value can only be attached (--cs=long); a bare
        // --cs means short.
        {OptionId::Cs, '\0', "--cs", kNone, kStable, A,
         "", SEC_IO,
         "output the cs tag; --cs=long for the long form"},
        {OptionId::Md, '\0', "--MD", kNone, kStable, A,
         "", SEC_IO,
         "output the MD tag"},
        // Applies to every realized CIGAR (SAM, BAM, SA:Z, cg:Z); plain PAF
        // has none.
        {OptionId::Eqx, '\0', "--eqx", kNone, kStable, A,
         "", SEC_IO,
         "write =/X CIGAR operators"},
        {OptionId::SoftClipSupp, 'Y', "--soft-clip-supp", kNone, kStable, A,
         "", SEC_IO, "use soft clipping for supplementary alignments"},
        {OptionId::ReadGroup, 'R', "--rg", kStr, kStable, A,
         "STR", SEC_IO,
         "SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []"},
        {OptionId::SamHitOnly, '\0', "--sam-hit-only", kNone, kStable, A,
         "", SEC_IO,
         "don't output unmapped reads in SAM/BAM"},
        // Also --secondary=yes, as minimap2 spells it. yes emits one record
        // per mapped alternative (FLAG 0x100, tp:A:S, MAPQ 0, no SEQ); either
        // way the alternatives appear on the primary's XA:Z, md:i and s2:i.
        {OptionId::Secondary, '\0', "--secondary", kStr, kStable, A,
         "yes|no", SEC_IO,
         "output secondary alignments [no]"},
        {OptionId::NoHeader, '\0', "--no-header", kNone, kStable, A,
         "", SEC_IO,
         "don't output the SAM header"},
        // PAF and SAM copy the comment verbatim; BAM requires SAM tag text.
        {OptionId::CopyComment, 'y', "", kNone, kStable, A,
         "", SEC_IO,
         "copy FASTA/Q comments to output"},
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
         "preset: lr, lr:hq, splice or splice:hq (see 'flashalign index') [lr]"},

        // align, not on the help screen; the manual page lists these.
        //
        // -a: minimap2's "output SAM". SAM is already the default; the row
        // lets minimap2 command lines parse.
        {OptionId::OutputSam, 'a', "", kNone, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // --batch-window INT [3]: batches in flight in the compute stage.
        // Output is identical for every value.
        {OptionId::BatchWindow, '\0', "--batch-window", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // --io-staging MIB [0]: ceiling on the parallel gzip decoder's
        // read-ahead; 0 scales it with -t.
        {OptionId::IoStaging, '\0', "--io-staging", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        // RNA presets only: --realize-max INT [5], the per-read realization
        // budget, and --rival-min-diff INT [30], the absolute retention band
        // beside -p.
        {OptionId::SpliceRealizeMax, '\0', "--realize-max", kInt, kDev, A,
         "", HIDDEN, "", STUDY_CLI_MANUAL},
        {OptionId::SpliceRivalMinDiff, '\0', "--rival-min-diff", kInt, kDev,
         A, "", HIDDEN, "", STUDY_CLI_MANUAL},

        // --tile-score hit=INT,block=INT,null=INT,miss=INT: the DNA 128-tile
        // objective, preset-owned; partial keys merge. hit: reward per
        // supported tile in a block (>0); block: cost per block opened; null:
        // cost per tile in no block; miss: penalty per unsupported tile in a
        // block.
        {OptionId::TileScore, '\0', "--tile-score", kStr, kDev, A, "",
         HIDDEN, "", STUDY_FLASH_NATIVE},

        // index, shown by --help
        {OptionId::Preset, 'x', "--preset", kStr, kStable, I,
         "STR", SEC_IOPT,
         // (k, s) per preset: resolve_preset_seeding() in options/presets.h.
         "preset [lr]\n"
         "  lr         noisy long reads (Nanopore) vs reference mapping                 [-k21 -s9]\n"
         "  lr:hq      accurate long reads (HiFi, error rate <1%) vs reference mapping  [-k21 -s5]\n"
         "  splice     spliced alignment for long RNA reads                             [-k15 -s10]\n"
         "  splice:hq  spliced alignment for accurate long RNA reads                    [-k15 -s10]"},
        {OptionId::K, 'k', "", kInt, kStable, I,
         "INT", SEC_IOPT, "k-mer size (no larger than 23) [21]"},
        {OptionId::SyncmerS, 's', "", kInt, kStable, I,
         "INT", SEC_IOPT, "closed-syncmer s-mer size [9]"},
        {OptionId::IndexBatchBp, 'I', "", kSize, kStable, I,
         "NUM", SEC_IOPT,
         "split index for every ~NUM reference bases [no split]"},
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

} // namespace

const std::vector<OptionSpec>& option_specs() { return specs_table(); }

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

} // namespace fa::cpu::cli
