#include "cli/option_registry.h"

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
constexpr LaneStage kPartition = LaneStage::Partition;

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
        // minimap2's -g. On a DNA preset it also bounds the dense chain's
        // reach, which never falls below the chain band; on a splice preset
        // it bounds base alignment only.
        {OptionId::DpMaxGap, 'g', "", kSize, kStable, A,
         "NUM", SEC_PLACEMENT,
         "stop alignment elongation if there are no seeds in NUM-bp [5000]"},
        // minimap2's -r. Base-level alignment only: NUM1 bands the read-end
        // extensions, NUM2 the gap fills.
        {OptionId::DpBw, 'r', "", kPair, kStable, A,
         "NUM[,NUM]", SEC_PLACEMENT,
         "alignment bandwidth and gap-fill bandwidth [500,20000]", "", kBaseAlign},
        // minimap2's -m. DNA presets only; the splice presets keep their own
        // chain floors.
        {OptionId::MinChainScore, 'm', "", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "minimal chaining score (matching bases minus log gap penalty) [40]"},
        {OptionId::PriRatio, 'p', "", kStr, kStable, A,
         "FLOAT", SEC_PLACEMENT,
         "min secondary-to-primary score ratio [0.8]", "", kSelection},
        // As in minimap2: -N 5 keeps the best locus and at most five rivals,
        // and -N 0 is --secondary no. Both realize up to INT rivals; the
        // DNA presets default to 1.
        {OptionId::SpliceMaxLoci, 'N', "", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "retain at most INT secondary alignments [1]", "", kSelection},
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
        // The Python Config.max_query_seeds. Above 128 a DNA preset keeps the
        // 128-seed selection and fills it from the larger one
        // (seeding/syncmer.h).
        {OptionId::VoteSeeds, '\0', "--vote-seeds", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "seeds per strand in voting; 0 for all [128]"},
        // The lane bound, DNA presets only, 1..64: the vote's peaks and the
        // catalogue's candidates per strand. Unset, 16 under ratio admission
        // and 4 under count admission (--vote-ratio 0).
        {OptionId::MaxCands, '\0', "--max-cands", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "max placement candidates per strand [16]"},
        // DNA presets only, 2..4096; the splice presets keep 128. The
        // --tile-score terms are per tile.
        {OptionId::Tiles, '\0', "--tiles", kInt, kStable, A,
         "INT", SEC_PLACEMENT,
         "query tiles in the read's placement partition [128]", "", kPartition},
        // DNA presets only. span: an accepted chain owns the tiles between
        // its first and last anchor, except those a same-locus rival keeps;
        // anchors: only its anchor tiles (dna/placement_chaining.h).
        {OptionId::TileOwner, '\0', "--tile-owner", kStr, kStable, A,
         "STR", SEC_PLACEMENT,
         "tiles a placed chain owns: span or anchors [span]", "", kPartition},

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
        // minimap2's -s. RNA admits a second family only when its DP maximum
        // reaches it; a primary below it keeps its placement. DNA uses it as
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
        // Unlike minimap2's -G, this does not also change -r.
        {OptionId::MaxIntron, 'G', "", kSize, kStable, A,
         "NUM", SEC_SPLICE,
         "max intron length (effective with -xsplice) [200k]"},
        {OptionId::SpliceStrand, 'u', "", kStr, kStable, A,
         "CHAR", SEC_SPLICE,
         "how to find GT-AG. f:transcript strand, b:both strands, r:reverse strand, n:don't match GT-AG [b]"},
        {OptionId::RnaJunctionBed, '\0', "--junc-bed", kStr,
         kStable, A, "FILE", SEC_SPLICE,
         "junctions to prefer during base alignment []"},
        {OptionId::RnaJunctionBonus, '\0', "--junc-bonus", kInt,
         kStable, A, "INT", SEC_SPLICE,
         "score bonus for a splice donor or acceptor found in annotation [9]"},

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
        // way the alternatives appear on the primary's XA:Z, md:i and s2:i,
        // output-only DNA secondaries (-N >= 2) aside.
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

        // --tile-score hit=INT,block=INT,null=INT,miss=INT: the DNA tile
        // objective (--tiles), preset-owned; partial keys merge. hit: reward per
        // supported tile in a block (>0); block: cost per block opened; null:
        // cost per tile in no block; miss: penalty per unsupported tile in a
        // block.
        {OptionId::TileScore, '\0', "--tile-score", kStr, kDev, A, "",
         HIDDEN, "", STUDY_FLASH_NATIVE, kPartition},
        // --dual yes|no, also --dual=yes|no: the overlap presets only. no
        // prints a pair once, from the read whose name sorts first, as
        // minimap2's --dual=no.
        {OptionId::Dual, '\0', "--dual", kStr, kDev, A, "", HIDDEN, "",
         STUDY_CLI_MANUAL, LaneStage::OverlapPairs},

        // index, shown by --help
        {OptionId::Preset, 'x', "--preset", kStr, kStable, I,
         "STR", SEC_IOPT,
         // (k, s) per preset: resolve_preset_seeding() in options/presets.h.
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
