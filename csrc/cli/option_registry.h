#pragma once

// The CLI option registry: one table drives both the parser (parse.cpp) and
// the help screens (help.cpp).
//
// Single letters keep minimap2's meanings, with one exception: -s is the
// closed-syncmer size, so minimap2's minimum DP score is -S here.
// Flash-specific tuning is preset-owned and reported by --show-config.

#include <cstddef>
#include <string_view>
#include <vector>

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
    // align: seed occurrence caps for the vote and the dense chain
    MaxVoteOcc, MaxChainOcc,
    // align: minimal chain score (minimap2 -m)
    MinChainScore,
    // align: vote peak admission
    VoteRatio,
    VoteDiagBinWidth,
    VoteDiagSlopeDen, VoteDiagWidthMax,
    MinSupport,
    // align: the vote's seeds per strand, its candidates per strand, the
    // partition's tiles and its first tile-ownership rule
    VoteSeeds, MaxCands, Tiles, TileOwner,
    // align: DP scoring (minimap2 -A/-B/-O/-E/-z, --score-N, --end-bonus).
    // DpMinScore is minimap2's -s, spelled -S here.
    DpMatch, DpMismatch, DpGapOpen, DpGapExtend, DpZdrop, DpScoreN, DpEndBonus,
    DpMinScore,
    // align: band geometry (minimap2 -r) and maximum gap (minimap2 -g).
    DpBw, DpMaxGap,
    // DNA tile objective, key=value list
    TileScore,
    // RNA intron bounds, strand and junction annotation
    MinIntron, MaxIntron, SpliceStrand, RnaJunctionBed, RnaJunctionBonus,
    // Secondary retention ratio (-p), the RNA retention band
    // (--rival-min-diff) and the catalogue depth and realization budget (-N)
    PriRatio, SpliceRivalMinDiff, SpliceMaxLoci,
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
};

// The registry, ordered align-Stable (grouped by section), align-Dev,
// index-Stable, index-Dev, so help can print each section header once.
const std::vector<OptionSpec>& option_specs();

// The spec a token selects in `mode`, or nullptr if none matches.
const OptionSpec* find_option_spec(std::string_view token, unsigned mode);

} // namespace fa::cpu::cli
