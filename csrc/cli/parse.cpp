#include "cli/parse.h"

#include "cli/errors.h"
#include "cli/help.h"
#include "cli/option_registry.h"
#include "cli/read_group.h"
#include "api/aligner.h"  // preset_is_valid

#include <cerrno>
#include <charconv>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace fa::cpu::cli {

namespace {

constexpr int kMaxTileObjectiveTerm =
    ::fa::cpu::voting::kMaxQueryTileObjectiveTerm;

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

void reject_composite_whitespace(
    const std::string& value, const std::string& option) {
    for (const unsigned char ch : value) {
        if (std::isspace(ch)) {
            throw UsageError(
                option + " rejects whitespace: " + value);
        }
    }
}

int parse_composite_int(
    const std::string& value, const std::string& option,
    const std::string& key) {
    int parsed = 0;
    const char* begin = value.data();
    const char* end = begin + value.size();
    const auto result = std::from_chars(begin, end, parsed, 10);
    if (value.empty() || result.ec != std::errc{} || result.ptr != end) {
        throw UsageError(
            "invalid integer for " + option + " key '" + key + "': " + value);
    }
    return parsed;
}

template <class Assign>
void parse_composite_entries(
    const std::string& value, const std::string& option, Assign assign) {
    if (value.empty()) {
        throw UsageError(option + " requires at least one key=value");
    }
    reject_composite_whitespace(value, option);
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const std::size_t comma = value.find(',', begin);
        const std::size_t end =
            comma == std::string::npos ? value.size() : comma;
        if (end == begin) {
            throw UsageError(option + " rejects empty entries");
        }
        const std::string entry = value.substr(begin, end - begin);
        const std::size_t equal = entry.find('=');
        if (equal == std::string::npos || equal == 0 ||
            equal + 1 == entry.size() ||
            entry.find('=', equal + 1) != std::string::npos) {
            throw UsageError(
                "malformed " + option + " entry '" + entry +
                "' (expected key=value)");
        }
        assign(entry.substr(0, equal), entry.substr(equal + 1));
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
}

void parse_tile_score(const std::string& value, AlignOptions& opt) {
    parse_composite_entries(
        value, "--tile-score",
        [&](const std::string& key, const std::string& raw) {
            const int parsed =
                parse_composite_int(raw, "--tile-score", key);
            if (key == "hit") {
                if (parsed <= 0) {
                    throw UsageError(
                        "--tile-score key 'hit' must be within [1," +
                        std::to_string(kMaxTileObjectiveTerm) + "]");
                }
                opt.tile_supported_reward = parsed;
            } else if (key == "block") {
                if (parsed < 0) {
                    throw UsageError(
                        "--tile-score key 'block' must be within [0," +
                        std::to_string(kMaxTileObjectiveTerm) + "]");
                }
                opt.tile_block_open_cost = parsed;
            } else if (key == "null") {
                if (parsed < 0) {
                    throw UsageError(
                        "--tile-score key 'null' must be within [0," +
                        std::to_string(kMaxTileObjectiveTerm) + "]");
                }
                opt.tile_null_cost = parsed;
            } else if (key == "miss") {
                if (parsed < 0) {
                    throw UsageError(
                        "--tile-score key 'miss' must be within [0," +
                        std::to_string(kMaxTileObjectiveTerm) + "]");
                }
                opt.tile_unsupported_cost = parsed;
            } else {
                throw UsageError(
                    "unknown --tile-score key '" + key +
                    "'; accepted keys: hit, block, null, miss");
            }
        });
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

// A size that must fit an int (-r, -g, -G).
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

// The refusal of an unknown -x, naming the preset here when minimap2's name
// has one.
std::string preset_refusal(const std::string& preset) {
    std::string message =
        "-x preset must be one of: " + fa::cpu::api::accepted_preset_names();
    std::string_view here;
    if (preset == "map-ont") here = "lr";
    else if (preset == "map-hifi" || preset == "map-ccs") here = "lr:hq";
    else if (preset == "cdna") here = "splice";
    if (!here.empty())
        message += "; minimap2's " + preset + " is -x " + std::string(here) +
                   " here";
    return message;
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
        // --secondary=VALUE, as minimap2 spells it.
        constexpr std::string_view kSecondaryPrefix = "--secondary=";
        const bool secondary_attached =
            arg.compare(0, kSecondaryPrefix.size(), kSecondaryPrefix) == 0;
        // --cs takes an optional attached value only, as in minimap2.
        constexpr std::string_view kCsPrefix = "--cs=";
        const bool cs_attached =
            arg.compare(0, kCsPrefix.size(), kCsPrefix) == 0;
        const OptionSpec* spec =
            find_option_spec(secondary_attached ? "--secondary"
                             : cs_attached      ? "--cs"
                                                : arg,
                             ModeAlign);
        if (spec == nullptr) {
            if (arg.size() > 1 && arg[0] == '-') {
                throw UsageError("unknown align option: " + arg +
                                 minimap2_hint(arg, ModeAlign));
            }
            // A bare "-" is the stdin marker, not an option -> positional.
            positional.push_back(arg);
            continue;
        }
        // --secondary accepts its value attached or separate.
        if (spec->id == OptionId::Secondary) {
            const std::string answer =
                secondary_attached
                    ? arg.substr(kSecondaryPrefix.size())
                    : require_value(i, ts, arg, ModeAlign);
            if (answer != "yes" && answer != "no")
                throw UsageError("--secondary must be yes or no");
            opt.output_secondary = answer == "yes";
            continue;
        }
        std::string val;
        if (spec->kind != ValueKind::None) {
            val = require_value(i, ts, arg, ModeAlign);
        }
        switch (spec->id) {
            case OptionId::Help:
                print_align_help(std::cout);
                std::exit(0);
            case OptionId::Output:     opt.output_path = val; break;
            case OptionId::PafCigar:   opt.paf_cigar = true; break;
            case OptionId::OutputSam:  opt.format = "sam"; break;
            case OptionId::Cs: {
              const std::string form =
                  cs_attached ? arg.substr(kCsPrefix.size()) : "short";
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
            case OptionId::TileScore: parse_tile_score(val, opt); break;
            case OptionId::MinSupport:    opt.min_support = parse_int(val, arg); break;
            case OptionId::VoteSeeds:
                opt.vote_seeds = parse_int(val, arg);
                if (*opt.vote_seeds < 0)
                    throw UsageError("--vote-seeds must be >= 0");
                break;
            case OptionId::MaxCands:
                opt.max_cands = parse_int(val, arg);
                if (*opt.max_cands < 1 || *opt.max_cands > 64)
                    throw UsageError("--max-cands must be within [1,64]");
                break;
            case OptionId::Tiles:
                opt.tiles = parse_int(val, arg);
                if (*opt.tiles < 2 || *opt.tiles > 4096)
                    throw UsageError("--tiles must be within [2,4096]");
                break;
            case OptionId::TileOwner:
                if (val != "span" && val != "anchors")
                    throw UsageError("--tile-owner must be span or anchors");
                opt.tile_owner = val;
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
                // it too, as in minimap2 (options/resolve.cpp).
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
            case OptionId::Secondary:
                // Handled above.
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
            // -G is a size ("200k"), as in minimap2.
            case OptionId::MaxIntron:
                opt.max_intron = parse_int_size(val, arg);
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
    // argument order does not matter. options/resolve.cpp repeats most of
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
    // The RNA-preset check for --max-chain-occ is in options/resolve.cpp.
    if (opt.max_chain_occ && *opt.max_chain_occ < 0)
        throw UsageError("--max-chain-occ must be >= 0");
    if (rna && opt.dna_vote_admission_ratio) {
        throw UsageError("--vote-ratio is valid only with a DNA preset");
    }
    if (rna && opt.min_chain_score)
        throw UsageError("-m is valid only with a DNA preset");
    if (rna && opt.max_cands)
        throw UsageError("--max-cands is valid only with a DNA preset");
    if (rna && opt.tiles)
        throw UsageError("--tiles is valid only with a DNA preset");
    if (rna && opt.tile_owner)
        throw UsageError("--tile-owner is valid only with a DNA preset");
    if (rna && (opt.tile_supported_reward || opt.tile_block_open_cost ||
                opt.tile_null_cost || opt.tile_unsupported_cost)) {
        throw UsageError(
            "--tile-score is incompatible with RNA preset '" + opt.preset + "'");
    }
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
        // -N n realizes up to n alternatives; the DNA default is 1. Ranks
        // 2..n print only as secondary records, so without secondary output
        // the count is left at the default.
        if (opt.rna_max_loci && opt.output_secondary)
          opt.dna_alternative_realize_max = *opt.rna_max_loci - 1;
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
