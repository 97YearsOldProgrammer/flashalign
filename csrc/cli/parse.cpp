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
                    ")");
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
    bool format_given = false;
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
                throw UsageError("unknown align option: " + arg);
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
            case OptionId::Format:     opt.format = val; format_given = true; break;
            case OptionId::PafCigar:   opt.paf_cigar = true; break;
            // minimap2 -a: the same as -f sam, so `-a -o out.paf` writes SAM.
            case OptionId::OutputSam:
                opt.format = "sam";
                format_given = true;
                break;
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
            case OptionId::SoftClipSupp: opt.soft_clip_supp = true; break;
            // -y: BAM's tag-text check is in the writer, which knows the read.
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
            case OptionId::DpZdrop:
                // minimap2's pair form is refused: there is no inversion Z-drop.
                if (val.find(',') != std::string::npos) {
                    throw UsageError(
                        "-z takes a single Z-drop score, not a pair");
                }
                opt.dp_zdrop = parse_int(val, arg);
                break;
            case OptionId::DpScoreN:   opt.dp_score_n = parse_int(val, arg); break;
            case OptionId::DpEndBonus: opt.dp_end_bonus = parse_int(val, arg); break;
            case OptionId::DpMinScore: opt.dp_min_score = parse_int(val, arg); break;
            case OptionId::DpBw: {
                // A lone value leaves the long-join bandwidth at its default,
                // as in minimap2.
                const auto comma = val.find(',');
                opt.dp_bw = parse_int(
                    comma == std::string::npos ? val : val.substr(0, comma), arg);
                if (comma != std::string::npos) {
                    opt.dp_bw_long = parse_int(val.substr(comma + 1), arg);
                }
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
            case OptionId::MaxIntron: {
                const int64_t bound = parse_size(val, arg);
                if (bound > std::numeric_limits<int>::max()) {
                    throw UsageError(
                        "integer out of range for " + arg + ": " + val);
                }
                opt.max_intron = static_cast<int>(bound);
                break;
            }
            case OptionId::RnaJunctionBed:
                opt.rna_junction_bed = val;
                break;
            case OptionId::RnaJunctionBonus:
                opt.rna_junction_bonus = parse_int(val, arg);
                break;
            case OptionId::SplicePriRatio:
                opt.rna_pri_ratio = parse_ratio(val, arg);
                break;
            case OptionId::SpliceRivalMinDiff:
                opt.rna_rival_min_diff = parse_int(val, arg);
                break;
            case OptionId::SpliceRealizeMax:
                opt.rna_realize_max = parse_int(val, arg);
                break;
            case OptionId::SpliceMaxLoci:
                opt.rna_max_loci = parse_int(val, arg);
                break;
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
    // Without -f, infer the format from the -o extension.
    if (!format_given && opt.output_path != "-") {
        const std::string& o = opt.output_path;
        const auto ends_with = [&](const char* suffix) {
            const std::string s(suffix);
            return o.size() >= s.size() &&
                   o.compare(o.size() - s.size(), s.size(), s) == 0;
        };
        if (ends_with(".bam")) opt.format = "bam";
        else if (ends_with(".paf")) opt.format = "paf";
        else if (ends_with(".sam")) opt.format = "sam";
    }
    if (opt.format != "sam" && opt.format != "paf" && opt.format != "bam") {
        throw UsageError("--format must be 'sam', 'bam', or 'paf'");
    }
    if (!fa::cpu::api::preset_is_valid(opt.preset)) {
        throw UsageError(
            "-x preset must be one of: " +
            fa::cpu::api::accepted_preset_names());
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
    // --cs and --MD imply a CIGAR, as minimap2's --cs does.
    if (!opt.cs.empty() || opt.emit_md)
      opt.paf_cigar = true;
    // The RNA-preset check for --max-chain-occ is in options/resolve.cpp.
    if (opt.max_chain_occ && *opt.max_chain_occ < 0)
        throw UsageError("--max-chain-occ must be >= 0");
    if (rna && opt.dna_vote_admission_ratio) {
        throw UsageError("--vote-ratio is valid only with a DNA preset");
    }
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
        if (opt.rna_pri_ratio || opt.rna_rival_min_diff ||
            opt.rna_realize_max || opt.rna_max_loci)
          throw UsageError(
              "-p, -N, --realize-max and --rival-min-diff are valid only "
              "with splice or splice:hq");
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
    if (opt.rna_realize_max && *opt.rna_realize_max < 0)
        throw UsageError("--realize-max must be >= 0");
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
                throw UsageError("unknown index option: " + arg);
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
        throw UsageError(
            "-x preset must be one of: " +
            fa::cpu::api::accepted_preset_names());
    }
    return opt;
}

}  // namespace fa::cpu::cli
