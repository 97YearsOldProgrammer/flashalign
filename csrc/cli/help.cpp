#include "cli/help.h"

#include "cli/option_registry.h"

#include <ostream>
#include <string>
#include <string_view>

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
    if (!s.metavar.empty()) {
        left += ' ';
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
