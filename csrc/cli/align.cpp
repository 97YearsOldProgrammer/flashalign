// The align command: the SAM and PAF output writer; what --stats and --show-config
// print; the run's plan (threads, the batch window, the decode loan) and run_align, last.
#include "cli/cli.h"

#include "api/api.h"
#include "core/core.h"
#include "fa_version.h"
#include "index/index.h"
#include "io/io.h"
#include "rna/rna.h"
#include "seed/seed.h"
#include "threading/threading.h"

#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

/**** output_writer ****/
// SAM and PAF records in input order: the output stream, the BAM inputs' header, the tags.

namespace fa::cpu::cli {

namespace {

void require_output(std::ostream& output, const char* operation) {
    if (!output) {
        throw std::runtime_error(
            std::string("failed to ") + operation + " alignment output");
    }
}

// Output sink for `-o -`: a 4 MiB buffer like `-o PATH` gets, instead of
// std::cout's small one. A failed write(2) (e.g. EPIPE) fails the stream, so
// require_output() reports it.
class FdSink final : public std::streambuf {
public:
    explicit FdSink(int fd) : fd_(fd), buffer_(1 << 22) {
        setp(buffer_.data(), buffer_.data() + buffer_.size());
    }
    ~FdSink() override { drain(); }

protected:
    int_type overflow(int_type ch) override {
        if (!drain()) return traits_type::eof();
        if (!traits_type::eq_int_type(ch, traits_type::eof())) {
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
        }
        return traits_type::not_eof(ch);
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        std::streamsize done = 0;
        while (done < n) {
            if (pptr() == epptr() && !drain()) break;
            const std::streamsize take = std::min(
                static_cast<std::streamsize>(epptr() - pptr()), n - done);
            std::memcpy(pptr(), s + done, static_cast<std::size_t>(take));
            pbump(static_cast<int>(take));
            done += take;
        }
        return done;
    }

    int sync() override { return drain() ? 0 : -1; }

private:
    // Writes every pending byte, resuming short writes and retrying EINTR. On
    // any other error the pending bytes are dropped and the stream fails.
    bool drain() {
        const char* cursor = pbase();
        std::size_t left = static_cast<std::size_t>(pptr() - pbase());
        bool ok = true;
        while (left > 0) {
            const ssize_t n = ::write(fd_, cursor, left);
            if (n < 0) {
                if (errno == EINTR) continue;
                ok = false;
                break;
            }
            cursor += n;
            left -= static_cast<std::size_t>(n);
        }
        setp(buffer_.data(), buffer_.data() + buffer_.size());
        return ok;
    }

    int fd_;
    std::vector<char> buffer_;
};

std::unique_ptr<std::ostream> open_output_owner(const std::string& path,
                                                std::streambuf* stdout_sink) {
    if (path == "-") return std::make_unique<std::ostream>(stdout_sink);
    static std::vector<char> output_buffer(1 << 22);
    auto output = std::make_unique<std::ofstream>();
    output->rdbuf()->pubsetbuf(
        output_buffer.data(),
        static_cast<std::streamsize>(output_buffer.size()));
    output->open(path, std::ios::binary);
    if (!*output)
        throw std::runtime_error("failed to open output: " + path);
    return output;
}

// A passthrough @RG with the same ID as -R/--rg would define one read group
// twice, so refuse it at startup.
void require_no_read_group_conflict(
    const std::vector<std::string>& passthrough,
    const std::string& read_group_id
) {
    if (read_group_id.empty()) return;
    for (const std::string& line : passthrough) {
        if (line.rfind("@RG\t", 0) != 0) continue;
        if (read_group_line_id(line) == read_group_id) {
            throw std::runtime_error(
                "-R/--rg: the input header already declares read group ID '" +
                read_group_id + "'; choose a different ID");
        }
    }
}

// True when a tag block already declares an RG:Z. Such a read keeps its group;
// -R only labels reads that have none.
bool sam_tags_declare_read_group(const std::string& tags) {
    return tags.rfind("RG:Z:", 0) == 0 ||
           tags.find("\tRG:Z:") != std::string::npos;
}

// MM, ML, MN and the draft spellings Mm, Ml: they count bases of SEQ.
bool is_base_modification_key(const char* key) {
    return key[0] == 'M' && (key[1] == 'M' || key[1] == 'L' || key[1] == 'N' ||
                             key[1] == 'm' || key[1] == 'l');
}

// A B array with one element per base of the read, such as PacBio's
// kinetics; a shorter array (sn:B:f) describes the whole read.
bool is_per_base_array(const char* tag, std::size_t read_length) {
    if (tag[2] != 'B') return false;
    std::uint32_t count;
    std::memcpy(&count, tag + 4, 4);
    return count == read_length;
}

// The read's input tags as SAM text, without those whose keys we write.
output::PassthroughTags passthrough_tags(const io::FastxRecord& read,
                                         const OwnTagKeys* own) {
    output::PassthroughTags tags;
    const std::string& blob = read.tag_bytes;
    for (std::size_t at = 0; io::bam_tag_at(blob, at);) {
        const std::size_t size = io::bam_tag_size(blob, at, read.name);
        const char* tag = blob.data() + at;
        at += size;
        if (own != nullptr && own->drop(tag)) continue;
        if (!tags.whole.empty()) tags.whole.push_back('\t');
        const std::size_t begin = tags.whole.size();
        io::append_bam_tag_sam_text(tags.whole, tag, size);
        if (!is_base_modification_key(tag) &&
            !is_per_base_array(tag, read.seq.size())) {
            if (!tags.clipped.empty()) tags.clipped.push_back('\t');
            tags.clipped.append(tags.whole, begin, std::string::npos);
        }
    }
    return tags;
}

}  // namespace

OwnTagKeys::OwnTagKeys(std::vector<std::string> keys)
    : keys_(std::move(keys)),
      warned_(std::make_unique<std::atomic<bool>[]>(keys_.size())) {}

bool OwnTagKeys::drop(const char* key) const {
    for (std::size_t k = 0; k < keys_.size(); ++k) {
        if (keys_[k][0] != key[0] || keys_[k][1] != key[1]) continue;
        if (!warned_[k].exchange(true)) {
            std::fprintf(stderr,
                         "flashalign: warning: input tag %s dropped: FlashAlign "
                         "writes %s\n",
                         keys_[k].c_str(), keys_[k].c_str());
        }
        return true;
    }
    return false;
}

namespace {

// Reads the header of every input that is a BAM file. Two inputs declaring
// one read group ID with different @RG lines is an error naming both.
InputHeader read_input_headers(const std::vector<std::string>& paths) {
    InputHeader merged;
    // Read group ID -> its line and the input that declared it first.
    std::unordered_map<std::string, std::pair<std::string, std::string>> groups;
    std::unordered_set<std::string> programs;
    for (const std::string& path : paths) {
        if (!io::path_is_bam(path)) continue;
        // The header only: a serial decoder, no thread team.
        io::IoPlan serial;
        serial.parallel_gz = false;
        serial.parallel_bgzf = false;
        const std::string text = io::BamReader(path, &serial).sam_header_text();
        std::string last_program;
        std::size_t begin = 0;
        while (begin < text.size()) {
            std::size_t end = text.find('\n', begin);
            if (end == std::string::npos) end = text.size();
            const std::string line = text.substr(begin, end - begin);
            begin = end + 1;
            if (line.rfind("@RG\t", 0) == 0) {
                const std::string id = read_group_line_id(line);
                const auto [it, added] = groups.emplace(id, std::make_pair(line, path));
                if (added) {
                    merged.passthrough.push_back(line);
                } else if (it->second.first != line) {
                    throw std::runtime_error(
                        "read group ID '" + id + "' is declared differently in " +
                        it->second.second + " and " + path);
                }
            } else if (line.rfind("@CO\t", 0) == 0) {
                merged.passthrough.push_back(line);
            } else if (line.rfind("@PG\t", 0) == 0) {
                last_program = read_group_line_id(line);
                if (programs.insert(last_program).second)
                    merged.program_lines.push_back(line);
            }
        }
        if (merged.previous_program.empty()) merged.previous_program = last_program;
    }
    return merged;
}

}  // namespace

AlignmentOutputWriter::AlignmentOutputWriter(
    const std::string& path,
    const std::string& format,
    const io::GenomeLoad& genome,
    const std::string& command_line,
    const std::string& read_group_line
) : format_(format),
    command_line_(command_line),
    read_group_line_(read_group_line),
    stdout_sink_(path == "-" ? std::make_unique<FdSink>(STDOUT_FILENO)
                             : nullptr),
    owner_(open_output_owner(path, stdout_sink_.get())) {
    read_group_id_ = read_group_line_.empty()
        ? std::string()
        : read_group_line_id(read_group_line_);
    text_ = owner_.get();
    references_.reserve(genome.references.size());
    for (const auto& reference : genome.references) {
        references_.push_back({reference.first, reference.second});
        reference_lengths_[reference.first] = reference.second;
    }
}

AlignmentOutputWriter::~AlignmentOutputWriter() = default;

void AlignmentOutputWriter::write_header(
    const InputHeader& input_header,
    bool suppress_sam_header
) {
    require_no_read_group_conflict(input_header.passthrough, read_group_id_);
    // Our @PG takes the first free ID, as htslib's sam_hdr_pg_id picks it.
    std::unordered_set<std::string> taken;
    for (const std::string& line : input_header.program_lines)
        taken.insert(read_group_line_id(line));
    std::string program_id = "flashalign";
    for (int k = 1; taken.count(program_id) != 0; ++k)
        program_id = "flashalign." + std::to_string(k);
    // The inputs' @PG lines go before ours; the declared @RG after it and
    // before the passthrough @RG/@CO lines.
    std::vector<std::string> header_lines =
        output::sam_header(references_, "flashalign", FA_VERSION,
                           command_line_, program_id,
                           input_header.previous_program);
    std::string program = std::move(header_lines.back());
    header_lines.pop_back();
    for (const std::string& line : input_header.program_lines)
        header_lines.push_back(line);
    header_lines.push_back(std::move(program));
    if (!read_group_line_.empty()) header_lines.push_back(read_group_line_);
    for (const std::string& line : input_header.passthrough)
        header_lines.push_back(line);
    if (format_ == "sam" && !suppress_sam_header) {
        for (const auto& line : header_lines) *text_ << line << '\n';
        require_output(*text_, "write");
    }
}

AlignmentOutputContext AlignmentOutputWriter::context(
    bool include_unmapped,
    bool paf_cigar,
    bool copy_comment,
    output::SamEmitOptions sam,
    std::vector<std::string> own_tag_keys
) {
    sam.read_group_id = read_group_id_;
    own_tag_keys_ = std::make_unique<OwnTagKeys>(std::move(own_tag_keys));
    // emit_secondary is already set by the caller from --secondary.
    return {
        *text_,
        &reference_lengths_,
        format_,
        include_unmapped,
        paf_cigar,
        copy_comment,
        sam,
        own_tag_keys_.get(),
    };
}

void AlignmentOutputWriter::close() {
    text_->flush();
    require_output(*text_, "close");
}

void AlignmentOutputWriter::flush() {
    text_->flush();
    require_output(*text_, "flush");
}

namespace {

// Writes one read's records in the requested format. Returns the number of
// records written.
int64_t write_alignment_records(
    const AlignmentOutputContext& context,
    const io::FastxRecord& read,
    const AlignResult& result
) {
    int64_t written = 0;
    output::validate_one_level_hypotheses(result);
    if (context.format == "sam") {
        if (result.mapped() || context.include_unmapped) {
            std::string line;
            const output::PassthroughTags passthrough =
                passthrough_tags(read, context.own_tag_keys);
            // -y: the comment goes on every record of the read, verbatim and
            // unchecked as in minimap2.
            std::string record_tags;
            if (context.copy_comment && !read.comment.empty())
              record_tags = read.comment;
            // Sized for the primary record; other records may still grow it.
            line.reserve(read.seq.size() + read.qual.size() +
                         result.cigar.size() + read.name.size() +
                         passthrough.whole.size() + record_tags.size() + 256);
            output::SamEmitOptions emit_options = context.sam;
            // A read that already carries an RG:Z keeps it.
            if (sam_tags_declare_read_group(passthrough.clipped) ||
                sam_tags_declare_read_group(record_tags))
                emit_options.read_group_id.clear();
            output::append_sam_records(
                line, read.name, result, read.seq, read.qual,
                passthrough, emit_options, record_tags);
            context.text << line;
            written += std::count(line.begin(), line.end(), '\n');
        }
    } else {
        if (result.mapped()) {
            output::write_paf_record(
                context.text, read, result, *context.reference_lengths,
                context.paf_cigar, output::EmittedRole::Primary,
                context.copy_comment);
            ++written;
            for (const auto& supplementary : result.supplementary) {
                if (!supplementary.mapped()) continue;
                output::write_paf_record(
                    context.text, read, supplementary,
                    *context.reference_lengths, context.paf_cigar,
                    output::EmittedRole::Supplementary,
                    context.copy_comment);
                ++written;
            }
            // With --secondary yes: each secondary, then its supplementaries,
            // in the same order as the SAM writer.
            if (context.sam.emit_secondary) {
                for (const auto& secondary : result.secondary) {
                    if (!secondary.mapped()) continue;
                    output::write_paf_record(
                        context.text, read, secondary,
                        *context.reference_lengths, context.paf_cigar,
                        output::EmittedRole::Secondary,
                        context.copy_comment);
                    ++written;
                    for (const auto& supplementary : secondary.supplementary) {
                        if (!supplementary.mapped()) continue;
                        output::write_paf_record(
                            context.text, read, supplementary,
                            *context.reference_lengths, context.paf_cigar,
                            output::EmittedRole::SecondarySupplementary,
                            context.copy_comment);
                        ++written;
                    }
                }
            }
        } else if (context.include_unmapped) {
            output::write_paf_no_hit_record(context.text, read);
            ++written;
        }
    }
    require_output(context.text, "write");
    return written;
}

}  // namespace

}  // namespace fa::cpu::cli

/**** report ****/
// What a run reports: --stats and the footer, the --show-config text.

namespace fa::cpu::cli {

namespace {

// Peak resident set size of this process, in bytes.
int64_t peak_rss_bytes() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0)
    return 0;
#ifdef __APPLE__
  return static_cast<int64_t>(ru.ru_maxrss);
#else
  return static_cast<int64_t>(ru.ru_maxrss) * 1024;
#endif
}

// User plus system time of every thread of the process.
double cpu_seconds() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0)
    return 0.0;
  return static_cast<double>(ru.ru_utime.tv_sec + ru.ru_stime.tv_sec) +
         1e-6 * static_cast<double>(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);
}

double gigabytes(int64_t bytes) {
  return static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
}

// --stats: the run summary block, written to stderr.
void write_stats(const RunSummary& summary, std::ostream& out) {
  const auto row = [&](const char* key, const std::string& value) {
    out << "  " << std::left << std::setw(22) << key << value << "\n";
  };
  const auto fmt = [](const char* format, double value) {
    char buf[64];
    std::snprintf(buf, sizeof buf, format, value);
    return std::string(buf);
  };
  const auto i2s = [](int64_t v) { return std::to_string(v); };
  const int64_t unmapped =
      summary.reads > summary.mapped_reads
          ? summary.reads - summary.mapped_reads
          : 0;
  const double mapped_pct =
      summary.reads > 0 ? 100.0 * static_cast<double>(summary.mapped_reads) /
                              static_cast<double>(summary.reads)
                        : 0.0;
  out << "[flashalign] run summary\n";
  row("preset:", summary.preset);
  row("output format:", summary.output_format);
  row("reads:", i2s(summary.reads));
  row("bases:", i2s(summary.bases));
  row("mapped reads:",
      i2s(summary.mapped_reads) + " (" + fmt("%.1f", mapped_pct) + "%)");
  row("unmapped reads:", i2s(unmapped));
  row("records written:", i2s(summary.records));
  row("supplementary:", i2s(summary.supplementary_records));
  row("secondary:", i2s(summary.secondary_records));
  if (summary.rna) {
    row("spliced reads:", i2s(summary.spliced_reads));
    row("junctions:", i2s(summary.junctions));
  }
  row("wall time:", fmt("%.2f s", summary.wall_seconds));
  const double map_s = summary.map_seconds;
  row("throughput:",
      fmt("%.1f reads/s, ",
          map_s > 0.0 ? static_cast<double>(summary.reads) / map_s : 0.0) +
          fmt("%.1f Mbp/s",
              map_s > 0.0 ? static_cast<double>(summary.bases) / map_s / 1e6
                          : 0.0));
  row("peak RSS:", fmt("%.1f GB", gigabytes(summary.peak_rss_bytes)));
}

} // namespace

void write_run_footer(int argc, char** argv, double real_seconds) {
  std::string cmd;
  for (int i = 0; i < argc; ++i) {
    if (i > 0)
      cmd += ' ';
    cmd += argv[i];
  }
  std::fprintf(stderr, "[flashalign] Version: %s\n", FA_VERSION);
  std::fprintf(stderr, "[flashalign] CMD: %s\n", cmd.c_str());
  std::fprintf(stderr,
               "[flashalign] Real time: %.1f sec; CPU: %.1f sec; "
               "Peak RSS: %.1f GB\n",
               real_seconds, cpu_seconds(), gigabytes(peak_rss_bytes()));
}

namespace {

// The resolved configuration as text, each value tagged with its source
// (builtin, preset, index, cli, ...). It needs no aligner, so no reference is
// opened or indexed.
std::string resolved_config_text(const AlignOptions& opt,
                                 const fa::cpu::options::ResolvedOptions& cfg,
                                 fa::cpu::api::BackendMode mode) {
  const auto& mapping = cfg.long_read();
  const auto* rna = cfg.rna();
  const bool is_rna = mode == fa::cpu::api::BackendMode::RnaSplice;
  // (k, s) come from a loaded index, otherwise from the preset.
  const bool seeding_from_index = cfg.index.from_index;

  std::ostringstream os;
  const auto row = [&](const char* key, const std::string& val,
                       const char* src) {
    os << "  " << std::left << std::setw(24) << key;
    // A key or value that fills its column still gets a separating space.
    if (std::string_view(key).size() >= 24)
      os << ' ';
    os << std::setw(16) << val;
    if (val.size() >= 16)
      os << ' ';
    os << "[" << src << "]\n";
    os << std::right; // reset
  };
  const auto i2s = [](long long v) { return std::to_string(v); };
  // Default precision, so 0.8 prints "0.8".
  const auto f2s = [](double value) {
    std::ostringstream text;
    text << value;
    return text.str();
  };
  const auto b2s = [](bool v) { return std::string(v ? "1" : "0"); };
  const auto src = [](bool given, const char* owner) {
    return given ? "explicit" : owner;
  };
  const auto dp_src = [&](const std::optional<int>& o) {
    return src(o.has_value(), "preset");
  };

  row("mode", is_rna ? "rna_long" : "dna_long", "resolved");
  row("preset", opt.preset, opt.preset_source.c_str());
  // asm5, asm10, asm20, ava-ont and ava-hifi are experimental.
  const bool assembly = fa::cpu::api::is_assembly_preset(opt.preset);
  const bool overlap = fa::cpu::api::is_overlap_preset(opt.preset);
  row("preset_status", assembly || overlap ? "experimental" : "supported",
      "preset");
  row("index_source", seeding_from_index ? "faix" : "built", "cli");
  row("k", i2s(cfg.index.k), seeding_from_index ? "index" : "preset");
  row("syncmer_s", i2s(cfg.index.syncmer_s),
      seeding_from_index ? "index" : "preset");
  row("syncmer_downsample", i2s(cfg.index.syncmer_downsample), "builtin");
  row("output_format", opt.format, "cli");
  row("threads", i2s(cfg.common.num_threads),
      src(opt.threads && *opt.threads > 0, "builtin"));

  row("min_support", i2s(cfg.common.min_support),
      src(opt.min_support.has_value(), "preset"));
  // With a vote spacing S a read of length L votes max(F, ceil(L / S)) seeds
  // per strand, F the preset's seeds per strand (rna/map.cpp).
  const int vote_spacing =
      is_rna && rna != nullptr ? rna->vote_seed_spacing : 0;
  if (vote_spacing > 0 && cfg.common.max_query_seeds_per_strand > 0)
    row("max_query_seeds",
        "max(" + i2s(cfg.common.max_query_seeds_per_strand) + ", ceil(L/" +
            i2s(vote_spacing) + "))",
        "rule");
  else
    row("max_query_seeds", i2s(cfg.common.max_query_seeds_per_strand),
        src(opt.vote_seeds.has_value(), "preset"));
  if (is_rna && rna != nullptr)
    row("vote_seed_spacing", i2s(rna->vote_seed_spacing),
        src(opt.vote_spacing.has_value(), "preset"));
  // Above kVoteSeedNestBase a DNA selection keeps the kVoteSeedNestBase one
  // and fills it (seed/seed.h).
  if (!is_rna)
    row("vote_seed_sampler",
        cfg.common.max_query_seeds_per_strand > fa::cpu::lr::kVoteSeedNestBase
            ? "nested"
            : "original",
        "derived");
  // max_seed_occ is the effective vote occurrence cap; 0 means no filtering.
  fa::cpu::lr::LongOccPolicyConfig occ_cfg;
  occ_cfg.policy = mapping.long_occ_policy;
  occ_cfg.primary_occ_cap = mapping.long_primary_occ_cap;
  occ_cfg.platform_occ_cap = mapping.long_occ_cap;
  const bool occ_explicit = opt.max_vote_occ.has_value();
  // The default DNA rule is max(preset floor, index quantile), capped by an
  // optional ceiling. Until an index has been scanned (e.g. --show-config),
  // the row prints the rule instead of a number.
  const auto& occ_index = mapping.long_occ_index_resolution;
  const bool occ_rule =
      !is_rna && !occ_explicit &&
      mapping.long_occ_policy == fa::cpu::lr::LongOccPolicy::Platform;
  if (occ_rule && !occ_index.resolved) {
    // As resolve_index_occ_cap_into: a ceiling applies when set and not
    // below the floor.
    const bool ceiling = mapping.long_occ_ceiling > 0 &&
                         mapping.long_occ_ceiling >= mapping.long_occ_cap;
    std::ostringstream rule;
    if (ceiling)
      rule << "min(" << mapping.long_occ_ceiling << ", ";
    rule << "max(" << mapping.long_occ_cap << ", index quantile f="
         << fa::cpu::lr::kLongOccQuantileF << ')';
    if (ceiling)
      rule << ')';
    row("max_seed_occ", rule.str(), "rule");
  } else {
    row("max_seed_occ",
        i2s(fa::cpu::lr::effective_long_primary_occ_cap(occ_cfg)),
        occ_index.resolved ? "index" : (occ_explicit ? "explicit" : "preset"));
  }
  row("occ_policy", fa::cpu::lr::long_occ_policy_name(mapping.long_occ_policy),
      occ_explicit ? "explicit" : "preset");
  row("long_occ_cap", i2s(mapping.long_occ_cap),
      occ_index.resolved ? "index" : "preset");
  row("long_primary_occ_cap", i2s(mapping.long_primary_occ_cap),
      occ_explicit ? "explicit" : "builtin");
  // How the index produced the cap.
  if (occ_index.resolved) {
    std::ostringstream detail;
    detail << "f=" << occ_index.quantile_f
           << " n=" << occ_index.distinct_keys
           << " raw=" << occ_index.raw_cap
           << " floor=" << occ_index.floor_cap;
    if (occ_index.ceiling_cap > 0)
      detail << " ceiling=" << occ_index.ceiling_cap;
    detail << " cap=" << occ_index.resolved_cap;
    row("seed_occ_quantile", detail.str(), "index");
  }
  // The vote's empty-tile rescue M; 0 on the RNA presets.
  row("tile_rescue_occ", i2s(mapping.dna_tile_rescue_occ), "preset");
  // The overlap presets' two mechanisms and --dual, printed only under them.
  if (overlap) {
    row("skip_self", b2s(mapping.skip_self), "preset");
    row("all_chains", b2s(mapping.all_chains), "preset");
    row("dual", mapping.dual ? "yes" : "no",
        src(opt.dual.has_value(), "preset"));
  }

  if (!is_rna)
    row("vote_ratio", f2s(mapping.vote_admission_ratio),
        src(opt.dna_vote_admission_ratio.has_value(), "preset"));
  row("vote_diag_bin_width", i2s(mapping.vote_diag_bin_width),
      src(opt.vote_diag_bin_width.has_value(), "preset"));
  // --dw turns the adaptive model off.
  row("vote_diag_adaptive", b2s(mapping.vote_diag_bin_width_adaptive),
      src(opt.vote_diag_bin_width.has_value(), "preset"));
  row("vote_diag_slope_den", i2s(mapping.vote_diag_slope_den),
      src(opt.vote_diag_slope_den.has_value(), "preset"));
  row("vote_diag_width_max", i2s(mapping.vote_diag_width_max),
      src(opt.vote_diag_width_max.has_value(), "preset"));
  row("mapq.output_range", "0..60", "fixed");

  // On a splice preset, the fine chain's query gap.
  row("chain_max_gap", i2s(mapping.cigar_local_interval_anchor_chain_max_gap),
      src(opt.dp_max_gap.has_value(), "preset"));
  row("interval_pad", i2s(mapping.cigar_local_interval_anchor_interval_pad),
      "builtin");
  // The bound the vote catalogue runs with.
  if (is_rna)
    row("chain_max_cands", i2s(mapping.chain_max_candidates_per_window),
        "builtin");
  else
    row("chain_max_cands", i2s(fa::cpu::lr::dna_chain_max_candidates(mapping)),
        src(opt.max_cands.has_value(), overlap ? "preset" : "derived"));
  // -b's band, only where its stage runs.
  if (!is_rna &&
      fa::cpu::api::lane_runs_stage(
          fa::cpu::api::LaneStage::CandidateChains, opt.preset))
    row("chain_band", i2s(mapping.chain_band),
        src(opt.chain_band.has_value(), "preset"));

  // lr, lr:hq: -A -B -O -E -z --score-N are the gap-fill row; the end row,
  // which prices every path, is the preset's (dp_end_row). The assembly
  // presets fill under their end row, so the two rows print as one.
  const bool fill = !is_rna;
  row("dp_match", i2s(fill ? mapping.fill_dp_match : mapping.cigar_dp_match),
      dp_src(opt.dp_match));
  row("dp_mismatch",
      i2s(fill ? mapping.fill_dp_mismatch : mapping.cigar_dp_mismatch),
      dp_src(opt.dp_mismatch));
  row("dp_score_n", i2s(fill ? mapping.fill_dp_ambi : mapping.cigar_dp_ambi),
      dp_src(opt.dp_score_n));
  row("dp_gap_open1",
      i2s(fill ? mapping.fill_dp_gap_open1 : mapping.cigar_dp_gap_open1),
      dp_src(opt.dp_gap_open1));
  row("dp_gap_extend1",
      i2s(fill ? mapping.fill_dp_gap_extend1 : mapping.cigar_dp_gap_extend1),
      dp_src(opt.dp_gap_extend1));
  row("dp_gap_open2",
      i2s(fill ? mapping.fill_dp_gap_open2 : mapping.cigar_dp_gap_open2),
      dp_src(opt.dp_gap_open2));
  row("dp_gap_extend2",
      i2s(fill ? mapping.fill_dp_gap_extend2 : mapping.cigar_dp_gap_extend2),
      dp_src(opt.dp_gap_extend2));
  row("dp_zdrop",
      i2s(fill ? mapping.fill_dp_tail_zdrop : mapping.cigar_dp_tail_zdrop),
      dp_src(opt.dp_zdrop));
  row("dp_end_bonus", i2s(mapping.cigar_dp_tail_end_bonus),
      dp_src(opt.dp_end_bonus));
  // On a splice preset dp_bw is the fine chain's band, and -G sets both.
  const bool intron_typed = is_rna && opt.max_intron.has_value();
  row("dp_bw", i2s(mapping.cigar_dp_bw),
      src(opt.dp_bw.has_value() || intron_typed, "preset"));
  row("dp_bw_long", i2s(mapping.cigar_dp_bw_long),
      src(opt.dp_bw_long.has_value() || intron_typed, "preset"));
  row("dp_max_gap", i2s(mapping.cigar_dp_max_gap), dp_src(opt.dp_max_gap));
  row("dp_min_score", i2s(mapping.cigar_dp_min_dp_max),
      dp_src(opt.dp_min_score));
  row("dp_split_min_anchors", i2s(mapping.cigar_dp_split_min_anchors),
      "builtin");
  // A splice preset's inversion Z-drop is splice_inv_zdrop.
  if (fill)
    row("dp_inversion_zdrop", i2s(mapping.fill_dp_inversion_zdrop),
        dp_src(opt.dp_zdrop));
  if (fill && !assembly)
    row("dp_end_row",
        "A" + i2s(mapping.cigar_dp_match) + " B" +
            i2s(mapping.cigar_dp_mismatch) + " N" +
            i2s(mapping.cigar_dp_ambi) + " O" +
            i2s(mapping.cigar_dp_gap_open1) + "," +
            i2s(mapping.cigar_dp_gap_open2) + " E" +
            i2s(mapping.cigar_dp_gap_extend1) + "," +
            i2s(mapping.cigar_dp_gap_extend2) + " z" +
            i2s(mapping.cigar_dp_tail_zdrop) + "," +
            i2s(mapping.cigar_dp_inversion_zdrop),
        "preset");
  if (fill) {
    // A value that is neither typed nor the struct's own is the preset row's.
    const int builtin =
        std::decay_t<decltype(mapping)>{}.min_chain_score;
    row("min_chain_score", i2s(mapping.min_chain_score),
        src(opt.min_chain_score.has_value(),
            mapping.min_chain_score != builtin ? "preset" : "builtin"));
  }

  if (!is_rna) {
    row("dense_diag_min_runs", i2s(mapping.dna_dense_diag_min_runs), "builtin");
    row("dp_min_ksw_len", i2s(mapping.cigar_dp_min_ksw_len), "builtin");
    // Unless --max-chain-occ sets it, the pool gate is max_seed_occ, printed
    // by name while that row is still a rule.
    if (mapping.dna_pool_gate_occ == fa::cpu::lr::kDnaPoolGateVoteCap) {
      row("pool_gate_occ",
          occ_rule && !occ_index.resolved
              ? std::string("max_seed_occ")
              : i2s(fa::cpu::lr::dna_chain_occ_thresholds(mapping)
                        .pool_gate_occ),
          "derived");
    } else {
      row("pool_gate_occ", i2s(mapping.dna_pool_gate_occ),
          src(opt.max_chain_occ.has_value(), "builtin"));
    }
    // The within-read repeat filter's fraction, where the filter runs.
    if (!overlap)
      row("query_occ_frac", f2s(mapping.dna_query_occ_frac),
          src(opt.q_occ_frac.has_value(), "builtin"));
  }

  if (is_rna && rna != nullptr) {
    row("rna_min_intron", i2s(rna->min_intron),
        src(opt.min_intron && *opt.min_intron > 0, "preset"));
    row("rna_max_intron", i2s(rna->max_intron),
        src(opt.max_intron && *opt.max_intron > 0, "preset"));
    row("rna_junction_bed",
        rna->junction_bed.empty() ? "off" : rna->junction_bed,
        src(opt.rna_junction_bed.has_value(), "builtin"));
    row("rna_junction_bonus", i2s(rna->splice_junction_bonus),
        src(opt.rna_junction_bonus.has_value(), "preset"));
    const char* strand = "unset";
    switch (rna->strand_mode) {
    case 0:
      strand = "auto";
      break;
    case 1:
      strand = "forward";
      break;
    case 2:
      strand = "reverse";
      break;
    case 3:
      strand = "none";
      break;
    default:
      break;
    }
    row("rna_strand_mode", strand,
        src(opt.splice_strand.has_value(), "preset"));
    row("splice_transition", i2s(rna->splice_transition), "preset");
    // -z's second value; as in minimap2, a scalar -z sets it too.
    row("splice_inv_zdrop", i2s(rna->splice_inversion_zdrop),
        dp_src(opt.dp_zdrop));
    row("rna_max_loci", i2s(rna->max_locus_chains),
        src(opt.rna_max_loci.has_value(), "builtin"));
  }

  if (!is_rna) {
    // HiFi only: the chain MAPQ margin rule (dna/dna.h).
    row("dna_chain_mapq.hifi_margin", b2s(mapping.chain_mapq_hifi_margin),
        "builtin");
    // HiFi only: the late inversion probe's local-chain gate
    // (dna/dna.h).
    row("dna_inv_probe.local_gate", b2s(mapping.inversion_probe_local_gate),
        "builtin");
    row("dna_rival.pri_ratio", f2s(mapping.pri_ratio),
        src(opt.pri_ratio.has_value(), "builtin"));
    // -N bounds both lanes; printed where the preset has a -c lane (not the
    // overlap presets, which keep every chain).
    if (fa::cpu::api::lane_runs_stage(fa::cpu::api::LaneStage::BaseAlignment,
                                      opt.preset))
      row("dna_rival.secondary_max", i2s(mapping.secondary_max),
          src(opt.dna_secondary_max.has_value(), "builtin"));
  } else if (rna != nullptr) {
    row("rna_rival.pri_ratio", f2s(rna->rival_pri_ratio),
        src(opt.pri_ratio.has_value(), "builtin"));
    row("rna_rival.min_diff", i2s(rna->rival_min_diff),
        src(opt.rna_rival_min_diff.has_value(), "builtin"));
    // Fixed thresholds are printed too, tagged [fixed]: the co-primary floor,
    // then the query a second-pass placement must explain.
    row("rna_rival.chim_score",
        i2s(fa::cpu::lr::rna::kRnaCoprimaryMinChainScore), "fixed");
    row("rna_rival.chim_anch", i2s(fa::cpu::lr::rna::kRnaCoprimaryMinAnchors),
        "fixed");
    row("rna_rival.chim_qbp", i2s(fa::cpu::lr::rna::kRnaCoprimaryMinQueryBases),
        "fixed");
    row("rna_rival.chimera_min_qbp",
        i2s(fa::cpu::lr::rna::kRnaChimeraMinQueryBases), "fixed");
    row("rna_rival.realize_max", i2s(rna->rival_realize_max),
        src(opt.rna_max_loci.has_value(), "builtin"));
    row("rna_mapq.qcov_tau", f2s(rna->rna_mapq_qcov_tau), "preset");
    row("rna_mapq.vote_damp", f2s(rna->rna_mapq_disjoint_vote_damp), "preset");
  }

  row("full_read_cigar", b2s(cfg.common.enable_full_read_cigar), "cli");
  const char* cs_mode = "none";
  if (cfg.common.cs == fa::cpu::options::CsMode::Short)
    cs_mode = "short";
  else if (cfg.common.cs == fa::cpu::options::CsMode::Long)
    cs_mode = "long";
  row("cs", cs_mode, src(!opt.cs.empty(), "builtin"));
  row("emit_md", b2s(cfg.common.emit_md), src(opt.emit_md, "builtin"));
  row("emit_eqx", b2s(cfg.common.emit_eqx), src(opt.emit_eqx, "builtin"));
  return os.str();
}

// --show-config: print the resolved configuration, each value tagged with its
// source; the caller then exits 0. `reference_sequences` is the sequence count
// when known (from an index header); nullopt omits the row.
void show_config(const AlignOptions& opt,
                 const fa::cpu::options::ResolvedOptions& cfg,
                 fa::cpu::api::BackendMode mode,
                 std::optional<std::size_t> reference_sequences,
                 std::ostream& out) {
  out << "[flashalign] resolved configuration"
         " (value [source]; source in"
         " {builtin,preset,index,explicit,cli,derived,fixed,resolved})\n";
  out << resolved_config_text(opt, cfg, mode);
  if (reference_sequences) {
    out << "  " << std::left << std::setw(24) << "reference_sequences"
        << *reference_sequences << "\n";
  }
}

// --show-config on a DNA preset with a .faix: the numbers that need the loaded
// index, or that the lane derives from the read length, each tagged
// [resolved]; printed after show_config's rows.
void show_index_resolved(const fa::cpu::api::LongReadAligner& aligner,
                         std::ostream& out) {
  // The read lengths the length-dependent rows are printed at.
  static constexpr int kReadLengths[] = {1000, 10000, 100000, 1000000};
  const fa::cpu::options::ResolvedOptions cfg = aligner.config();
  const auto& mapping = cfg.long_read();
  std::ostringstream os;
  const auto row = [&](const char* key, const std::string& val) {
    os << "  " << std::left << std::setw(24) << key;
    if (std::string_view(key).size() >= 24)
      os << ' ';
    os << std::setw(16) << val;
    if (val.size() >= 16)
      os << ' ';
    os << "[resolved]\n" << std::right;
  };
  const auto cap_text = [](long long cap) {
    return cap > 0 ? std::to_string(cap) : std::string("none");
  };

  // N, the vote's occurrence cap; 0 is none.
  fa::cpu::lr::LongOccPolicyConfig occ_cfg;
  occ_cfg.policy = mapping.long_occ_policy;
  occ_cfg.primary_occ_cap = mapping.long_primary_occ_cap;
  occ_cfg.platform_occ_cap = mapping.long_occ_cap;
  const int n_cap = fa::cpu::lr::effective_long_primary_occ_cap(occ_cfg);
  row("occ_cap.N", cap_text(n_cap));
  const auto& occ = mapping.long_occ_index_resolution;
  if (occ.resolved) {
    std::ostringstream parts;
    parts << "f=" << occ.quantile_f << " raw=" << occ.raw_cap
          << " floor=" << occ.floor_cap
          << " ceiling=" << cap_text(occ.ceiling_cap);
    row("occ_cap.parts", parts.str());
  }
  const fa::cpu::lr::DnaChainOccThresholds chain_occ =
      fa::cpu::lr::dna_chain_occ_thresholds(mapping);
  const auto chain_cap_text = [&](int cap) {
    return cap <= 0 || cap == INT_MAX ? std::string("none")
                                      : std::to_string(cap);
  };
  row("chain_occ.pool_gate", chain_cap_text(chain_occ.pool_gate_occ));

  // M, the empty-tile rescue: a tile's rarest seed over N and at most M
  // (seed/seed.h VoteTileRescue).
  const int m_cap = mapping.dna_tile_rescue_occ;
  row("tile_rescue.M", cap_text(m_cap));
  if (m_cap > 0 && n_cap > 0 && n_cap < m_cap) {
    std::ostringstream window;
    window << "occ in (" << n_cap << ", " << m_cap
           << "], L >= " << fa::cpu::lr::kDnaTileRescueMinLen
           << ", budget " << fa::cpu::lr::kDnaTileRescueBudget;
    row("tile_rescue.window", window.str());
  } else {
    row("tile_rescue.window", "none");
  }

  // W(L) at each of kReadLengths; the chain's band/reach, one for every L.
  std::ostringstream rule;
  if (mapping.vote_diag_bin_width_adaptive && mapping.vote_diag_slope_den > 0)
    rule << "clamp(L / " << mapping.vote_diag_slope_den << ", "
         << std::max(1, mapping.vote_diag_bin_width) << ", "
         << mapping.vote_diag_width_max << ')';
  else
    rule << std::max(1, mapping.vote_diag_bin_width);
  row("vote_diag_width.rule", rule.str());
  std::ostringstream width;
  for (const int length : kReadLengths) {
    const fa::cpu::api::DnaReadGeometry g = aligner.dna_read_geometry(length);
    const char* sep = length == kReadLengths[0] ? "" : " ";
    width << sep << "L=" << length << ':' << g.vote_diag_bin_width;
  }
  row("vote_diag_width.at", width.str());
  const fa::cpu::api::DnaChainGeometry chain = aligner.dna_chain_geometry();
  row("chain.band/reach",
      std::to_string(chain.band) + '/' + std::to_string(chain.reach));
  out << os.str();
}

} // namespace

} // namespace fa::cpu::cli

/**** align_command ****/
// The run's plan and inputs, the runs that end early, and run_align last.

namespace fa::cpu::cli {

namespace {
// Lends idle decode-team threads to the compute window. The read source (and
// its decode team) starts before the compute window exists, so the team calls
// through this hub, which answers "no work" until install() and after
// revoke(). revoke() waits for loaned threads to leave the target, so the
// window can be destroyed while the team still runs; it is called only after
// the pipeline has joined, so no new work can keep it waiting.
class HelperHub {
public:
  bool call(int tid) {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return target_ ? target_(tid) : false;
  }
  void install(std::function<bool(int)> target) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    target_ = std::move(target);
  }
  void revoke() {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    target_ = nullptr;
  }

private:
  std::shared_mutex mutex_;
  std::function<bool(int)> target_;
};

// The smallest -t at which gzip decoding that is not BGZF borrows mapping
// threads. BGZF borrows them in every windowed run whose inputs are all
// eligible.
constexpr int kDecodeLoanMinThreads = 32;

// Threads the decode team borrows from the -t budget: -t/8, clamped to [1, 8].
int decode_loan_threads(int num_threads) {
  const int wanted = num_threads / 8;
  if (wanted < 1)
    return 1;
  return wanted > 8 ? 8 : wanted;
}

// A loan needs every input to be eligible for a parallel source the run
// allows (`gz`: member-parallel gzip; `bgzf`: block-parallel BGZF, BAM
// included); otherwise the loaned threads would sit idle while the serial
// reader runs.
bool decode_loan_can_be_fed(const std::vector<std::string>& reads_paths,
                            bool gz, bool bgzf) {
  if (reads_paths.empty())
    return false;
  for (const std::string& path : reads_paths) {
    if (bgzf && ::fa::cpu::io::bgzf_eligible(path))
      continue;
    // Generic-gzip BAM uses strict serial decoding, so it cannot feed a loan.
    if (!(gz && ::fa::cpu::io::gz_member_parallel_eligible(path) &&
          !::fa::cpu::io::path_is_bam(path)))
      return false;
  }
  return true;
}

// Decompressed bytes the decoder may hold ahead of the reader: 64 MiB per
// thread, clamped to [512 MiB, 4 GiB], unless --io-staging sets it.
uint64_t decode_staging_budget(int num_threads, int io_staging_mib) {
  if (io_staging_mib > 0)
    return static_cast<uint64_t>(io_staging_mib) << 20;
  const uint64_t scaled =
      static_cast<uint64_t>(num_threads < 1 ? 1 : num_threads) * (64ull << 20);
  if (scaled < (512ull << 20))
    return 512ull << 20;
  return scaled > (4ull << 30) ? (4ull << 30) : scaled;
}

// Mapping threads: -t, or all cores but two (for the reader and writer).
int resolve_num_threads(const AlignOptions& opt) {
  return (opt.threads && *opt.threads > 0)
             ? *opt.threads
             : std::max(1, ::fa::cpu::threading::default_thread_count() - 2);
}

// Batches in flight in the compute stage. Forced to 1 at -t 1, so that run
// uses a single OS thread (run_pipeline_serial).
int resolve_batch_window(const AlignOptions& opt, int num_threads) {
  int window = opt.batch_window < 1 ? 1 : opt.batch_window;
  if (num_threads == 1)
    window = 1;
  return window;
}

// The input decode plan, minus the loan callback. The parallel decoders
// borrow threads only with a compute window and inputs that can feed them.
// The member-parallel gzip decoder also needs -t >= kDecodeLoanMinThreads:
// below that, serial decoding on the reader thread hides behind compute and
// its extra work makes lending threads a net loss; its loan is then at least
// four threads, so one thread stuck on a long read cannot starve the reader.
// The BGZF decoder costs no extra work, and its reader inflates the run it
// needs itself when no worker has claimed it, so it borrows from -t 2.
fa::cpu::io::IoPlan build_io_plan(const AlignOptions& opt, int num_threads,
                                  int batch_window) {
  fa::cpu::io::IoPlan plan;
  const bool gz = batch_window >= 2 && num_threads >= kDecodeLoanMinThreads;
  const bool bgzf = batch_window >= 2;
  const bool fed = decode_loan_can_be_fed(opt.reads_paths, gz, bgzf);
  plan.parallel_gz = gz && fed;
  plan.parallel_bgzf = bgzf && fed;
  if (!plan.parallel_gz && !plan.parallel_bgzf)
    return plan;
  plan.gz.workers = static_cast<unsigned>(decode_loan_threads(num_threads));
  plan.gz.staging_budget =
      decode_staging_budget(num_threads, opt.io_staging_mib);
  // The source raises a budget below twice this cap, so keep the cap at half
  // the budget.
  plan.gz.segment_stage_cap =
      std::min<uint64_t>(512ull << 20, plan.gz.staging_budget / 2);
  return plan;
}

// SAM always carries a CIGAR; PAF only with -c (or --cs/--MD).
bool run_needs_cigar(const AlignOptions& opt) {
  return opt.format == "sam" || opt.paf_cigar;
}

// DNA map-only reads no reference bases, so its index load skips unpacking
// the reference. CIGAR realization and RNA need them.
bool run_needs_reference_bases(const AlignOptions& opt) {
  return run_needs_cigar(opt) || fa::cpu::api::is_rna_preset(opt.preset);
}

// The SAM tag keys this run writes; a uBAM's own tags with these keys give
// way to ours.
std::vector<std::string> own_sam_tag_keys(const AlignOptions& opt) {
  std::vector<std::string> keys = {"AS", "ms", "md", "NM", "nn", "tp", "cm",
                                   "s1", "s2", "de", "SA", "XA"};
  if (!opt.cs.empty())
    keys.push_back("cs");
  if (opt.emit_md)
    keys.push_back("MD");
  if (fa::cpu::api::is_rna_preset(opt.preset))
    keys.push_back("ts");
  return keys;
}

// Parsed CLI values as typed overrides, plus the CLI's preset-dependent
// runtime defaults.
fa::cpu::api::UserOverrides build_user_overrides(const AlignOptions& opt) {
  fa::cpu::api::UserOverrides u;
  const bool need_cigar = run_needs_cigar(opt);

  // (k, s) is not set here: it comes from the .faix or the preset.
  u.num_threads = resolve_num_threads(opt);
  u.enable_full_read_cigar = need_cigar;
  if (!opt.cs.empty()) {
    u.cs = opt.cs == "long" ? fa::cpu::options::CsMode::Long
                            : fa::cpu::options::CsMode::Short;
  }
  if (opt.emit_md)
    u.emit_md = true;
  if (opt.emit_eqx)
    u.emit_eqx = true;
  if (opt.dna_vote_admission_ratio)
    u.dna_vote_admission_ratio = *opt.dna_vote_admission_ratio;
  if (opt.max_chain_occ)
    u.dna_pool_gate_occ = *opt.max_chain_occ;
  if (opt.q_occ_frac)
    u.dna_query_occ_frac = *opt.q_occ_frac;
  if (opt.min_support)
    u.min_support = *opt.min_support;
  if (opt.vote_seeds)
    u.max_query_seeds_per_strand = *opt.vote_seeds;
  if (opt.max_cands)
    u.max_cands = *opt.max_cands;
  if (opt.min_chain_score)
    u.min_chain_score = *opt.min_chain_score;
  if (opt.dual)
    u.dual = *opt.dual;
  // --max-vote-occ: 0 turns occurrence filtering off, >0 fixes the cap.
  if (opt.max_vote_occ) {
    if (*opt.max_vote_occ <= 0) {
      u.occ_policy = std::string("off");
    } else {
      u.occ_policy = std::string("fixed");
      u.primary_occ_cap = *opt.max_vote_occ;
    }
  }
  if (opt.vote_diag_bin_width)
    u.vote_diag_bin_width = *opt.vote_diag_bin_width;
  if (opt.vote_diag_slope_den)
    u.vote_diag_slope_den = *opt.vote_diag_slope_den;
  if (opt.vote_diag_width_max)
    u.vote_diag_width_max = *opt.vote_diag_width_max;
  if (opt.dp_match)
    u.dp_match = *opt.dp_match;
  if (opt.dp_mismatch)
    u.dp_mismatch = *opt.dp_mismatch;
  if (opt.dp_score_n)
    u.dp_ambi = *opt.dp_score_n;
  if (opt.dp_gap_open1)
    u.dp_gap_open1 = *opt.dp_gap_open1;
  if (opt.dp_gap_open2)
    u.dp_gap_open2 = *opt.dp_gap_open2;
  if (opt.dp_gap_extend1)
    u.dp_gap_extend1 = *opt.dp_gap_extend1;
  if (opt.dp_gap_extend2)
    u.dp_gap_extend2 = *opt.dp_gap_extend2;
  if (opt.dp_zdrop)
    u.dp_tail_zdrop = *opt.dp_zdrop;
  if (opt.dp_zdrop_inv)
    u.dp_inversion_zdrop = *opt.dp_zdrop_inv;
  if (opt.dp_end_bonus)
    u.dp_tail_end_bonus = *opt.dp_end_bonus;
  if (opt.dp_min_score)
    u.dp_min_dp_max = *opt.dp_min_score;
  if (opt.dp_bw)
    u.dp_bw = *opt.dp_bw;
  if (opt.dp_bw_long)
    u.dp_bw_long = *opt.dp_bw_long;
  if (opt.dp_max_gap)
    u.dp_max_gap = *opt.dp_max_gap;
  if (opt.chain_band)
    u.chain_band = *opt.chain_band;
  // Passed as given; the resolver validates the intron bounds.
  if (opt.min_intron)
    u.rna_min_intron = *opt.min_intron;
  if (opt.max_intron)
    u.rna_max_intron = *opt.max_intron;
  if (opt.rna_junction_bed)
    u.rna_junction_bed = *opt.rna_junction_bed;
  if (opt.rna_junction_bonus)
    u.rna_junction_bonus = *opt.rna_junction_bonus;
  if (opt.pri_ratio)
    u.pri_ratio = *opt.pri_ratio;
  if (opt.rna_rival_min_diff)
    u.rna_rival_min_diff = *opt.rna_rival_min_diff;
  if (opt.rna_max_loci)
    u.rna_max_loci = *opt.rna_max_loci;
  if (opt.vote_spacing)
    u.rna_vote_seed_spacing = *opt.vote_spacing;
  if (opt.dna_secondary_max)
    u.dna_secondary_max = *opt.dna_secondary_max;
  // "auto" (-u b) maps to Unknown, the preset default.
  if (opt.splice_strand) {
    u.rna_strand = static_cast<int>(fa::cpu::lr::rna::parse_strand_mode(
        opt.splice_strand->c_str(), fa::cpu::lr::rna::StrandMode::Unknown));
  }
  return u;
}

// The align target, told apart by its content (api::reference_kind). Of an index only the
// header is read here.
struct ResolvedInputs {
  std::string ref_path;   // FASTA to load; empty when the target is an index
  std::string index_path; // .faix to load
  bool use_index = false; // load index_path instead of building from FASTA
  bool embedded_ref = false; // the reference comes from the index: no FASTA
  // False for an --idx-no-seq index, which supports plain PAF only.
  bool index_has_reference = true;
  // The preset the .faix records, or empty. Used when no -x is given.
  std::string index_preset;
  // The index's seeding and sequence count, from its header.
  int index_k = 0;
  int index_s = 0;
  std::size_t reference_sequences = 0;
  // A multi-part (-I) index, which run_align maps one part at a time.
  bool multipart = false;
};

ResolvedInputs resolve_inputs(const AlignOptions& opt) {
  ResolvedInputs in;
  // Only --show-config runs without a target.
  if (opt.target_path.empty())
    return in;
  if (fa::cpu::api::reference_kind(opt.target_path) ==
      fa::cpu::api::ReferenceKind::Sequences) {
    in.ref_path = opt.target_path;
    return in;
  }
  const fa::cpu::api::IndexProbe probe =
      fa::cpu::api::LongReadAligner::probe_index_header(opt.target_path);
  if (!probe.error.empty())
    throw std::runtime_error(opt.target_path + ": " + probe.error);
  // An --idx-no-seq index still serves plain PAF, which needs only contig
  // names and lengths.
  if (!probe.has_reference && run_needs_cigar(opt)) {
    throw std::runtime_error(
        opt.target_path +
        ": this .faix does not embed its reference (built with --idx-no-seq), "
        "so it cannot realize a CIGAR; map to plain PAF (without "
        "-a/-c/--cs/--MD), or rebuild it with 'flashalign index'");
  }
  in.index_path = opt.target_path;
  in.use_index = true;
  in.embedded_ref = true;
  in.index_has_reference = probe.has_reference;
  in.index_preset = probe.preset;
  in.index_k = probe.k;
  in.index_s = probe.syncmer_s;
  in.reference_sequences = static_cast<std::size_t>(probe.reference_sequences);
  in.multipart = probe.multipart;
  return in;
}

// Refuses SAM text given as a reads file, which kseq would otherwise parse as
// nonsense FASTQ. Only regular uncompressed files are probed: stdin, FIFOs and
// gzip inputs are left alone.
void refuse_sam_text_reads(const std::vector<std::string>& paths) {
  for (const std::string& path : paths) {
    if (path == "-")
      continue;
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode))
      continue; // left for the reader to report
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr)
      continue;
    char head[4] = {0, 0, 0, 0};
    const std::size_t got = std::fread(head, 1, sizeof(head), file);
    std::fclose(file);
    if (got != sizeof(head) || head[0] != '@' || head[3] != '\t')
      continue;
    const std::string tag(head + 1, 2);
    if (tag != "HD" && tag != "SQ" && tag != "RG" && tag != "PG" && tag != "CO")
      continue;
    throw std::runtime_error(
        "reads file '" + path +
        "' is SAM text; flashalign reads FASTA/FASTQ (plain, gzip, bgzf) and "
        "unaligned BAM -- convert with 'samtools fastq'");
  }
}

// --show-config resolves options only: no read is mapped and no reference base
// is read; an index is known by its header, and on a DNA preset by its postings
// too (below).
void run_show_config(const AlignOptions& opt, const ResolvedInputs& in) {
  fa::cpu::api::ResolveRequest req;
  req.preset = opt.preset;
  req.user = build_user_overrides(opt);
  std::optional<size_t> reference_sequences;
  if (in.use_index) {
    req.index.has_index = true;
    req.index.k = in.index_k;
    req.index.syncmer_s = in.index_s;
    reference_sequences = in.reference_sequences;
  }
  const fa::cpu::api::ResolvedMapOptions resolved =
      fa::cpu::api::resolve_options(req);
  show_config(opt, resolved.resolved(), resolved.mode(), reference_sequences,
              std::cout);
  // A DNA preset's occurrence cap needs the index's postings, so the index
  // is loaded, without its bases. A multi-part index resolves per part and
  // prints nothing here.
  if (in.use_index && !in.multipart &&
      resolved.mode() == fa::cpu::api::BackendMode::DnaLong) {
    fa::cpu::api::LongReadAligner aligner =
        fa::cpu::api::LongReadAligner::from_index(
            {}, in.index_path, opt.min_support.value_or(3),
            /*chain_syncmer_downsample=*/1, req.user.num_threads.value_or(0),
            /*reference_bases_needed=*/false);
    aligner.reconfigure(resolved.resolved());
    show_index_resolved(aligner, std::cout);
  }
}

} // namespace

int run_align(const AlignOptions& cli_options) {
  // Only the preset is rewritten, when a .faix supplies it.
  AlignOptions opt = cli_options;
  refuse_sam_text_reads(opt.reads_paths);
  const ResolvedInputs in = resolve_inputs(opt);
  // A multi-part index re-reads the reads once per part, so stdin cannot
  // be used.
  if (in.multipart) {
    for (const std::string& path : opt.reads_paths) {
      if (path == "-") {
        throw std::runtime_error(
            "reads from stdin ('-') cannot be mapped against a multi-part "
            "index, which reads the input once per part; give the reads as "
            "files");
      }
    }
  }
  // The preset is an explicit -x, else the one the .faix records, else lr.
  // An unknown recorded preset is an error rather than a silent lr.
  if (opt.preset_source == "builtin" && !in.index_preset.empty()) {
    if (!fa::cpu::api::preset_is_valid(in.index_preset)) {
      throw std::runtime_error(
          in.index_path + ": this .faix records preset '" + in.index_preset +
          "', which this build does not know; pass an explicit -x, or rebuild "
          "the index with 'flashalign index'");
    }
    opt.preset = in.index_preset;
    opt.preset_source = "index";
  }
  // parse.cpp's checks that read the preset, for the preset an index records.
  check_max_cands_range(opt);
  if (const auto refusal = option_stage_refusal(opt.given_options, opt.preset))
    throw UsageError(*refusal);
  // RNA mapping needs the reference bases even without a CIGAR. Checked
  // after the preset is settled and before --show-config, so a dry run
  // refuses what a real run would.
  if (in.use_index && !in.index_has_reference &&
      fa::cpu::api::is_rna_preset(opt.preset)) {
    throw std::runtime_error(
        in.index_path +
        ": this .faix does not embed its reference (built with --idx-no-seq) "
        "and the spliced presets map from it; rebuild it with 'flashalign "
        "index -x " + opt.preset + "'");
  }
  if (opt.show_config) {
    run_show_config(opt, in);
    return 0;
  }
  // The BAM inputs' header lines, merged before the index loads so that
  // conflicting read groups fail fast.
  const InputHeader input_header = read_input_headers(opt.reads_paths);
  // Taken before the reference or index is loaded, for end-to-end wall time.
  const auto program_t0 = std::chrono::steady_clock::now();
  // A multi-part index: its table names every contig and locates each part's image.
  fa::cpu::FaixMultipartTable table;
  if (in.multipart) {
    fa::cpu::FaixLoadStatus status;
    table = fa::cpu::load_faix_multipart_table(in.index_path, &status);
    if (!table.loaded()) {
      throw std::runtime_error(in.index_path + ": " +
                               fa::cpu::faix_load_error_message(status));
    }
  }
  // A .faix target needs no FASTA: the reference comes from the index.
  const bool use_index = in.use_index;
  const bool embedded_ref = in.embedded_ref;
  fa::cpu::io::GenomeLoad genome;
  if (!embedded_ref) {
    genome = fa::cpu::io::load_fasta_genome(in.ref_path);
  }
  // align_batch maps a batch's sequences, given their read-name hashes and,
  // where the preset leaves a read out of its own vote, their own contigs, and
  // where it prints a pair once, their names' ranks among the contigs.
  std::unique_ptr<fa::cpu::api::LongReadAligner> dna_aligner;
  std::function<std::vector<fa::cpu::AlignResult>(
      const std::vector<std::string>&, const std::vector<std::uint32_t>&,
      const std::vector<int>&, const std::vector<int>&)>
      align_batch;

  // Fills the output's reference table when no FASTA was loaded.
  auto fill_refs_from = [&](const std::vector<std::string>& names,
                            const std::vector<int64_t>& lens) {
    if (!genome.references.empty())
      return;
    if (names.size() != lens.size()) {
      throw std::runtime_error(
          "flashalign: index reference names/lengths shape mismatch");
    }
    genome.references.reserve(names.size());
    for (size_t i = 0; i < names.size(); ++i) {
      genome.references.emplace_back(names[i], lens[i]);
    }
  };

  fa::cpu::api::BackendMode backend_mode = fa::cpu::api::BackendMode::DnaLong;

  // Index parts to map against: one, or each part of a multi-part index in
  // turn, as minimap2 does without --split-prefix. Reads are re-read per part
  // and output is part-major; each part is its own engine run, with its own
  // MAPQ rivals and occurrence cap. The @SQ table covers the whole reference.
  std::vector<fa::cpu::api::IndexImageSpan> part_spans(1);
  if (in.multipart) {
    part_spans.clear();
    for (const fa::cpu::FaixPartEntry& part : table.parts)
      part_spans.push_back({part.image_offset, part.image_bytes});
  }

  // Resolved once, on the first attach. With an index, (k, s) is read from
  // it; otherwise the resolver takes it from the preset. Every part of a
  // multi-part index has the same seeding.
  std::optional<fa::cpu::api::ResolvedMapOptions> resolved;
  // Loads part `part` (or builds the in-memory index for a FASTA target),
  // installs the resolved options and rebinds align_batch. The caller frees
  // the previous part first.
  auto attach_part = [&](size_t part) {
    fa::cpu::api::ResolveRequest req;
    req.preset = opt.preset;
    req.user = build_user_overrides(opt);
    if (use_index) {
      // The reference unpacks on -t threads, not on every core.
      dna_aligner = std::make_unique<fa::cpu::api::LongReadAligner>(
          fa::cpu::api::LongReadAligner::from_index(
              std::move(genome.sequences), in.index_path,
              opt.min_support.value_or(3),
              /*chain_syncmer_downsample=*/1,
              req.user.num_threads.value_or(0),
              run_needs_reference_bases(opt), part_spans[part]));
      if (!resolved) {
        req.index.has_index = true;
        req.index.k = dna_aligner->config().index.k;
        req.index.syncmer_s = dna_aligner->get_index_syncmer_s();
        resolved = fa::cpu::api::resolve_options(req);
        // As minimap2, note when the index's (k, s) differ from the preset's.
        if (!opt.quiet) {
          const fa::cpu::api::PresetSeeding seeding =
              fa::cpu::api::resolve_preset_seeding(opt.preset);
          const int got_k = dna_aligner->config().index.k;
          const int got_s = dna_aligner->get_index_syncmer_s();
          if (seeding.k != got_k) {
            std::fprintf(
                stderr,
                "[flashalign] note: prebuilt index k=%d overrides preset "
                "%s (k=%d); rebuild the index to change k\n",
                got_k, opt.preset.c_str(), seeding.k);
          }
          if (seeding.syncmer_s != got_s) {
            std::fprintf(stderr,
                         "[flashalign] note: prebuilt index s=%d overrides "
                         "preset %s (s=%d); rebuild the index to change s\n",
                         got_s, opt.preset.c_str(), seeding.syncmer_s);
          }
        }
      }
    } else {
      req.index.has_index = false;
      resolved = fa::cpu::api::resolve_options(req);
      dna_aligner = std::make_unique<fa::cpu::api::LongReadAligner>(
          std::move(genome.sequences), resolved->k(),
          resolved->min_support(), resolved->syncmer_s(),
          resolved->syncmer_downsample(), req.user.num_threads.value_or(0));
    }
    dna_aligner->reconfigure(resolved->resolved());
    backend_mode = resolved->mode();
    auto* p = dna_aligner.get();
    align_batch = [p](const std::vector<std::string>& seqs,
                      const std::vector<std::uint32_t>& name_hashes,
                      const std::vector<int>& self_contigs,
                      const std::vector<int>& name_ranks) {
      return p->align_batch(
          seqs, &name_hashes, self_contigs.empty() ? nullptr : &self_contigs,
          name_ranks.empty() ? nullptr : &name_ranks);
    };
  };
  attach_part(0);
  // The output's reference table is always the whole reference's.
  if (in.multipart) {
    std::vector<int64_t> lengths;
    lengths.reserve(table.chrom_names.size());
    for (size_t c = 0; c < table.chrom_names.size(); ++c) {
      lengths.push_back(static_cast<int64_t>(table.chr_offsets[c + 1] -
                                             table.chr_offsets[c]));
    }
    fill_refs_from(table.chrom_names, lengths);
  } else {
    fill_refs_from(dna_aligner->chromosome_names(),
                   dna_aligner->chromosome_lengths());
  }

  AlignmentOutputWriter output(opt.output_path, opt.format, genome,
                               opt.command_line, opt.read_group_line);

  // With a window of 1, each batch is aligned and joined in turn, so it costs
  // its slowest read. With two or more, a worker that runs out of reads in
  // batch k starts on k+1; output order is unchanged. Resolved before the
  // reader, since only a window can borrow decode threads.
  const int num_threads = dna_aligner->config().common.num_threads;
  const int batch_window = resolve_batch_window(opt, num_threads);

  // State shared across passes over the reads (one pass per index part).
  // The reader is rebound per pass; the output context is bound once, after
  // the header. batch_bp bounds buffered input, not alignment cost.
  const int64_t batch_bp = opt.batch_bp;
  std::unique_ptr<fa::cpu::io::BatchReader> batch_reader;
  std::optional<AlignmentOutputContext> output_context;
  size_t current_part = 0;
  // Reads seen in the current pass: written by the reader stage, read by the
  // writer stage.
  std::atomic<int64_t> seen{0};
  // The --stats counters, written by the writer stage only.
  RunSummary summary;
  const auto align_t0 = std::chrono::steady_clock::now();

  struct ReadBatch {
    std::vector<fa::cpu::io::FastxRecord> records;
  };
  struct AlignedBatch {
    std::vector<fa::cpu::io::FastxRecord> records;
    std::vector<fa::cpu::AlignResult> results;
  };

  // The records of each in-flight batch, whose sequences were moved into the
  // window. Between collect_windowed() calls its size equals in_flight(),
  // including after a batch throws. Caller thread only.
  std::deque<std::vector<fa::cpu::io::FastxRecord>> windowed_records;
  // Declared last so its workers are joined first on unwind. One per pass.
  std::unique_ptr<fa::cpu::api::WindowedAlignSession> session;

  auto produce = [&]() -> std::optional<ReadBatch> {
    auto records = batch_reader->next(seen);
    if (!records)
      return std::nullopt;
    ReadBatch batch;
    batch.records = std::move(*records);
    return batch;
  };

  // The compute stage runs on the caller thread while the reader and writer
  // stages overlap it. The batch's sequences are moved into the aligner and
  // back rather than copied. Workers see a read's name only as its X31 hash
  // (seed/seed.h), which breaks exact ties as minimap2's qname does.
  const auto read_name_hashes =
      [](const std::vector<fa::cpu::io::FastxRecord>& records) {
        std::vector<std::uint32_t> hashes;
        hashes.reserve(records.size());
        for (const fa::cpu::io::FastxRecord& record : records)
          hashes.push_back(fa::cpu::lr::tie_name_hash(record.name));
        return hashes;
      };
  // Where the preset leaves a read out of its own vote: each read's own contig
  // in the index part being mapped, by exact name, -1 when the read is not in
  // it. Empty otherwise.
  const bool skip_self = resolved->resolved().long_read().skip_self;
  const auto read_self_contigs =
      [&](const std::vector<fa::cpu::io::FastxRecord>& records) {
        std::vector<int> contigs;
        if (!skip_self)
          return contigs;
        contigs.reserve(records.size());
        for (const fa::cpu::io::FastxRecord& record : records)
          contigs.push_back(dna_aligner->contig_id(record.name));
        return contigs;
      };
  // Where the preset prints a pair once (--dual=no): each read's rank in the
  // pair order among the contigs of the index part being mapped, from its name
  // and the length of its sequence in `seqs`. Empty otherwise.
  const bool skip_dual = resolved->resolved().long_read().all_chains &&
                         !resolved->resolved().long_read().dual;
  const auto read_name_ranks =
      [&](const std::vector<fa::cpu::io::FastxRecord>& records,
          const std::vector<std::string>& seqs) {
        std::vector<int> ranks;
        if (!skip_dual)
          return ranks;
        ranks.reserve(records.size());
        for (std::size_t i = 0; i < records.size(); ++i)
          ranks.push_back(dna_aligner->contig_dual_rank(
              records[i].name, static_cast<std::int64_t>(seqs[i].size())));
        return ranks;
      };

  // Window of 1: one barrier per batch.
  auto process_barrier = [&](ReadBatch in) -> std::optional<AlignedBatch> {
    std::vector<std::string> seqs;
    seqs.reserve(in.records.size());
    for (auto& rec : in.records)
      seqs.push_back(std::move(rec.seq));
    AlignedBatch done;
    done.results = align_batch(seqs, read_name_hashes(in.records),
                               read_self_contigs(in.records),
                               read_name_ranks(in.records, seqs));
    for (size_t i = 0; i < seqs.size(); ++i)
      in.records[i].seq = std::move(seqs[i]);
    done.records = std::move(in.records);
    return std::make_optional(std::move(done));
  };

  // Window of two or more: takes the oldest batch from the window and puts
  // its sequences back into its records. Caller thread only.
  auto collect_windowed = [&]() -> AlignedBatch {
    // Pop before collect(), which drops its batch before rethrowing that
    // batch's exception, so the two queues stay in step.
    AlignedBatch out;
    out.records = std::move(windowed_records.front());
    windowed_records.pop_front();
    fa::cpu::api::AlignedReadBatch done = session->collect();
    for (size_t i = 0; i < done.reads.size(); ++i)
      out.records[i].seq = std::move(done.reads[i]);
    out.results = std::move(done.results);
    return out;
  };

  // Submits a batch and collects the oldest once the window is full, so
  // in_flight() never exceeds the window.
  auto process_windowed = [&](ReadBatch in) -> std::optional<AlignedBatch> {
    std::vector<std::string> seqs;
    seqs.reserve(in.records.size());
    for (auto& rec : in.records)
      seqs.push_back(std::move(rec.seq));
    std::vector<std::uint32_t> name_hashes = read_name_hashes(in.records);
    std::vector<int> self_contigs = read_self_contigs(in.records);
    std::vector<int> name_ranks = read_name_ranks(in.records, seqs);
    windowed_records.push_back(std::move(in.records));
    session->submit(std::move(seqs), std::move(name_hashes),
                    std::move(self_contigs), std::move(name_ranks));
    if (session->in_flight() < static_cast<size_t>(batch_window))
      return std::nullopt;
    return std::make_optional(collect_windowed());
  };

  // At end of input, drains the window oldest-first.
  auto flush_windowed = [&]() -> std::optional<AlignedBatch> {
    if (session->in_flight() == 0)
      return std::nullopt;
    return std::make_optional(collect_windowed());
  };

  fa::cpu::output::SamEmitOptions emit_opts;
  emit_opts.hard_clip_supp = !opt.soft_clip_supp;
  emit_opts.emit_secondary = opt.output_secondary;

  // Writer stage, on its own thread: counters and record output, in input
  // order. Over a multi-part index the counters sum over passes, except input
  // bases, which are counted on the first pass only.
  auto consume = [&](AlignedBatch batch) {
    const auto& records = batch.records;
    const auto& results = batch.results;
    for (size_t i = 0; i < records.size(); ++i) {
      const auto& result = results[i];
      if (current_part == 0)
        summary.bases += static_cast<int64_t>(records[i].seq.size());
      if (result.mapped()) {
        ++summary.mapped_reads;
        const int64_t introns =
            std::count(result.cigar.begin(), result.cigar.end(), 'N');
        if (introns > 0) {
          ++summary.spliced_reads;
          summary.junctions += introns;
        }
      }
      for (const auto& s : result.supplementary)
        if (s.mapped())
          ++summary.supplementary_records;
      // Secondary records and their supplementaries, written only with
      // --secondary yes.
      if (opt.output_secondary && result.mapped()) {
        for (const auto& secondary : result.secondary) {
          if (!secondary.mapped())
            continue;
          ++summary.secondary_records;
          for (const auto& s : secondary.supplementary)
            if (s.mapped())
              ++summary.secondary_records;
        }
      }
      summary.records +=
          write_alignment_records(*output_context, records[i], result);
    }
    if (opt.progress && !opt.quiet) {
      if (in.multipart) {
        std::cerr << "[flashalign] part " << (current_part + 1) << "/"
                  << part_spans.size() << ": aligned " << seen
                  << " reads, wrote " << summary.records << " records\n";
      } else {
        std::cerr << "[flashalign] aligned " << seen << " reads, wrote "
                  << summary.records << " records\n";
      }
    }
  };

  // One pass over the reads per index part. The previous part is freed
  // before the next is loaded, and the header is written on the first pass.
  for (size_t part = 0; part < part_spans.size(); ++part) {
    current_part = part;
    if (part > 0) {
      dna_aligner.reset();
      attach_part(part);
    }
    if (in.multipart && !opt.quiet) {
      const fa::cpu::FaixPartEntry& entry = table.parts[part];
      std::fprintf(stderr,
                   "[flashalign] part %zu/%zu: %llu sequences, %llu bp\n",
                   part + 1, part_spans.size(),
                   static_cast<unsigned long long>(entry.contig_count),
                   static_cast<unsigned long long>(table.part_bp(part)));
    }
    seen.store(0);

    // A parallel decoder borrows some of the -t mapping threads when they are
    // idle rather than adding its own. Its team starts with the reader, before
    // the window exists, so it reaches the window through the hub, which
    // must outlive the team.
    HelperHub helper_hub;
    fa::cpu::io::IoPlan io_plan = build_io_plan(opt, num_threads, batch_window);
    int decode_loan = io_plan.parallel_gz || io_plan.parallel_bgzf
                          ? static_cast<int>(io_plan.gz.workers)
                          : 0;
    if (decode_loan > 0) {
      // Helper tids are the top of the compute width, which the window
      // reserves for them.
      const int helper_tid_base = num_threads - decode_loan;
      HelperHub* hub = &helper_hub;
      io_plan.gz.helper = [hub, helper_tid_base](unsigned ordinal) {
        return hub->call(helper_tid_base + static_cast<int>(ordinal));
      };
    }

    fa::cpu::io::SequenceReader reader(opt.reads_paths, &io_plan,
                                       /*keep_comment=*/opt.copy_comment);
    // The source may still have fallen back to serial decoding; then drop
    // the loan so the window is not left short of workers. Only the first
    // input is checked.
    if (decode_loan > 0 && !reader.decode_parallel_active())
      decode_loan = 0;
    if (part == 0) {
      output.write_header(input_header, opt.no_header);
      output_context.emplace(output.context(
          /*include_unmapped=*/opt.format == "paf" ? opt.paf_no_hit
                                                   : !opt.sam_hit_only,
          opt.paf_cigar,
          opt.copy_comment, emit_opts, own_sam_tag_keys(opt)));
    }
    batch_reader = std::make_unique<fa::cpu::io::BatchReader>(reader, batch_bp);

    if (batch_window >= 2)
      session = dna_aligner->make_windowed_session(batch_window, decode_loan);
    // From here the decode team's idle threads map reads.
    if (session != nullptr && decode_loan > 0) {
      fa::cpu::api::WindowedAlignSession* target = session.get();
      helper_hub.install([target](int tid) { return target->try_help(tid); });
    }
    // Revokes the loan before the window is torn down, on unwind too.
    struct HelperLoanGuard {
      HelperHub& hub;
      ~HelperLoanGuard() { hub.revoke(); }
    } helper_loan_guard{helper_hub};

    // Release the persistent parallel_for team left over from unpacking the
    // index. The windowed path never uses it; the barrier path respawns it.
    fa::cpu::threading::join_persistent_pool();

    // Reader, compute and writer run concurrently with two batches per
    // hand-off, plus up to the window in the compute stage. Output order is
    // independent of batching, threads and window. -t 1 runs all three stages
    // on this thread.
    if (batch_window >= 2) {
      fa::cpu::threading::run_pipeline<ReadBatch, AlignedBatch>(
          produce, process_windowed, consume, /*depth=*/2, flush_windowed);
    } else if (num_threads == 1) {
      fa::cpu::threading::run_pipeline_serial<ReadBatch, AlignedBatch>(
          produce, process_barrier, consume);
    } else {
      fa::cpu::threading::run_pipeline<ReadBatch, AlignedBatch>(
          produce, process_barrier, consume, /*depth=*/2);
    }
    // Revoke the loan before the window goes.
    helper_hub.revoke();
    session.reset();
    batch_reader.reset();
  }
  output.close();
  // Zero reads usually means an unrecognized input format.
  if (seen == 0) {
    std::cerr << "flashalign: warning: 0 reads read from input"
                 " (empty file, or an unrecognized format?)\n";
  }
  if (opt.progress && !opt.quiet) {
    std::cerr << "[flashalign] done: aligned " << seen << " reads, wrote "
              << summary.records << " records\n";
  }

  if (opt.stats) {
    const auto done = std::chrono::steady_clock::now();
    summary.preset = opt.preset;
    summary.output_format = opt.format;
    summary.rna = backend_mode == fa::cpu::api::BackendMode::RnaSplice;
    summary.reads = seen.load();
    summary.wall_seconds =
        std::chrono::duration<double>(done - program_t0).count();
    summary.map_seconds =
        std::chrono::duration<double>(done - align_t0).count();
    summary.peak_rss_bytes = peak_rss_bytes();
    write_stats(summary, std::cerr);
  }
  return 0;
}

} // namespace fa::cpu::cli
