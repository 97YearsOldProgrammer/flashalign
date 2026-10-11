#include "output_writer.h"

#include "../io/bam_reader.h"
#include "../io/paf.h"
#include "read_group.h"
#include "fa_version.h"

#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

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

}  // namespace fa::cpu::cli
