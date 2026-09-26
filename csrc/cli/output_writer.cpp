#include "output_writer.h"

#include "../io/bam.h"
#include "../io/bam_reader.h"
#include "../io/bgzf.h"
#include "../io/paf.h"
#include "../io/sam_tags.h"
#include "read_group.h"
#include "fa_version.h"

#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
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

std::vector<std::string> passthrough_header_lines(
    const std::string& source
) {
    std::vector<std::string> lines;
    std::size_t begin = 0;
    while (begin < source.size()) {
        std::size_t end = source.find('\n', begin);
        if (end == std::string::npos) end = source.size();
        const std::string line = source.substr(begin, end - begin);
        if (line.rfind("@RG\t", 0) == 0 ||
            line.rfind("@CO\t", 0) == 0) {
            lines.push_back(line);
        }
        begin = end + 1;
    }
    return lines;
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

// -y for BAM output. SAM and PAF copy the comment verbatim, as minimap2 does,
// but BAM aux fields are typed, so a comment that is not SAM tag text has no
// encoding and stops the run.
void require_sam_tag_text(const io::FastxRecord& read) {
    std::string offending;
    if (io::sam_tag_text_valid(read.comment, &offending)) return;
    throw std::runtime_error("-y: the FASTA/Q comment of read " + read.name +
                             " is not SAM tag text (" + offending +
                             "); BAM aux fields are typed");
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
    const std::string& source_header_text,
    bool suppress_sam_header
) {
    const std::vector<std::string> source_lines =
        passthrough_header_lines(source_header_text);
    require_no_read_group_conflict(source_lines, read_group_id_);
    // The declared @RG goes after @PG and before the passthrough @RG/@CO lines.
    std::vector<std::string> header_lines =
        output::sam_header(references_, "flashalign", FA_VERSION,
                           command_line_);
    if (!read_group_line_.empty()) header_lines.push_back(read_group_line_);
    for (const std::string& line : source_lines) header_lines.push_back(line);
    if (format_ == "sam" && !suppress_sam_header) {
        for (const auto& line : header_lines) *text_ << line << '\n';
        require_output(*text_, "write");
    }
    if (format_ != "bam") return;

    std::string header_text;
    for (const auto& line : header_lines) {
        header_text += line;
        header_text.push_back('\n');
    }
    bam_writer_ = std::make_unique<io::BgzfWriter>(*text_);
    bam_writer_->write(output::encode_bam_header(references_, header_text));
    for (std::size_t index = 0; index < references_.size(); ++index) {
        bam_reference_ids_.emplace(
            references_[index].name, static_cast<int>(index));
    }
}

AlignmentOutputContext AlignmentOutputWriter::context(
    bool include_unmapped,
    bool paf_cigar,
    bool copy_comment,
    output::SamEmitOptions sam
) {
    sam.read_group_id = read_group_id_;
    // emit_secondary is already set by the caller from --secondary.
    return {
        *text_,
        bam_writer_.get(),
        &bam_reference_ids_,
        &reference_lengths_,
        format_,
        include_unmapped,
        paf_cigar,
        copy_comment,
        sam,
    };
}

void AlignmentOutputWriter::close() {
    if (bam_writer_) bam_writer_->close();
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
            const std::string passthrough_tags =
                io::bam_tags_to_sam_text(read.tag_bytes);
            // -y: the comment goes on every record of the read, verbatim and
            // unchecked as in minimap2; the passthrough tags stay primary-only.
            std::string record_tags;
            if (context.copy_comment && !read.comment.empty())
              record_tags = read.comment;
            // Sized for the primary record; other records may still grow it.
            line.reserve(read.seq.size() + read.qual.size() +
                         result.cigar.size() + read.name.size() +
                         passthrough_tags.size() + record_tags.size() + 256);
            output::SamEmitOptions emit_options = context.sam;
            // A read that already carries an RG:Z keeps it.
            if (sam_tags_declare_read_group(passthrough_tags) ||
                sam_tags_declare_read_group(record_tags))
                emit_options.read_group_id.clear();
            output::append_sam_records(
                line, read.name, result, read.seq, read.qual,
                passthrough_tags, emit_options, record_tags);
            context.text << line;
            written += std::count(line.begin(), line.end(), '\n');
        }
    } else if (context.format == "bam") {
        if (result.mapped() || context.include_unmapped) {
            output::SamEmitOptions emit_options = context.sam;
            // -y: the comment encoded as aux bytes on every record of the read.
            std::string record_extra_tags;
            bool comment_declares_read_group = false;
            if (context.copy_comment && !read.comment.empty()) {
                require_sam_tag_text(read);
                record_extra_tags = io::sam_tag_text_to_bam(read.comment);
                comment_declares_read_group =
                    sam_tags_declare_read_group(read.comment);
            }
            if (!emit_options.read_group_id.empty() &&
                (comment_declares_read_group ||
                 sam_tags_declare_read_group(
                     io::bam_tags_to_sam_text(read.tag_bytes))))
                emit_options.read_group_id.clear();
            const std::vector<std::string> bam_records =
                output::encode_bam_records(
                    read.name, result, read.seq, read.qual,
                    *context.bam_reference_ids, read.tag_bytes, emit_options,
                    record_extra_tags);
            for (const std::string& record : bam_records) {
                context.bam_writer->write(record);
                ++written;
            }
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
            // in the same order as the SAM and BAM writers.
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
        }
    }
    if (context.format != "bam") require_output(context.text, "write");
    return written;
}

}  // namespace fa::cpu::cli
