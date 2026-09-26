// Reads query records from one or more paths in turn. A BAM/uBAM input goes to
// BamReader, anything else (FASTA/FASTQ, plain or gzip, or "-" for stdin) to
// FastxReader. uBAM tags (MM/ML/RG...) are kept in FastxRecord::tag_bytes.
#pragma once

#include "bam_reader.h"
#include "fastx.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace io {

class SequenceReader {
public:
    // `plan` is borrowed and must outlive the reader; null takes the defaults.
    // BAM inputs ignore it. `keep_comment` (-y) applies to FASTA/FASTQ only.
    explicit SequenceReader(std::string path, const IoPlan* plan = nullptr,
                            bool keep_comment = false)
        : SequenceReader(std::vector<std::string>{std::move(path)}, plan,
                         keep_comment) {}

    explicit SequenceReader(std::vector<std::string> paths,
                            const IoPlan* plan = nullptr,
                            bool keep_comment = false)
        : paths_(std::move(paths)), plan_(plan), keep_comment_(keep_comment)
    {
        if (paths_.empty()) throw std::runtime_error("no input path given");
        open_current();
    }

    bool next(FastxRecord& rec) {
        for (;;) {
            const bool got = bam_ ? bam_->next(rec) : fx_->next(rec);
            if (got) return true;
            ++idx_;
            if (idx_ >= paths_.size()) return false;
            open_current();
        }
    }

    // SAM header text of the current input if it is BAM, else empty.
    std::string source_header_text() const {
        return bam_ ? bam_->sam_header_text() : std::string();
    }

    // Whether the current input is decoded by a thread team (never for BAM).
    bool decode_parallel_active() const {
        return fx_ && fx_->decode_parallel_active();
    }

private:
    void open_current() {
        fx_.reset();
        bam_.reset();
        const std::string& p = paths_[idx_];
        if (path_is_cram(p)) {
            throw std::runtime_error(
                "CRAM is not supported natively (it needs the reference/refget "
                "machinery only samtools/htslib provides); decode it first, e.g.: "
                "samtools fastq " + p + " | flashalign align ref.fa -");
        }
        if (path_is_bam(p)) bam_ = std::make_unique<BamReader>(p);
        else fx_ = std::make_unique<FastxReader>(p, plan_, keep_comment_);
    }

    std::vector<std::string> paths_;
    size_t idx_ = 0;
    const IoPlan* plan_ = nullptr;
    bool keep_comment_ = false;
    std::unique_ptr<FastxReader> fx_;
    std::unique_ptr<BamReader> bam_;
};

}}}  // namespace fa::cpu::io
