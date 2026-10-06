// Streaming FASTA/FASTQ reader over kseq.h, reading from a ByteSource
// (io/byte_source.h): plain, gzip or bgzf input, or "-" for stdin. Several
// paths are read in sequence.
#pragma once

#include "byte_source.h"
#include "gz_member_source.h"
#include "vendor/kseq.h"

#include <cctype>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace io {

KSEQ_INIT(ByteSource*, bytesource_read)

// kseq_read(), except that when `src` is a gzip file cut inside its data the
// record the cut interrupts is dropped and the input ends before it, so only
// whole records are read. A record is whole when it stopped before the end of
// input (the next header ended it), or when it is FASTQ with its full quality
// string. When the data is intact, parse errors are returned as kseq_read()
// returns them.
inline int64_t read_whole_record(kseq_t* ks, const ByteSource& src) {
    const int64_t len = kseq_read(ks);
    if (len == -1 || !src.cut_in_data()) return len;
    const bool at_end = ks_eof(ks->f);
    if (len >= 0) {
        const bool whole = !at_end || (ks->is_fastq && ks->qual.l == ks->seq.l);
        return whole ? len : -1;
    }
    // A truncated FASTQ record (-2) is the cut only when it ran into the end of input.
    return len == -2 && at_end ? -1 : len;
}

struct FastxRecord {
    std::string name;
    std::string seq;
    std::string qual;
    // Raw BAM tag block from uBAM input (e.g. MM/ML/RG), passed through to
    // the output; empty for FASTA/FASTQ.
    std::string tag_bytes;
    // Header text after the first space; filled only for -y.
    std::string comment;
};

struct GenomeLoad {
    std::unordered_map<std::string, std::string> sequences;
    std::vector<std::pair<std::string, int64_t>> references;
};

inline void uppercase_inplace(std::string& seq) {
    // Soft-masked (lowercase) bases compare equal to uppercase ones.
    for (char& ch : seq) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
}

// Opens one input. A large gzip file that is not BGZF gets the parallel
// source, which returns the same bytes; everything else, and any failure to
// start it, gets ZlibByteSource. `plan` (or null for defaults) sets threads and
// memory only.
inline std::unique_ptr<ByteSource>
make_byte_source(const std::string& path, const IoPlan* plan = nullptr) {
  if ((plan == nullptr || plan->parallel_gz) &&
      gz_member_parallel_eligible(path)) {
    try {
      return std::make_unique<MemberParallelGzByteSource>(
          path, plan == nullptr ? GzMemberSourceOptions() : plan->gz);
    } catch (const std::exception&) {
      // Fall back to the serial source, which reports any open error.
    }
  }
  return std::make_unique<ZlibByteSource>(path);
}

class FastxReader {
public:
    // `plan` is borrowed and must outlive the reader. `keep_comment` fills
    // FastxRecord::comment.
    explicit FastxReader(std::string path, const IoPlan* plan = nullptr,
                         bool keep_comment = false)
        : FastxReader(std::vector<std::string>{std::move(path)}, plan,
                      keep_comment) {}

    explicit FastxReader(std::vector<std::string> paths,
                         const IoPlan* plan = nullptr,
                         bool keep_comment = false)
        : paths_(std::move(paths)), plan_(plan), keep_comment_(keep_comment)
    {
        if (paths_.empty()) {
            throw std::runtime_error("no FASTA/FASTQ input path given");
        }
        open_current();
    }

    ~FastxReader() { close_current(); }

    FastxReader(const FastxReader&) = delete;
    FastxReader& operator=(const FastxReader&) = delete;

    bool next(FastxRecord& rec) {
        for (;;) {
            if (!ks_) return false;
            const int64_t len = read_whole_record(ks_, *src_);
            if (len >= 0) {
                rec.name.assign(ks_->name.s, ks_->name.l);
                rec.seq.assign(ks_->seq.s, ks_->seq.l);
                if (ks_->qual.l) rec.qual.assign(ks_->qual.s, ks_->qual.l);
                else rec.qual.clear();
                if (keep_comment_ && ks_->comment.l)
                    rec.comment.assign(ks_->comment.s, ks_->comment.l);
                else rec.comment.clear();
                uppercase_inplace(rec.seq);
                return true;
            }
            if (len == -3) {
                throw std::runtime_error("error reading FASTA/FASTQ: " + current_path());
            }
            if (len == -2) {
                throw std::runtime_error(
                    "truncated FASTQ or qual/seq length mismatch in " + current_path());
            }
            // len == -1: end of this file.
            close_current();
            ++idx_;
            if (idx_ >= paths_.size()) return false;
            open_current();
        }
    }

    // Whether the current input is decoded by a thread team; the parallel
    // source may have fallen back to serial.
    bool decode_parallel_active() const {
        return src_ && src_->is_parallel();
    }

private:
    const std::string& current_path() const { return paths_[idx_]; }

    void open_current() {
        src_ = make_byte_source(paths_[idx_], plan_);
        ks_ = kseq_init(src_.get());
    }

    void close_current() {
        if (ks_) { kseq_destroy(ks_); ks_ = nullptr; }
        src_.reset();  // kseq does not own the ByteSource
    }

    std::vector<std::string> paths_;
    size_t idx_ = 0;
    const IoPlan* plan_ = nullptr;
    bool keep_comment_ = false;
    std::unique_ptr<ByteSource> src_;
    kseq_t* ks_ = nullptr;
};

inline GenomeLoad load_fasta_genome(const std::string& path) {
    FastxReader reader(path);
    GenomeLoad out;
    FastxRecord rec;
    while (reader.next(rec)) {
        if (rec.seq.empty()) continue;
        if (out.sequences.count(rec.name))
            throw std::runtime_error(
                "duplicate sequence name in " + path + ": " + rec.name);
        out.references.push_back({rec.name, static_cast<int64_t>(rec.seq.size())});
        out.sequences.emplace(std::move(rec.name), std::move(rec.seq));
    }
    if (out.sequences.empty()) {
        throw std::runtime_error("reference FASTA has no sequences: " + path);
    }
    return out;
}

}}}  // namespace fa::cpu::io
