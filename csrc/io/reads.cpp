// Query records from FASTA, FASTQ or BAM paths, in input order: the byte-source chain, kseq, the BAM dispatch, bp-budgeted batches.
#include "io.h"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <set>
#include <stdexcept>

/**** byte_source ****/

namespace fa { namespace cpu { namespace io {

void warn_truncated_gzip(const std::string& path) {
    // A file may be read more than once (flashalign index reads the reference once per
    // part), so the warning is kept to one per path.
    static std::mutex mu;
    static std::set<std::string> warned;
    {
        std::lock_guard<std::mutex> lock(mu);
        if (!warned.insert(path).second) return;
    }
    std::fprintf(stderr,
                 "flashalign: warning: %s: truncated gzip file (it ends before the gzip "
                 "trailer); the sequences after that point are missing\n",
                 path == "-" ? "stdin" : path.c_str());
}

ZlibByteSource::ZlibByteSource(const std::string& path, InputPolicy policy)
    : path_(path), policy_(policy), in_(new unsigned char[kInputBytes]),
      out_(new unsigned char[kOutputBytes]) {
    next_ = in_.get();
    if (path == "-") {
        fd_ = STDIN_FILENO;
        return;
    }
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) throw std::runtime_error("failed to open FASTA/FASTQ: " + path);
    owns_fd_ = true;
#ifdef POSIX_FADV_SEQUENTIAL
    (void)posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
}

ZlibByteSource::ZlibByteSource(const std::string& path, uint64_t start_offset,
                               bool after_member, InputPolicy policy, bool last_bgzf_empty)
    : ZlibByteSource(path, policy) {
    after_member_ = after_member;
    status_.compressed_offset = start_offset;
    status_.last_bgzf_empty = last_bgzf_empty;
    if (start_offset == 0) return;
    if (!owns_fd_ || ::lseek(fd_, static_cast<off_t>(start_offset), SEEK_SET) < 0)
        throw std::runtime_error("failed to seek in FASTA/FASTQ: " + path);
    bytes_read_ = start_offset;
}

ZlibByteSource::~ZlibByteSource() {
    if (zs_ready_) inflateEnd(&zs_);
    if (owns_fd_) ::close(fd_);
}

// Appends the next read(2) to the unread input. False on a read error.
bool ZlibByteSource::fill() {
    if (next_ != in_.get()) {
        std::memmove(in_.get(), next_, avail_);
        next_ = in_.get();
    }
    for (;;) {
        const ssize_t n = ::read(fd_, in_.get() + avail_, kInputBytes - avail_);
        if (n < 0) {
            if (errno == EINTR) continue;
            return fail(std::string("input read failed: ") + std::strerror(errno), false);
        }
        if (n == 0) eof_ = true;
        avail_ += static_cast<size_t>(n);
        bytes_read_ += static_cast<uint64_t>(n);
        return true;
    }
}

int ZlibByteSource::deliver(size_t n) {
    delivered_ += n;
    if (owns_fd_ && delivered_ >= next_advice_) advise();
    return static_cast<int>(n);
}

int ZlibByteSource::read(void* buf, int len) {
    if (len <= 0) return 0;
    const size_t want = static_cast<size_t>(len);
    for (;;) {
        if (out_pos_ < out_len_) {
            const size_t n = std::min(want, out_len_ - out_pos_);
            std::memcpy(buf, out_.get() + out_pos_, n);
            out_pos_ += n;
            return deliver(n);
        }
        switch (mode_) {
        case Mode::kLook:
            // As gz_look(): a member starts only where the next two bytes are 1f 8b.
            while (avail_ < 2 && !eof_) {
                if (!fill()) {
                    mode_ = Mode::kError;
                    return -1;
                }
            }
            if (avail_ >= 2 && next_[0] == 0x1f && next_[1] == 0x8b) {
                member_start_ = bytes_read_ - avail_;
                status_.compressed_offset = member_start_;
                if (policy_ == InputPolicy::bam_integrity && !prepare_member()) {
                    mode_ = Mode::kError;
                    return -1;
                }
                const int ret = zs_ready_ ? inflateReset(&zs_) : inflateInit2(&zs_, 15 + 16);
                if (ret != Z_OK) {
                    fail("cannot initialize gzip decoder");
                    mode_ = Mode::kError;
                    return -1;
                }
                zs_ready_ = true;
                final_block_ended_ = false;
                mode_ = Mode::kGzip;
            } else {
                // BAM must consume the complete stream, including after an empty block.
                if (policy_ == InputPolicy::bam_integrity && avail_ > 0 &&
                    (zs_ready_ || after_member_)) {
                    member_start_ = bytes_read_ - avail_;
                    fail(avail_ == 1 ? "incomplete next gzip header" :
                         "invalid next gzip header or trailing bytes");
                    mode_ = Mode::kError;
                    return -1;
                }
                status_.compressed_offset = bytes_read_ - avail_;
                // Plain input, or data after the last member, which gzread ignores.
                mode_ = zs_ready_ || after_member_ ? Mode::kEnd : Mode::kCopy;
            }
            continue;
        case Mode::kCopy: {
            if (avail_ == 0 && !eof_ && !fill()) {
                mode_ = Mode::kError;
                return -1;
            }
            if (avail_ == 0) {
                mode_ = Mode::kEnd;
                continue;
            }
            const size_t n = std::min(want, avail_);
            std::memcpy(buf, next_, n);
            next_ += n;
            avail_ -= n;
            return deliver(n);
        }
        case Mode::kGzip:
            if (!inflate_some()) {
                mode_ = Mode::kError;
                return -1;
            }
            continue;
        case Mode::kEnd:
            return 0;
        case Mode::kError:
            return -1;
        }
    }
}

// Inflates into out_ until it is full, the member ends (mode_ becomes kLook) or
// the input ends inside it (kEnd). False on corrupt data or a read error.
bool ZlibByteSource::inflate_some() {
    zs_.next_out = out_.get();
    zs_.avail_out = static_cast<uInt>(kOutputBytes);
    int stalled = 0;
    while (zs_.avail_out > 0) {
        // The last kTailBytes of input wait until the end of input is known.
        if (avail_ <= kTailBytes && !eof_) {
            if (!fill()) return false;
            continue;
        }
        if (avail_ == 0) {
            // gzread's "unexpected end of file": what was decoded is delivered,
            // then EOF.
            cut_in_data_ = !final_block_ended_;
            if (policy_ == InputPolicy::bam_integrity)
                return fail(final_block_ended_ ? "incomplete gzip trailer" :
                            "incomplete gzip member");
            if (policy_ != InputPolicy::gzip_probe) warn_truncated_gzip(path_);
            mode_ = Mode::kEnd;
            break;
        }
        // A member cut inside its trailer ends its final block within the last
        // few bytes of input. Those go through Z_BLOCK, which returns at each
        // block end, so that end is seen (data_type bits 64 and 128 together);
        // the rest goes through Z_NO_FLUSH, which is faster.
        const size_t feed = eof_ ? avail_ : avail_ - kTailBytes;
        zs_.next_in = next_;
        zs_.avail_in = static_cast<uInt>(feed);
        const uInt out_before = zs_.avail_out;
        const int ret = inflate(&zs_, eof_ ? Z_BLOCK : Z_NO_FLUSH);
        const size_t used = feed - zs_.avail_in;
        next_ += used;
        avail_ -= used;
        if (eof_ && (zs_.data_type & 192) == 192) final_block_ended_ = true;
        if (policy_ == InputPolicy::bam_integrity && member_bgzf_bytes_ != 0 &&
            (bytes_read_ - avail_ - member_start_ > member_bgzf_bytes_ ||
             zs_.total_out > kBgzfMaxBytes))
            return fail("BGZF BSIZE or uncompressed block size mismatch");
        if (ret == Z_STREAM_END) {
            if (policy_ == InputPolicy::bam_integrity) {
                if (member_bgzf_bytes_ != 0 &&
                    bytes_read_ - avail_ - member_start_ != member_bgzf_bytes_)
                    return fail("BGZF BSIZE does not match the complete gzip member");
                status_.last_bgzf_empty = member_bgzf_bytes_ != 0 && zs_.total_out == 0;
                status_.compressed_offset = bytes_read_ - avail_;
            }
            mode_ = Mode::kLook;
            break;
        }
        if (ret != Z_OK && ret != Z_BUF_ERROR) {
            if (zs_.msg && std::strcmp(zs_.msg, "incorrect data check") == 0)
                return fail("gzip CRC32 mismatch");
            if (zs_.msg && std::strcmp(zs_.msg, "incorrect length check") == 0)
                return fail("gzip ISIZE mismatch");
            return fail(std::string("corrupt gzip member: ") +
                        (zs_.msg ? zs_.msg : "inflate failed"));
        }
        // A Z_BLOCK return at a block end may consume and produce nothing
        // (empty blocks decode from bits zlib already holds); anything else
        // that makes no progress twice is a zlib fault.
        const bool progress =
            used > 0 || zs_.avail_out != out_before || (zs_.data_type & 128) != 0;
        stalled = progress ? 0 : stalled + 1;
        if (stalled >= 2) return fail("gzip decoder made no progress");
    }
    out_pos_ = 0;
    out_len_ = kOutputBytes - zs_.avail_out;
    return true;
}

bool ZlibByteSource::fail(const std::string& reason, bool at_member_start) {
    status_.compressed_offset = at_member_start ? member_start_ : bytes_read_ - avail_;
    status_.at_member_start = at_member_start;
    status_.error = reason;
    return false;
}

bool ZlibByteSource::prepare_member() {
    member_bgzf_bytes_ = 0;
    for (;;) {
        const BgzfHeader header = parse_bgzf_header(next_, avail_);
        if (header.kind == BgzfHeaderKind::need_more) {
            if (eof_) return fail("incomplete gzip header");
            if (!fill()) return false;
            continue;
        }
        if (header.kind == BgzfHeaderKind::malformed) return fail(header.error);
        if (header.kind == BgzfHeaderKind::bgzf)
            member_bgzf_bytes_ = header.block_bytes;
        return true;
    }
}

void ZlibByteSource::advise() {
    next_advice_ = delivered_ + kAdviceEvery;
#if defined(POSIX_FADV_WILLNEED) && defined(POSIX_FADV_DONTNEED)
    const uint64_t pos = bytes_read_ - avail_;
    (void)posix_fadvise(fd_, static_cast<off_t>(pos), static_cast<off_t>(kReadahead),
                        POSIX_FADV_WILLNEED);
    if (pos > kDropLag) {
        const uint64_t drop_end = pos - kDropLag;
        if (drop_end >= dropped_to_ + kDropStep) {
            (void)posix_fadvise(fd_, static_cast<off_t>(dropped_to_),
                                static_cast<off_t>(drop_end - dropped_to_),
                                POSIX_FADV_DONTNEED);
            dropped_to_ = drop_end;
        }
    }
#endif
}

}}}  // namespace fa::cpu::io

/**** fastx ****/

namespace fa { namespace cpu { namespace io {

// The kseq.h functions io.h declares with KSEQ_DECLARE.
__KS_BASIC(ByteSource*, 16384)
__KS_GETC(bytesource_read, 16384)
__KS_GETUNTIL(bytesource_read, 16384)
__KSEQ_BASIC(, ByteSource*)
__KSEQ_READ()

int64_t read_whole_record(kseq_t* ks, const ByteSource& src) {
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

void uppercase_inplace(std::string& seq) {
    // Soft-masked (lowercase) bases compare equal to uppercase ones.
    for (char& ch : seq) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
}

std::unique_ptr<ByteSource>
make_byte_source(const std::string& path, const IoPlan* plan,
                 InputPolicy policy) {
  // Generic-gzip BAM is supported serially, where every member is validated.
  if (policy == InputPolicy::gzip_compatible &&
      (plan == nullptr || plan->parallel_gz) &&
      gz_member_parallel_eligible(path)) {
    return std::make_unique<MemberParallelGzByteSource>(
        path, plan == nullptr ? GzMemberSourceOptions() : plan->gz);
  }
  if ((plan == nullptr || plan->parallel_bgzf) &&
      bgzf_eligible(path)) {
    return std::make_unique<BgzfByteSource>(
        path, plan == nullptr ? GzMemberSourceOptions() : plan->gz, policy);
  }
  return std::make_unique<ZlibByteSource>(path, policy);
}

FastxReader::FastxReader(std::string path, const IoPlan* plan,
                         bool keep_comment)
    : FastxReader(std::vector<std::string>{std::move(path)}, plan,
                  keep_comment) {}

FastxReader::FastxReader(std::vector<std::string> paths,
                         const IoPlan* plan,
                         bool keep_comment)
    : paths_(std::move(paths)), plan_(plan), keep_comment_(keep_comment)
{
    if (paths_.empty()) {
        throw std::runtime_error("no FASTA/FASTQ input path given");
    }
    open_current();
}

FastxReader::~FastxReader() { close_current(); }

bool FastxReader::next(FastxRecord& rec) {
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
        // len == -1: end of this file. After the last one the source is
        // kept until the reader is destroyed, so a decode team lent by
        // the caller keeps running the caller's work until then.
        if (idx_ + 1 >= paths_.size()) {
            kseq_destroy(ks_);
            ks_ = nullptr;
            return false;
        }
        close_current();
        ++idx_;
        open_current();
    }
}

const std::string& FastxReader::current_path() const { return paths_[idx_]; }

void FastxReader::open_current() {
    src_ = make_byte_source(paths_[idx_], plan_);
    ks_ = kseq_init(src_.get());
}

void FastxReader::close_current() {
    if (ks_) { kseq_destroy(ks_); ks_ = nullptr; }
    src_.reset();  // kseq does not own the ByteSource
}

GenomeLoad load_fasta_genome(const std::string& path) {
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

/**** reads ****/

namespace fa { namespace cpu { namespace io {

SequenceReader::SequenceReader(std::string path, const IoPlan* plan,
                               bool keep_comment)
    : SequenceReader(std::vector<std::string>{std::move(path)}, plan,
                     keep_comment) {}

SequenceReader::SequenceReader(std::vector<std::string> paths,
                               const IoPlan* plan,
                               bool keep_comment)
    : paths_(std::move(paths)), plan_(plan), keep_comment_(keep_comment)
{
    if (paths_.empty()) throw std::runtime_error("no input path given");
    open_current();
}

bool SequenceReader::next(FastxRecord& rec) {
    for (;;) {
        const bool got = bam_ ? bam_->next(rec) : fx_->next(rec);
        if (got) return true;
        ++idx_;
        if (idx_ >= paths_.size()) return false;
        open_current();
    }
}

void SequenceReader::open_current() {
    fx_.reset();
    bam_.reset();
    const std::string& p = paths_[idx_];
    if (path_is_cram(p)) {
        throw std::runtime_error(
            "CRAM is not supported natively (it needs the reference/refget "
            "machinery only samtools/htslib provides); decode it first, e.g.: "
            "samtools fastq " + p + " | flashalign align ref.fa -");
    }
    if (path_is_bam(p)) bam_ = std::make_unique<BamReader>(p, plan_);
    else fx_ = std::make_unique<FastxReader>(p, plan_, keep_comment_);
}

}}}  // namespace fa::cpu::io

/**** batch_reader ****/

namespace fa { namespace cpu { namespace io {

BatchReader::BatchReader(SequenceReader& reader, int64_t batch_bp)
    : reader_(reader), batch_bp_(batch_bp) {}

std::optional<std::vector<FastxRecord>>
BatchReader::next(std::atomic<int64_t>& seen) {
    if (done_) return std::nullopt;
    std::vector<FastxRecord> batch;
    int64_t batch_bp_seen = 0;
    FastxRecord rec;
    for (;;) {
        if (!reader_.next(rec)) { done_ = true; break; }
        batch_bp_seen += static_cast<int64_t>(rec.seq.size());
        batch.push_back(std::move(rec));
        ++seen;
        if (batch_bp_seen >= batch_bp_) break;
    }
    if (batch.empty()) return std::nullopt;
    return batch;
}

}}}  // namespace fa::cpu::io
