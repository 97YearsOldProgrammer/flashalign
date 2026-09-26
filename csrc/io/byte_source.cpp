#include "byte_source.h"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <set>
#include <stdexcept>

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

ZlibByteSource::ZlibByteSource(const std::string& path)
    : path_(path), in_(new unsigned char[kInputBytes]),
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
            return false;
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
                const int ret = zs_ready_ ? inflateReset(&zs_) : inflateInit2(&zs_, 15 + 16);
                if (ret != Z_OK) {
                    mode_ = Mode::kError;
                    return -1;
                }
                zs_ready_ = true;
                final_block_ended_ = false;
                mode_ = Mode::kGzip;
            } else {
                // Plain input, or data after the last member, which gzread ignores.
                mode_ = zs_ready_ ? Mode::kEnd : Mode::kCopy;
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
            warn_truncated_gzip(path_);
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
        if (ret == Z_STREAM_END) {
            mode_ = Mode::kLook;
            break;
        }
        if (ret != Z_OK && ret != Z_BUF_ERROR) return false;
        // A Z_BLOCK return at a block end may consume and produce nothing
        // (empty blocks decode from bits zlib already holds); anything else
        // that makes no progress twice is a zlib fault.
        const bool progress =
            used > 0 || zs_.avail_out != out_before || (zs_.data_type & 128) != 0;
        stalled = progress ? 0 : stalled + 1;
        if (stalled >= 2) return false;
    }
    out_pos_ = 0;
    out_len_ = kOutputBytes - zs_.avail_out;
    return true;
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
