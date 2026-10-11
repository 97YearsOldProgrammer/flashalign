// The decompressed-byte source kseq.h reads from, so the decompressor can be
// swapped without modifying kseq. The default policy delivers the bytes zlib's
// gzread would, in order; BAM policy instead fails on confirmed damage. A source
// is read from one thread only. The
// implementations are ZlibByteSource below, MemberParallelGzByteSource
// (io/gz_member_source.h) and BgzfByteSource
// (io/bgzf_source.h); io/fastx.h make_byte_source() chooses. A gzip file
// cut before a member's trailer reads as ending at the cut, as with gzread;
// gzip-compatible sources call warn_truncated_gzip() when they reach the cut.
#pragma once

#include <zlib.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace fa { namespace cpu { namespace io {

// A format probe uses gzip-compatible bytes but leaves warnings to the reader.
enum class InputPolicy { gzip_compatible, gzip_probe, bam_integrity };

struct ByteSourceStatus {
    // A known compressed position; a source error uses its member's start.
    uint64_t compressed_offset = ~uint64_t{0};
    bool last_bgzf_empty = false;
    bool at_member_start = false;
    std::string error;
};

// read() returns the byte count, 0 at EOF, or < 0 on error, as gzread does.
struct ByteSource {
    virtual ~ByteSource() = default;
    virtual int read(void* buf, int len) = 0;
    // Whether this source decodes with a thread team. The bytes are the same
    // either way; a caller that lent threads needs to know if they are used.
    virtual bool is_parallel() const { return false; }
    // Whether the input is a gzip file cut inside its compressed data, rather
    // than only inside the last member's trailer, so the record the cut
    // interrupts is incomplete. Set by the time read() returns 0.
    virtual bool cut_in_data() const { return false; }
    virtual ByteSourceStatus status() const { return {}; }
};

// The read function passed to KSEQ_INIT.
inline int bytesource_read(ByteSource* s, void* buf, int len) {
    return s->read(buf, len);
}

// Warns on stderr that the gzip input `path` ("-" for stdin) ends before a
// member's trailer, so the sequences after the cut are missing. Printed once per
// path per process.
void warn_truncated_gzip(const std::string& path);

// Plain, gzip and bgzf input, including "-" for stdin, decoded with zlib's
// inflate as gzread decodes it: a member follows another while the next two
// bytes are the gzip magic, anything else after a member is ignored, and input
// that does not start with the magic is passed through as is.
class ZlibByteSource : public ByteSource {
public:
    // Throws std::runtime_error if `path` cannot be opened.
    explicit ZlibByteSource(const std::string& path,
                            InputPolicy policy = InputPolicy::gzip_compatible);
    // Reads the file `path` from byte `start_offset`. With `after_member` the
    // bytes there follow a gzip member, so anything but the gzip magic ends
    // the input under gzip-compatible policy, rather than being passed through.
    // BAM policy refuses such trailing bytes and retains the last marker state.
    ZlibByteSource(const std::string& path, uint64_t start_offset,
                   bool after_member,
                   InputPolicy policy = InputPolicy::gzip_compatible,
                   bool last_bgzf_empty = false);
    ~ZlibByteSource() override;
    ZlibByteSource(const ZlibByteSource&) = delete;
    ZlibByteSource& operator=(const ZlibByteSource&) = delete;

    int read(void* buf, int len) override;
    bool cut_in_data() const override { return cut_in_data_; }
    ByteSourceStatus status() const override { return status_; }

private:
    enum class Mode { kLook, kCopy, kGzip, kEnd, kError };

    bool fill();
    bool inflate_some();
    bool prepare_member();
    bool fail(const std::string& reason, bool at_member_start = true);
    int deliver(size_t n);
    void advise();

    // Compressed bytes per read(2), and decompressed bytes per refill of out_,
    // gzread's sizes under a 1 MiB gzbuffer().
    static constexpr size_t kInputBytes = 1u << 20;
    static constexpr size_t kOutputBytes = 2u << 20;
    // Input held back from inflate until the end of input is known: more than
    // a trailer plus the end-of-block code before it.
    static constexpr size_t kTailBytes = 32;
    // Every kAdviceEvery output bytes, request readahead past the compressed
    // position and drop the pages well behind it. The input is read once, and
    // without the drop a large FASTQ would evict the memory-mapped index from
    // the page cache.
    static constexpr uint64_t kAdviceEvery = 32ull << 20;
    static constexpr uint64_t kReadahead = 32ull << 20;
    static constexpr uint64_t kDropLag = 64ull << 20;
    static constexpr uint64_t kDropStep = 64ull << 20;

    std::string path_;
    InputPolicy policy_;
    ByteSourceStatus status_;
    uint64_t member_start_ = 0;
    size_t member_bgzf_bytes_ = 0;
    int fd_ = -1;
    bool owns_fd_ = false; // false for stdin
    std::unique_ptr<unsigned char[]> in_;
    unsigned char* next_ = nullptr; // unread input in in_
    size_t avail_ = 0;
    bool eof_ = false;
    std::unique_ptr<unsigned char[]> out_;
    size_t out_pos_ = 0; // out_[out_pos_, out_len_) is decoded, not yet read
    size_t out_len_ = 0;
    uint64_t bytes_read_ = 0; // compressed bytes read from fd_
    Mode mode_ = Mode::kLook;
    z_stream zs_{};
    bool zs_ready_ = false;
    bool after_member_ = false;
    // The current member's final deflate block has ended; only its trailer is
    // left. Tracked in the last kTailBytes of input only.
    bool final_block_ended_ = false;
    bool cut_in_data_ = false;
    uint64_t delivered_ = 0;
    uint64_t next_advice_ = kAdviceEvery;
    uint64_t dropped_to_ = 0;
};

}}}  // namespace fa::cpu::io
