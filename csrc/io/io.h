// Byte sources, the readers, the IO plan, and the SAM and PAF writers' declarations.
#pragma once

#include "../core/core.h"
#include "fa_version.h"
#include "vendor/kseq.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**** byte_source ****/
// The decompressed-byte source kseq.h reads from, so the decompressor can be
// swapped without modifying kseq. The default policy delivers the bytes zlib's
// gzread would, in order; BAM policy instead fails on confirmed damage. A source
// is read from one thread only. The
// implementations are ZlibByteSource, MemberParallelGzByteSource and
// BgzfByteSource, below; make_byte_source() chooses. A gzip file
// cut before a member's trailer reads as ending at the cut, as with gzread;
// gzip-compatible sources call warn_truncated_gzip() when they reach the cut.

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

inline uint32_t read_le32(const unsigned char* p) {
  return static_cast<uint32_t>(p[0]) | static_cast<uint32_t>(p[1]) << 8 |
         static_cast<uint32_t>(p[2]) << 16 | static_cast<uint32_t>(p[3]) << 24;
}

// Reads the whole request, retrying EINTR. Returns the byte count (short only
// at end of file) or -1 on error.
inline int64_t pread_exact(int fd, void* buf, size_t count, uint64_t offset) {
  char* out = static_cast<char*>(buf);
  size_t done = 0;
  while (done < count) {
    const ssize_t n = ::pread(fd, out + done, count - done,
                              static_cast<off_t>(offset + done));
    if (n < 0) {
      if (errno == EINTR)
        continue;
      return -1;
    }
    if (n == 0)
      break;
    done += static_cast<size_t>(n);
  }
  return static_cast<int64_t>(done);
}

}}}  // namespace fa::cpu::io

/**** gz_member_source ****/
// A ByteSource that inflates a gzip file on several threads and delivers the
// bytes in order. The members of a multi-member file (e.g. `cat a.fq.gz
// b.fq.gz`) are inflated in parallel, and so are pieces of one member, cut at
// deflate block boundaries (deflate_split.cpp). A piece is inflated against a
// stand-in window and repaired once its real history is known.
//
// The output is exactly what zlib's gzread would deliver:
//   - A member is accepted only when zlib parsed its header, inflate reached
//     the end of its final block, the CRC-32 and ISIZE of its output (combined
//     across its pieces) match its trailer, checked in zlib's order, and it
//     starts where the previous accepted member ended; the chain starts at
//     offset 0. Speculative segments that fail either test are discarded.
//   - A piece's bytes are used only once the inflate that reached it from the
//     member's start stopped at exactly the bit where the piece starts, which
//     a position that is not a block boundary of the stream cannot do. A piece
//     inflated without its history is delivered only after its prefix up to
//     the last stand-in byte has been inflated again against the real history.
//   - As in zlib's gz_look(), a member is followed by another only if the next
//     two bytes are the gzip magic; anything else is ignored trailing data.
//   - A member truncated by end of file yields its decodable prefix, then EOF,
//     and warn_truncated_gzip() is called.
//   - A corrupt member yields -1 after every byte before the error. gzread may
//     deliver less before its error (how much depends on its buffering), so on
//     error gzread's output is a prefix of ours.
//
// Memory is bounded by GzMemberSourceOptions::staging_budget, charged per
// allocated chunk, plus about workers * chunk_bytes for chunks in progress, a
// recycled-buffer pool of 2 * workers * chunk_bytes, and per worker 1 MiB of
// input scratch, the block search's 1.25 MiB window and 64 KiB of inflate
// scratch. Each piece opened ahead of delivery (at most 2 * workers) adds a
// zlib state and up to two 32 KiB windows. A producer over budget suspends or
// waits for the reader.
//
// Pieces cost CPU a straight inflate does not (the stand-in pass and the
// repair), so they are opened only while the reader is short of input: it is
// waiting, or less than a quarter of the staging budget is ready for it.
// Otherwise one worker inflates the member front to back, as gzread would.
//
// With a `helper`, the decode threads are borrowed from the caller: a worker
// with nothing to inflate runs one unit of the caller's work through it.
// Decoding keeps priority, but a unit cannot be preempted, so decoding can
// stall while every worker is inside a long unit; the caller should lend
// several threads.

namespace fa {
namespace cpu {
namespace io {

// Smaller files are read serially.
inline constexpr uint64_t kGzMemberParallelMinBytes = 64ull << 20;

struct GzMemberSourceOptions {
  // Inflate threads; 0 means min(8, threading::default_thread_count()). The
  // member scan adds one thread, or with a `helper` runs on worker 0 so the
  // team is exactly `workers` (a one-thread loaned team does no scan).
  unsigned workers = 0;
  // Runs one unit of the caller's work and returns false if there is none.
  // Called with no lock held by an idle worker, with its ordinal in
  // [0, workers). Must not call back into this source. Empty: idle workers
  // sleep.
  std::function<bool(unsigned)> helper;
  // Size of one output buffer, and so the granularity of back-pressure
  // checks. Must fit a zlib uInt.
  size_t chunk_bytes = 8u << 20;
  // Ceiling on staged-but-unconsumed bytes of any single member.
  uint64_t segment_stage_cap = 512ull << 20;
  // Ceiling on staged-but-unconsumed bytes across all members. Must exceed
  // segment_stage_cap, so the member being delivered always has room.
  uint64_t staging_budget = 4ull << 30;
};

// How a run decodes its inputs. A null plan means the defaults.
struct IoPlan {
  // False forces the serial ZlibByteSource for gzip that is not BGZF.
  bool parallel_gz = true;
  // False forces it for BGZF (io/bgzf.cpp).
  bool parallel_bgzf = true;
  // The team of either parallel source.
  GzMemberSourceOptions gz;
};

// True for a regular file of at least `min_bytes` that starts with a deflate
// gzip member (1f 8b 08) and is not BGZF, which has its own source
// (io/bgzf.cpp). Never throws; a failed probe answers false.
bool gz_member_parallel_eligible(
    const std::string& path, uint64_t min_bytes = kGzMemberParallelMinBytes);

// read() may be called from one thread only. The constructor starts the
// thread team and the destructor joins it, at any point in the stream.
class MemberParallelGzByteSource : public ByteSource {
public:
  // Throws std::runtime_error if `path` cannot be opened or is not a regular
  // file.
  explicit MemberParallelGzByteSource(const std::string& path,
                                      GzMemberSourceOptions options = {});
  ~MemberParallelGzByteSource() override;

  MemberParallelGzByteSource(const MemberParallelGzByteSource&) = delete;
  MemberParallelGzByteSource&
  operator=(const MemberParallelGzByteSource&) = delete;

  // As gzread: bytes produced, 0 at EOF, or exactly -1 on error (kseq tests
  // for -1). Never throws; internal errors surface as -1.
  int read(void* buf, int len) override;

  bool is_parallel() const override { return true; }
  bool cut_in_data() const override;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace io
} // namespace cpu
} // namespace fa

/**** bgzf_source ****/
// A ByteSource that inflates a BGZF file (bgzip output, BAM) on several
// threads and delivers the bytes in order. A BGZF block states its own length
// (the 'BC' extra subfield) and its output's length and CRC-32 (the trailer),
// so the file is cut into runs of whole blocks, about 4 MiB of input each,
// without inflating anything, and every run inflates on its own into an
// output buffer of known size, one raw inflate per block.
//
// The output is exactly what ZlibByteSource delivers:
//   - A block is accepted only when its raw inflate ends exactly at its
//     trailer and its output matches the trailer's CRC-32 and ISIZE.
//   - The cut stops at the first block that is not BGZF-shaped (gzip magic,
//     deflate, FLG = FEXTRA alone, a 'BC' subfield of length 2, ISIZE at
//     most 64 KiB) or that runs past the end of the file, and delivery stops
//     at the first block that fails a check. From there the rest of the file
//     is read by ZlibByteSource(path, offset, after_member), so truncation,
//     trailing data, a gzip member that is not BGZF and a corrupt block give
//     the serial source's bytes, warning, cut_in_data() and -1.
//   - BAM integrity policy follows the same fallback with strict framing and
//     fatal truncation/trailing-data checks; it never treats damage as EOF.
//
// Memory: a ring of 2 * (workers + 1) runs, each a 32 MiB output buffer of
// which only the pages its largest run wrote are resident (a FASTQ run's
// output is about 3-5 times its input), and 4 MiB of input per thread that has
// cut a run.
//
// Demand: a worker starts a run only while the reader is waiting or fewer
// than two runs are cut and not yet consumed, the one the reader is on
// included; the ring is the hard cap. Otherwise it runs one unit of the
// caller's work through `helper`, as in MemberParallelGzByteSource. The
// reader inflates a run itself only when the run it needs is unclaimed, and
// otherwise waits for that run's filler, so with every worker busy elsewhere
// the file decodes serially on the reader's thread.

namespace fa {
namespace cpu {
namespace io {

// Smaller files are read serially.
inline constexpr uint64_t kBgzfMinBytes = 16ull << 20;

// True for a regular file of at least `min_bytes` whose first block is
// BGZF-shaped (above). Never throws; a failed probe answers false.
bool bgzf_eligible(const std::string& path,
                   uint64_t min_bytes = kBgzfMinBytes);

// read() may be called from one thread only. The constructor starts the
// thread team and the destructor joins it, at any point in the stream.
class BgzfByteSource : public ByteSource {
public:
  // `options` as for MemberParallelGzByteSource; chunk_bytes,
  // segment_stage_cap and staging_budget are unused. Throws
  // std::runtime_error if `path` cannot be opened or is not a regular file.
  explicit BgzfByteSource(const std::string& path,
                          GzMemberSourceOptions options = {},
                          InputPolicy policy = InputPolicy::gzip_compatible);
  ~BgzfByteSource() override;

  BgzfByteSource(const BgzfByteSource&) = delete;
  BgzfByteSource& operator=(const BgzfByteSource&) = delete;

  // As gzread: bytes produced, 0 at EOF, or -1 on error. Never throws.
  int read(void* buf, int len) override;

  bool is_parallel() const override { return true; }
  bool cut_in_data() const override;
  ByteSourceStatus status() const override;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace io
} // namespace cpu
} // namespace fa

/**** deflate_split ****/
// Where one deflate stream can be cut so that its pieces inflate in parallel.
// A block boundary is a place inflate can start again, except that the
// block's back-references reach up to 32 KiB into output it has not seen; the
// caller inflates such a piece against a stand-in window and repairs it once
// the real window is known. Nothing here does I/O or keeps state.
//
// find_block_start() returns a candidate, not a boundary: a position that
// parses as a dynamic block header zlib would accept and whose block zlib
// then inflates whole. It becomes a piece boundary only when the inflate that
// runs from the true start of the stream stops at exactly that bit; a
// candidate that inflate steps over is absorbed into the piece before it.

namespace fa {
namespace cpu {
namespace io {
namespace deflate_split {

// The deflate window: the furthest a back-reference can reach, and so all the
// history a piece needs from before its start.
inline constexpr size_t kWindowBytes = 32768;

// Fill byte of the stand-in window. Every output byte is a literal, a copy of
// earlier output or a copy of the stand-in, so a byte that is not this value
// is already correct, and only the output up to the last one needs inflating
// again. 0xFF does not occur in FASTA or FASTQ text; input that contains it
// still decodes exactly, with a longer repair.
inline constexpr unsigned char kSentinelByte = 0xFF;

// kWindowBytes of kSentinelByte, for inflateSetDictionary.
const unsigned char* sentinel_window();

// Starts raw inflate at bit `bit_in_byte` (0-7, LSB first) of `first_byte`:
// its remaining 8 - bit_in_byte bits go in through inflatePrime, and the
// caller hands zlib the stream from the next byte. With bit_in_byte == 0
// nothing is primed and the caller starts at `first_byte` itself. False if
// zlib refuses the prime.
bool prime_at(z_stream& zs, unsigned char first_byte, unsigned bit_in_byte);

inline constexpr uint64_t kNotFound = ~0ull;

// First candidate block start among the byte-aligned positions (aligned_only)
// or all bit positions in buf[scan_from, scan_to), skipping the byte-aligned
// positions below skip_aligned_below (an aligned pass already tried them).
// Returns a bit offset into buf, or kNotFound.
//
// A candidate is a dynamic-block header with BFINAL = 0 whose three codes pass
// zlib's inflate_table() rules, confirmed by zlib inflating the whole block
// against the stand-in window. The final block has nothing after it to run in
// parallel, and BFINAL = 0 keeps a gzip magic byte (1f, BFINAL = 1) from
// looking like a block start. A block that runs past buf[nbytes) is not
// confirmed, which costs a split point, never correctness. `out` is scratch
// for the confirming inflate.
uint64_t find_block_start(const unsigned char* buf, size_t nbytes,
                          size_t scan_from, size_t scan_to,
                          size_t skip_aligned_below, bool aligned_only,
                          std::vector<unsigned char>& out);

} // namespace deflate_split
} // namespace io
} // namespace cpu
} // namespace fa

/**** bgzf_header ****/
// Bounded gzip extra-field scan shared by BGZF probes and decoders.

namespace fa { namespace cpu { namespace io {

inline constexpr size_t kGzipExtraHeaderBytes = 12 + 65535;
inline constexpr size_t kBgzfMaxBytes = 65536;

enum class BgzfHeaderKind { need_more, gzip, bgzf, malformed };
struct BgzfHeader {
    BgzfHeaderKind kind = BgzfHeaderKind::need_more;
    size_t header_bytes = 0;
    size_t block_bytes = 0;
    const char* error = nullptr;
};

BgzfHeader parse_bgzf_header(const unsigned char* p, size_t n);

}}} // namespace fa::cpu::io

/**** fastx ****/
// Streaming FASTA/FASTQ reader over kseq.h, reading from a ByteSource
// (above): plain, gzip or bgzf input, or "-" for stdin. Several
// paths are read in sequence.

namespace fa { namespace cpu { namespace io {

KSEQ_DECLARE(ByteSource*)

// kseq_read(), except that when `src` is a gzip file cut inside its data the
// record the cut interrupts is dropped and the input ends before it, so only
// whole records are read. A record is whole when it stopped before the end of
// input (the next header ended it), or when it is FASTQ with its full quality
// string. When the data is intact, parse errors are returned as kseq_read()
// returns them.
int64_t read_whole_record(kseq_t* ks, const ByteSource& src);

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

void uppercase_inplace(std::string& seq);

// Opens one input. A large gzip file that is not BGZF gets the member-parallel
// source and a large BGZF file the block-parallel one, which return the same
// bytes; everything else, and any failure to start either, gets
// ZlibByteSource. `plan` (or null for defaults) sets threads and memory only.
std::unique_ptr<ByteSource>
make_byte_source(const std::string& path, const IoPlan* plan = nullptr,
                 InputPolicy policy = InputPolicy::gzip_compatible);

class FastxReader {
public:
    // `plan` is borrowed and must outlive the reader. `keep_comment` fills
    // FastxRecord::comment.
    explicit FastxReader(std::string path, const IoPlan* plan = nullptr,
                         bool keep_comment = false);

    explicit FastxReader(std::vector<std::string> paths,
                         const IoPlan* plan = nullptr,
                         bool keep_comment = false);

    ~FastxReader();

    FastxReader(const FastxReader&) = delete;
    FastxReader& operator=(const FastxReader&) = delete;

    bool next(FastxRecord& rec);

    // Whether the current input is decoded by a thread team; the parallel
    // source may have fallen back to serial.
    bool decode_parallel_active() const;

private:
    const std::string& current_path() const;

    void open_current();

    void close_current();

    std::vector<std::string> paths_;
    size_t idx_ = 0;
    const IoPlan* plan_ = nullptr;
    bool keep_comment_ = false;
    std::unique_ptr<ByteSource> src_;
    kseq_t* ks_ = nullptr;
};

GenomeLoad load_fasta_genome(const std::string& path);

}}}  // namespace fa::cpu::io

/**** bgzf ****/
// Flat decoded stream adapter for BAM, with source error context.

namespace fa { namespace cpu { namespace io {

class BgzfReader {
public:
    BgzfReader(std::string path, std::unique_ptr<ByteSource> src);

    // Read up to `len` bytes; returns less only at source EOF.
    size_t read(void* dst, size_t len);

    void read_exact(void* dst, size_t len, const char* field);

    bool decode_parallel_active() const;

    std::runtime_error error(const std::string& reason) const;

    void check_eof_marker();

private:
    std::string path_;
    std::unique_ptr<ByteSource> src_;
    bool checked_eof_ = false;
};

}}}  // namespace fa::cpu::io

/**** bam_reader ****/
// BAM/uBAM input. Records become FastxRecords with the raw tag block kept, so
// base-modification tags (MM/ML) and read groups carry through to the output.
// An aligned BAM reads like a uBAM: secondary and supplementary records are
// skipped, SEQ is returned in the original read orientation, and a mapped
// record loses the alignment tags of the aligner that wrote it.

namespace fa { namespace cpu { namespace io {

// True for a gzip/BGZF file whose first decompressed bytes are "BAM\1".
bool path_is_bam(const std::string& path);

// CRAM is detected only to refuse it with a clear message; otherwise the
// FASTA/FASTQ reader would silently find no records.
bool path_is_cram(const std::string& path);

class BamReader {
public:
    // The blocks are decoded by make_byte_source(path, plan); `plan` is
    // borrowed, and null takes the defaults.
    explicit BamReader(const std::string& path, const IoPlan* plan = nullptr);

    // Whether a thread team decodes the blocks.
    bool decode_parallel_active() const;

    const std::string& sam_header_text() const;

    bool next(FastxRecord& rec);

private:
    bool read_record(FastxRecord& rec, uint16_t& flag_out);

    void warn_aligned_input();

    void drop_alignment_tags(FastxRecord& rec);

    void parse_header();

    int32_t read_i32();
    void skip(size_t n, const char* what);

    std::string path_;
    std::unique_ptr<BgzfReader> bgzf_;
    std::string header_text_;
    bool warned_aligned_ = false;
};

// Whether a tag starts at blob[at]. Fewer than 3 bytes after the last tag
// hold no key and type and are ignored, as htslib's bam_aux_next ignores them;
// a tag with its key and type but a short value is not.
bool bam_tag_at(const std::string& blob, size_t at);

// The byte size of the tag at blob[at], key and type included. Throws when
// the tag runs past the block or its type is unknown.
size_t bam_tag_size(const std::string& blob, size_t at,
                    const std::string& read_name);

namespace detail {

template <typename T>
inline T load_le(const char* p) {
    T v;
    std::memcpy(&v, p, sizeof v);
    return v;
}

template <typename T>
inline void append_decimal(std::string& out, T v) {
    char buf[24];
    const auto r = std::to_chars(buf, buf + sizeof buf, v);
    out.append(buf, r.ptr);
}

}  // namespace detail

// One tag of `size` bytes (as bam_tag_size gives) as SAM text, spelled as
// htslib's sam_format_aux1 spells it: every integer type as i, f with %g,
// B:f elements with kputd.
void append_bam_tag_sam_text(std::string& out, const char* tag, size_t size);

}}}  // namespace fa::cpu::io

/**** reads ****/
// Reads query records from one or more paths in turn. A BAM/uBAM input goes to
// BamReader, anything else (FASTA/FASTQ, plain or gzip, or "-" for stdin) to
// FastxReader. uBAM tags (MM/ML/RG...) are kept in FastxRecord::tag_bytes.

namespace fa { namespace cpu { namespace io {

class SequenceReader {
public:
    // `plan` is borrowed and must outlive the reader; null takes the defaults.
    // `keep_comment` (-y) applies to FASTA/FASTQ only.
    explicit SequenceReader(std::string path, const IoPlan* plan = nullptr,
                            bool keep_comment = false);

    explicit SequenceReader(std::vector<std::string> paths,
                            const IoPlan* plan = nullptr,
                            bool keep_comment = false);

    bool next(FastxRecord& rec);

    // Whether the current input is decoded by a thread team.
    bool decode_parallel_active() const;

private:
    void open_current();

    std::vector<std::string> paths_;
    size_t idx_ = 0;
    const IoPlan* plan_ = nullptr;
    bool keep_comment_ = false;
    std::unique_ptr<FastxReader> fx_;
    std::unique_ptr<BamReader> bam_;
};

}}}  // namespace fa::cpu::io

/**** batch_reader ****/
// Splits the input into batches of about batch_bp bases, in input order. A
// batch ends with the record that reaches the budget.

namespace fa { namespace cpu { namespace io {

class BatchReader {
public:
    // batch_bp: bp budget per batch (must be > 0).
    BatchReader(SequenceReader& reader, int64_t batch_bp);

    // Next batch, or std::nullopt at EOF. `seen` counts the records returned.
    std::optional<std::vector<FastxRecord>> next(std::atomic<int64_t>& seen);

private:
    SequenceReader& reader_;
    int64_t batch_bp_;
    bool done_ = false;
};

}}}  // namespace fa::cpu::io

/**** emitted_role ****/

namespace fa::cpu::output {

// A record's role in the output, from its place in the result tree.
enum class EmittedRole {
  Primary,
  Supplementary,
  Secondary,
  SecondarySupplementary,
};

void validate_one_level_hypotheses(const AlignResult& result);

}  // namespace fa::cpu::output

/**** alignment_aux_tags ****/
// Auxiliary tag values shared by the PAF and SAM writers.

namespace fa::cpu::output {

struct AlignmentAuxiliaryTags {
  std::optional<int> alignment_score;
  std::optional<int> edit_distance;
  std::optional<int> ambiguities;
  // Absent when the engine did not compute the value; no tag is written.
  std::optional<int> max_segment_score;     // ms:i
  std::optional<int> max_score_margin;      // md:i
  std::optional<int> chain_anchors;         // cm:i
  std::optional<int> chain_score;           // s1:i
  std::optional<int> secondary_chain_score; // s2:i
};

}  // namespace fa::cpu::output

/**** block_divergence ****/
// Aligned-block accounting shared by the PAF and SAM writers: the query
// interval, block length, matches, edit distance and gap counts behind PAF
// columns 3, 4, 10 and 11 and the NM:i and de:f tags.

namespace fa {
namespace cpu {
namespace output {

struct BlockAccounting {
  int qstart = 0;
  int qend = 0;
  int block_len = 0;
  int matches = 0;
  int nm = -1; // -1 => not computable (the record carries no CIGAR)
  int ambiguities = 0; // N bases the replay excluded from block_len/matches
  int gap_bases = 0;
  int gap_opens = 0;
};

} // namespace output
} // namespace cpu
} // namespace fa

/**** sam ****/
// SAM record and header output.

namespace fa { namespace cpu { namespace output {

struct SamReference {
    std::string name;
    int64_t length = 0;
};

int reference_span(const AlignResult& result);

// Appends the reverse complement of `seq`, for reverse-strand SEQ.
void append_reverse_complement(std::string& out, const std::string& seq);

void append_sam_record_line(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::string& extra_tags = {},
    int flag_extra = 0,  // OR-ed into FLAG when mapped (e.g. 0x800 supplementary)
    // For hard-clipped supplementaries: `cigar_override` replaces
    // result.cigar, and `seq_preoriented` means seq/qual are already
    // reference-oriented and clipped.
    const std::string* cigar_override = nullptr,
    bool seq_preoriented = false,
    // RG:Z value; empty emits no tag.
    const std::string& read_group_id = {},
    // For tp:A; consistent with flag_extra.
    EmittedRole role = EmittedRole::Primary,
    // Tag text written last, after cs/MD, on every record of the read (the -y
    // comment). extra_tags, by contrast, belongs to this record only.
    const std::string& record_tags = {}
);

// A supplementary segment hard-clipped as minimap2 does.
struct HardClippedRecord {
    std::string cigar;  // terminal soft clips rewritten to hard clips ('S' -> 'H')
    std::string seq;    // reference-oriented read, aligned portion only (clips dropped)
    std::string qual;   // matching QUAL slice; empty when the read has no QUAL
};

// A read's input tags (the uBAM passthrough) as SAM text, in input order.
// `whole` holds them all; `clipped` leaves out the per-base ones (the
// base-modification tags and a B array with one element per base of the
// read), for a record whose SEQ is not the whole read.
struct PassthroughTags {
    std::string whole;
    std::string clipped;
};

struct SamEmitOptions {
    bool hard_clip_supp = false;  // the CLI sets this unless -Y (minimap2's default)
    // One record per mapped alternative (FLAG 0x100), as --secondary yes.
    bool emit_secondary = false;
    // RG:Z for every record of the read; empty for none, including a read
    // that already carries its own RG:Z.
    std::string read_group_id;
};

// All SAM records of one read: the primary, then its supplementaries (FLAG
// 0x800), each with an SA:Z listing the others, then secondaries if enabled.
// Only supplementaries are hard-clipped; SA:Z lists the soft CIGARs, as
// minimap2 does. The passthrough goes on every record: `whole` where SEQ is
// the whole read (the primary, an unmapped read, a soft-clipped
// supplementary), `clipped` on a hard-clipped supplementary and on
// secondaries, whose SEQ is '*'. record_tags goes on every record.
void append_sam_records(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const PassthroughTags& passthrough = {},
    SamEmitOptions opts = {},
    const std::string& record_tags = {}
);

// @HD, @SQ and @PG lines. A non-empty `command_line` becomes @PG CL:, with
// tabs and newlines neutralized so an argument cannot forge a header line.
// The @PG ID is `program_id`, or `program_name` when empty; a non-empty
// `previous_program` becomes PP:.
std::vector<std::string> sam_header(
    const std::vector<SamReference>& references,
    const std::string& program_name = "flashalign",
    const std::string& version = FA_VERSION,
    const std::string& command_line = {},
    const std::string& program_id = {},
    const std::string& previous_program = {}
);

}}} // namespace fa::cpu::output

/**** paf ****/
// PAF record output.

namespace fa { namespace cpu { namespace output {

// `copy_comment` (-y) appends the read's FASTA/Q comment verbatim as the last
// field, as minimap2 does.
void write_paf_record(
    std::ostream& out,
    const ::fa::cpu::io::FastxRecord& read,
    const AlignResult& result,
    const std::unordered_map<std::string, int64_t>& ref_lengths,
    bool with_cigar,
    EmittedRole role = EmittedRole::Primary,
    bool copy_comment = false
);

// minimap2's --paf-no-hit row for a read with no alignment.
void write_paf_no_hit_record(
    std::ostream& out, const ::fa::cpu::io::FastxRecord& read);

}}}  // namespace fa::cpu::output
