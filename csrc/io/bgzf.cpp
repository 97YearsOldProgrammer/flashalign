// Block-parallel BGZF source for bgzip output and BAM; BAM and uBAM records become FastxRecords with the raw tag block kept.
//
// Threads:
//   reader   the one caller of read(); owns `cur`, `tail`, `finished` and its
//            scratch. When the run it needs is not ready it cuts and inflates
//            that run itself if no worker has claimed it, and otherwise waits.
//   workers  each cuts and inflates one run at a time while the demand gate
//            (header) is open; otherwise it runs one unit of the caller's
//            work through `helper`, or sleeps.
//
// A run is cut by the holder of the cut token (`cutting`): one pread of up to
// kRunInputBytes at `cut_pos` and a walk over the block headers. The token is
// handed back when the walk ends, before the run is inflated, so cuts are
// serial and inflates are not. Run `seq` lives in ring[seq % ring.size()];
// seq order is file order, and only the reader frees slots, in seq order, so
// the slot for next_seq is free whenever fewer than ring.size() runs are
// outstanding.
//
// One mutex, mu, guards all shared state except a slot's output while its
// filler owns it (kRunning), and the page-cache advice, which only the token
// holder touches.
#include "io.h"

#include "threading/threading.h"    // default_thread_count

#include <zlib.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

/**** bgzf_source ****/

namespace fa {
namespace cpu {
namespace io {
namespace bgzf_detail {

constexpr size_t kTrailerBytes = 8;
constexpr uint32_t kMaxBlockOutput = 1u << 16;
// Compressed input per run, and a ceiling on a run's output.
constexpr size_t kRunInputBytes = 4u << 20;
constexpr uint64_t kRunMaxOutput = 32u << 20;
// Readahead past each cut and drop-behind, as in the member-parallel source.
constexpr uint64_t kInputReadaheadBytes = 32ull << 20;
constexpr uint64_t kInputDropLagBytes = 64ull << 20;
constexpr uint64_t kInputDropStepBytes = 64ull << 20;
// Workers keep at most this many runs cut and not yet consumed, the reader's
// own included, unless the reader is waiting.
constexpr uint64_t kRunsAhead = 2;

struct Block {
  uint32_t in_off;  // in the run's input
  uint32_t total;   // BSIZE + 1
  uint32_t header_bytes;
  uint64_t out_off; // in the run's output
  uint32_t isize;
  uint32_t crc;
};

// What a thread needs to cut and inflate a run.
struct Scratch {
  std::unique_ptr<unsigned char[]> in;
  std::vector<Block> blocks;
  z_stream zs{};
  bool zs_ready = false;
  ~Scratch() {
    if (zs_ready)
      inflateEnd(&zs);
  }
};

enum class SlotState { kFree, kRunning, kReady };

struct Slot {
  SlotState state = SlotState::kFree;
  uint64_t seq = 0;
  std::unique_ptr<unsigned char[]> out; // kRunMaxOutput bytes
  // Deliverable bytes: the whole run, or the blocks before the failed one.
  uint64_t out_len = 0;
  uint64_t out_pos = 0; // reader-owned once kReady
  bool failed = false;
  uint64_t end_offset = 0;
  bool last_bgzf_empty = false;
  uint64_t fail_offset = 0; // file offset of the first block that failed
};

enum class Stop { kBoundary, kEnd, kTail };

} // namespace bgzf_detail

using namespace bgzf_detail;

struct BgzfByteSource::Impl {
  std::string path;
  InputPolicy policy = InputPolicy::gzip_compatible;
  ByteSourceStatus status;
  int fd = -1;
  uint64_t size = 0;
  unsigned workers = 1;
  std::function<bool(unsigned)> helper;

  // ---- token holder only ----
  uint64_t input_dropped_to = 0;

  // ---- reader-owned ----
  Slot* cur = nullptr;
  std::unique_ptr<ZlibByteSource> tail;
  bool finished = false;
  bool broken = false; // an internal error was reported; every read fails
  Scratch reader_scratch;

  // ---- mu-guarded ----
  std::mutex mu;
  std::condition_variable work_cv; // workers wait for a run to claim
  std::condition_variable data_cv; // the reader waits for a run or the token
  bool abort = false;
  bool stopped = false;   // a run failed: nothing further is cut
  uint64_t fail_seq = ~uint64_t{0}; // the first run that failed
  bool cutting = false;   // the cut token is held
  bool walk_done = false; // no run follows the last one cut
  bool has_tail = false;  // the walk stopped at tail_offset before EOF
  uint64_t tail_offset = 0;
  uint64_t cut_pos = 0;
  uint64_t next_seq = 0; // the next run to be cut
  uint64_t read_seq = 0; // the run the reader is on
  bool reader_waiting = false;
  std::vector<Slot> ring;

  std::atomic<bool> stopping{false};
  std::vector<std::thread> threads;

  ~Impl() {
    shutdown();
    if (fd >= 0)
      ::close(fd);
  }

  // ---- helpers, called with mu held ----

  bool can_cut_locked() const {
    return !abort && !stopped && !cutting && !walk_done &&
           next_seq - read_seq < ring.size();
  }
  // `pending` counts a run whose cut is in progress.
  bool worker_wanted_locked(uint64_t pending = 0) const {
    return reader_waiting || next_seq + pending - read_seq < kRunsAhead;
  }

  // ---- the job: cut one run, then inflate it ----

  // Called with the token held and without mu. Reads the run at `pos` into
  // sc.in and lists its blocks; returns where the walk stopped and why.
  Stop cut(uint64_t pos, Scratch& sc, uint64_t& end) {
    sc.blocks.clear();
    end = pos;
    uint64_t out_size = 0;
    if (pos >= size)
      return pos == size ? Stop::kEnd : Stop::kTail;
    if (!sc.in)
      sc.in.reset(new unsigned char[kRunInputBytes]);
    const size_t want =
        static_cast<size_t>(std::min<uint64_t>(kRunInputBytes, size - pos));
    const int64_t got = pread_exact(fd, sc.in.get(), want, pos);
    if (got < 0)
      return Stop::kTail; // the serial tail reports the read error
    advise_input(pos + static_cast<uint64_t>(got));
    const size_t have = static_cast<size_t>(got);
    // A short read means the file shrank; the tail sees what is there.
    const uint64_t file_end = have < want ? pos + have : size;
    const unsigned char* in = sc.in.get();
    size_t p = 0;
    Stop stop = Stop::kBoundary;
    for (;;) {
      const uint64_t at = pos + p;
      if (at == file_end) {
        stop = file_end == size ? Stop::kEnd : Stop::kTail;
        break;
      }
      const unsigned char* h = in + p;
      const BgzfHeader header = parse_bgzf_header(h, have - p);
      if (header.kind == BgzfHeaderKind::need_more) {
        if (pos + have == file_end) stop = Stop::kTail;
        break;
      }
      if (header.kind != BgzfHeaderKind::bgzf) {
        stop = Stop::kTail;
        break;
      }
      const size_t total = header.block_bytes;
      if (at + total > file_end) {
        stop = Stop::kTail;
        break;
      }
      if (have - p < total)
        break;
      const uint32_t isize = read_le32(h + total - 4);
      if (isize > kMaxBlockOutput) {
        stop = Stop::kTail;
        break;
      }
      if (!sc.blocks.empty() && out_size + isize > kRunMaxOutput)
        break;
      sc.blocks.push_back(Block{static_cast<uint32_t>(p),
                                static_cast<uint32_t>(total),
                                static_cast<uint32_t>(header.header_bytes), out_size, isize,
                                read_le32(h + total - kTrailerBytes)});
      out_size += isize;
      p += total;
    }
    end = pos + p;
    // A whole block always fits in a run's input.
    assert(!(stop == Stop::kBoundary && sc.blocks.empty()));
    return stop;
  }

  // Inflates sc's blocks into `slot`, without mu: the caller set it kRunning.
  // Stops at the first block that fails a check.
  void inflate_run(Scratch& sc, Slot& slot, uint64_t run_start) {
    slot.failed = false;
    slot.out_pos = 0;
    slot.end_offset = run_start;
    slot.last_bgzf_empty = false;
    try {
      // Sized once for any run; only the pages a run writes become resident.
      if (!slot.out)
        slot.out.reset(new unsigned char[kRunMaxOutput]);
      if (!sc.zs_ready) {
        sc.zs = z_stream{};
        if (inflateInit2(&sc.zs, -15) != Z_OK)
          throw std::runtime_error("inflateInit2");
        sc.zs_ready = true;
      }
    } catch (...) {
      slot.failed = true;
      slot.fail_offset = run_start;
      slot.out_len = 0;
      return;
    }
    unsigned char* out = slot.out.get();
    // Written to the slot once, at the end: out_len shares a cache line with
    // the reader's out_pos.
    uint64_t out_len = 0;
    for (const Block& b : sc.blocks) {
      if (stopping.load(std::memory_order_relaxed)) {
        slot.failed = true;
        slot.fail_offset = run_start + b.in_off;
        break;
      }
      unsigned char* dst = out + b.out_off;
      z_stream& zs = sc.zs;
      inflateReset(&zs);
      zs.next_in = sc.in.get() + b.in_off + b.header_bytes;
      zs.avail_in = static_cast<uInt>(b.total - b.header_bytes - kTrailerBytes);
      zs.next_out = dst;
      zs.avail_out = static_cast<uInt>(b.isize);
      const int ret = inflate(&zs, Z_FINISH);
      const bool ok = ret == Z_STREAM_END && zs.avail_in == 0 &&
                      zs.avail_out == 0 &&
                      static_cast<uint32_t>(crc32_z(0, dst, b.isize)) == b.crc;
      if (!ok) {
        slot.failed = true;
        slot.fail_offset = run_start + b.in_off;
        break;
      }
      out_len = b.out_off + b.isize;
      slot.end_offset = run_start + b.in_off + b.total;
      slot.last_bgzf_empty = b.isize == 0;
    }
    slot.out_len = out_len;
  }

  // Takes the token (held on entry, under `lock`), cuts and inflates one run,
  // and publishes it. Returns with `lock` held.
  void run_job(std::unique_lock<std::mutex>& lock, Scratch& sc) {
    const uint64_t pos = cut_pos;
    lock.unlock();
    uint64_t end = pos;
    Stop stop = Stop::kTail;
    try {
      stop = cut(pos, sc, end);
    } catch (...) {
      sc.blocks.clear();
      end = pos;
      stop = Stop::kTail;
    }
    lock.lock();
    cutting = false;
    cut_pos = end;
    if (stop != Stop::kBoundary) {
      walk_done = true;
      has_tail = stop == Stop::kTail;
      tail_offset = end;
    }
    if (sc.blocks.empty() || abort || stopped) {
      work_cv.notify_all();
      data_cv.notify_all();
      return;
    }
    const uint64_t seq = next_seq++;
    Slot& slot = ring[seq % ring.size()];
    slot.state = SlotState::kRunning;
    slot.seq = seq;
    // The token is free again.
    work_cv.notify_all();
    data_cv.notify_all();
    lock.unlock();
    inflate_run(sc, slot, pos);
    lock.lock();
    slot.state = SlotState::kReady;
    if (slot.failed) {
      stopped = true;
      fail_seq = std::min(fail_seq, seq);
    }
    // A run after a failed one is never read.
    if (seq > fail_seq)
      slot.out.reset();
    data_cv.notify_all();
  }

  void advise_input(uint64_t pos) {
#if defined(POSIX_FADV_WILLNEED) && defined(POSIX_FADV_DONTNEED)
    if (pos < size) {
      const uint64_t ahead =
          std::min<uint64_t>(kInputReadaheadBytes, size - pos);
      (void)posix_fadvise(fd, static_cast<off_t>(pos),
                          static_cast<off_t>(ahead), POSIX_FADV_WILLNEED);
    }
    if (pos > kInputDropLagBytes) {
      const uint64_t drop_end = pos - kInputDropLagBytes;
      if (drop_end >= input_dropped_to + kInputDropStepBytes) {
        (void)posix_fadvise(fd, static_cast<off_t>(input_dropped_to),
                            static_cast<off_t>(drop_end - input_dropped_to),
                            POSIX_FADV_DONTNEED);
        input_dropped_to = drop_end;
      }
    }
#else
    (void)pos;
#endif
  }

  void worker_main(unsigned ordinal) {
    Scratch sc;
    std::unique_lock<std::mutex> lock(mu);
    // Set when the helper just found no work, so the loop does not spin.
    bool helper_was_empty = false;
    for (;;) {
      if (abort)
        return;
      if (can_cut_locked()) {
        if (worker_wanted_locked()) {
          cutting = true;
          run_job(lock, sc);
          helper_was_empty = false;
          continue;
        }
      }
      // Nothing is cut after the walk ends, so the input buffer goes.
      if (sc.in && (walk_done || stopped))
        sc.in.reset();
      // A cut in progress ends within one pread: wait for the token rather
      // than leave for a unit of the caller's work.
      const bool token_soon = cutting && !stopped && !walk_done &&
                              next_seq + 1 - read_seq < ring.size() &&
                              worker_wanted_locked(1);
      if (!token_soon && helper && !helper_was_empty) {
        lock.unlock();
        const bool ran = helper(ordinal);
        lock.lock();
        helper_was_empty = !ran;
        continue;
      }
      if (helper && !token_soon) {
        // Poll: the caller's work does not signal work_cv.
        work_cv.wait_for(lock, std::chrono::milliseconds(2));
      } else {
        work_cv.wait(lock);
      }
      helper_was_empty = false;
    }
  }

  // Ends the parallel part: frees the ring's idle buffers and the reader's
  // input buffer.
  void release_ring_locked() {
    reader_scratch.in.reset();
    for (Slot& slot : ring) {
      if (slot.state == SlotState::kRunning)
        continue;
      slot.state = SlotState::kFree;
      slot.out.reset();
    }
  }

  int read(char* dst, size_t want) {
    for (;;) {
      if (broken)
        return -1;
      if (tail)
        return tail->read(dst, static_cast<int>(want));
      if (finished)
        return 0;
      if (cur != nullptr) {
        if (cur->out_pos < cur->out_len) {
          const size_t n = static_cast<size_t>(
              std::min<uint64_t>(want, cur->out_len - cur->out_pos));
          std::memcpy(dst, cur->out.get() + cur->out_pos, n);
          cur->out_pos += n;
          status.compressed_offset = cur->end_offset;
          return static_cast<int>(n);
        }
        std::unique_lock<std::mutex> lock(mu);
        Slot* done = cur;
        cur = nullptr;
        status.compressed_offset = done->end_offset;
        status.last_bgzf_empty = done->last_bgzf_empty;
        if (done->failed) {
          // Everything before the failed block is delivered; the serial
          // source takes it from there.
          stopped = true;
          const uint64_t at = done->fail_offset;
          release_ring_locked();
          work_cv.notify_all();
          lock.unlock();
          tail.reset(new ZlibByteSource(path, at, at > 0, policy,
                                        status.last_bgzf_empty));
          continue;
        }
        done->state = SlotState::kFree;
        ++read_seq;
        work_cv.notify_all();
        continue;
      }
      std::unique_lock<std::mutex> lock(mu);
      for (;;) {
        Slot& next = ring[read_seq % ring.size()];
        if (next.state == SlotState::kReady && next.seq == read_seq) {
          cur = &next;
          break;
        }
        if (read_seq == next_seq && walk_done && !cutting) {
          release_ring_locked();
          work_cv.notify_all();
          if (has_tail) {
            const uint64_t at = tail_offset;
            lock.unlock();
            tail.reset(new ZlibByteSource(path, at, at > 0, policy,
                                        status.last_bgzf_empty));
          } else {
            finished = true;
          }
          break;
        }
        // Only the run the reader needs, and only when no worker has it.
        if (read_seq == next_seq && can_cut_locked()) {
          cutting = true;
          run_job(lock, reader_scratch);
          continue;
        }
        reader_waiting = true;
        work_cv.notify_all();
        data_cv.wait(lock);
        reader_waiting = false;
      }
    }
  }

  void shutdown() {
    stopping.store(true, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> guard(mu);
      abort = true;
    }
    work_cv.notify_all();
    data_cv.notify_all();
    for (std::thread& thread : threads) {
      if (thread.joinable())
        thread.join();
    }
    threads.clear();
  }
};

bool bgzf_eligible(const std::string& path, uint64_t min_bytes) {
  if (path.empty() || path == "-")
    return false;
  struct stat info;
  if (::stat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode))
    return false;
  if (info.st_size < 0 || static_cast<uint64_t>(info.st_size) < min_bytes)
    return false;
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    return false;
  unsigned char head[kGzipExtraHeaderBytes];
  const int64_t got = pread_exact(fd, head, sizeof(head), 0);
  ::close(fd);
  return got >= 0 && parse_bgzf_header(head, static_cast<size_t>(got)).kind ==
                         BgzfHeaderKind::bgzf;
}

BgzfByteSource::BgzfByteSource(const std::string& path,
                               GzMemberSourceOptions options, InputPolicy policy)
    : impl_(new Impl()) {
  Impl* self = impl_.get();
  self->path = path;
  self->policy = policy;
  self->status.compressed_offset = 0;
  self->fd = ::open(path.c_str(), O_RDONLY);
  if (self->fd < 0)
    throw std::runtime_error("failed to open FASTA/FASTQ: " + path);
  struct stat info;
  if (::fstat(self->fd, &info) != 0 || !S_ISREG(info.st_mode))
    throw std::runtime_error("not a regular file: " + path);
  self->size = static_cast<uint64_t>(info.st_size);
  // advise_input() adds explicit readahead and drop-behind.
#ifdef POSIX_FADV_SEQUENTIAL
  (void)posix_fadvise(self->fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

  unsigned wanted = options.workers;
  if (wanted == 0) {
    const int available = ::fa::cpu::threading::default_thread_count();
    wanted =
        std::min(8u, available < 1 ? 1u : static_cast<unsigned>(available));
  }
  self->workers = std::max(1u, wanted);
  self->helper = options.helper;
  self->ring.resize(2 * (static_cast<size_t>(self->workers) + 1));

  self->threads.reserve(self->workers);
  for (unsigned i = 0; i < self->workers; ++i) {
    self->threads.emplace_back([self, i] {
      fa::cpu::threading::name_current_thread("fa-gz", static_cast<int>(i));
      self->worker_main(i);
    });
  }
}

BgzfByteSource::~BgzfByteSource() = default;

int BgzfByteSource::read(void* buf, int len) {
  if (buf == nullptr || len <= 0)
    return 0;
  try {
    return impl_->read(static_cast<char*>(buf), static_cast<size_t>(len));
  } catch (const std::exception& e) {
    // Nothing may propagate into kseq.
    impl_->status.error = e.what();
    impl_->broken = true;
    return -1;
  } catch (...) {
    impl_->status.error = "BGZF decoder failed";
    impl_->broken = true;
    return -1;
  }
}

bool BgzfByteSource::cut_in_data() const {
  return impl_->tail && impl_->tail->cut_in_data();
}

ByteSourceStatus BgzfByteSource::status() const {
  return impl_->tail ? impl_->tail->status() : impl_->status;
}

} // namespace io
} // namespace cpu
} // namespace fa

/**** bgzf_header ****/

namespace fa { namespace cpu { namespace io {

static unsigned bgzf_le16(const unsigned char* p) {
    return static_cast<unsigned>(p[0]) | static_cast<unsigned>(p[1]) << 8;
}

BgzfHeader parse_bgzf_header(const unsigned char* p, size_t n) {
    if (n < 10) return {};
    if (p[0] != 0x1f || p[1] != 0x8b || p[2] != 8 || (p[3] & 0xe0))
        return {BgzfHeaderKind::malformed, 0, 0, "invalid gzip header"};
    if (!(p[3] & 4)) return {BgzfHeaderKind::gzip};
    if (n < 12) return {};
    const size_t end = 12 + bgzf_le16(p + 10);
    if (n < end) return {};
    size_t block_bytes = 0;
    for (size_t at = 12; at < end;) {
        if (end - at < 4)
            return {BgzfHeaderKind::malformed, 0, 0, "incomplete gzip extra subfield"};
        const size_t len = bgzf_le16(p + at + 2);
        if (len > end - at - 4)
            return {BgzfHeaderKind::malformed, 0, 0, "gzip extra subfield exceeds XLEN"};
        if (p[at] == 'B' && p[at + 1] == 'C') {
            if (len != 2 || block_bytes != 0)
                return {BgzfHeaderKind::malformed, 0, 0, "invalid or duplicate BGZF BC subfield"};
            block_bytes = 1 + bgzf_le16(p + at + 4);
        }
        at += 4 + len;
    }
    if (block_bytes == 0) return {BgzfHeaderKind::gzip};
    if (p[3] != 4)
        return {BgzfHeaderKind::malformed, 0, 0, "invalid BGZF flags"};
    if (block_bytes < end + 2 + 8 || block_bytes > kBgzfMaxBytes)
        return {BgzfHeaderKind::malformed, 0, 0, "BGZF BSIZE is smaller than its header and trailer"};
    return {BgzfHeaderKind::bgzf, end, block_bytes};
}

}}} // namespace fa::cpu::io

/**** bgzf ****/

namespace fa { namespace cpu { namespace io {

BgzfReader::BgzfReader(std::string path, std::unique_ptr<ByteSource> src)
    : path_(std::move(path)), src_(std::move(src)) {}

size_t BgzfReader::read(void* dst, size_t len) {
    char* out = static_cast<char*>(dst);
    size_t got = 0;
    while (got < len) {
        const size_t ask = std::min<size_t>(len - got, INT_MAX);
        const int n = src_->read(out + got, static_cast<int>(ask));
        if (n < 0) {
            const auto status = src_->status();
            throw error(status.error.empty() ? "corrupt or unreadable BAM input" :
                        status.error);
        }
        if (n == 0) break;
        got += static_cast<size_t>(n);
    }
    return got;
}

void BgzfReader::read_exact(void* dst, size_t len, const char* field) {
    if (read(dst, len) != len) throw error(std::string("truncated ") + field);
}

bool BgzfReader::decode_parallel_active() const { return src_->is_parallel(); }

std::runtime_error BgzfReader::error(const std::string& reason) const {
    const auto status = src_->status();
    std::string text = path_ + ": " + reason;
    if (status.compressed_offset != ~uint64_t{0}) {
        text += status.at_member_start ? " (member starts at compressed offset " :
                                        " (decoder compressed offset ";
        text += std::to_string(status.compressed_offset) + ")";
    }
    return std::runtime_error(text);
}

void BgzfReader::check_eof_marker() {
    if (checked_eof_) return;
    checked_eof_ = true;
    if (!src_->status().last_bgzf_empty) {
        const auto warning = error("BGZF EOF marker absent; input may be incomplete");
        std::fprintf(stderr, "flashalign: warning: %s\n", warning.what());
    }
}

}}}  // namespace fa::cpu::io

/**** bam_reader ****/

namespace fa { namespace cpu { namespace io {

bool path_is_bam(const std::string& path) {
    if (path == "-") return false;  // stdin BAM detection deferred
    // Only regular files are probed: probing a pipe would consume its bytes.
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    unsigned char m[2] = {0, 0};
    f.read(reinterpret_cast<char*>(m), 2);
    if (f.gcount() < 2 || m[0] != 0x1f || m[1] != 0x8b) return false;  // not gzip/BGZF
    // A probe decodes only enough bytes to recognize BAM, across member boundaries.
    // Its gzip-compatible policy leaves FASTQ prefix recovery unchanged.
    BgzfReader r(path, std::make_unique<ZlibByteSource>(path, InputPolicy::gzip_probe));
    char magic[4] = {0, 0, 0, 0};
    const size_t got = r.read(magic, 4);
    if (got > 0 && got < 4 && std::memcmp(magic, "BAM\1", got) == 0)
        throw r.error("truncated BAM magic");
    if (got == 4 && std::memcmp(magic, "BAM", 3) == 0) {
        if (magic[3] != 1) throw r.error("unsupported BAM version");
        return true;
    }
    return false;
}

bool path_is_cram(const std::string& path) {
    if (path == "-") return false;
    struct stat st;  // non-rewindable inputs: skip the probe (see path_is_bam)
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char m[4] = {0, 0, 0, 0};
    f.read(m, 4);
    return f.gcount() == 4 && std::string(m, 4) == std::string("CRAM", 4);
}

BamReader::BamReader(const std::string& path, const IoPlan* plan)
    : path_(path),
      bgzf_(std::make_unique<BgzfReader>(
          path, make_byte_source(path, plan, InputPolicy::bam_integrity))) {
    parse_header();
}

bool BamReader::decode_parallel_active() const { return bgzf_->decode_parallel_active(); }

const std::string& BamReader::sam_header_text() const { return header_text_; }

bool BamReader::next(FastxRecord& rec) {
    uint16_t flag = 0;
    while (read_record(rec, flag)) {
        if (!(flag & 0x4) || (flag & 0x900)) warn_aligned_input();
        if (flag & 0x900) continue;
        if (!(flag & 0x4)) drop_alignment_tags(rec);
        return true;
    }
    return false;
}

bool BamReader::read_record(FastxRecord& rec, uint16_t& flag_out) {
    int32_t block_size = 0;
    const size_t got = bgzf_->read(&block_size, 4);
    if (got == 0) {
        bgzf_->check_eof_marker();
        return false;
    }
    if (got != 4) throw bgzf_->error("truncated BAM record size");
    // Every field is bounds-checked against the block (`need`); the fixed
    // part of a record is 32 bytes.
    if (block_size < 32) {
        throw bgzf_->error("malformed BAM: record block_size too small");
    }
    std::string buf(static_cast<size_t>(block_size), '\0');
    bgzf_->read_exact(&buf[0], buf.size(), "BAM record");

    const size_t n = buf.size();
    size_t c = 0;
    auto need = [&](size_t k) {
        if (k > n || c > n - k) {  // c + k > n, overflow-safe
            throw bgzf_->error("malformed BAM: record field exceeds block");
        }
    };
    auto i32 = [&] { need(4); int32_t v; std::memcpy(&v, buf.data() + c, 4); c += 4; return v; };
    auto u16 = [&] { need(2); uint16_t v; std::memcpy(&v, buf.data() + c, 2); c += 2; return v; };
    auto u8  = [&] { need(1); return static_cast<uint8_t>(buf[c++]); };
    auto skip = [&](size_t k) { need(k); c += k; };

    skip(4);                      // refID (we re-align, ignore)
    skip(4);                      // pos
    const uint8_t l_read_name = u8();
    skip(1);                      // mapq
    skip(2);                      // bin
    const uint16_t n_cigar = u16();
    const uint16_t flag = u16();
    flag_out = flag;
    const int32_t l_seq = i32();
    skip(4 + 4 + 4);              // next_refID, next_pos, tlen
    if (l_seq < 0) throw bgzf_->error("malformed BAM: negative l_seq");

    need(l_read_name);
    rec.name.assign(buf.data() + c, l_read_name ? l_read_name - 1 : 0);  // drop NUL
    c += l_read_name;
    skip(static_cast<size_t>(n_cigar) * 4);  // CIGAR (unused for re-alignment)

    // Check the packed bytes exist before allocating rec.seq.
    static const char NT16[17] = "=ACMGRSVTWYHKDBN";
    const size_t seq_bytes = (static_cast<size_t>(l_seq) + 1) / 2;
    need(seq_bytes);
    rec.seq.assign(static_cast<size_t>(l_seq), 'N');
    for (int32_t k = 0; k < l_seq; ++k) {
        const uint8_t byte = static_cast<uint8_t>(buf[c + (k >> 1)]);
        const uint8_t nib = (k & 1) ? (byte & 0x0f) : (byte >> 4);
        rec.seq[k] = NT16[nib];
    }
    c += seq_bytes;

    // QUAL: Phred + 33; a leading 0xff means no qualities.
    rec.qual.clear();
    if (l_seq > 0) {
        need(static_cast<size_t>(l_seq));
        if (static_cast<uint8_t>(buf[c]) != 0xff) {
            rec.qual.resize(static_cast<size_t>(l_seq));
            for (int32_t k = 0; k < l_seq; ++k) {
                rec.qual[k] = static_cast<char>(static_cast<uint8_t>(buf[c + k]) + 33);
            }
        }
        c += static_cast<size_t>(l_seq);
    }

    // The rest is the tag block.
    rec.tag_bytes.assign(buf.data() + c, n - c);
    rec.comment.clear();

    // Undo the reverse complement of a reverse-mapped record.
    if (flag & 0x10) {
        std::string rcseq;
        output::append_reverse_complement(rcseq, rec.seq);
        rec.seq.swap(rcseq);
        if (!rec.qual.empty()) std::reverse(rec.qual.begin(), rec.qual.end());
    }
    return true;
}

void BamReader::warn_aligned_input() {
    if (warned_aligned_) return;
    warned_aligned_ = true;
    // One warning per path, though a multi-part index reads it per part.
    static std::mutex mu;
    static std::set<std::string> warned;
    {
        std::lock_guard<std::mutex> lock(mu);
        if (!warned.insert(path_).second) return;
    }
    std::fprintf(stderr,
                 "flashalign: warning: %s: aligned BAM input; secondary and "
                 "supplementary records skipped, alignment tags dropped\n",
                 path_.c_str());
}

void BamReader::parse_header() {
    char magic[4];
    if (bgzf_->read(magic, 4) != 4 || std::string(magic, 4) != std::string("BAM\1", 4)) {
        throw bgzf_->error("invalid BAM magic");
    }
    int32_t l_text = read_i32();
    if (l_text < 0) throw bgzf_->error("malformed BAM: negative l_text");
    // Kept so @RG/@CO lines can be carried into the output header.
    header_text_.resize(static_cast<size_t>(l_text));
    if (l_text) bgzf_->read_exact(&header_text_[0], header_text_.size(), "BAM header text");
    // The reference dictionary is skipped: the reads are re-aligned.
    const int32_t n_ref = read_i32();
    if (n_ref < 0) throw bgzf_->error("malformed BAM: negative n_ref");
    for (int32_t i = 0; i < n_ref; ++i) {
        const int32_t l_name = read_i32();
        if (l_name < 0) throw bgzf_->error("malformed BAM: negative l_name");
        skip(static_cast<size_t>(l_name), "BAM reference name");
        read_i32();  // l_ref
    }
}

int32_t BamReader::read_i32() { int32_t v = 0; bgzf_->read_exact(&v, 4, "BAM header"); return v; }

void BamReader::skip(size_t n, const char* what) {
    std::string tmp(n, '\0');
    if (n) bgzf_->read_exact(&tmp[0], n, what);
}

[[noreturn]] static void throw_malformed_tag_block(const std::string& read_name,
                                                  const char* what) {
    throw std::runtime_error("malformed BAM: tag block of read " + read_name +
                             ": " + what);
}

bool bam_tag_at(const std::string& blob, size_t at) {
    return blob.size() - at >= 3;
}

size_t bam_tag_size(const std::string& blob, size_t at,
                    const std::string& read_name) {
    const size_t left = blob.size() - at;
    if (left < 4) throw_malformed_tag_block(read_name, "a tag runs past the block");
    auto fits = [&](size_t size) {
        if (size > left) throw_malformed_tag_block(read_name, "a tag runs past the block");
        return size;
    };
    switch (blob[at + 2]) {
        case 'A': case 'c': case 'C': return 4;
        case 's': case 'S': return fits(5);
        case 'i': case 'I': case 'f': return fits(7);
        case 'd': return fits(11);
        case 'Z': case 'H': {
            const void* nul = std::memchr(blob.data() + at + 3, '\0', left - 3);
            if (nul == nullptr) {
                throw_malformed_tag_block(read_name, "a Z/H tag has no terminating NUL");
            }
            return static_cast<size_t>(static_cast<const char*>(nul) - (blob.data() + at)) + 1;
        }
        case 'B': {
            fits(8);
            size_t width = 0;
            switch (blob[at + 3]) {
                case 'c': case 'C': width = 1; break;
                case 's': case 'S': width = 2; break;
                case 'i': case 'I': case 'f': width = 4; break;
                default: throw_malformed_tag_block(read_name, "a B tag has an unknown subtype");
            }
            uint32_t count;
            std::memcpy(&count, blob.data() + at + 4, 4);
            if ((left - 8) / width < count) {
                throw_malformed_tag_block(read_name, "a tag runs past the block");
            }
            return 8 + static_cast<size_t>(count) * width;
        }
        default: throw_malformed_tag_block(read_name, "a tag has an unknown type");
    }
}

// The alignment tags of the aligner that wrote a mapped record: dorado's list
// (minimap2's tags and bh, TS) with bwa's XS and FlashAlign's md and XA.
void BamReader::drop_alignment_tags(FastxRecord& rec) {
    static const char kKeys[][3] = {"NM", "MD", "AS", "XS", "SA", "cs", "ms",
                                    "md", "nn", "de", "dv", "tp", "cm", "s1",
                                    "s2", "zd", "rl", "bh", "TS", "XA"};
    const std::string& blob = rec.tag_bytes;
    std::string kept;
    kept.reserve(blob.size());
    for (size_t at = 0; bam_tag_at(blob, at);) {
        const size_t size = bam_tag_size(blob, at, rec.name);
        bool drop = false;
        for (const char* key : kKeys) {
            if (blob[at] == key[0] && blob[at + 1] == key[1]) drop = true;
        }
        if (!drop) kept.append(blob, at, size);
        at += size;
    }
    rec.tag_bytes.swap(kept);
}

namespace detail {

static void append_printf_g(std::string& out, double d) {
    char buf[32];
    const int n = std::snprintf(buf, sizeof buf, "%g", d);
    out.append(buf, static_cast<size_t>(n));
}

// htslib's kputd (kstring.c), which spells B:f elements: six significant
// digits, trailing zeros cut, printf's %g outside [1e-4, 999999].
static void append_kputd(std::string& out, double d) {
    if (d == 0) {
        out += std::signbit(d) ? "-0" : "0";
        return;
    }
    if (d < 0) {
        out.push_back('-');
        d = -d;
    }
    if (!(d >= 0.0001 && d <= 999999)) {
        append_printf_g(out, d);
        return;
    }
    char buf[21];
    char* cp = buf + 20;
    uint32_t i;
    if (d < 0.001)       i = static_cast<uint32_t>(std::rint(d * 1000000000)), cp -= 1;
    else if (d < 0.01)   i = static_cast<uint32_t>(std::rint(d * 100000000)),  cp -= 2;
    else if (d < 0.1)    i = static_cast<uint32_t>(std::rint(d * 10000000)),   cp -= 3;
    else if (d < 1)      i = static_cast<uint32_t>(std::rint(d * 1000000)),    cp -= 4;
    else if (d < 10)     i = static_cast<uint32_t>(std::rint(d * 100000)),     cp -= 5;
    else if (d < 100)    i = static_cast<uint32_t>(std::rint(d * 10000)),      cp -= 6;
    else if (d < 1000)   i = static_cast<uint32_t>(std::rint(d * 1000)),       cp -= 7;
    else if (d < 10000)  i = static_cast<uint32_t>(std::rint(d * 100)),        cp -= 8;
    else if (d < 100000) i = static_cast<uint32_t>(std::rint(d * 10)),         cp -= 9;
    else                 i = static_cast<uint32_t>(std::rint(d)),              cp -= 10;
    static const char dig2[] =
        "00010203040506070809" "10111213141516171819" "20212223242526272829"
        "30313233343536373839" "40414243444546474849" "50515253545556575859"
        "60616263646566676869" "70717273747576777879" "80818283848586878889"
        "90919293949596979899";
    std::memcpy(cp -= 2, &dig2[2 * (i % 100)], 2); i /= 100;
    std::memcpy(cp -= 2, &dig2[2 * (i % 100)], 2); i /= 100;
    std::memcpy(cp -= 2, &dig2[2 * (i % 100)], 2);
    if (i >= 100) *--cp = static_cast<char>('0' + i / 100);  // rounded up to 7 digits
    char* ep;
    int p = static_cast<int>(buf + 20 - cp);
    if (p <= 10) {  // d < 1
        ep = cp + 5;
        while (p < 10) {
            *--cp = '0';
            ++p;
        }
        *--cp = '.';
        *--cp = '0';
    } else {
        char* xp = --cp;
        ep = cp + 6;
        while (p > 10) {
            xp[0] = xp[1];
            ++xp;
            --p;
        }
        xp[0] = '.';
    }
    while (*ep == '0' && ep > cp) --ep;
    if (*ep != '.') ++ep;
    out.append(cp, static_cast<size_t>(ep - cp));
}

}  // namespace detail

void append_bam_tag_sam_text(std::string& out, const char* tag, size_t size) {
    out.append(tag, 2);
    out.push_back(':');
    const char* v = tag + 3;
    switch (tag[2]) {
        case 'A': out += "A:"; out.push_back(v[0]); break;
        case 'c': out += "i:"; detail::append_decimal(out, static_cast<int>(detail::load_le<int8_t>(v))); break;
        case 'C': out += "i:"; detail::append_decimal(out, static_cast<int>(detail::load_le<uint8_t>(v))); break;
        case 's': out += "i:"; detail::append_decimal(out, detail::load_le<int16_t>(v)); break;
        case 'S': out += "i:"; detail::append_decimal(out, detail::load_le<uint16_t>(v)); break;
        case 'i': out += "i:"; detail::append_decimal(out, detail::load_le<int32_t>(v)); break;
        case 'I': out += "i:"; detail::append_decimal(out, detail::load_le<uint32_t>(v)); break;
        case 'f': out += "f:"; detail::append_printf_g(out, detail::load_le<float>(v)); break;
        case 'd': out += "d:"; detail::append_printf_g(out, detail::load_le<double>(v)); break;
        case 'Z': case 'H':
            out.push_back(tag[2]);
            out.push_back(':');
            out.append(v, size - 4);  // without the NUL
            break;
        case 'B': {
            const char sub = v[0];
            const uint32_t count = detail::load_le<uint32_t>(v + 1);
            const char* e = v + 5;
            out += "B:";
            out.push_back(sub);
            for (uint32_t k = 0; k < count; ++k) {
                out.push_back(',');
                switch (sub) {
                    case 'c': detail::append_decimal(out, static_cast<int>(detail::load_le<int8_t>(e))); e += 1; break;
                    case 'C': detail::append_decimal(out, static_cast<int>(detail::load_le<uint8_t>(e))); e += 1; break;
                    case 's': detail::append_decimal(out, detail::load_le<int16_t>(e)); e += 2; break;
                    case 'S': detail::append_decimal(out, detail::load_le<uint16_t>(e)); e += 2; break;
                    case 'i': detail::append_decimal(out, detail::load_le<int32_t>(e)); e += 4; break;
                    case 'I': detail::append_decimal(out, detail::load_le<uint32_t>(e)); e += 4; break;
                    default: detail::append_kputd(out, detail::load_le<float>(e)); e += 4; break;  // 'f'
                }
            }
            break;
        }
    }
}

}}}  // namespace fa::cpu::io
