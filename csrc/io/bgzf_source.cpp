// Block-parallel BGZF source; the output contract is in the header.
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
#include "bgzf_source.h"

#include "threading/parallel_for.h" // default_thread_count
#include "threading/thread_name.h"

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
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace fa {
namespace cpu {
namespace io {
namespace bgzf_detail {

constexpr size_t kHeaderBytes = 18;
constexpr size_t kTrailerBytes = 8;
// Header, the shortest deflate stream (an empty fixed block) and trailer.
constexpr size_t kMinBlockBytes = 28;
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

inline uint32_t read_le32(const unsigned char* p) {
  return static_cast<uint32_t>(p[0]) | static_cast<uint32_t>(p[1]) << 8 |
         static_cast<uint32_t>(p[2]) << 16 | static_cast<uint32_t>(p[3]) << 24;
}

// Gzip magic, CM = deflate, FLG = FEXTRA only, XLEN = 6 holding one 'BC'
// subfield with SLEN = 2; the deflate data then starts at byte 18.
inline bool bgzf_shaped_header(const unsigned char* h) {
  return h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && h[3] == 4 &&
         h[10] == 6 && h[11] == 0 && h[12] == 'B' && h[13] == 'C' &&
         h[14] == 2 && h[15] == 0;
}

// Reads the whole request, retrying EINTR. Returns the byte count (short only
// at end of file) or -1 on error.
int64_t pread_exact(int fd, void* buf, size_t count, uint64_t offset) {
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

struct Block {
  uint32_t in_off;  // in the run's input
  uint32_t total;   // BSIZE + 1
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
  uint64_t fail_offset = 0; // file offset of the first block that failed
};

enum class Stop { kBoundary, kEnd, kTail };

} // namespace bgzf_detail

using namespace bgzf_detail;

struct BgzfByteSource::Impl {
  std::string path;
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
      if (have - p < kHeaderBytes) {
        if (at + kHeaderBytes > file_end)
          stop = Stop::kTail;
        break;
      }
      const unsigned char* h = in + p;
      if (!bgzf_shaped_header(h)) {
        stop = Stop::kTail;
        break;
      }
      const size_t total = (static_cast<size_t>(h[16]) |
                            static_cast<size_t>(h[17]) << 8) + 1;
      if (total < kMinBlockBytes || at + total > file_end) {
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
                                static_cast<uint32_t>(total), out_size, isize,
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
      zs.next_in = sc.in.get() + b.in_off + kHeaderBytes;
      zs.avail_in = static_cast<uInt>(b.total - kHeaderBytes - kTrailerBytes);
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
          return static_cast<int>(n);
        }
        std::unique_lock<std::mutex> lock(mu);
        Slot* done = cur;
        cur = nullptr;
        if (done->failed) {
          // Everything before the failed block is delivered; the serial
          // source takes it from there.
          stopped = true;
          const uint64_t at = done->fail_offset;
          release_ring_locked();
          work_cv.notify_all();
          lock.unlock();
          tail.reset(new ZlibByteSource(path, at, at > 0));
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
            tail.reset(new ZlibByteSource(path, at, at > 0));
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
  unsigned char head[kHeaderBytes] = {0};
  const int64_t got = pread_exact(fd, head, sizeof(head), 0);
  ::close(fd);
  return got == static_cast<int64_t>(sizeof(head)) && bgzf_shaped_header(head);
}

BgzfByteSource::BgzfByteSource(const std::string& path,
                               GzMemberSourceOptions options)
    : impl_(new Impl()) {
  Impl* self = impl_.get();
  self->path = path;
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
  } catch (...) {
    // Nothing may propagate into kseq.
    impl_->broken = true;
    return -1;
  }
}

bool BgzfByteSource::cut_in_data() const {
  return impl_->tail && impl_->tail->cut_in_data();
}

} // namespace io
} // namespace cpu
} // namespace fa
