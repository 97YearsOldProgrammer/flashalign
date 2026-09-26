// Member-parallel gzip source; the output contract is in the header.
//
// Threads:
//   reader   the one caller of read(); owns `held`, `held_off` and `cut_in_data`.
//   scanner  preads the file and publishes candidate member starts (offsets of
//            1f 8b 08). Speculation only: correctness never depends on it.
//   workers  each runs one job at a time: inflating a Segment, repairing a
//            finished piece against its real history (a fix-up), or searching
//            a split point for a block start and opening a piece there. A
//            worker with nothing to do runs one unit of the caller's work
//            through `helper` if one is set, and otherwise sleeps on work_cv.
// Without a helper the scanner has its own thread (team = workers + 1). With
// one, worker 0 scans first and then inflates (folded_scan), so the team is
// exactly `workers`; a one-thread loaned team does not scan.
//
// Pieces:
//   A Segment is a piece of one member's deflate stream, keyed by the bit
//   where it starts. A member head starts at a gzip header. A block piece
//   starts at a candidate block boundary the finder found near a split point
//   (one every kPieceBytes of input). A piece runs until it lands (stops at
//   exactly the bit where another piece starts) or reaches its member's end,
//   an error or end of file. The chain is `current`, the piece being
//   delivered, and every piece reached from it by landings; only the chain is
//   known to have started at a true boundary.
//
//   A block piece started before the chain reached it has no history, so it
//   runs against the 0xFF stand-in window and is sentinel: its bytes are not
//   deliverable until a fix-up re-inflates its prefix up to the last 0xFF
//   against the real window. Member heads, and pieces started once the chain
//   has reached them, are real from the start.
//
//   Pieces cost CPU a straight inflate does not, so they are opened only while
//   the reader is short of bytes (want_speculation_locked), and only ahead of
//   every real producer, since a real producer stops only where a piece is
//   open. Otherwise one real producer inflates the chain front to back.
//
// One mutex, mu, guards all shared state except the page-cache advice, which
// has its own (advise_mu) because chain jobs on several threads issue it. The
// only lock-free fields are three hints: chain_input_pos, piece_epoch and
// stopping. These Segment fields are not guarded but owned by the worker that
// set the segment kRunning or `fixing`: zs, zs_ready, last_boundary, out_len,
// last_ff, crc, trailer and trailer_have. Ownership changes only under mu,
// which orders each owner's writes before the next owner and the reader
// (which reads out_len, crc and trailer only of a piece that has ended).
// in_offset is the owner's too, but it is written back under mu at every
// hand-off, so others may read it under mu as a lower bound. A fix-up moves
// the piece's chunks out of the Segment under mu and back under mu; in
// between only the fixing worker touches them, so it rewrites them in place
// without mu.
//
// Condition variables, all on mu:
//   work_cv   workers wait for something to claim.
//   data_cv   the reader waits for the current piece to gain deliverable bytes
//             or change state.
//   space_cv  the current piece's producer waits for the reader to free room.
//
// Liveness:
//   1. The chain is claimed before any speculation: first the head (to inflate
//      or fix it up), then the rest of the chain.
//   2. A speculative producer never blocks: at its cap or the speculative
//      budget it suspends, keeping its inflate state in the Segment. A piece
//      that is not the head also suspends at every chunk where the head needs
//      a worker. So every worker is back in claim_locked within one chunk or
//      one short job.
//   3. spec_cap_locked() keeps at least one live inflate thread out of member
//      speculation, so one is always running the chain or free to take it;
//      the reader wakes work_cv before it waits. In the folded team the
//      scanning thread does not inflate until scan_done, so the cap is one
//      lower until then.
//   4. Without a helper the head's producer blocks at its cap only on the
//      reader, which is draining that piece. With a helper it parks instead,
//      and claim_locked holds the head back until read() frees room and wakes
//      work_cv. A sentinel head is claimable at any size: its owner's first
//      act is the fix-up, which allocates nothing and makes every byte
//      deliverable.
//   5. The head's history is always known: it is the tail of the piece the
//      reader just finished, captured when that piece landed or rebuilt by its
//      fix-up. So a sentinel head's fix-up is always claimable.
//   6. A thread inside `helper` holds no lock and owns no Segment, so decoding
//      (and teardown) waits at most for one unit of borrowed work.
//   7. Teardown sets abort under mu and wakes all three condition variables;
//      every wait checks abort, inflate checks it at every chunk, and a fix-up
//      and the finder check `stopping` at every read.
#include "gz_member_source.h"

#include "deflate_split.h"
#include "threading/parallel_for.h" // default_thread_count
#include "threading/thread_name.h"

#include <zlib.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <new>
#include <stdexcept>
#include <thread>
#include <vector>

namespace fa {
namespace cpu {
namespace io {
// Named rather than anonymous because Impl's fields use these types
// (-Wsubobject-linkage).
namespace gz_member_detail {

namespace ds = ::fa::cpu::io::deflate_split;

// Compressed bytes per pread, in a per-worker buffer. Input left unconsumed at
// suspension is re-read, since a resume starts from Segment::in_offset.
constexpr size_t kInputBufferBytes = 1u << 20;

// The chain producer's readahead past each pread, and how far behind it
// consumed pages are dropped, in steps. Without the drop a large input evicts
// the memory-mapped index from the page cache. The lag covers a parked
// producer resuming from its in_offset, within one input buffer of its last
// read, and every piece the chain has reached, within the piece lookahead.
constexpr uint64_t kInputReadaheadBytes = 32ull << 20;
constexpr uint64_t kInputDropLagBytes = 64ull << 20;
constexpr uint64_t kInputDropStepBytes = 64ull << 20;
// How far the member scan may run ahead of the chain producer. Unbounded, it
// would read a large file a second time and flood the page cache; speculation
// only needs candidates within about the staging budget.
constexpr uint64_t kScanAheadBytes = 2ull << 30;
// Candidate-scan window. The last two bytes are carried into the next window
// so a 1f 8b 08 straddling the boundary is still found exactly once.
constexpr size_t kScanBufferBytes = 4u << 20;
// A chunk goes to zlib as a uInt, and abort and back-pressure are checked per
// chunk.
constexpr size_t kMaxChunkBytes = 64u << 20;
// Cap on the candidate list, against inputs full of false 1f 8b 08 hits.
// Dropping candidates costs only parallelism.
constexpr size_t kMaxCandidates = 1u << 20;

// Compressed bytes between split points. The unknown window stops propagating
// after a roughly fixed amount of output, so smaller pieces repair a larger
// share of their output, and larger ones hold more memory per piece in flight.
constexpr uint64_t kPieceBytes = 4ull << 20;
// A piece inflates with Z_NO_FLUSH until it has handed zlib all but this many
// bytes before its target's byte, so every boundary it crosses meanwhile lies
// before the target. Only the rest steps block by block with Z_BLOCK, since
// each Z_BLOCK return costs zlib a 32 KiB window copy. Any guard of one byte
// or more is exact.
constexpr uint64_t kPieceGuardBytes = 16;
// The finder tries the byte-aligned positions of a split's first stretch
// before it scans every bit: pigz ends a block on a byte boundary every
// 128 KiB of output, so an aligned candidate usually comes early.
constexpr size_t kFinderAlignedBytes = 256u << 10;
// The finder reads a split in windows of this size, plus a slack holding the
// header and block of a candidate near the window's end, which zlib inflates
// to confirm it.
constexpr size_t kFinderWindowBytes = 1u << 20;
constexpr size_t kFinderSlackBytes = 256u << 10;
// Pieces opened ahead of the chain's input position, per worker: two keep
// every worker busy while the chain fixes up and delivers the one before.
constexpr unsigned kPiecesPerWorker = 2;
// Pieces are opened only while the reader's lead is below staging_budget / 4.
// The mark must cover the time speculation takes to deliver its first piece,
// and sit well below what the sequential producer may stage (the chain cap is
// at least half the budget), or the lead could never climb past it and
// speculation would never stop.
constexpr uint64_t kSpeculationWatermarkDivisor = 4;

constexpr uint64_t kNoTarget = ~0ull;

inline uint64_t guard_byte(uint64_t bit) {
  const uint64_t byte = bit >> 3;
  return byte > kPieceGuardBytes ? byte - kPieceGuardBytes : 0;
}

inline uint32_t read_le32(const unsigned char* p) {
  return static_cast<uint32_t>(p[0]) | static_cast<uint32_t>(p[1]) << 8 |
         static_cast<uint32_t>(p[2]) << 16 | static_cast<uint32_t>(p[3]) << 24;
}

// The last `value` in p[0, n), or null.
inline const unsigned char* find_last(const unsigned char* p, size_t n,
                                      unsigned char value) {
#if defined(__GLIBC__)
  return static_cast<const unsigned char*>(memrchr(p, value, n));
#else
  for (size_t i = n; i-- > 0;) {
    if (p[i] == value)
      return p + i;
  }
  return nullptr;
#endif
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

// A BGZF member is a gzip member with an FEXTRA 'BC' subfield.
bool looks_like_bgzf(const unsigned char* head, size_t len) {
  if (len < 14)
    return false;
  if ((head[3] & 0x04) == 0)
    return false; // FLG.FEXTRA
  return head[12] == 'B' && head[13] == 'C';
}

enum class SegState {
  kPending,   // unowned: never started, or suspended with inflate state kept
  kRunning,   // exactly one worker owns zs and the owner-only fields
  kLanded,    // stopped at exactly `next`, where another piece starts
  kDone,      // its member's deflate stream ended; `trailer` holds what follows
  kTruncated, // input ran out mid-member: gzread's "unexpected end of file"
  kFailed,    // zlib data error, allocation failure, or a read error
};

inline bool has_ended(SegState state) {
  return state == SegState::kLanded || state == SegState::kDone ||
         state == SegState::kTruncated || state == SegState::kFailed;
}

struct Chunk {
  std::unique_ptr<char[]> data;
  size_t size = 0;
};

// Up to 32 KiB of real output: the history a piece needs from before its
// start, or the tail a landed piece hands to the piece it landed on. Shared,
// because a fix-up reads its piece's history without mu while the reader may
// be retiring the piece that produced it.
using Window = std::vector<unsigned char>;
using WindowPtr = std::shared_ptr<const Window>;

// One piece of one member's deflate stream.
struct Segment {
  Segment(uint64_t start_bit, bool is_member_head)
      : start(start_bit), member_head(is_member_head) {}
  ~Segment() {
    if (zs_ready)
      inflateEnd(&zs);
  }
  Segment(const Segment&) = delete;
  Segment& operator=(const Segment&) = delete;

  const uint64_t start;    // file bit offset this piece starts at
  const bool member_head;  // starts at a gzip header, byte aligned

  // Owned by the worker that set kRunning or `fixing`.
  z_stream zs{};
  bool zs_ready = false;
  uint64_t in_offset = 0;     // next compressed byte zlib takes, once zs_ready
  uint64_t last_boundary = 0; // last block boundary this piece stepped onto
  uint64_t out_len = 0;       // bytes of output produced
  int64_t last_ff = -1;       // sentinel only: last 0xFF in the output
  // CRC-32 of the output after last_ff (all of it when last_ff < 0), so a
  // fix-up only checksums the prefix it rewrites.
  uint32_t crc = 0;
  unsigned char trailer[8] = {0}; // kDone: the bytes after the deflate stream
  unsigned trailer_have = 0;      // ...how many of them the file holds

  // mu-guarded.
  SegState state = SegState::kPending;
  bool counted_speculative = false; // this segment holds a spec_running slot
  bool discarded = false;           // contradicted; its worker must stop
  bool in_chain = false;            // reached from `current` by landings
  bool sentinel = false; // output so far used the stand-in window
  bool fixing = false;   // a worker holds the chunks for a fix-up
  uint64_t next = 0;     // kLanded: the bit where the next piece starts
  uint64_t end = 0;      // kDone: first byte after the trailer; else file size
  WindowPtr history;     // the real window before `start`, once known
  WindowPtr tail;        // kLanded: the real window at `next`, once known
  std::deque<Chunk> chunks;
  uint64_t staged = 0; // sum of chunks[].size, i.e. undelivered bytes
};

using SegPtr = std::shared_ptr<Segment>;

// One split point: the candidate block start the finder found in
// [j * kPieceBytes, (j + 1) * kPieceBytes). Searched at most once, by the
// first worker that needs it; others wait in call_once.
struct Split {
  std::once_flag once;
  std::atomic<uint64_t> bit{ds::kNotFound};
};

// A pread-filled window onto the compressed file, addressed by file offset.
struct InputCursor {
  explicit InputCursor(std::vector<unsigned char>& storage)
      : buf(storage.data()), cap(storage.size()) {}
  unsigned char* buf;
  size_t cap;
  uint64_t base = 0; // file offset of buf[0]
  size_t have = 0;   // valid bytes in buf
  size_t pos = 0;    // next byte to hand to zlib
  uint64_t offset() const { return base + pos; }
  void seek(uint64_t at) {
    base = at;
    have = 0;
    pos = 0;
  }
};

struct WorkerScratch {
  std::vector<unsigned char> input;   // the running piece's compressed input
  std::vector<unsigned char> aux;     // finder windows, fix-up input
  std::vector<unsigned char> confirm; // the finder's one-block inflate output
};

enum class JobKind { kNone, kRun, kFix, kOpen };

struct Job {
  JobKind kind = JobKind::kNone;
  SegPtr seg;
  size_t split = 0;
};

// How far a fix-up's re-inflate got.
struct Reinflate {
  bool ok = true;
  uint64_t produced = 0; // bytes rewritten; all of them correct even if !ok
  uint32_t crc = 0;      // CRC-32 of those bytes
  uint64_t in_offset = 0;
};

} // namespace gz_member_detail

using namespace gz_member_detail;

struct MemberParallelGzByteSource::Impl {
  ~Impl() {
    shutdown();
    if (fd >= 0)
      ::close(fd);
  }

  // ---- immutable after construction ----
  std::string path; // for the truncation warning
  int fd = -1;
  uint64_t size = 0;
  unsigned workers = 1;
  unsigned spec_worker_cap = 0; // workers - 1, or 0 when single-threaded
  size_t chunk_bytes = 0;
  uint64_t segment_cap = 0;
  uint64_t spec_budget = 0;
  uint64_t staging_budget = 0; // segment_cap + spec_budget, the whole allowance
  // Split points 1..split_count; 0 means the file is never cut (one worker,
  // or no longer than one piece). splits[0] is unused.
  size_t split_count = 0;
  std::unique_ptr<Split[]> splits;
  size_t piece_lookahead = 0; // kPiecesPerWorker * workers
  uint64_t spec_watermark = 0; // staging_budget / kSpeculationWatermarkDivisor
  // Worker 0 scans before inflating (loaned team only).
  bool folded_scan = false;
  // The caller's work, run by idle workers; empty for a private team.
  std::function<bool(unsigned)> helper;

  // ---- page-cache advice (advise_mu) ----
  std::mutex advise_mu;
  // Page cache below this offset has been dropped.
  uint64_t input_dropped_to = 0;
  // The furthest input position a chain job has reported. Written under
  // advise_mu and read without it: the scanner's throttle and the piece
  // lookahead need it only roughly.
  std::atomic<uint64_t> chain_input_pos{0};
  // Bumped whenever a piece is opened, so a producer re-reads the open pieces
  // (under mu) only when there is something new.
  std::atomic<uint64_t> piece_epoch{0};
  // Set by shutdown(); checked by work that runs without mu for longer than a
  // chunk (a fix-up, the finder).
  std::atomic<bool> stopping{false};

  // ---- reader-owned ----
  Chunk held;          // one chunk detached from the current piece
  size_t held_off = 0; // bytes of `held` already delivered
  bool cut_in_data = false; // EOF came inside a member's deflate data or header

  // ---- mu-guarded ----
  std::mutex mu;
  std::condition_variable work_cv;
  std::condition_variable data_cv;
  std::condition_variable space_cv;

  bool abort = false;
  bool eof = false;
  bool failed = false;
  uint64_t frontier = 0; // start byte of the member being delivered
  SegPtr current;        // the piece being delivered, once magic is confirmed
  std::map<uint64_t, SegPtr> segments; // by start bit, ascending
  std::vector<uint64_t> candidates;    // ascending 1f 8b 08 offsets
  size_t next_candidate = 0;
  bool scan_done = false; // no further candidates will be published
  uint64_t total_staged = 0;
  unsigned spec_running = 0;
  std::vector<std::unique_ptr<char[]>> recycle; // spare chunk buffers
  size_t open_cursor = 1;      // next split point to consider opening
  bool reader_waiting = false; // read() is blocked on data_cv
  // CRC-32 and length of the current member's pieces the reader has
  // finished, for the trailer check.
  uint32_t member_crc = 0;
  uint64_t member_len = 0;

  std::vector<std::thread> threads; // the whole team; [0] scans first

  // ---- helpers, all called with mu held ----

  // Staged bytes outside the current piece. It may over-report, which only
  // stops speculation sooner: `held` belongs to no segment, and `current` is
  // briefly null while attach_current_locked reads the magic.
  uint64_t speculative_staged() const {
    const uint64_t in_current = current ? current->staged : 0;
    return total_staged >= in_current ? total_staged - in_current : 0;
  }

  // Stage cap of the current piece: whatever of staging_budget speculation is
  // not using, and at least segment_cap. The larger runway matters under a
  // loan, where a parked producer may be away for a whole unit of borrowed
  // work.
  uint64_t chain_cap_locked() const {
    const uint64_t spec = speculative_staged();
    const uint64_t room = staging_budget > spec ? staging_budget - spec : 0;
    return std::max(segment_cap, room);
  }

  // Called by a chain job after a pread ending at `pos`. Speculative producers
  // do not call it: the chain still needs what they read. Several chain jobs
  // can run at once (the head's producer and a fix-up further on), so the
  // position only moves forward.
  void advise_chain_input(uint64_t pos) {
    std::lock_guard<std::mutex> guard(advise_mu);
    if (pos <= chain_input_pos.load(std::memory_order_relaxed))
      return;
    chain_input_pos.store(pos, std::memory_order_relaxed);
#if defined(POSIX_FADV_WILLNEED) && defined(POSIX_FADV_DONTNEED)
    if (pos < size) {
      const uint64_t ahead = std::min<uint64_t>(kInputReadaheadBytes, size - pos);
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
#endif
  }

  // Threads allowed in member speculation at once: one fewer than the live
  // inflate team (liveness point 3 above).
  unsigned spec_cap_locked() const {
    if (!folded_scan || scan_done)
      return spec_worker_cap;
    return workers > 2 ? workers - 2 : 0;
  }

  // Chunk buffers are all chunk_bytes, so released ones are reused rather than
  // reallocated. The pool holds at most 2 * workers buffers.
  std::unique_ptr<char[]> take_buffer_locked() {
    if (recycle.empty())
      return nullptr;
    std::unique_ptr<char[]> buffer = std::move(recycle.back());
    recycle.pop_back();
    return buffer;
  }
  void put_buffer_locked(std::unique_ptr<char[]> buffer) {
    if (!buffer)
      return;
    // Nothing inflates after a terminal state.
    if (eof || failed) {
      release_spares_locked();
      return;
    }
    if (recycle.size() < static_cast<size_t>(workers) * 2)
      recycle.push_back(std::move(buffer));
  }

  // Frees the pool, so a drained source left open holds no spare buffers.
  void release_spares_locked() {
    recycle.clear();
    recycle.shrink_to_fit();
  }

  // Hands a segment back. kPending means suspended and keeps the inflate
  // state for a later worker; any other state is terminal and frees it. A
  // failed speculation's chunks stay staged until the chain passes it.
  void release_locked(const SegPtr& seg, SegState state) {
    seg->state = state;
    if (state != SegState::kPending && seg->zs_ready) {
      inflateEnd(&seg->zs);
      seg->zs_ready = false;
    }
    if (seg->counted_speculative) {
      seg->counted_speculative = false;
      --spec_running;
    }
    data_cv.notify_all();
    work_cv.notify_all();
  }

  // Drops a contradicted speculation's output and tells its worker to stop;
  // the worker gives up its spec_running slot when it releases. The caller
  // must erase the segment from `segments` at once: its in_offset may count
  // input whose output was dropped, so it must never become `current`.
  void discard_locked(const SegPtr& seg) {
    seg->discarded = true;
    total_staged -= seg->staged;
    seg->staged = 0;
    for (Chunk& chunk : seg->chunks)
      put_buffer_locked(std::move(chunk.data));
    seg->chunks.clear();
  }

  // Discards every segment strictly between two bits: what the chain stepped
  // over.
  void discard_range_locked(uint64_t after_bit, uint64_t before_bit) {
    for (auto it = segments.upper_bound(after_bit);
         it != segments.end() && it->first < before_bit;) {
      discard_locked(it->second);
      it = segments.erase(it);
    }
  }

  // The piece `p` landed on, if it has landed and that piece exists.
  SegPtr chain_next_locked(const Segment& p) const {
    if (p.state != SegState::kLanded || p.discarded)
      return nullptr;
    auto it = segments.find(p.next);
    return it == segments.end() ? nullptr : it->second;
  }

  // A finished sentinel piece whose history is known: ready for a fix-up.
  static bool fix_ready(const Segment& s) {
    return s.sentinel && !s.fixing && !s.discarded && s.history &&
           has_ended(s.state);
  }

  // Whether the head has work nobody is doing. Other pieces yield to it at
  // their next chunk.
  bool head_needs_worker_locked() const {
    if (!current || current->discarded)
      return false;
    const Segment& head = *current;
    if (head.state == SegState::kPending)
      return head.sentinel || !helper || head.staged < chain_cap_locked();
    return fix_ready(head);
  }

  // Deliverable bytes staged ahead of the reader: the head and every piece
  // after it, while each is final and the one before it landed on it.
  uint64_t chain_lead_locked() const {
    uint64_t lead = 0;
    const Segment* p = current.get();
    while (p != nullptr && !p->sentinel && !p->discarded) {
      lead += p->staged;
      const SegPtr q = chain_next_locked(*p);
      p = q.get();
    }
    return lead;
  }

  // Pieces are opened (and speculative ones resumed) only while the reader is
  // short of bytes: waiting, or with less than spec_watermark staged ahead.
  // Otherwise one worker inflating the chain keeps up, and a piece would spend
  // the stand-in pass and the fix-up on bytes nobody is waiting for.
  bool want_speculation_locked() const {
    if (split_count == 0)
      return false;
    return reader_waiting || chain_lead_locked() < spec_watermark;
  }

  // The last bit the chain has proven anything about: a block piece that
  // starts at or before it is either in the chain or was stepped over.
  uint64_t chain_floor_locked() const {
    uint64_t floor = frontier * 8;
    for (const Segment* p = current.get(); p != nullptr;) {
      floor = std::max(floor, p->start);
      if (p->state == SegState::kLanded)
        floor = std::max(floor, p->next);
      const SegPtr q = chain_next_locked(*p);
      p = q.get();
    }
    return std::max(floor,
                    chain_input_pos.load(std::memory_order_relaxed) * 8);
  }

  // Follows the landings forward from the head. The piece a chained piece
  // landed on joins the chain (created unstarted if nobody opened it); every
  // segment strictly between the two was a false start or was stepped over;
  // and a landed piece's tail is the next one's history. Idempotent; called
  // whenever a piece lands, a fix-up finishes or the head moves.
  void propagate_locked() {
    bool changed = false;
    for (Segment* p = current.get();
         p != nullptr && p->state == SegState::kLanded && !p->discarded;) {
      auto it = segments.find(p->next);
      if (it == segments.end()) {
        it = segments
                 .emplace(p->next, std::make_shared<Segment>(p->next, false))
                 .first;
        changed = true;
      }
      Segment& q = *it->second;
      if (!q.in_chain) {
        q.in_chain = true;
        discard_range_locked(p->start, q.start);
        changed = true;
      }
      if (!q.history && p->tail) {
        q.history = p->tail;
        changed = true;
      }
      p = &q;
    }
    if (changed) {
      work_cv.notify_all();
      data_cv.notify_all();
    }
  }

  // Whether a real-window producer has inflated past compressed byte `at`, or
  // will within one read, without stopping there. Such a producer stops only
  // where a piece is open, so a piece opened behind it would never be reached.
  bool covered_by_real_producer_locked(uint64_t at) const {
    for (const auto& entry : segments) {
      const Segment& p = *entry.second;
      if (p.discarded || p.sentinel || (p.start >> 3) >= at)
        continue;
      const uint64_t reach = p.state == SegState::kLanded
                                 ? p.next >> 3
                                 : p.in_offset + kInputBufferBytes;
      if (at < reach)
        return true;
    }
    return false;
  }

  // The next split point worth opening: at least one read beyond the chain's
  // input position (closer, the chain producer may already be past it),
  // within the lookahead of it, and not covered by a real-window producer.
  bool next_split_to_open_locked(size_t& split) {
    if (split_count == 0)
      return false;
    const uint64_t base =
        std::max(frontier, chain_input_pos.load(std::memory_order_relaxed));
    const size_t first = static_cast<size_t>(
        (base + kInputBufferBytes + kPieceBytes - 1) / kPieceBytes);
    size_t j = std::max<size_t>(open_cursor, std::max<size_t>(first, 1));
    const size_t limit = std::min(
        split_count, static_cast<size_t>(base / kPieceBytes) + piece_lookahead);
    while (j <= limit && covered_by_real_producer_locked(j * kPieceBytes))
      ++j;
    if (j > limit)
      return false;
    open_cursor = j + 1;
    split = j;
    return true;
  }

  // Highest priority first: the head (inflate it, or fix it up), then the rest
  // of the chain, then, only while the reader is short of bytes, pieces ahead
  // of it, then speculation on member starts.
  Job claim_locked() {
    Job job;
    if (abort)
      return job;
    if (current && !current->discarded) {
      Segment& head = *current;
      if (head.state == SegState::kPending) {
        // Under a loan the head may be parked at its cap; hold it back until
        // the reader frees room, unless it is sentinel (liveness point 4).
        if (head.sentinel || !helper || head.staged < chain_cap_locked()) {
          head.state = SegState::kRunning;
          job.kind = JobKind::kRun;
          job.seg = current;
          return job;
        }
      } else if (fix_ready(head)) {
        head.fixing = true;
        job.kind = JobKind::kFix;
        job.seg = current;
        return job;
      }
      // A chained piece with no history yet can only start against the
      // stand-in, so it waits for its history unless the reader is short.
      const bool speculate = want_speculation_locked();
      for (SegPtr q = chain_next_locked(head); q; q = chain_next_locked(*q)) {
        if (fix_ready(*q)) {
          q->fixing = true;
          job.kind = JobKind::kFix;
          job.seg = q;
          return job;
        }
        if (q->state == SegState::kPending && (q->history || speculate) &&
            q->staged < segment_cap && speculative_staged() < spec_budget) {
          q->state = SegState::kRunning;
          job.kind = JobKind::kRun;
          job.seg = q;
          return job;
        }
      }
    }
    if (eof || failed)
      return job;
    if (want_speculation_locked() && speculative_staged() < spec_budget) {
      const uint64_t floor = chain_floor_locked();
      for (auto& entry : segments) {
        Segment& q = *entry.second;
        if (q.member_head || q.in_chain || q.discarded ||
            q.state != SegState::kPending || entry.first <= floor ||
            q.staged >= segment_cap)
          continue;
        q.state = SegState::kRunning;
        job.kind = JobKind::kRun;
        job.seg = entry.second;
        return job;
      }
      size_t split = 0;
      if (next_split_to_open_locked(split)) {
        job.kind = JobKind::kOpen;
        job.split = split;
        return job;
      }
    }
    if (spec_running >= spec_cap_locked())
      return job;
    if (speculative_staged() >= spec_budget)
      return job;
    for (auto& entry : segments) {
      if (entry.first <= frontier * 8)
        continue;
      const SegPtr& seg = entry.second;
      if (!seg->member_head || seg->state != SegState::kPending ||
          seg->discarded)
        continue;
      if (seg->staged >= segment_cap)
        continue;
      seg->state = SegState::kRunning;
      seg->counted_speculative = true;
      ++spec_running;
      job.kind = JobKind::kRun;
      job.seg = seg;
      return job;
    }
    while (next_candidate < candidates.size()) {
      const uint64_t offset = candidates[next_candidate++];
      if (offset <= frontier || segments.count(offset * 8) != 0)
        continue;
      SegPtr seg = std::make_shared<Segment>(offset * 8, true);
      seg->state = SegState::kRunning;
      seg->counted_speculative = true;
      ++spec_running;
      segments.emplace(offset * 8, seg);
      job.kind = JobKind::kRun;
      job.seg = seg;
      return job;
    }
    return job;
  }

  // Attaches the member at `frontier`, or sets EOF: as zlib's gz_look(), the
  // stream continues only if the next two bytes are the gzip magic.
  void attach_current_locked(std::unique_lock<std::mutex>& lock) {
    if (frontier + 2 > size) {
      eof = true;
      return;
    }
    unsigned char magic[2] = {0, 0};
    lock.unlock();
    const int64_t got = pread_exact(fd, magic, sizeof(magic), frontier);
    lock.lock();
    if (got < 0) {
      failed = true;
      return;
    }
    if (got < 2 || magic[0] != 0x1f || magic[1] != 0x8b) {
      eof = true;
      return;
    }
    // No block piece starts at a member start: the finder takes only
    // BFINAL = 0 headers, and the magic's 1f reads as BFINAL = 1.
    auto entry = segments.find(frontier * 8);
    if (entry == segments.end()) {
      entry = segments
                  .emplace(frontier * 8,
                           std::make_shared<Segment>(frontier * 8, true))
                  .first;
    }
    current = entry->second;
    current->in_chain = true;
    // An adopted speculation leaves the speculative pool.
    if (current->counted_speculative) {
      current->counted_speculative = false;
      --spec_running;
    }
    propagate_locked();
    work_cv.notify_all();
    space_cv.notify_all();
  }

  // The head landed and the reader has taken all of it: the piece it landed
  // on becomes the head.
  void advance_piece_locked() {
    const SegPtr done = current;
    member_crc = static_cast<uint32_t>(
        crc32_combine(member_crc, done->crc, static_cast<z_off_t>(done->out_len)));
    member_len += done->out_len;
    propagate_locked(); // guarantees the piece at done->next exists
    const SegPtr next = segments.find(done->next)->second;
    segments.erase(done->start);
    current = next;
    work_cv.notify_all();
    space_cv.notify_all();
  }

  enum class MemberVerdict { kIntact, kTruncated, kCorrupt };

  // The head reached its member's end. zlib's checks in zlib's order: the
  // CRC-32 of the member's whole output (all its pieces combined), then ISIZE.
  // A trailer cut short by end of file is "unexpected end of file", which
  // gzread treats as EOF, unless the part that is there is already wrong.
  MemberVerdict member_verdict_locked() const {
    const Segment& last = *current;
    const uint32_t crc = static_cast<uint32_t>(crc32_combine(
        member_crc, last.crc, static_cast<z_off_t>(last.out_len)));
    const uint64_t len = member_len + last.out_len;
    if (last.trailer_have < 4)
      return MemberVerdict::kTruncated;
    if (read_le32(last.trailer) != crc)
      return MemberVerdict::kCorrupt;
    if (last.trailer_have < 8)
      return MemberVerdict::kTruncated;
    if (read_le32(last.trailer + 4) != static_cast<uint32_t>(len))
      return MemberVerdict::kCorrupt;
    return MemberVerdict::kIntact;
  }

  // The current member ended at a validated offset; every segment that
  // started inside it was a false candidate.
  void advance_member_locked(std::unique_lock<std::mutex>& lock) {
    const uint64_t member_end = current->end;
    segments.erase(current->start);
    current.reset();
    frontier = member_end;
    member_crc = 0;
    member_len = 0;
    for (auto it = segments.begin(); it != segments.end();) {
      if (it->first >= frontier * 8)
        break;
      discard_locked(it->second);
      it = segments.erase(it);
    }
    work_cv.notify_all();
    space_cv.notify_all();
    attach_current_locked(lock);
  }

  int read(char* out, size_t want) {
    // The detached chunk is reader-private, so serving it needs no lock.
    if (held_off < held.size) {
      const size_t take = std::min(want, held.size - held_off);
      std::memcpy(out, held.data.get() + held_off, take);
      held_off += take;
      return static_cast<int>(take);
    }
    std::unique_lock<std::mutex> lock(mu);
    for (;;) {
      if (failed)
        return -1;
      if (!current && !eof)
        attach_current_locked(lock);
      // A sentinel piece's bytes are not final until its fix-up.
      if (current && !current->sentinel && !current->chunks.empty()) {
        // The held chunk stays charged until it is replaced.
        if (held.data) {
          put_buffer_locked(std::move(held.data));
          total_staged -= chunk_bytes;
        }
        held = std::move(current->chunks.front());
        current->chunks.pop_front();
        held_off = 0;
        current->staged -= chunk_bytes;
        space_cv.notify_all();
        // Under a loan a parked head producer waits on work_cv instead.
        if (helper)
          work_cv.notify_all();
        lock.unlock();
        const size_t take = std::min(want, held.size);
        std::memcpy(out, held.data.get(), take);
        held_off = take;
        return static_cast<int>(take);
      }
      if (failed) {
        release_spares_locked();
        return -1;
      }
      if (eof) {
        release_spares_locked();
        return 0;
      }
      // attach_current_locked() sets current, eof or failed; a broken
      // invariant becomes a read error rather than a crash.
      if (!current) {
        failed = true;
        release_spares_locked();
        return -1;
      }
      if (!current->sentinel) {
        switch (current->state) {
        case SegState::kLanded:
          advance_piece_locked();
          continue;
        case SegState::kDone:
          switch (member_verdict_locked()) {
          case MemberVerdict::kIntact:
            advance_member_locked(lock);
            continue;
          case MemberVerdict::kTruncated:
            eof = true;
            warn_truncated_gzip(path);
            continue;
          case MemberVerdict::kCorrupt:
            failed = true;
            release_spares_locked();
            return -1;
          }
          continue;
        case SegState::kTruncated:
          // As gzread: the decodable prefix was delivered; now EOF.
          eof = true;
          cut_in_data = true;
          warn_truncated_gzip(path);
          continue;
        case SegState::kFailed:
          failed = true;
          release_spares_locked();
          return -1;
        default:
          break;
        }
      }
      // Unowned, running with nothing staged yet, or a finished sentinel piece
      // awaiting its fix-up. Waiting also turns speculation on.
      reader_waiting = true;
      work_cv.notify_all();
      data_cv.wait(lock);
      reader_waiting = false;
    }
  }

  // Returns false once the candidate list is full, which stops the scanner.
  bool publish_candidates(const std::vector<uint64_t>& found) {
    std::lock_guard<std::mutex> guard(mu);
    candidates.insert(candidates.end(), found.begin(), found.end());
    work_cv.notify_all();
    return candidates.size() < kMaxCandidates;
  }

  void scanner_main(unsigned ordinal) {
    std::vector<unsigned char> buffer(kScanBufferBytes);
    std::vector<uint64_t> found;
    uint64_t next = 0;   // next file offset to read
    uint64_t window = 0; // absolute offset of buffer[0]
    size_t carry = 0;    // bytes carried from the previous window
    while (next < size) {
      {
        // Stay within kScanAheadBytes of the chain producer; a loaned scanner
        // runs borrowed work meanwhile. The position is published without mu,
        // so the wait is timed.
        std::unique_lock<std::mutex> lock(mu);
        for (;;) {
          if (abort)
            return;
          if (next <= chain_input_pos.load(std::memory_order_relaxed) +
                          kScanAheadBytes)
            break;
          if (helper) {
            lock.unlock();
            const bool ran = helper(ordinal);
            lock.lock();
            if (ran)
              continue;
          }
          space_cv.wait_for(lock, std::chrono::milliseconds(20));
        }
      }
      const size_t want = static_cast<size_t>(
          std::min<uint64_t>(size - next, buffer.size() - carry));
      const int64_t got = pread_exact(fd, buffer.data() + carry, want, next);
      if (got <= 0)
        break; // read error or EOF: speculation simply stops
      const size_t have = carry + static_cast<size_t>(got);
      next += static_cast<uint64_t>(got);
      found.clear();
      for (size_t i = 0; i + 3 <= have;) {
        const void* hit = std::memchr(buffer.data() + i, 0x1f, have - i - 2);
        if (hit == nullptr)
          break;
        const size_t at = static_cast<size_t>(
            static_cast<const unsigned char*>(hit) - buffer.data());
        if (buffer[at + 1] == 0x8b && buffer[at + 2] == 0x08) {
          found.push_back(window + at);
        }
        i = at + 1;
      }
      if (!found.empty() && !publish_candidates(found))
        return;
      carry = have >= 2 ? 2 : have;
      std::memmove(buffer.data(), buffer.data() + have - carry, carry);
      window = next - carry;
    }
  }

  void search_split(size_t j, Split& split, WorkerScratch& sc) {
    const uint64_t lo = static_cast<uint64_t>(j) * kPieceBytes;
    const uint64_t hi = std::min<uint64_t>(size, lo + kPieceBytes);
    const uint64_t aligned_end = std::min<uint64_t>(hi, lo + kFinderAlignedBytes);
    uint64_t found = ds::kNotFound;
    try {
      sc.aux.resize(kFinderWindowBytes + kFinderSlackBytes);
      // Pass 0 tries the byte-aligned positions of the first stretch; pass 1
      // scans every bit of the split, skipping what pass 0 tried.
      for (int pass = 0; pass < 2 && found == ds::kNotFound; ++pass) {
        const uint64_t pass_end = pass == 0 ? aligned_end : hi;
        for (uint64_t w = lo; w < pass_end && found == ds::kNotFound;
             w += kFinderWindowBytes) {
          if (stopping.load(std::memory_order_relaxed))
            break;
          const uint64_t scan_end = std::min<uint64_t>(pass_end, w + kFinderWindowBytes);
          const uint64_t read_end = std::min<uint64_t>(size, scan_end + kFinderSlackBytes);
          const int64_t got = pread_exact(fd, sc.aux.data(),
                                          static_cast<size_t>(read_end - w), w);
          if (got < static_cast<int64_t>(scan_end - w))
            break; // read error: no split point here
          const size_t skip =
              pass == 1 && aligned_end > w ? static_cast<size_t>(aligned_end - w) : 0;
          const uint64_t rel = ds::find_block_start(
              sc.aux.data(), static_cast<size_t>(got), 0,
              static_cast<size_t>(scan_end - w), skip, pass == 0, sc.confirm);
          if (rel != ds::kNotFound)
            found = w * 8 + rel;
        }
      }
    } catch (const std::bad_alloc&) {
      found = ds::kNotFound; // no memory for the search: no split point here
    }
    split.bit.store(found, std::memory_order_release);
  }

  uint64_t find_split(size_t j, WorkerScratch& sc) {
    Split& split = splits[j];
    std::call_once(split.once, [&] { search_split(j, split, sc); });
    return split.bit.load(std::memory_order_acquire);
  }

  // The first candidate block start after `bit`, searching split points as
  // needed: where a sentinel piece the chain has not reached stops, whether
  // or not anyone has opened a piece there.
  uint64_t next_split_start_after(uint64_t bit, WorkerScratch& sc) {
    size_t j = std::max<size_t>(
        1, static_cast<size_t>((bit >> 3) / kPieceBytes));
    for (; j <= split_count; ++j) {
      if (stopping.load(std::memory_order_relaxed))
        return kNoTarget;
      const uint64_t found = find_split(j, sc);
      if (found != ds::kNotFound && found > bit)
        return found;
    }
    return kNoTarget;
  }

  // Where a real producer, or a piece the chain has reached, stops: the first
  // open block piece ahead of it, so with none open it inflates straight on.
  // Ahead means past its last boundary and past what zlib may already have
  // decoded from the input it was handed (well under kHeldBits); a piece
  // behind that can no longer be stopped on. For a chain producer such a
  // piece is dead, so it is dropped now rather than left running.
  static constexpr uint64_t kHeldBits = 128;
  uint64_t open_target_locked(const Segment& s, uint64_t consumed_byte,
                              bool chain_role) {
    const uint64_t consumed_bit = consumed_byte * 8;
    const uint64_t after = std::max(
        s.last_boundary, consumed_bit > kHeldBits ? consumed_bit - kHeldBits : 0);
    for (auto it = segments.upper_bound(s.start); it != segments.end();) {
      const Segment& q = *it->second;
      if (q.member_head || q.discarded) {
        ++it;
        continue;
      }
      if (it->first > after)
        return it->first;
      if (chain_role && !q.in_chain) {
        discard_locked(it->second);
        it = segments.erase(it);
        continue;
      }
      ++it;
    }
    return kNoTarget;
  }

  // Refills `in` at its offset. 1: bytes available; 0: end of file; -1: a read
  // error.
  int refill(InputCursor& in) const {
    const uint64_t at = in.offset();
    if (at >= size)
      return 0;
    const size_t want = static_cast<size_t>(std::min<uint64_t>(size - at, in.cap));
    const int64_t got = pread_exact(fd, in.buf, want, at);
    if (got < 0)
      return -1;
    if (got == 0)
      return 0;
    in.base = at;
    in.have = static_cast<size_t>(got);
    in.pos = 0;
    return 1;
  }

  // Points raw `zs` at bit `bit` of the file: the sub-byte bits go in through
  // inflatePrime and `in` is left at the next byte.
  bool prime_raw_at(z_stream& zs, InputCursor& in, uint64_t bit) {
    in.seek(bit >> 3);
    if ((bit & 7) == 0)
      return true;
    if (refill(in) <= 0)
      return false;
    if (!ds::prime_at(zs, in.buf[0], static_cast<unsigned>(bit & 7)))
      return false;
    in.pos = 1;
    return true;
  }

  // A member head: zlib parses the gzip header (windowBits 15 + 16; Z_BLOCK
  // returns right after it), then the deflate data is inflated raw like any
  // other piece, so the member's CRC-32 and ISIZE can be checked across its
  // pieces. Returns kRunning with `s.zs` open and `in` at the data, kTruncated
  // if the file ends inside the header, or kFailed if zlib refuses it.
  SegState open_member_head(Segment& s, InputCursor& in) {
    z_stream header;
    std::memset(&header, 0, sizeof(header));
    if (inflateInit2(&header, 15 + 16) != Z_OK)
      return SegState::kFailed;
    in.seek(s.start >> 3);
    unsigned char sink[16];
    SegState verdict = SegState::kFailed;
    uint64_t data_bit = 0;
    int stalled = 0;
    for (;;) {
      if (in.pos == in.have) {
        const int r = refill(in);
        if (r <= 0) {
          verdict = r < 0 ? SegState::kFailed : SegState::kTruncated;
          break;
        }
      }
      header.next_in = in.buf + in.pos;
      header.avail_in = static_cast<uInt>(in.have - in.pos);
      header.next_out = sink;
      header.avail_out = sizeof(sink);
      const uInt before = header.avail_in;
      const int status = inflate(&header, Z_BLOCK);
      in.pos += before - header.avail_in;
      if (status == Z_OK && (header.data_type & 128)) {
        data_bit = in.offset() * 8 - static_cast<uint64_t>(header.data_type & 7);
        verdict = SegState::kRunning;
        break;
      }
      if (status != Z_OK && status != Z_BUF_ERROR)
        break; // kFailed: zlib refused the header
      if (before == header.avail_in && ++stalled >= 2)
        break;
      if (before != header.avail_in)
        stalled = 0;
    }
    inflateEnd(&header);
    if (verdict != SegState::kRunning)
      return verdict;
    std::memset(&s.zs, 0, sizeof(s.zs));
    if (inflateInit2(&s.zs, -15) != Z_OK)
      return SegState::kFailed;
    s.zs_ready = true;
    if (!prime_raw_at(s.zs, in, data_bit))
      return SegState::kFailed;
    s.last_boundary = data_bit;
    return SegState::kRunning;
  }

  // A block piece: raw inflate from its start bit against its real history,
  // or against the stand-in when `history` is null.
  SegState open_block_piece(Segment& s, InputCursor& in,
                            const Window* history) {
    std::memset(&s.zs, 0, sizeof(s.zs));
    if (inflateInit2(&s.zs, -15) != Z_OK)
      return SegState::kFailed;
    s.zs_ready = true;
    const unsigned char* dict =
        history != nullptr ? history->data() : ds::sentinel_window();
    const size_t dict_len =
        history != nullptr ? history->size() : ds::kWindowBytes;
    if (dict_len > 0 &&
        inflateSetDictionary(&s.zs, dict, static_cast<uInt>(dict_len)) != Z_OK)
      return SegState::kFailed;
    if (!prime_raw_at(s.zs, in, s.start))
      return SegState::kFailed;
    s.last_boundary = s.start;
    return SegState::kRunning;
  }

  // Folds n new output bytes at p into the piece's CRC. While the piece is
  // sentinel the CRC restarts after every 0xFF, so it covers exactly the bytes
  // a fix-up will not rewrite.
  static void account_output(Segment& s, bool sentinel, const unsigned char* p,
                             size_t n) {
    if (sentinel) {
      const unsigned char* ff = find_last(p, n, ds::kSentinelByte);
      if (ff != nullptr) {
        const size_t i = static_cast<size_t>(ff - p);
        s.last_ff = static_cast<int64_t>(s.out_len + i);
        s.crc = static_cast<uint32_t>(
            crc32_z(0, ff + 1, static_cast<z_size_t>(n - i - 1)));
        s.out_len += n;
        return;
      }
    }
    s.crc = static_cast<uint32_t>(crc32_z(s.crc, p, static_cast<z_size_t>(n)));
    s.out_len += n;
  }

  // The real window at the end of what `zs` has inflated: its last 32 KiB of
  // output, preceded by its dictionary while the output is shorter. Null if
  // zlib cannot say.
  static WindowPtr window_of(z_stream& zs) {
    auto window = std::make_shared<Window>(ds::kWindowBytes);
    uInt len = 0;
    if (inflateGetDictionary(&zs, window->data(), &len) != Z_OK)
      return nullptr;
    window->resize(len);
    return window;
  }

  // The real window at the end of a piece whose bytes are all final: its last
  // 32 KiB, preceded by the end of its history when it is shorter.
  static WindowPtr window_from_chunks(const std::deque<Chunk>& chunks,
                                      uint64_t out_len,
                                      const Window& history) {
    const size_t from_output =
        static_cast<size_t>(std::min<uint64_t>(ds::kWindowBytes, out_len));
    const size_t from_history = std::min<size_t>(
        history.size(), ds::kWindowBytes - from_output);
    auto window = std::make_shared<Window>(from_history + from_output);
    std::memcpy(window->data(), history.data() + history.size() - from_history,
                from_history);
    size_t need = from_output;
    for (auto it = chunks.rbegin(); it != chunks.rend() && need > 0; ++it) {
      const size_t take = std::min(need, it->size);
      std::memcpy(window->data() + from_history + need - take,
                  it->data.get() + it->size - take, take);
      need -= take;
    }
    return window;
  }

  // Inflates `s` again from its start against `history`, writing the first
  // `upto` bytes over `chunks` in place. That is safe: zlib reads back only
  // output it wrote earlier in the same call, or its own window, never bytes
  // ahead of next_out. `fix` is left open iff `keep` and the re-inflate
  // succeeded. A stop short of `upto` is a data error the stand-in hid (a
  // back-reference beyond the member's real history) or a broken invariant;
  // either way the bytes before it are correct and nothing after it is.
  Reinflate reinflate(const Segment& s, std::deque<Chunk>& chunks,
                      uint64_t upto, const Window& history, z_stream& fix,
                      bool keep, std::vector<unsigned char>& storage) {
    Reinflate r;
    std::memset(&fix, 0, sizeof(fix));
    if (inflateInit2(&fix, -15) != Z_OK) {
      r.ok = false;
      return r;
    }
    if (storage.size() < kInputBufferBytes)
      storage.resize(kInputBufferBytes);
    InputCursor in(storage);
    if ((!history.empty() &&
         inflateSetDictionary(&fix, history.data(),
                              static_cast<uInt>(history.size())) != Z_OK) ||
        !prime_raw_at(fix, in, s.start)) {
      inflateEnd(&fix);
      r.ok = false;
      return r;
    }
    size_t index = 0;
    size_t offset = 0;
    int stalled = 0;
    while (r.produced < upto) {
      if (stopping.load(std::memory_order_relaxed)) {
        r.ok = false;
        break;
      }
      if (in.pos == in.have && refill(in) <= 0) {
        r.ok = false;
        break;
      }
      Chunk& chunk = chunks[index];
      unsigned char* dst =
          reinterpret_cast<unsigned char*>(chunk.data.get()) + offset;
      fix.next_out = dst;
      fix.avail_out = static_cast<uInt>(
          std::min<uint64_t>(chunk.size - offset, upto - r.produced));
      fix.next_in = in.buf + in.pos;
      fix.avail_in = static_cast<uInt>(in.have - in.pos);
      const uInt before = fix.avail_in;
      const int status = inflate(&fix, Z_NO_FLUSH);
      in.pos += before - fix.avail_in;
      const size_t made = static_cast<size_t>(fix.next_out - dst);
      r.crc = static_cast<uint32_t>(
          crc32_z(r.crc, dst, static_cast<z_size_t>(made)));
      r.produced += made;
      offset += made;
      if (offset == chunk.size) {
        ++index;
        offset = 0;
      }
      if (status == Z_STREAM_END) {
        r.ok = r.produced == upto;
        break;
      }
      if (status != Z_OK && status != Z_BUF_ERROR) {
        r.ok = false;
        break;
      }
      stalled = (made == 0 && before == fix.avail_in) ? stalled + 1 : 0;
      if (stalled >= 2) {
        r.ok = false;
        break;
      }
    }
    r.in_offset = in.offset();
    if (!keep || !r.ok)
      inflateEnd(&fix);
    return r;
  }

  // Cuts `chunks` to their first `keep` bytes, returning the rest to the pool
  // and its charge to the budget. Only for a piece that is still charged, i.e.
  // not discarded.
  void truncate_chunks_locked(Segment& s, std::deque<Chunk>& chunks,
                              uint64_t keep) {
    uint64_t kept = 0;
    std::deque<Chunk> out;
    for (Chunk& chunk : chunks) {
      if (kept >= keep) {
        put_buffer_locked(std::move(chunk.data));
        s.staged -= chunk_bytes;
        total_staged -= chunk_bytes;
        continue;
      }
      if (kept + chunk.size > keep)
        chunk.size = static_cast<size_t>(keep - kept);
      kept += chunk.size;
      out.push_back(std::move(chunk));
    }
    chunks.swap(out);
  }

  // Hands chunks a worker held outside the Segment back to it, or to the pool
  // if the piece was discarded meanwhile (its charge went with it).
  void restore_chunks_locked(Segment& s, std::deque<Chunk>& chunks) {
    if (s.discarded) {
      for (Chunk& chunk : chunks)
        put_buffer_locked(std::move(chunk.data));
      chunks.clear();
      return;
    }
    s.chunks.swap(chunks);
  }

  // The owner of a sentinel piece the chain has reached makes it final while
  // it is still inflating, by re-inflating in place the prefix that may hold
  // stand-in bytes. If the stand-in may still be in zlib's live window (a
  // short piece, or a 0xFF in its last 32 KiB), the re-inflate runs to the end
  // of the output and its stream carries on instead; otherwise the sentinel
  // stream's window holds only real bytes and simply continues. Returns false
  // if the piece was released (its prefix failed to decode, or it was
  // discarded).
  bool make_real(const SegPtr& seg, InputCursor& in, WorkerScratch& sc,
                 const WindowPtr& history) {
    Segment& s = *seg;
    std::deque<Chunk> chunks;
    {
      std::lock_guard<std::mutex> guard(mu);
      chunks.swap(s.chunks);
      s.fixing = true;
    }
    advise_chain_input(s.start >> 3);
    const bool dirty_tail =
        s.out_len < ds::kWindowBytes ||
        s.last_ff >= static_cast<int64_t>(s.out_len - ds::kWindowBytes);
    const uint64_t upto =
        dirty_tail ? s.out_len : static_cast<uint64_t>(s.last_ff + 1);
    z_stream fix;
    Reinflate r;
    if (dirty_tail || upto > 0)
      r = reinflate(s, chunks, upto, *history, fix, dirty_tail, sc.aux);
    bool swapped = true;
    if (r.ok && dirty_tail) {
      inflateEnd(&s.zs);
      s.zs_ready = false;
      swapped = inflateCopy(&s.zs, &fix) == Z_OK;
      inflateEnd(&fix);
      s.zs_ready = swapped;
      in.seek(r.in_offset);
      s.crc = r.crc;
    } else if (r.ok && upto > 0) {
      s.crc = static_cast<uint32_t>(crc32_combine(
          r.crc, s.crc, static_cast<z_off_t>(s.out_len - upto)));
    }
    std::lock_guard<std::mutex> guard(mu);
    s.fixing = false;
    s.sentinel = false;
    s.last_ff = -1;
    if (s.discarded) {
      restore_chunks_locked(s, chunks);
      release_locked(seg, SegState::kPending);
      return false;
    }
    if (!r.ok || !swapped) {
      // The bytes before the failure are the member's true output; the stream
      // is corrupt right after them, where gzread would report it.
      if (!r.ok) {
        truncate_chunks_locked(s, chunks, r.produced);
        s.out_len = r.produced;
        s.crc = r.crc;
      }
      restore_chunks_locked(s, chunks);
      release_locked(seg, SegState::kFailed);
      return false;
    }
    restore_chunks_locked(s, chunks);
    data_cv.notify_all();
    return true;
  }

  // A fix-up job: the piece finished its stand-in pass and the chain has since
  // reached it. Re-inflates the prefix up to its last 0xFF against the real
  // history, in place, rebuilds its tail if the landing could not take one,
  // and makes it deliverable.
  void run_fix(const SegPtr& seg, WorkerScratch& sc) {
    Segment& s = *seg;
    std::deque<Chunk> chunks;
    WindowPtr history;
    bool needs_tail = false;
    {
      std::lock_guard<std::mutex> guard(mu);
      chunks.swap(s.chunks);
      history = s.history;
      needs_tail = s.state == SegState::kLanded && !s.tail;
    }
    advise_chain_input(s.start >> 3);
    const uint64_t upto = s.last_ff >= 0 ? static_cast<uint64_t>(s.last_ff + 1) : 0;
    Reinflate r;
    r.produced = upto;
    if (upto > 0) {
      z_stream fix;
      r = reinflate(s, chunks, upto, *history, fix, false, sc.aux);
    }
    WindowPtr tail;
    if (r.ok) {
      if (upto > 0)
        s.crc = static_cast<uint32_t>(crc32_combine(
            r.crc, s.crc, static_cast<z_off_t>(s.out_len - upto)));
      if (needs_tail)
        tail = window_from_chunks(chunks, s.out_len, *history);
    } else {
      s.out_len = r.produced;
      s.crc = r.crc;
    }
    s.last_ff = -1;
    std::lock_guard<std::mutex> guard(mu);
    if (!r.ok && !s.discarded) {
      truncate_chunks_locked(s, chunks, r.produced);
      s.state = SegState::kFailed;
    }
    restore_chunks_locked(s, chunks);
    if (tail)
      s.tail = tail;
    s.sentinel = false;
    s.fixing = false;
    propagate_locked();
    data_cv.notify_all();
    work_cv.notify_all();
  }

  // An open job: searches split point j and, if it holds a candidate the chain
  // has not passed, opens a sentinel piece there and runs it.
  void run_open(size_t j, WorkerScratch& sc) {
    const uint64_t bit = find_split(j, sc);
    if (bit == ds::kNotFound)
      return;
    SegPtr seg;
    {
      std::lock_guard<std::mutex> guard(mu);
      if (abort || eof || failed || segments.count(bit) != 0 ||
          bit <= chain_floor_locked())
        return;
      seg = std::make_shared<Segment>(bit, false);
      seg->state = SegState::kRunning;
      segments.emplace(bit, seg);
      piece_epoch.fetch_add(1, std::memory_order_relaxed);
    }
    run_segment(seg, sc);
  }

  // Inflates `seg` until it lands, ends, fails, suspends or is cancelled.
  // Called without mu.
  void run_segment(const SegPtr& seg, WorkerScratch& sc) {
    Segment& s = *seg;
    InputCursor in(sc.input);
    // Roles are re-read under mu at every hand-off: the chain may reach a
    // piece while it runs, and the reader may make it the head.
    bool chain_role = false;
    bool sentinel = false;
    bool fresh = false;
    WindowPtr history;
    {
      std::lock_guard<std::mutex> guard(mu);
      chain_role = s.in_chain || seg == current;
      history = s.history;
      if (!s.zs_ready) {
        fresh = true;
        s.sentinel = !s.member_head && !history;
      }
      sentinel = s.sentinel;
    }
    if (fresh) {
      const SegState opened =
          s.member_head ? open_member_head(s, in)
                        : open_block_piece(s, in, sentinel ? nullptr : history.get());
      if (opened != SegState::kRunning) {
        std::lock_guard<std::mutex> guard(mu);
        s.end = size;
        release_locked(seg, opened);
        return;
      }
    } else {
      in.seek(s.in_offset);
      if (sentinel && chain_role && history) {
        if (!make_real(seg, in, sc, history))
          return;
        sentinel = false;
      }
    }

    // Where to stop. A sentinel piece the chain has not reached stops at the
    // next split point's candidate, so the piece after it can run in
    // parallel; everything else stops only at a piece someone opened.
    uint64_t target = kNoTarget;
    bool member_ends_here = false; // the final block has begun: no more stops
    uint64_t epoch_seen = piece_epoch.load(std::memory_order_relaxed);
    auto choose_target = [&] {
      if (member_ends_here) {
        target = kNoTarget;
        return;
      }
      if (sentinel && !chain_role) {
        target = next_split_start_after(s.last_boundary, sc);
      } else {
        std::lock_guard<std::mutex> guard(mu);
        target = open_target_locked(s, in.offset(), chain_role);
      }
    };
    choose_target();

    for (;;) {
      Chunk chunk;
      {
        std::lock_guard<std::mutex> guard(mu);
        chunk.data = take_buffer_locked();
      }
      if (!chunk.data) {
        // Allocate outside mu.
        try {
          chunk.data.reset(new char[chunk_bytes]);
        } catch (const std::bad_alloc&) {
          std::lock_guard<std::mutex> guard(mu);
          s.in_offset = in.offset();
          release_locked(seg, SegState::kFailed);
          return;
        }
      }

      s.zs.next_out = reinterpret_cast<Bytef*>(chunk.data.get());
      s.zs.avail_out = static_cast<uInt>(chunk_bytes);
      bool stream_end = false;
      bool truncated = false;
      bool broken = false;
      bool landed = false;
      // inflate() always progresses with input and output space; this only
      // guards against hanging on a zlib bug. A Z_BLOCK return at a boundary
      // is progress even without new input: pigz pads with empty fixed blocks
      // that decode from bits zlib already holds.
      int stalled = 0;
      while (s.zs.avail_out > 0) {
        if (in.pos == in.have) {
          const int r = refill(in);
          if (r < 0) {
            broken = true;
            break;
          }
          if (r == 0) {
            truncated = true;
            break;
          }
          if (chain_role)
            advise_chain_input(in.base + in.have);
          if (!(sentinel && !chain_role)) {
            const uint64_t epoch = piece_epoch.load(std::memory_order_relaxed);
            if (epoch != epoch_seen) {
              epoch_seen = epoch;
              choose_target();
            }
          }
        }
        const uint64_t at = in.offset();
        const bool stepping = target != kNoTarget && at >= guard_byte(target);
        size_t avail = in.have - in.pos;
        if (target != kNoTarget && !stepping)
          avail = static_cast<size_t>(
              std::min<uint64_t>(avail, guard_byte(target) - at));
        s.zs.next_in = in.buf + in.pos;
        s.zs.avail_in = static_cast<uInt>(avail);
        unsigned char* out_before = s.zs.next_out;
        const int status = inflate(&s.zs, stepping ? Z_BLOCK : Z_NO_FLUSH);
        const size_t used = avail - s.zs.avail_in;
        in.pos += used;
        const size_t made = static_cast<size_t>(s.zs.next_out - out_before);
        if (made > 0)
          account_output(s, sentinel, out_before, made);
        if (status == Z_STREAM_END) {
          stream_end = true;
          break;
        }
        const bool at_boundary =
            stepping && status == Z_OK && (s.zs.data_type & 128) != 0;
        if (used == 0 && made == 0 && !at_boundary) {
          if (++stalled >= 2) {
            broken = true;
            break;
          }
        } else {
          stalled = 0;
        }
        if (status == Z_BUF_ERROR)
          continue; // input ran dry: loop and refill (or declare truncation)
        if (status != Z_OK) {
          broken = true; // Z_DATA_ERROR / Z_MEM_ERROR / Z_NEED_DICT / ...
          break;
        }
        if (!at_boundary)
          continue;
        if (s.zs.data_type & 64) {
          // The block that ended was the member's last: what follows is its
          // trailer, and any candidate beyond it belongs to another member.
          member_ends_here = true;
          target = kNoTarget;
          continue;
        }
        const uint64_t boundary =
            in.offset() * 8 - static_cast<uint64_t>(s.zs.data_type & 7);
        s.last_boundary = boundary;
        if (boundary < target)
          continue;
        if (boundary == target) {
          landed = true;
          break;
        }
        // Stepped over the target: it was not a block boundary of this
        // stream, or was adopted when inflate was already past it.
        if (chain_role)
          drop_stepped_over(target);
        choose_target();
      }
      const size_t produced = chunk_bytes - static_cast<size_t>(s.zs.avail_out);

      // Before the hand-off: the tail a landing leaves for the next piece, and
      // the trailer a member end leaves for the reader.
      WindowPtr tail;
      if (landed) {
        const bool tail_real =
            !sentinel ||
            (s.out_len >= ds::kWindowBytes &&
             s.last_ff < static_cast<int64_t>(s.out_len - ds::kWindowBytes));
        if (tail_real) {
          tail = window_of(s.zs);
          if (!tail)
            broken = true;
        }
      }
      if (stream_end) {
        const uint64_t trailer_at = in.offset();
        const int64_t got = pread_exact(fd, s.trailer, sizeof(s.trailer), trailer_at);
        if (got < 0) {
          broken = true;
        } else {
          s.trailer_have = static_cast<unsigned>(got);
          s.end = trailer_at + sizeof(s.trailer);
        }
      }

      std::unique_lock<std::mutex> lock(mu);
      s.in_offset = in.offset();
      if (abort || s.discarded) {
        put_buffer_locked(std::move(chunk.data));
        release_locked(seg, SegState::kPending);
        return;
      }
      if (produced > 0) {
        chunk.size = produced;
        s.chunks.push_back(std::move(chunk));
        // Charged at chunk_bytes however full, so the caps bound memory.
        s.staged += chunk_bytes;
        total_staged += chunk_bytes;
        if (seg == current && !s.sentinel)
          data_cv.notify_all();
      } else {
        put_buffer_locked(std::move(chunk.data));
      }
      if (broken) {
        release_locked(seg, SegState::kFailed);
        return;
      }
      if (landed) {
        s.next = target;
        s.tail = tail;
        release_locked(seg, SegState::kLanded);
        propagate_locked();
        return;
      }
      if (stream_end) {
        release_locked(seg, SegState::kDone);
        return;
      }
      if (truncated) {
        s.end = size;
        release_locked(seg, SegState::kTruncated);
        return;
      }
      chain_role = s.in_chain || seg == current;
      history = s.history;
      if (s.sentinel && chain_role && history) {
        lock.unlock();
        if (!make_real(seg, in, sc, history))
          return;
        sentinel = false;
        choose_target();
        lock.lock();
      } else if (!(sentinel && !chain_role) && !member_ends_here) {
        target = open_target_locked(s, in.offset(), chain_role);
      }
      if (seg == current) {
        if (helper && s.staged >= chain_cap_locked()) {
          // A loaned thread must not wait on the reader: park the head like a
          // suspended speculation and look for other work.
          release_locked(seg, SegState::kPending);
          return;
        }
        // The reader is draining this piece, so waiting is safe.
        while (!abort && seg == current && s.staged >= chain_cap_locked()) {
          space_cv.wait(lock);
        }
        if (abort) {
          release_locked(seg, SegState::kPending);
          return;
        }
        continue;
      }
      if (s.staged >= segment_cap || speculative_staged() >= spec_budget) {
        // Speculation never blocks: suspend and free the worker.
        release_locked(seg, SegState::kPending);
        return;
      }
      if (!s.member_head && head_needs_worker_locked()) {
        // A piece ahead of the head yields to it and resumes from here.
        release_locked(seg, SegState::kPending);
        return;
      }
    }
  }

  // A chain producer stepped over an open piece: nothing the chain reaches
  // starts there, so it is stopped now.
  void drop_stepped_over(uint64_t bit) {
    std::lock_guard<std::mutex> guard(mu);
    auto it = segments.find(bit);
    if (it == segments.end() || it->second->in_chain ||
        it->second->member_head)
      return;
    discard_locked(it->second);
    segments.erase(it);
  }

  void run_job(const Job& job, WorkerScratch& sc) {
    switch (job.kind) {
    case JobKind::kRun:
      run_segment(job.seg, sc);
      break;
    case JobKind::kFix:
      run_fix(job.seg, sc);
      break;
    case JobKind::kOpen:
      run_open(job.split, sc);
      break;
    case JobKind::kNone:
      break;
    }
  }

  // No more candidates. In the folded team the scanning thread now inflates.
  void finish_scan() {
    std::lock_guard<std::mutex> guard(mu);
    scan_done = true;
    work_cv.notify_all();
  }

  void worker_main(unsigned ordinal) {
    WorkerScratch sc;
    sc.input.resize(kInputBufferBytes);
    std::unique_lock<std::mutex> lock(mu);
    // Set when the helper just found no work, so the loop does not spin;
    // cleared whenever the state may have changed.
    bool helper_was_empty = false;
    for (;;) {
      if (abort)
        return;
      Job job = claim_locked();
      if (job.kind != JobKind::kNone) {
        lock.unlock();
        run_job(job, sc);
        job = Job(); // a retired piece is freed here, outside mu
        lock.lock();
        helper_was_empty = false;
        continue;
      }
      if (helper && !helper_was_empty) {
        // Reached only with nothing to claim. Called without mu.
        lock.unlock();
        const bool ran = helper(ordinal);
        lock.lock();
        // Re-claim rather than sleep: a notify may have fired meanwhile.
        helper_was_empty = !ran;
        continue;
      }
      if (helper) {
        // Poll: the caller's work does not signal work_cv, and decoding
        // stops signalling once the reader is done.
        work_cv.wait_for(lock, std::chrono::milliseconds(2));
      } else {
        work_cv.wait(lock);
      }
      helper_was_empty = false;
    }
  }

  // Stops and joins the team; a worker inside the helper finishes its unit
  // first.
  void shutdown() {
    stopping.store(true, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> guard(mu);
      abort = true;
    }
    work_cv.notify_all();
    data_cv.notify_all();
    space_cv.notify_all();
    for (std::thread& thread : threads) {
      if (thread.joinable())
        thread.join();
    }
    threads.clear();
  }
};

bool gz_member_parallel_eligible(const std::string& path, uint64_t min_bytes) {
  if (path.empty() || path == "-")
    return false;
  struct stat info;
  if (::stat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode))
    return false;
  if (info.st_size < 0)
    return false;
  if (static_cast<uint64_t>(info.st_size) < min_bytes)
    return false;
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    return false;
  unsigned char head[18] = {0};
  const int64_t got = pread_exact(fd, head, sizeof(head), 0);
  ::close(fd);
  if (got < 3)
    return false;
  // Beyond gzread's magic test, require CM = deflate, so a malformed header
  // takes the serial path, which reports the error.
  if (head[0] != 0x1f || head[1] != 0x8b || head[2] != 0x08)
    return false;
  return !looks_like_bgzf(head, static_cast<size_t>(got));
}

MemberParallelGzByteSource::MemberParallelGzByteSource(
    const std::string& path, GzMemberSourceOptions options)
    : impl_(new Impl()) {
  Impl* self = impl_.get();
  self->path = path;
  self->fd = ::open(path.c_str(), O_RDONLY);
  if (self->fd < 0) {
    throw std::runtime_error("failed to open FASTA/FASTQ: " + path);
  }
  struct stat info;
  if (::fstat(self->fd, &info) != 0 || !S_ISREG(info.st_mode)) {
    throw std::runtime_error("not a regular file: " + path);
  }
  self->size = static_cast<uint64_t>(info.st_size);
  // advise_chain_input() adds explicit readahead and drop-behind.
#ifdef POSIX_FADV_SEQUENTIAL
  (void)posix_fadvise(self->fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

  unsigned wanted = options.workers;
  if (wanted == 0) {
    // Cgroup-aware, so a container with a CPU quota gets a small team.
    const int available = ::fa::cpu::threading::default_thread_count();
    wanted =
        std::min(8u, available < 1 ? 1u : static_cast<unsigned>(available));
  }
  self->workers = std::max(1u, wanted);
  self->spec_worker_cap = self->workers > 1 ? self->workers - 1 : 0;
  self->helper = options.helper;
  self->chunk_bytes =
      std::min(kMaxChunkBytes, std::max<size_t>(1024, options.chunk_bytes));
  self->segment_cap =
      std::max<uint64_t>(self->chunk_bytes, options.segment_stage_cap);
  // The budget must leave the current member a full segment_cap beyond
  // speculation's share, so it is raised to at least twice the cap
  // (saturating).
  const uint64_t budget = std::max<uint64_t>(
      options.staging_budget, self->segment_cap > (UINT64_MAX / 2)
                                  ? UINT64_MAX
                                  : 2 * self->segment_cap);
  self->spec_budget = budget - self->segment_cap;
  self->staging_budget = budget;
  // Cutting a member into pieces needs a second worker and a file longer than
  // one piece.
  if (self->workers > 1 && self->size > kPieceBytes) {
    self->split_count = static_cast<size_t>((self->size - 1) / kPieceBytes);
    self->splits.reset(new Split[self->split_count + 1]);
  }
  self->piece_lookahead = static_cast<size_t>(kPiecesPerWorker) * self->workers;
  self->spec_watermark = budget / kSpeculationWatermarkDivisor;

  // Team shape: see the top of this file.
  self->folded_scan = static_cast<bool>(options.helper) && self->workers > 1;
  const bool dedicated_scanner = !options.helper;
  // Nothing will publish candidates, so the reserve never needs widening.
  self->scan_done = !dedicated_scanner && !self->folded_scan;
  self->threads.reserve(self->workers + (dedicated_scanner ? 1u : 0u));
  if (dedicated_scanner) {
    self->threads.emplace_back([self] {
      fa::cpu::threading::name_current_thread("fa-gz-scan");
      self->scanner_main(0);
      self->finish_scan();
    });
  }
  for (unsigned i = 0; i < self->workers; ++i) {
    const bool scan_first = i == 0 && self->folded_scan;
    self->threads.emplace_back([self, scan_first, i] {
      fa::cpu::threading::name_current_thread("fa-gz", static_cast<int>(i));
      if (scan_first) {
        self->scanner_main(i);
        self->finish_scan();
      }
      self->worker_main(i);
    });
  }
}

MemberParallelGzByteSource::~MemberParallelGzByteSource() = default;

int MemberParallelGzByteSource::read(void* buf, int len) {
  if (buf == nullptr || len <= 0)
    return 0;
  try {
    return impl_->read(static_cast<char*>(buf), static_cast<size_t>(len));
  } catch (...) {
    // Nothing may propagate into kseq.
    return -1;
  }
}

bool MemberParallelGzByteSource::cut_in_data() const {
  return impl_->cut_in_data;
}

} // namespace io
} // namespace cpu
} // namespace fa
