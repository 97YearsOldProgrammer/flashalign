// A ByteSource that inflates a gzip file on several threads and delivers the
// bytes in order. The members of a multi-member file (e.g. `cat a.fq.gz
// b.fq.gz`) are inflated in parallel, and so are pieces of one member, cut at
// deflate block boundaries (deflate_split.h). A piece is inflated against a
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
#pragma once

#include "byte_source.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>

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
  // False forces the serial ZlibByteSource.
  bool parallel_gz = true;
  GzMemberSourceOptions gz;
};

// True for a regular file of at least `min_bytes` that starts with a deflate
// gzip member (1f 8b 08) and is not BGZF, whose small members would make
// parallel decoding slower. Never throws; a failed probe answers false.
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
