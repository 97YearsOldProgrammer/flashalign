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
#pragma once

#include "byte_source.h"
#include "gz_member_source.h"

#include <cstdint>
#include <memory>
#include <string>

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
