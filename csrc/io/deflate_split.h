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
#pragma once

#include <zlib.h>

#include <cstddef>
#include <cstdint>
#include <vector>

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
