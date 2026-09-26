// The header check here is a fast filter in front of zlib's one-block
// inflate, which is what makes a position a candidate.
#include "deflate_split.h"

#include <algorithm>
#include <cstring>

namespace fa {
namespace cpu {
namespace io {
namespace deflate_split {
namespace {

// A dynamic header is at most 17 + 19 * 3 bits of counts and code-length code
// plus 286 + 30 code lengths of at most 7 + 7 bits: under 600 bytes. A
// position with less than this left in the buffer is not examined.
constexpr size_t kHeaderLookaheadBytes = 640;
// The confirming inflate discards its output in slices of this size; a block
// that produces more than the cap is treated as not confirmed.
constexpr size_t kConfirmSliceBytes = 64u << 10;
constexpr size_t kConfirmOutputCap = 64u << 20;

// Code-length code order (RFC 1951 3.2.7).
constexpr int kClOrder[19] = {16, 17, 18, 0, 8,  7, 9,  6, 10, 5,
                              11, 4,  12, 3, 13, 2, 14, 1, 15};

// n <= 57 bits from bit position `pos`, LSB first.
inline uint64_t getbits(const unsigned char* buf, size_t nbytes, uint64_t pos,
                        int n) {
  const size_t b = static_cast<size_t>(pos >> 3);
  const int shift = static_cast<int>(pos & 7);
  uint64_t v = 0;
  if (b + 8 <= nbytes) {
    std::memcpy(&v, buf + b, 8); // little-endian host
  } else {
    for (size_t i = 0; b + i < nbytes && i < 8; ++i)
      v |= static_cast<uint64_t>(buf[b + i]) << (8 * i);
  }
  return (v >> shift) & ((1ull << n) - 1);
}

// zlib inflate_table()'s rule for one code: over-subscribed is refused;
// incomplete is refused unless it is a single length-1 code (never for the
// code-length code); no symbols at all is accepted for the literal/length and
// distance codes (the caller requires end-of-block) and refused for the
// code-length code.
bool code_acceptable(const unsigned char* lens, int n, bool code_length_code) {
  int count[16] = {0};
  int max = 0;
  for (int i = 0; i < n; ++i) {
    ++count[lens[i]];
    max = std::max<int>(max, lens[i]);
  }
  if (max == 0)
    return !code_length_code;
  int left = 1;
  for (int len = 1; len <= 15; ++len) {
    left <<= 1;
    left -= count[len];
    if (left < 0)
      return false;
  }
  return !(left > 0 && (code_length_code || max != 1));
}

struct CodeLengthHuffman {
  short count[8];
  short symbol[19];
};

int decode_code_length(const CodeLengthHuffman& h, const unsigned char* buf,
                       size_t nbytes, uint64_t& pos) {
  int code = 0, first = 0, index = 0;
  for (int len = 1; len <= 7; ++len) {
    code |= static_cast<int>(getbits(buf, nbytes, pos++, 1));
    const int count = h.count[len];
    if (code - count < first)
      return h.symbol[index + (code - first)];
    index += count;
    first += count;
    first <<= 1;
    code <<= 1;
  }
  return -1;
}

// A dynamic block header, BFINAL = 0, that zlib would accept.
bool dynamic_header_acceptable(const unsigned char* buf, size_t nbytes,
                               uint64_t pos) {
  if ((pos >> 3) + kHeaderLookaheadBytes > nbytes)
    return false;
  const uint64_t head = getbits(buf, nbytes, pos, 17);
  if ((head & 7) != 4) // BFINAL = 0, BTYPE = 10
    return false;
  const int hlit = static_cast<int>((head >> 3) & 31);
  const int hdist = static_cast<int>((head >> 8) & 31);
  const int hclen = static_cast<int>((head >> 13) & 15);
  if (hlit > 29 || hdist > 29)
    return false;
  pos += 17;
  unsigned char cl[19] = {0};
  const int ncl = hclen + 4;
  const uint64_t packed = getbits(buf, nbytes, pos, 57);
  for (int i = 0; i < ncl; ++i)
    cl[kClOrder[i]] = static_cast<unsigned char>((packed >> (3 * i)) & 7);
  pos += 3 * static_cast<uint64_t>(ncl);
  if (!code_acceptable(cl, 19, true))
    return false;
  CodeLengthHuffman h;
  std::memset(&h, 0, sizeof(h));
  for (int i = 0; i < 19; ++i)
    ++h.count[cl[i]];
  h.count[0] = 0;
  short offs[9];
  offs[1] = 0;
  for (int len = 1; len < 8; ++len)
    offs[len + 1] = static_cast<short>(offs[len] + h.count[len]);
  for (int i = 0; i < 19; ++i)
    if (cl[i])
      h.symbol[offs[cl[i]]++] = static_cast<short>(i);
  const int nlen = hlit + 257;
  const int n = nlen + hdist + 1;
  unsigned char lens[320];
  int i = 0;
  while (i < n) {
    const int sym = decode_code_length(h, buf, nbytes, pos);
    if (sym < 0)
      return false;
    if (sym < 16) {
      lens[i++] = static_cast<unsigned char>(sym);
      continue;
    }
    int repeat = 0;
    unsigned char value = 0;
    if (sym == 16) {
      if (i == 0) // nothing to repeat
        return false;
      value = lens[i - 1];
      repeat = 3 + static_cast<int>(getbits(buf, nbytes, pos, 2));
      pos += 2;
    } else if (sym == 17) {
      repeat = 3 + static_cast<int>(getbits(buf, nbytes, pos, 3));
      pos += 3;
    } else {
      repeat = 11 + static_cast<int>(getbits(buf, nbytes, pos, 7));
      pos += 7;
    }
    if (i + repeat > n)
      return false;
    while (repeat-- > 0)
      lens[i++] = value;
  }
  if (lens[256] == 0) // no end-of-block code
    return false;
  return code_acceptable(lens, nlen, false) &&
         code_acceptable(lens + nlen, n - nlen, false);
}

// Raw inflate from `bit` against the stand-in window must finish one whole
// block (Z_BLOCK returns with data_type bit 128) without an error. The
// stand-in makes every back-reference legal, so this tests the block's
// structure.
bool zlib_accepts_block(const unsigned char* buf, size_t nbytes, uint64_t bit,
                        std::vector<unsigned char>& out) {
  const size_t first = static_cast<size_t>(bit >> 3);
  if (first >= nbytes)
    return false;
  z_stream zs;
  std::memset(&zs, 0, sizeof(zs));
  if (inflateInit2(&zs, -15) != Z_OK)
    return false;
  bool ok = false;
  if (inflateSetDictionary(&zs, sentinel_window(),
                           static_cast<uInt>(kWindowBytes)) == Z_OK &&
      prime_at(zs, buf[first], static_cast<unsigned>(bit & 7))) {
    const size_t start = first + ((bit & 7) != 0 ? 1 : 0);
    zs.next_in = const_cast<Bytef*>(buf + start);
    zs.avail_in = static_cast<uInt>(nbytes - start);
    if (out.size() < kConfirmSliceBytes)
      out.resize(kConfirmSliceBytes);
    size_t produced = 0;
    for (;;) {
      zs.next_out = out.data();
      zs.avail_out = static_cast<uInt>(kConfirmSliceBytes);
      const int status = inflate(&zs, Z_BLOCK);
      produced += kConfirmSliceBytes - zs.avail_out;
      if (status != Z_OK)
        break;
      if (zs.data_type & 128) {
        ok = true;
        break;
      }
      if (zs.avail_in == 0 || produced > kConfirmOutputCap)
        break; // the block runs past what we were given
    }
  }
  inflateEnd(&zs);
  return ok;
}

} // namespace

const unsigned char* sentinel_window() {
  static const std::vector<unsigned char> window(kWindowBytes, kSentinelByte);
  return window.data();
}

bool prime_at(z_stream& zs, unsigned char first_byte, unsigned bit_in_byte) {
  if (bit_in_byte == 0)
    return true;
  return inflatePrime(&zs, static_cast<int>(8 - bit_in_byte),
                      static_cast<int>(first_byte >> bit_in_byte)) == Z_OK;
}

uint64_t find_block_start(const unsigned char* buf, size_t nbytes,
                          size_t scan_from, size_t scan_to,
                          size_t skip_aligned_below, bool aligned_only,
                          std::vector<unsigned char>& out) {
  scan_to = std::min(scan_to, nbytes);
  if (aligned_only) {
    for (size_t b = scan_from; b < scan_to; ++b) {
      if ((buf[b] & 7) != 4) // BFINAL = 0, BTYPE = 10 in the low three bits
        continue;
      const uint64_t bit = static_cast<uint64_t>(b) * 8;
      if (!dynamic_header_acceptable(buf, nbytes, bit))
        continue;
      if (zlib_accepts_block(buf, nbytes, bit, out))
        return bit;
    }
    return kNotFound;
  }
  const uint64_t end = static_cast<uint64_t>(scan_to) * 8;
  for (uint64_t bit = static_cast<uint64_t>(scan_from) * 8; bit < end; ++bit) {
    if ((bit & 7) == 0 && (bit >> 3) < skip_aligned_below)
      continue;
    if (getbits(buf, nbytes, bit, 3) != 4)
      continue;
    if (!dynamic_header_acceptable(buf, nbytes, bit))
      continue;
    if (zlib_accepts_block(buf, nbytes, bit, out))
      return bit;
  }
  return kNotFound;
}

} // namespace deflate_split
} // namespace io
} // namespace cpu
} // namespace fa
