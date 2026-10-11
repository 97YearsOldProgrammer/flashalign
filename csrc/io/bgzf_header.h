// Bounded gzip extra-field scan shared by BGZF probes and decoders.
#pragma once

#include <cstddef>
#include <cstdint>

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

inline unsigned bgzf_le16(const unsigned char* p) {
    return static_cast<unsigned>(p[0]) | static_cast<unsigned>(p[1]) << 8;
}

inline BgzfHeader parse_bgzf_header(const unsigned char* p, size_t n) {
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
