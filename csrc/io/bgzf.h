// BGZF reader over zlib raw inflate, for uBAM input. BGZF, the container BAM
// uses, is a series of gzip members of at most 64 KiB uncompressed, each with a
// 'BC' extra subfield giving its size, ended by a 28-byte empty block.
#pragma once

#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <istream>
#include <stdexcept>
#include <string>

namespace fa { namespace cpu { namespace io {

// Streaming BGZF reader: inflates blocks on demand and serves a flat byte
// stream. Skips empty blocks (incl. the EOF marker) transparently.
class BgzfReader {
public:
    explicit BgzfReader(std::istream& in) : in_(in) {}

    // Read up to `len` bytes; returns the count actually read (< len only at EOF).
    size_t read(void* dst, size_t len) {
        char* out = static_cast<char*>(dst);
        size_t got = 0;
        while (got < len) {
            if (upos_ >= ubuf_.size()) {
                if (!fill_block()) break;
            }
            const size_t n = std::min(len - got, ubuf_.size() - upos_);
            std::memcpy(out + got, ubuf_.data() + upos_, n);
            upos_ += n;
            got += n;
        }
        return got;
    }

    // Read exactly `len` bytes or throw (for fixed-width fields).
    void read_exact(void* dst, size_t len) {
        if (read(dst, len) != len) throw std::runtime_error("bgzf: unexpected EOF");
    }

private:
    bool fill_block() {
        for (;;) {
            unsigned char hdr[18];
            in_.read(reinterpret_cast<char*>(hdr), 18);
            const std::streamsize n = in_.gcount();
            if (n == 0) return false;  // clean EOF at a block boundary
            if (n != 18 || hdr[0] != 0x1f || hdr[1] != 0x8b) {
                throw std::runtime_error("bgzf: bad block header");
            }
            const uint16_t bsize = static_cast<uint16_t>(hdr[16] | (hdr[17] << 8));
            const size_t total = static_cast<size_t>(bsize) + 1;
            if (total < 26) throw std::runtime_error("bgzf: block too small");
            const size_t clen = total - 18 - 8;

            std::string cdata;
            cdata.resize(clen);
            if (clen) {
                in_.read(&cdata[0], static_cast<std::streamsize>(clen));
                if (static_cast<size_t>(in_.gcount()) != clen) {
                    throw std::runtime_error("bgzf: truncated block payload");
                }
            }
            unsigned char tail[8];
            in_.read(reinterpret_cast<char*>(tail), 8);
            if (in_.gcount() != 8) throw std::runtime_error("bgzf: truncated block footer");
            const uint32_t isize = tail[4] | (tail[5] << 8) | (tail[6] << 16) |
                                   (static_cast<uint32_t>(tail[7]) << 24);
            if (isize == 0) continue;  // empty block (e.g. EOF marker) -> next

            // ISIZE above 64 KiB means this is not BGZF (e.g. path_is_bam
            // probing a plain .fq.gz). Refuse before allocating.
            if (isize > 0x10000) {
                throw std::runtime_error(
                    "bgzf: uncompressed block size exceeds 64 KiB (not a BGZF stream)");
            }

            ubuf_.assign(isize, '\0');
            upos_ = 0;
            z_stream zs;
            std::memset(&zs, 0, sizeof(zs));
            if (inflateInit2(&zs, -15) != Z_OK) throw std::runtime_error("bgzf: inflateInit2");
            zs.next_in = reinterpret_cast<Bytef*>(cdata.data());
            zs.avail_in = static_cast<uInt>(clen);
            zs.next_out = reinterpret_cast<Bytef*>(&ubuf_[0]);
            zs.avail_out = static_cast<uInt>(isize);
            const int r = inflate(&zs, Z_FINISH);
            inflateEnd(&zs);
            if (r != Z_STREAM_END) throw std::runtime_error("bgzf: inflate failed");
            return true;
        }
    }

    std::istream& in_;
    std::string ubuf_;
    size_t upos_ = 0;
};

}}}  // namespace fa::cpu::io
