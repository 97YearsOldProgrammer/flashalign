// BGZF writer and reader over zlib raw deflate. BGZF, the container BAM uses,
// is a series of gzip members of at most 64 KiB uncompressed, each with a 'BC'
// extra subfield giving its size, ended by a 28-byte empty block.
#pragma once

#include <zlib.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace fa { namespace cpu { namespace io {

class BgzfWriter {
public:
    explicit BgzfWriter(std::ostream& out, int level = 6)
        : out_(out), level_(level) { ubuf_.reserve(MAX_BLOCK_UNCOMP); }

    ~BgzfWriter() { try { close(); } catch (...) {} }

    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;

    void write(const void* data, size_t len) {
        const char* p = static_cast<const char*>(data);
        while (len > 0) {
            const size_t take = std::min(len, MAX_BLOCK_UNCOMP - ubuf_.size());
            ubuf_.append(p, take);
            p += take;
            len -= take;
            if (ubuf_.size() == MAX_BLOCK_UNCOMP) flush_block();
        }
    }

    void write(const std::string& s) { write(s.data(), s.size()); }

    // Virtual offset of the next byte: (block's compressed offset << 16) |
    // offset within the pending block.
    uint64_t virtual_offset() const {
        return (coffset_ << 16) | static_cast<uint64_t>(ubuf_.size());
    }

    // Flushes the pending block and writes the EOF block. Idempotent.
    void close() {
        if (closed_) return;
        if (!ubuf_.empty()) flush_block();
        write_eof_block();
        out_.flush();
        require_output("flush");
        closed_ = true;
    }

private:
    void require_output(const char* operation) const {
        if (!out_) {
            throw std::runtime_error(
                std::string("bgzf: output ") + operation + " failed");
        }
    }

    static void put16(std::string& b, uint16_t v) {
        b.push_back(static_cast<char>(v & 0xff));
        b.push_back(static_cast<char>((v >> 8) & 0xff));
    }
    static void put32(std::string& b, uint32_t v) {
        for (int i = 0; i < 4; ++i) b.push_back(static_cast<char>((v >> (8 * i)) & 0xff));
    }

    void flush_block() {
        const size_t n = ubuf_.size();

        z_stream zs;
        std::memset(&zs, 0, sizeof(zs));
        if (deflateInit2(&zs, level_, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            throw std::runtime_error("bgzf: deflateInit2 failed");
        }
        std::string comp;
        comp.resize(deflateBound(&zs, static_cast<uLong>(n)));
        zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(ubuf_.data()));
        zs.avail_in = static_cast<uInt>(n);
        zs.next_out = reinterpret_cast<Bytef*>(&comp[0]);
        zs.avail_out = static_cast<uInt>(comp.size());
        const int r = deflate(&zs, Z_FINISH);
        if (r != Z_STREAM_END) {
            deflateEnd(&zs);
            throw std::runtime_error("bgzf: deflate did not finish in one block");
        }
        const size_t clen = comp.size() - zs.avail_out;
        deflateEnd(&zs);

        const uint32_t crc = crc32(
            crc32(0L, Z_NULL, 0),
            reinterpret_cast<const Bytef*>(ubuf_.data()), static_cast<uInt>(n));

        const size_t block_size = 18 + clen + 8;  // header(18) + payload + crc(4) + isize(4)
        if (block_size - 1 > 0xffff) {
            throw std::runtime_error("bgzf: compressed block exceeds 64 KiB");
        }

        // 18-byte BGZF header.
        std::string hdr;
        hdr.reserve(18);
        const unsigned char gz[10] = {0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff};
        hdr.append(reinterpret_cast<const char*>(gz), 10);  // ID1 ID2 CM FLG MTIME XFL OS
        put16(hdr, 6);                                      // XLEN
        hdr.push_back(66);                                  // SI1 = 'B'
        hdr.push_back(67);                                  // SI2 = 'C'
        put16(hdr, 2);                                      // SLEN
        put16(hdr, static_cast<uint16_t>(block_size - 1));  // BSIZE = total - 1

        std::string ftr;
        ftr.reserve(8);
        put32(ftr, crc);
        put32(ftr, static_cast<uint32_t>(n));               // ISIZE

        out_.write(hdr.data(), static_cast<std::streamsize>(hdr.size()));
        out_.write(comp.data(), static_cast<std::streamsize>(clen));
        out_.write(ftr.data(), static_cast<std::streamsize>(ftr.size()));
        require_output("write");

        coffset_ += block_size;
        ubuf_.clear();
    }

    void write_eof_block() {
        // The canonical 28-byte empty BGZF block htslib appends so readers can
        // detect a complete (non-truncated) file.
        static const unsigned char eof[28] = {
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
            0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00};
        out_.write(reinterpret_cast<const char*>(eof), 28);
        require_output("write");
        coffset_ += 28;
    }

    std::ostream& out_;
    int level_;
    std::string ubuf_;        // pending uncompressed bytes (<= MAX_BLOCK_UNCOMP)
    uint64_t coffset_ = 0;    // compressed bytes written so far (block-aligned)
    bool closed_ = false;
    static constexpr size_t MAX_BLOCK_UNCOMP = 0xff00;  // 65280, htslib's safe cap
};

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
