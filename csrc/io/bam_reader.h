// BAM/uBAM input. Records become FastxRecords with the raw tag block kept, so
// base-modification tags (MM/ML) and read groups carry through to the output.
// An aligned BAM reads like a uBAM: secondary and supplementary records are
// skipped, SEQ is returned in the original read orientation, and a mapped
// record loses the alignment tags of the aligner that wrote it.
#pragma once

#include "bgzf.h"
#include "fastx.h"
#include "sam.h"

#include <sys/stat.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fa { namespace cpu { namespace io {

// True for a gzip/BGZF file whose first decompressed bytes are "BAM\1".
inline bool path_is_bam(const std::string& path) {
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

// CRAM is detected only to refuse it with a clear message; otherwise the
// FASTA/FASTQ reader would silently find no records.
inline bool path_is_cram(const std::string& path) {
    if (path == "-") return false;
    struct stat st;  // non-rewindable inputs: skip the probe (see path_is_bam)
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char m[4] = {0, 0, 0, 0};
    f.read(m, 4);
    return f.gcount() == 4 && std::string(m, 4) == std::string("CRAM", 4);
}

class BamReader {
public:
    // The blocks are decoded by make_byte_source(path, plan); `plan` is
    // borrowed, and null takes the defaults.
    explicit BamReader(const std::string& path, const IoPlan* plan = nullptr)
        : path_(path),
          bgzf_(std::make_unique<BgzfReader>(
              path, make_byte_source(path, plan, InputPolicy::bam_integrity))) {
        parse_header();
    }

    // Whether a thread team decodes the blocks.
    bool decode_parallel_active() const { return bgzf_->decode_parallel_active(); }

    const std::string& sam_header_text() const { return header_text_; }

    bool next(FastxRecord& rec) {
        uint16_t flag = 0;
        while (read_record(rec, flag)) {
            if (!(flag & 0x4) || (flag & 0x900)) warn_aligned_input();
            if (flag & 0x900) continue;
            if (!(flag & 0x4)) drop_alignment_tags(rec);
            return true;
        }
        return false;
    }

private:
    bool read_record(FastxRecord& rec, uint16_t& flag_out) {
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

    void warn_aligned_input() {
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

    void drop_alignment_tags(FastxRecord& rec);

    void parse_header() {
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

    int32_t read_i32() { int32_t v = 0; bgzf_->read_exact(&v, 4, "BAM header"); return v; }
    void skip(size_t n, const char* what) {
        std::string tmp(n, '\0');
        if (n) bgzf_->read_exact(&tmp[0], n, what);
    }

    std::string path_;
    std::unique_ptr<BgzfReader> bgzf_;
    std::string header_text_;
    bool warned_aligned_ = false;
};

[[noreturn]] inline void throw_malformed_tag_block(const std::string& read_name,
                                                  const char* what) {
    throw std::runtime_error("malformed BAM: tag block of read " + read_name +
                             ": " + what);
}

// Whether a tag starts at blob[at]. Fewer than 3 bytes after the last tag
// hold no key and type and are ignored, as htslib's bam_aux_next ignores them;
// a tag with its key and type but a short value is not.
inline bool bam_tag_at(const std::string& blob, size_t at) {
    return blob.size() - at >= 3;
}

// The byte size of the tag at blob[at], key and type included. Throws when
// the tag runs past the block or its type is unknown.
inline size_t bam_tag_size(const std::string& blob, size_t at,
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
inline void BamReader::drop_alignment_tags(FastxRecord& rec) {
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

inline void append_printf_g(std::string& out, double d) {
    char buf[32];
    const int n = std::snprintf(buf, sizeof buf, "%g", d);
    out.append(buf, static_cast<size_t>(n));
}

// htslib's kputd (kstring.c), which spells B:f elements: six significant
// digits, trailing zeros cut, printf's %g outside [1e-4, 999999].
inline void append_kputd(std::string& out, double d) {
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

template <typename T>
inline T load_le(const char* p) {
    T v;
    std::memcpy(&v, p, sizeof v);
    return v;
}

template <typename T>
inline void append_decimal(std::string& out, T v) {
    char buf[24];
    const auto r = std::to_chars(buf, buf + sizeof buf, v);
    out.append(buf, r.ptr);
}

}  // namespace detail

// One tag of `size` bytes (as bam_tag_size gives) as SAM text, spelled as
// htslib's sam_format_aux1 spells it: every integer type as i, f with %g,
// B:f elements with kputd.
inline void append_bam_tag_sam_text(std::string& out, const char* tag, size_t size) {
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
