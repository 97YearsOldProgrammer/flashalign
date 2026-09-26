// BAM/uBAM input. Records become FastxRecords with the raw tag block kept, so
// base-modification tags (MM/ML) and read groups carry through to the output.
// SEQ is returned in the original read orientation, so an aligned BAM reads
// like a uBAM.
#pragma once

#include "bgzf.h"
#include "fastx.h"
#include "sam.h"

#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fa { namespace cpu { namespace io {

// True for a BGZF file whose first decompressed bytes are "BAM\1".
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
    f.clear();
    f.seekg(0);
    try {
        BgzfReader r(f);
        char magic[4] = {0, 0, 0, 0};
        return r.read(magic, 4) == 4 && std::string(magic, 4) == std::string("BAM\1", 4);
    } catch (...) {
        return false;
    }
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
    explicit BamReader(const std::string& path)
        : in_(path, std::ios::binary), bgzf_(in_)
    {
        if (!in_) throw std::runtime_error("failed to open BAM: " + path);
        parse_header(path);
    }

    const std::vector<output::SamReference>& references() const { return refs_; }
    const std::string& sam_header_text() const { return header_text_; }

    bool next(FastxRecord& rec) {
        int32_t block_size = 0;
        if (bgzf_.read(&block_size, 4) != 4) return false;  // EOF
        // Every field is bounds-checked against the block (`need`); the fixed
        // part of a record is 32 bytes.
        if (block_size < 32) {
            throw std::runtime_error("malformed BAM: record block_size too small");
        }
        std::string buf(static_cast<size_t>(block_size), '\0');
        bgzf_.read_exact(&buf[0], buf.size());

        const size_t n = buf.size();
        size_t c = 0;
        auto need = [&](size_t k) {
            if (k > n || c > n - k) {  // c + k > n, overflow-safe
                throw std::runtime_error("malformed BAM: record field exceeds block");
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
        const int32_t l_seq = i32();
        skip(4 + 4 + 4);              // next_refID, next_pos, tlen
        if (l_seq < 0) throw std::runtime_error("malformed BAM: negative l_seq");

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

private:
    void parse_header(const std::string& path) {
        char magic[4];
        if (bgzf_.read(magic, 4) != 4 || std::string(magic, 4) != std::string("BAM\1", 4)) {
            throw std::runtime_error("not a BAM file: " + path);
        }
        int32_t l_text = read_i32();
        if (l_text < 0) throw std::runtime_error("malformed BAM: negative l_text");
        // Kept so @RG/@CO lines can be carried into the output header.
        header_text_.resize(static_cast<size_t>(l_text));
        if (l_text) bgzf_.read_exact(&header_text_[0], header_text_.size());
        const int32_t n_ref = read_i32();
        if (n_ref < 0) throw std::runtime_error("malformed BAM: negative n_ref");
        refs_.reserve(static_cast<size_t>(std::max(0, n_ref)));
        for (int32_t i = 0; i < n_ref; ++i) {
            const int32_t l_name = read_i32();
            if (l_name < 0) throw std::runtime_error("malformed BAM: negative l_name");
            std::string name(static_cast<size_t>(l_name ? l_name - 1 : 0), '\0');
            if (l_name) {
                std::string tmp(static_cast<size_t>(l_name), '\0');
                bgzf_.read_exact(&tmp[0], tmp.size());
                name.assign(tmp.data(), tmp.size() ? tmp.size() - 1 : 0);
            }
            const int32_t l_ref = read_i32();
            refs_.push_back({name, l_ref});
        }
    }

    int32_t read_i32() { int32_t v = 0; bgzf_.read_exact(&v, 4); return v; }
    void skip(size_t n) {
        std::string tmp(n, '\0');
        if (n) bgzf_.read_exact(&tmp[0], n);
    }

    std::ifstream in_;
    BgzfReader bgzf_;
    std::vector<output::SamReference> refs_;
    std::string header_text_;
};

// A raw BAM tag block as tab-separated SAM tag text. Integer types become
// 'i'; 'B' keeps its subtype.
inline std::string bam_tags_to_sam_text(const std::string& blob) {
    std::string out;
    size_t c = 0;
    auto emit_tab = [&] { if (!out.empty()) out.push_back('\t'); };
    while (c + 3 <= blob.size()) {
        const char t0 = blob[c], t1 = blob[c + 1], type = blob[c + 2];
        c += 3;
        emit_tab();
        out.push_back(t0); out.push_back(t1); out.push_back(':');
        auto rd_i = [&](int width, bool sign) -> long long {
            long long v = 0;
            for (int k = 0; k < width; ++k) v |= static_cast<long long>(static_cast<uint8_t>(blob[c + k])) << (8 * k);
            if (sign && width < 8) {  // sign-extend
                const long long m = 1LL << (8 * width - 1);
                v = (v ^ m) - m;
            }
            c += width;
            return v;
        };
        switch (type) {
            case 'A': out += "A:"; out.push_back(blob[c++]); break;
            case 'c': out += "i:"; out += std::to_string(rd_i(1, true)); break;
            case 'C': out += "i:"; out += std::to_string(rd_i(1, false)); break;
            case 's': out += "i:"; out += std::to_string(rd_i(2, true)); break;
            case 'S': out += "i:"; out += std::to_string(rd_i(2, false)); break;
            case 'i': out += "i:"; out += std::to_string(rd_i(4, true)); break;
            case 'I': out += "i:"; out += std::to_string(rd_i(4, false)); break;
            case 'f': {
                float fv; std::memcpy(&fv, blob.data() + c, 4); c += 4;
                out += "f:"; out += std::to_string(fv); break;
            }
            case 'Z': case 'H': {
                out.push_back(type); out.push_back(':');
                while (c < blob.size() && blob[c] != '\0') out.push_back(blob[c++]);
                ++c;  // skip NUL
                break;
            }
            case 'B': {
                const char sub = blob[c++];
                int32_t count; std::memcpy(&count, blob.data() + c, 4); c += 4;
                out += "B:"; out.push_back(sub);
                const int w = (sub == 'c' || sub == 'C') ? 1 : (sub == 's' || sub == 'S') ? 2 : 4;
                const bool sgn = (sub == 'c' || sub == 's' || sub == 'i');
                for (int32_t k = 0; k < count; ++k) {
                    out.push_back(',');
                    if (sub == 'f') { float fv; std::memcpy(&fv, blob.data() + c, 4); c += 4; out += std::to_string(fv); }
                    else { out += std::to_string(rd_i(w, sgn)); }
                }
                break;
            }
            default: return out;  // unknown type: stop
        }
    }
    return out;
}

}}}  // namespace fa::cpu::io
