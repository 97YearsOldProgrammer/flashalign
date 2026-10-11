// Flat decoded stream adapter for BAM, with source error context.
#pragma once

#include "byte_source.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace fa { namespace cpu { namespace io {

class BgzfReader {
public:
    BgzfReader(std::string path, std::unique_ptr<ByteSource> src)
        : path_(std::move(path)), src_(std::move(src)) {}

    // Read up to `len` bytes; returns less only at source EOF.
    size_t read(void* dst, size_t len) {
        char* out = static_cast<char*>(dst);
        size_t got = 0;
        while (got < len) {
            const size_t ask = std::min<size_t>(len - got, INT_MAX);
            const int n = src_->read(out + got, static_cast<int>(ask));
            if (n < 0) {
                const auto status = src_->status();
                throw error(status.error.empty() ? "corrupt or unreadable BAM input" :
                            status.error);
            }
            if (n == 0) break;
            got += static_cast<size_t>(n);
        }
        return got;
    }

    void read_exact(void* dst, size_t len, const char* field) {
        if (read(dst, len) != len) throw error(std::string("truncated ") + field);
    }

    bool decode_parallel_active() const { return src_->is_parallel(); }

    std::runtime_error error(const std::string& reason) const {
        const auto status = src_->status();
        std::string text = path_ + ": " + reason;
        if (status.compressed_offset != ~uint64_t{0}) {
            text += status.at_member_start ? " (member starts at compressed offset " :
                                            " (decoder compressed offset ";
            text += std::to_string(status.compressed_offset) + ")";
        }
        return std::runtime_error(text);
    }

    void check_eof_marker() {
        if (checked_eof_) return;
        checked_eof_ = true;
        if (!src_->status().last_bgzf_empty) {
            const auto warning = error("BGZF EOF marker absent; input may be incomplete");
            std::fprintf(stderr, "flashalign: warning: %s\n", warning.what());
        }
    }

private:
    std::string path_;
    std::unique_ptr<ByteSource> src_;
    bool checked_eof_ = false;
};

}}}  // namespace fa::cpu::io
