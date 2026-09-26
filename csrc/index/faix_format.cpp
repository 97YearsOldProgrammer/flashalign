// The header decoder, the load errors, the 4-bit reference codec and the name table.

#include "format.h"

#include "faix.h"
#include "../threading/parallel_for.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

namespace fa { namespace cpu {

namespace {

// One unpack task: [lo, hi) in global base coordinates within contig `chrom`. Large enough
// to amortize task overhead, small enough that one long contig spreads over all threads.
constexpr uint64_t kUnpackSliceBases = 32ULL << 20;

struct UnpackSlice {
    uint64_t chrom;
    uint64_t lo;
    uint64_t hi;
};

}  // namespace

FaixLoadError faix_header_decode(const char* p, FaixHeader& h) {
    if (std::memcmp(p, kFaixDevelopmentMagic, 8) == 0)
        return FaixLoadError::DevelopmentBuild;
    if (std::memcmp(p, kFaixMagic, 8) != 0)
        return FaixLoadError::NotIndex;
    std::memcpy(&h.version, p + 8, 4);
    if (h.version == 0)
        return FaixLoadError::InvalidHeader;
    if (h.version > kFaixFormatVersion)
        return FaixLoadError::NewerFormat;
    std::memcpy(&h.flags, p + 12, 4);
    std::memcpy(&h.k, p + 16, 4);
    std::memcpy(&h.syncmer_s, p + 20, 4);
    std::memcpy(&h.shard_bits, p + 24, 4);
    std::memcpy(&h.rank_block_bits, p + 28, 4);
    std::memcpy(&h.contig_count, p + 32, 8);
    std::memcpy(&h.total_bp, p + 40, 8);
    std::memcpy(&h.name_bytes, p + 48, 8);
    std::memcpy(&h.key_count, p + 56, 8);
    std::memcpy(&h.slot_count, p + 64, 8);
    std::memcpy(&h.multi_posting_count, p + 72, 8);
    std::memcpy(&h.part_count, p + 80, 8);
    std::memcpy(h.preset, p + kFaixPresetOffset, kFaixPresetBytes);
    // A later writer may use a flag or the zero bytes without a new version.
    if ((h.flags & ~kFaixKnownFlags) != 0)
        return FaixLoadError::NewerFormat;
    for (size_t i = 88; i < kFaixPresetOffset; ++i) {
        if (p[i] != 0)
            return FaixLoadError::NewerFormat;
    }

    const int k = static_cast<int>(h.k);
    if (h.k < 1 || h.k > static_cast<uint32_t>(kFaixMaxK) || h.syncmer_s < 1 ||
        h.syncmer_s > h.k || h.contig_count == 0 ||
        h.contig_count > kFaixMaxContigCount) {
        return FaixLoadError::InvalidHeader;
    }
    if ((h.flags & kFaixFlagMultipart) != 0) {
        if (h.part_count == 0 || h.part_count > h.contig_count ||
            h.shard_bits != 0 || h.rank_block_bits != 0 || h.key_count != 0 ||
            h.slot_count != 0 || h.multi_posting_count != 0) {
            return FaixLoadError::InvalidHeader;
        }
    } else if (h.part_count != 0 || h.shard_bits < 1 || h.shard_bits > 24 ||
               2 * k - static_cast<int>(h.shard_bits) > 31 ||
               h.rank_block_bits < 1 || h.rank_block_bits > 24 ||
               h.key_count > h.slot_count) {
        return FaixLoadError::InvalidHeader;
    }
    return FaixLoadError::None;
}

std::string faix_load_error_message(const FaixLoadStatus& status) {
    switch (status.error) {
      case FaixLoadError::None: return "no error";
      case FaixLoadError::OpenFailed: return "cannot open the file";
      case FaixLoadError::ReadFailed: return "cannot read the file";
      case FaixLoadError::OutOfMemory: return "not enough memory";
      case FaixLoadError::NotIndex: return "not a FlashAlign index";
      case FaixLoadError::DevelopmentBuild:
        return "an index from a FlashAlign development build; rebuild it with "
               "`flashalign index`";
      case FaixLoadError::NewerFormat:
        if (status.observed_version > kFaixFormatVersion) {
          return "written by a newer FlashAlign (index format " +
                 std::to_string(status.observed_version) +
                 "; this build reads format " +
                 std::to_string(kFaixFormatVersion) + ")";
        }
        return "written by a newer FlashAlign (index format " +
               std::to_string(kFaixFormatVersion) +
               " with features this build does not know)";
      case FaixLoadError::Multipart:
        return "a multi-part index (flashalign index -I), which only "
               "`flashalign align` reads";
      case FaixLoadError::Truncated: return "the file is truncated";
      case FaixLoadError::InvalidHeader: return "the header is damaged";
      case FaixLoadError::InvalidLayout: return "the file is damaged";
      case FaixLoadError::ContigTooLong:
        return "a reference sequence is longer than 2147483647 bp, the most SAM "
               "and BAM can address";
    }
    return "unknown error";
}

void pack_reference_4bit_into(
    const uint8_t* reference,
    uint64_t decoded_bytes,
    uint64_t decoded_offset,
    std::vector<uint8_t>& out
) {
    if (decoded_bytes == 0) return;
    if (!reference) return;
    const uint64_t needed = packed_reference_bytes(decoded_offset + decoded_bytes);
    if (out.size() < needed) out.resize(static_cast<size_t>(needed), 0);
    for (uint64_t i = 0; i < decoded_bytes; ++i) {
        const uint64_t pos = decoded_offset + i;
        const uint8_t base = reference[i] & 0x0fU;
        uint8_t& byte = out[static_cast<size_t>(pos >> 1)];
        if ((pos & 1ULL) == 0) {
            byte = static_cast<uint8_t>((byte & 0xf0U) | base);
        } else {
            byte = static_cast<uint8_t>((byte & 0x0fU) | static_cast<uint8_t>(base << 4));
        }
    }
}

void unpack_reference_4bit_to(
    const uint8_t* packed,
    uint64_t lo,
    uint64_t hi,
    uint8_t* out
) {
    if (hi <= lo) return;
    uint64_t i = lo;
    // An odd first base is the high nibble of its byte; peel it so the loop starts even.
    if ((i & 1ULL) != 0) {
        out[0] = static_cast<uint8_t>(
            (packed[static_cast<size_t>(i >> 1)] >> 4) & 0x0fU);
        ++i;
    }
    // Each source byte is two output bases, low nibble first; this loop vectorizes.
    const uint64_t pairs = (hi - i) >> 1;
    const uint8_t* src = packed + static_cast<size_t>(i >> 1);
    uint8_t* dst = out + static_cast<size_t>(i - lo);
    for (uint64_t p = 0; p < pairs; ++p) {
        const uint8_t byte = src[static_cast<size_t>(p)];
        dst[static_cast<size_t>(2 * p)] = static_cast<uint8_t>(byte & 0x0fU);
        dst[static_cast<size_t>(2 * p + 1)] = static_cast<uint8_t>(byte >> 4);
    }
    i += pairs * 2;
    // A leftover last base, at an even index, is the low nibble of its byte.
    if (i < hi) {
        out[static_cast<size_t>(i - lo)] = static_cast<uint8_t>(
            packed[static_cast<size_t>(i >> 1)] & 0x0fU);
    }
}

// Contigs are independent and a contig's output bases split into disjoint ranges, so the
// unpack runs in parallel: one pass sizes the contigs, one unpacks the slices.
bool FaixIndex::reference_chromosomes_u8(
    std::vector<std::vector<uint8_t>>& out,
    int n_threads
) const {
    out.clear();
    if (!has_reference_payload() || !chr_offsets_) return false;
    const size_t n_chrom = static_cast<size_t>(header_.contig_count);

    // Validate the whole offset table first, so a refusal leaves `out` empty.
    for (size_t i = 0; i < n_chrom; ++i) {
        if (chr_offsets_[i + 1] < chr_offsets_[i] ||
            chr_offsets_[i + 1] > header_.total_bp) {
            return false;
        }
    }

    std::vector<UnpackSlice> slices;
    for (size_t i = 0; i < n_chrom; ++i) {
        const uint64_t end = chr_offsets_[i + 1];
        for (uint64_t lo = chr_offsets_[i]; lo < end;
             lo += kUnpackSliceBases) {
            slices.push_back({static_cast<uint64_t>(i), lo,
                              std::min(lo + kUnpackSliceBases, end)});
        }
    }

    int threads = n_threads > 0 ? n_threads
                                : threading::default_thread_count();
    if (threads < 1) threads = 1;

    out.resize(n_chrom);
    // vector::resize zero-fills every base; spread that over the workers, one per contig.
    threading::parallel_for(
        threads, static_cast<int64_t>(n_chrom), [&](int64_t i, int) {
            const size_t c = static_cast<size_t>(i);
            out[c].resize(
                static_cast<size_t>(chr_offsets_[c + 1] - chr_offsets_[c]));
        });

    // Slices partition the output, so workers never write the same byte; neighbours may only
    // read the same packed byte.
    threading::parallel_for(
        threads, static_cast<int64_t>(slices.size()), [&](int64_t s, int) {
            const UnpackSlice& slice = slices[static_cast<size_t>(s)];
            const size_t c = static_cast<size_t>(slice.chrom);
            uint8_t* dst = out[c].data()
                + static_cast<size_t>(slice.lo - chr_offsets_[c]);
            unpack_reference_4bit_to(packed_reference_, slice.lo, slice.hi, dst);
        });
    return true;
}

bool parse_chrom_names(
    const char* data,
    uint64_t bytes,
    uint64_t count,
    std::vector<std::string>& out
) {
    if (!data && bytes > 0) return false;
    out.clear();
    out.reserve(static_cast<size_t>(count));
    uint64_t offset = 0;
    for (uint64_t i = 0; i < count; ++i) {
        if (offset + sizeof(uint32_t) > bytes) return false;
        uint32_t len = 0;
        std::memcpy(&len, data + offset, sizeof(len));
        offset += sizeof(uint32_t);
        if (offset + static_cast<uint64_t>(len) > bytes) return false;
        out.emplace_back(data + offset, data + offset + len);
        offset += static_cast<uint64_t>(len);
    }
    return offset == bytes;
}

bool validate_chrom_offsets(
    const std::vector<uint64_t>& offsets,
    uint64_t chrom_count,
    uint64_t total_bp
) {
    if (offsets.size() != static_cast<size_t>(chrom_count + 1)) return false;
    if (offsets.empty() || offsets.front() != 0) return false;
    for (size_t i = 1; i < offsets.size(); ++i) {
        if (offsets[i] < offsets[i - 1]) return false;
    }
    return offsets.back() == total_bp;
}

}}  // namespace fa::cpu
