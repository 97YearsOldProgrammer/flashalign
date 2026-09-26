// Non-owning view of a loaded reference. The aligner owns every pointed-to object; a
// ReferenceContext is valid only while the aligner lives and must not escape a mapping call.
#pragma once

#include "../index/index.h" // fa::cpu::SeedIndex (FaixIndex alias)

#include <cstdint>
#include <string>
#include <vector>

namespace fa {
namespace cpu {
namespace engine {

struct ReferenceContext {
    const SeedIndex* index = nullptr;
    const std::vector<std::string>* names = nullptr;
    // One length per contig, in contig-id order, aligned with `names`. Read from the
    // index's offsets table, so it is present whether or not the bases are.
    const std::vector<int64_t>* lengths = nullptr;
    // The unpacked reference, one byte per base. Empty when the mode reads no bases (DNA
    // map-only); check has_bases() first, and take contig counts and lengths from
    // `lengths`.
    const std::vector<std::vector<uint8_t>>* encoded = nullptr;

    bool valid() const noexcept {
        return index != nullptr && names != nullptr && encoded != nullptr;
    }
    // The reference bases are materialized (CIGAR output, RNA, a FASTA target).
    bool has_bases() const noexcept {
        return encoded != nullptr && !encoded->empty();
    }
    int contig_count() const noexcept {
        return lengths == nullptr ? 0 : static_cast<int>(lengths->size());
    }
    // 0 for an id outside the table (a zero-length contig admits no seed span).
    int64_t contig_length(int contig) const noexcept {
        if (lengths == nullptr || contig < 0 || contig >= contig_count())
            return 0;
        return (*lengths)[static_cast<size_t>(contig)];
    }
};

} // namespace engine
} // namespace cpu
} // namespace fa
