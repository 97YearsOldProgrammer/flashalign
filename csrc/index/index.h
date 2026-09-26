// Seed-index facade: the header DNA, RNA and seeding code include. The runtime (FaixIndex)
// is in faix.h, the file format in format.h and the builder in build.h.
#pragma once

#include "build.h"

namespace fa { namespace cpu {

using SeedIndex = FaixIndex;
using IndexView = SeedIndex;

// `image_offset` / `image_bytes` name one part of a multi-part container
// (FaixPartEntry); both 0 is the whole file, i.e. a single-part .faix.
inline SeedIndex load_seed_index(
    const std::string& path, FaixLoadStatus* status = nullptr,
    int n_threads = 0, uint64_t image_offset = 0, uint64_t image_bytes = 0) {
    return SeedIndex::load_image(path, image_offset, image_bytes, status,
                                 n_threads);
}

// Builds an index in memory. `status`, when supplied, names the refusal behind an empty
// result (see FaixBuildError).
inline SeedIndex build_seed_index(
    const std::vector<std::vector<uint8_t>>& chr_encs,
    const FaixBuildConfig& cfg,
    FaixBuildStatus* status = nullptr
) {
    return build_faix_index(chr_encs, cfg, status);
}

template <class EmitEncodedSequences>
inline SeedIndex build_seed_index_streamed(
    const std::vector<uint64_t>& chr_offsets,
    const FaixBuildConfig& cfg,
    EmitEncodedSequences&& emit_encoded_sequences,
    FaixBuildStatus* status = nullptr
) {
    return build_faix_index_streamed(
        chr_offsets,
        cfg,
        std::forward<EmitEncodedSequences>(emit_encoded_sequences), status);
}

}}  // namespace fa::cpu
