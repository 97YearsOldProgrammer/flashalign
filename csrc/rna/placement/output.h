#pragma once

#include "types.h"
#include "../../core/types.h"

#include <cstdint>
#include <string>
#include <vector>

namespace fa { namespace cpu { namespace lr { namespace rna { namespace placement {

// Fills `out` from the selected locus and its exact anchor path. Returns false when the
// path geometry or the MAPQ is out of range; there is no coarse fallback projection.
bool project_fine_path_placement(
    const CoarseLocus& selected_locus,
    std::int32_t path_reference_begin, std::int32_t path_reference_end,
    std::int32_t path_query_begin, std::int32_t path_query_end,
    std::uint64_t selected_query_covered_bases,
    int mapq, int read_len, int reference_count,
    const std::vector<std::string>& reference_names,
    AlignResult& out);

}}}}}  // namespace fa::cpu::lr::rna::placement
