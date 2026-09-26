#include "output.h"

#include <algorithm>
#include <climits>

namespace fa { namespace cpu { namespace lr { namespace rna { namespace placement {

bool project_fine_path_placement(
    const CoarseLocus& selected_locus,
    std::int32_t path_reference_begin, std::int32_t path_reference_end,
    std::int32_t path_query_begin, std::int32_t path_query_end,
    std::uint64_t selected_query_covered_bases,
    int mapq, int read_len, int reference_count,
    const std::vector<std::string>& reference_names,
    AlignResult& out) {
    const int ref_id = selected_locus.reference_id;
    // MAPQ 0 is a legal mapped row, as in minimap2; only a negative MAPQ is invalid.
    if (ref_id < 0 || ref_id >= reference_count || mapq < 0 ||
        path_query_begin < 0 || path_query_end <= path_query_begin ||
        path_query_end > read_len || path_reference_begin < 0 ||
        path_reference_end <= path_reference_begin ||
        selected_query_covered_bases == 0 ||
        selected_query_covered_bases >
            static_cast<std::uint64_t>(path_query_end - path_query_begin) ||
        selected_query_covered_bases > static_cast<std::uint64_t>(INT_MAX))
      return false;

    out.chromosome = reference_names[static_cast<size_t>(ref_id)];
    out.pos = path_reference_begin;
    out.target_end = path_reference_end;
    out.is_reverse = selected_locus.reverse;
    if (!out.is_reverse) {
      out.query_start = path_query_begin;
      out.query_end = path_query_end;
    } else {
      out.query_start = read_len - path_query_end;
      out.query_end = read_len - path_query_begin;
    }
    out.score = static_cast<int>(std::max<int64_t>(
        1, std::min<int64_t>(selected_locus.score, INT_MAX)));
    out.matches = static_cast<int>(selected_query_covered_bases);
    out.block_len = path_query_end - path_query_begin;
    out.mapq = mapq;
    return true;
}

}}}}}  // namespace fa::cpu::lr::rna::placement
