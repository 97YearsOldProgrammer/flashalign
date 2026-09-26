// Value types produced by RNA coarse locus placement.
#pragma once

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

struct CoarseDiagonalPeak {
  int reference_id = -1;
  bool reverse = false;
  std::uint32_t oriented_query_begin = 0;
  std::uint32_t oriented_query_end = 0;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  std::int64_t median_diagonal = 0;
  std::uint32_t query_covered_bases = 0;
  std::uint32_t distinct_seed_support = 0;
};

struct CoarseLocus {
  std::vector<std::uint32_t> peak_indices;
  // Bounded intron-aware colinear score over the envelope's peaks; the RNA MAPQ reads it
  // as coarse confidence.
  std::int64_t score = 0;
  // Catalogue order key: how much of the read the envelope explains with seed matches.
  std::int64_t rank_score = 0;
  int reference_id = -1;
  bool reverse = false;
  std::uint32_t oriented_query_begin = 0;
  std::uint32_t oriented_query_end = 0;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  std::uint32_t intron_count = 0;
  // Union of the peaks' oriented query intervals, in read bases; a selector tie-break.
  std::int64_t union_query_coverage_bases = 0;
};

}  // namespace placement
}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
