// Coarse locus-envelope selection over RNA diagonal peaks.
#pragma once

#include "types.h"
#include "../../seeding/tie_hash.h" // tie_locus_hash (the exact-tie break)

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace placement {

struct CoarseLocusOptions {
  std::uint32_t diagonal_band = 24;
  std::uint32_t max_query_gap = 200;
  std::uint32_t max_query_overlap = 64;
  std::uint32_t min_intron = 20;
  std::uint32_t max_intron = 200000;
  // Reference gap across which two envelopes may still be coalesced: several
  // consecutive exons below the vote threshold can push the next peak beyond one legal
  // intron. Defaults to 4 * max_intron; 0 disables coalescing.
  std::uint32_t coalesce_gap = 800000;
  // Ceiling on a coalesced envelope's reference span, so the envelope plus the
  // +/- max_intron harvest pad fits skeleton_regions' 4 Mb locus bound. Defaults to
  // 16 * max_intron.
  std::uint32_t max_coalesced_span = 3200000;
  // At least one of two coalesced envelopes must carry this much distinct-seed support,
  // so two threshold-level peaks cannot fabricate a transcript. Set by the caller from
  // the vote threshold.
  std::uint32_t coalesce_min_support = 8;
  // Minimum oriented query extent of either piece of a coalesce (default 2 seed spans):
  // a narrower scrap is not an exon, and its few novel bases could lift a rival's
  // rank_score past the true locus.
  std::uint32_t coalesce_min_piece_span = 30;
  std::uint32_t max_predecessors = 64;
  std::uint32_t max_locus_chains = 6;
  // How many ranked envelopes beyond the cut the selector may additionally report in
  // `overflow_out`; 0 leaves the selection unchanged.
  std::uint32_t max_overflow_loci = 0;
  std::int32_t splice_open = 12;
  std::int32_t splice_log_scale = 2;
  std::int32_t query_gap_scale = 1;
  std::int32_t diagonal_drift_scale = 1;
  std::int32_t overlap_scale = 1;
};

CoarseLocusOptions make_coarse_locus_options(std::uint32_t seed_span,
                                             std::uint32_t min_intron,
                                             std::uint32_t max_intron);

// A locus's rank in an exact tie: minimap2's read-seeded hash (tie_locus_hash) of its
// contig slot, strand and reference start, lower first. The rival lifecycle breaks ties
// by catalogue index, so this hash decides every RNA tie between windows equal on every
// vote key. The start enters by its low 32 bits, as the DNA bin does.
inline std::uint64_t coarse_locus_tie_hash(std::uint32_t tie_seed,
                                           const CoarseLocus& locus) noexcept {
  return tie_locus_hash(tie_seed, locus.reference_id, locus.reverse,
                        static_cast<std::int64_t>(locus.reference_begin));
}

// Caller-owned reusable buffers for one worker's coarse locus selection: cleared per
// call, never shrunk, holding no state any decision reads.
struct CoarseLocusScratch {
  // Per-call.
  std::vector<std::uint32_t> order;
  std::vector<CoarseLocus> envelopes;
  std::vector<std::uint64_t> support;
  std::vector<std::int64_t> coverage;
  // One flag per surviving envelope slot: did it absorb another envelope in the
  // coalescing pass? Only those slots pay the evidence recompute.
  std::vector<std::uint8_t> coalesce_absorbed;
  std::vector<std::uint32_t> ranked;
  // Per-envelope temporaries; distinct buffers, so no helper can alias another.
  std::vector<std::uint32_t> query_order;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> spans;
  std::vector<std::int64_t> chain_score;
  std::vector<std::uint32_t> by_reference;
};

// Selects candidate locus envelopes for the downstream harvest. A reference-ordered sweep
// merges same-contig, same-strand peaks separated by at most one legal intron (covering
// one unpeaked exon); a second pass coalesces envelopes that are colinear pieces of one
// transcript within `coalesce_gap` (covering a run of unpeaked exons). The catalogue is
// ranked by seed evidence, then query coverage, support and peak count, then
// coarse_locus_tie_hash under `tie_seed`; contig slot and start settle only a hash
// collision. At most one envelope is kept per reference window, strand ignored.
//
// `coalesced`, when non-null, receives how many envelopes the second pass absorbed.
// `overflow_out`, when non-null and max_overflow_loci > 0, receives the next ranked
// distinct envelopes past the cut, unfinished (no colinear score, no intron estimate,
// peaks in reference order); the kept prefix is unaffected.
std::vector<CoarseLocus> select_coarse_loci(
    const std::vector<CoarseDiagonalPeak>& peaks,
    const CoarseLocusOptions& options, CoarseLocusScratch& scratch,
    std::uint32_t tie_seed, std::size_t* coalesced = nullptr,
    std::vector<CoarseLocus>* overflow_out = nullptr);

} // namespace placement
} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
