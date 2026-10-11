// The re-vote over stretches of the query, the unexplained pieces of a read, the
// leftover harvest over the query no printed record covers, and the segment geometry the
// shadow rule and the pieces read. The backend's second coarse pass re-votes the whole
// read.
#pragma once

#include "placement/fused_capture.h" // FusedCaptureScratch (the re-vote's logs)
#include "placement/types.h"
#include "anchoring/query_span.h" // RnaQuerySpan (header-only arithmetic)
#include "anchoring/skeleton_harvest.h" // CandidatePool (the leftover's pools)

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Posting budget per read of the re-vote and the leftover harvest: past it they decline
// the read.
inline constexpr std::int64_t kRnaExplainMaxPostings = 65536;

// The geometry every colinearity decision reads.
struct RnaSegmentGeometry {
  int reference_id = -1;
  bool reverse = false;
  std::uint64_t reference_begin = 0;
  std::uint64_t reference_end = 0;
  RnaQuerySpan forward_span; // forward-frame query hull
};

// Same contig and strand, query order agreeing with reference order at both ends
// (mirrored on the reverse strand), and reference spans overlapping or within
// max_intron. Overlap is allowed because an exon the realizer skipped lies inside the
// primary's span; a piece that steps back on the reference (a duplicated exon, a
// back-splice) is a second place. Query overlap is not consulted; a censored span is
// never a continuation.
bool rna_segment_continuation(const RnaSegmentGeometry& a,
                              const RnaSegmentGeometry& b,
                              int max_intron) noexcept;

// Same contig and strand, and reference spans that overlap or continue each other: two
// descriptions of one segment, never competing (not a rival or a demoted vote).
bool rna_same_place(const RnaSegmentGeometry& a, const RnaSegmentGeometry& b,
                    int max_intron) noexcept;

// Caller-owned reusable buffers: cleared per call, never shrunk.
struct RnaExplainScratch {
  // The re-vote's own capture area (the two strand logs and bin counters the
  // fused capture's drain reads; its view/selection members stay empty).
  placement::FusedCaptureScratch revote_capture;
  // A realized family's segments, unexplained_pieces's sorted copy, and the
  // unexplained pieces of the query.
  std::vector<RnaQuerySpan> segments, sorted, pieces;
};

// The re-vote's own accounting: how many extracted seeds were re-voted, how
// many postings they hold, and whether that exceeded kRnaExplainMaxPostings
// (in which case nothing at all was voted and both counters are 0).
struct RnaRevoteStats {
  int seeds_revoted = 0;
  std::int64_t postings = 0;
  bool over_budget = false;
  // Records logged per strand, read off the logs. `postings_inspected` counts postings
  // that reached the strand test, the same for both strands since each posting is
  // tested for both.
  std::int64_t emitted_forward = 0;
  std::int64_t emitted_reverse = 0;
  std::int64_t postings_inspected = 0;
  std::int64_t rejected_forward = 0;
  std::int64_t rejected_reverse = 0;
};

// Re-votes every extracted seed whose span lies inside one of `intervals`, selected or
// not, with the fused capture's admission (a held, non-empty posting list within the
// occurrence policy) and bin key, on both strand frames. `views` is the capture's
// posting view per forward seed, so nothing probes the index. The posting budget is
// checked before any posting is walked; over it the read is declined. Returns the peaks
// drained at `min_support`; `capture` is cleared and used as scratch.
std::vector<placement::CoarseDiagonalPeak> rna_revote_intervals(
    const LongReadSeedContext& seed_ctx,
    const std::vector<QuerySeed>& forward_seeds,
    const std::vector<KmerPostingView>& views,
    const std::vector<RnaQuerySpan>& intervals, int read_len, int k,
    int min_support, int n_chr, placement::FusedCaptureScratch& capture,
    RnaRevoteStats& stats);

// The unexplained pieces of the query: the hull's terminal clips, then the gaps
// between `segments` in query order (a running reach, so nested or overlapping
// segments leave no false gap), each at least `min_bp` wide (0 keeps every non-empty
// piece). `sorted` is scratch.
void unexplained_pieces(RnaQuerySpan hull,
                        const std::vector<RnaQuerySpan>& segments, int read_len,
                        int min_bp, std::vector<RnaQuerySpan>& sorted,
                        std::vector<RnaQuerySpan>& out);

// The leftover's intervals: every forward-frame interval of at least `min_bp` that
// neither `primary_segments` nor `explained` covers, the clips outside their joint hull
// with `primary_hull` included.
void rna_leftover_intervals(RnaQuerySpan primary_hull,
                            const std::vector<RnaQuerySpan>& primary_segments,
                            const std::vector<RnaQuerySpan>& explained,
                            int read_len, int min_bp,
                            std::vector<RnaQuerySpan>& out);

// minimap2's chain floor (min_cnt, min_chain_score), on the leftover's chains only.
inline constexpr int kRnaLeftoverChainMinCount = 3;
inline constexpr int kRnaLeftoverChainMinScore = 40;

// One contig and lane's anchors from rna_leftover_pools, in that lane's query frame.
struct RnaLeftoverPool {
  int reference_id = -1;
  bool reverse = false;
  CandidatePool pool;
};

// The leftover harvest: rna_revote_intervals's walk over `intervals` (every extracted
// seed whose k-mer lies inside one, under the vote's occurrence gate, each posting in the
// lane its orientation bit allows), each posting an anchor in its lane's frame instead of
// a vote. The anchors are pooled per contig and lane in reference, then query order,
// duplicates removed. Over kRnaExplainMaxPostings postings nothing is pooled and false is
// returned.
bool rna_leftover_pools(const LongReadSeedContext& seed_ctx,
                        const std::vector<QuerySeed>& forward_seeds,
                        const std::vector<KmerPostingView>& views,
                        const std::vector<RnaQuerySpan>& intervals, int read_len,
                        int k, int n_chr, std::vector<RnaLeftoverPool>& pools);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
