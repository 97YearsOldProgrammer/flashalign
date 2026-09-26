// The query partition, run after the primary's MAPQ is fixed. The primary is its
// realized segments and everything else on the read is evidence. The decision is an
// explanation of the query: an assignment of the read's 128 forward-frame tiles to
// segments (one contig, strand and spliced chain each). The output is a list of windows,
// stretches of the unexplained query one candidate owns, which the caller harvests,
// chains, realizes, prices and emits as second families.
//
//   E0  unexplained intervals: the complement of the primary's realized segments,
//       pieces >= kRnaChimeraMinQueryBases, longest first, at most
//       kRnaExplainMaxIntervals
//   E1  re-vote: every extracted seed inside those intervals walks its posting view
//       into the vote's bins; a peak needs kRnaCoprimaryMinAnchors distinct seeds
//   E2  units: colinear staircases (the selector's recurrence with a backtrack) over
//       the read peaks query-disjoint from the realized hull and outside the re-voted
//       intervals, plus the re-vote peaks
//   E3  candidates: the primary (realized-segment tiles, pinned), every non-elected
//       chained hypothesis (anchor tiles) and every unit (peak tiles); every
//       non-primary mask has the primary's tiles subtracted
//   E3' contenders: a unit whose tiles compete with another candidate's is chained by
//       the caller over its unexplained piece before the solve, so a contested
//       interval is decided chain against chain
//   E4  solve: solve_query_partition, once
//   E5  windows: per owner and unexplained piece, the hull of the owner's blocks
//       clipped to the piece; continuing owners merge; a window continuing the primary
//       is refused; one narrower than kRnaChimeraMinQueryBases is not harvested
//   E6  emission (backend.cpp): every nominee, strongest chain first
#pragma once

#include "placement/coarse_chain.h" // CoarseLocusOptions (the staircase DP recipe)
#include "placement/fused_capture.h" // FusedCaptureScratch (the re-vote's logs)
#include "placement/types.h"
#include "anchoring/exact_anchor_path.h" // ExactAnchorPath (chain masks)
#include "anchoring/query_span.h" // RnaQuerySpan (header-only arithmetic)
#include "anchoring/selected_locus_anchoring.h" // Rank1HarvestResult
#include "../voting/query_partition.h"

#include <cstdint>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// The solver's own ceiling; voting/query_partition.cpp refuses a larger catalogue.
inline constexpr int kRnaPartitionCatalogueBound =
    2 * ::fa::cpu::voting::kCatalogueLaneBound; // 32
// The solver's own block floor: a candidate with fewer supported tiles left
// after the primary's are subtracted can never own a block and is not offered.
inline constexpr int kRnaPartitionMinBlockTiles = 2;
// Unexplained intervals re-voted per read, longest first.
inline constexpr int kRnaExplainMaxIntervals = 2;
// Re-vote posting budget per read: past it the re-vote declines the read.
inline constexpr std::int64_t kRnaExplainMaxPostings = 65536;
// Windows harvested per read (families cap at kRnaChimeraMaxFamilies).
inline constexpr int kRnaExplainMaxWindows = 4;
// Units chained before the solve because they compete for tiles with another candidate
// (rna_explain_contender), strongest vote first, one per place.
inline constexpr int kRnaExplainMaxContenders = 2;
// kRnaChimeraMinQueryBases (rival_lifecycle.h) is the segment floor throughout: a
// narrower unexplained piece is not re-voted, a narrower window is not harvested, and a
// realized family with fewer exclusive bases is not emitted.

enum class RnaSegmentKind : std::uint8_t {
  Primary = 0,
  Chained = 1,
  Unit = 2,
  Revote = 3
};

// The geometry every colinearity decision reads. One predicate, used for
// merging adjacent windows and for refusing a window that continues the
// primary.
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
// descriptions of one segment, never competing (not a contender, a rival or a demoted
// vote).
bool rna_same_place(const RnaSegmentGeometry& a, const RnaSegmentGeometry& b,
                    int max_intron) noexcept;

struct RnaSegmentCandidate {
  RnaSegmentKind kind = RnaSegmentKind::Unit;
  RnaSegmentGeometry geometry;
  // Forward-frame tiles, the primary's tiles already subtracted (the primary
  // itself carries its own).
  ::fa::cpu::voting::QueryTileMask support;
  // Query bases explained by vote evidence: a hypothesis's locus rank_score,
  // a unit's union of its peaks' forward spans.
  std::int64_t vote_evidence = 0;
  // Anchors when chained, summed distinct seed support otherwise.
  std::int64_t chain_evidence = 0;
  int hypothesis_index = -1; // Chained (and a contender after the caller
                             // chained it)
  int staircase_index = -1;  // Unit / Revote: index into
                             // RnaExplainResult::staircases
  // Chained: the hypothesis's locus, copied into any window it owns.
  const placement::CoarseLocus* locus = nullptr;
  bool offered = false;    // reached the solver (after pruning)
  int catalogue_slot = -1; // solver slot when offered
};

struct RnaExplainWindow {
  RnaQuerySpan span; // forward frame, clipped to its unexplained interval,
                     // >= kRnaChimeraMinQueryBases
  RnaQuerySpan piece; // the unexplained piece of the query the window lies
                      // in: the harvest's query window (the chain decides
                      // the segment's extent; the solve's tiles only say
                      // where to look)
  int tile_begin = 0; // the hull of the owner's blocks, in tiles
  int tile_end = 0;
  std::vector<int> owners; // the owner, then any candidate merged into the
                           // window as its continuation
  RnaSegmentKind kind = RnaSegmentKind::Unit; // the first owner's kind
  int candidate_index = -1;  // owner in RnaExplainResult::candidates
  int hypothesis_index = -1; // the owner's hypothesis when it has one (a
                             // chained hypothesis, or a unit the caller
                             // chained as a contender): the harvest's twin
  // What to harvest: a copy of the hypothesis's locus (reference window
  // intact, peak_indices as is), or built from the window's peaks in query
  // order with the union of their forward spans as rank_score.
  placement::CoarseLocus locus;
  // Harvest skeleton = RnaExplainResult::unit_peaks (Unit / Revote owners)
  // else the read's own peaks (a Chained owner's locus indexes those).
  bool unit_peaks = false;
  int supporting_tiles = 0;
  // Largest vote_evidence among the other candidates sharing more than half of the
  // shorter tile mask with the window, excluding the primary, the window's owners and
  // candidates at the window's own place (rna_same_place).
  std::int64_t competing_vote = 0;
};

// One colinear chain of peaks = one candidate segment.
struct RnaStaircase {
  int reference_id = -1;
  bool reverse = false;
  std::uint64_t reference_begin = 0; // hull over its peaks
  std::uint64_t reference_end = 0;
  std::uint32_t oriented_query_begin = 0;
  std::uint32_t oriented_query_end = 0;
  std::uint32_t forward_query_begin = 0;
  std::uint32_t forward_query_end = 0;
  // The colinear DP score of the chain's own peaks (query reward minus
  // penalties, from the peak the walk stopped at); a post-merged staircase
  // carries the sum of its pieces, the gap the DP refused uncharged.
  std::int64_t score = 0;
  std::uint64_t seed_support = 0;          // summed distinct_seed_support
  std::vector<std::uint32_t> peak_indices; // in reference order
};

// The staircase DP's per-read buffers: the peak subset in group order, each
// position's colinear score and best predecessor, the extraction order and the
// assigned flags, the post-merge order.
struct RnaStaircaseScratch {
  std::vector<std::uint32_t> order;
  std::vector<std::int64_t> score;
  std::vector<std::int32_t> prev;
  std::vector<std::uint32_t> extract;
  std::vector<std::uint8_t> assigned;
  std::vector<std::uint32_t> merge_order;
};

// Caller-owned reusable buffers: cleared per call, never shrunk.
struct RnaExplainScratch {
  ::fa::cpu::voting::QueryPartitionProblem problem;
  // The re-vote's own capture area (the two strand logs and bin counters the
  // fused capture's drain reads; its view/selection members stay empty).
  placement::FusedCaptureScratch revote_capture;
  RnaStaircaseScratch staircase;
  // The caller's realized segments, the interval rule's sorted copy, and the
  // unexplained pieces of the query (the clip intervals of E5).
  std::vector<RnaQuerySpan> segments, sorted, pieces;
  std::vector<std::uint32_t> window_peaks;
  std::vector<std::size_t> offer_order, offered_index;
  std::vector<::fa::cpu::voting::QueryTileMask> masks;
};

struct RnaExplainResult {
  // E0 / E1
  std::vector<RnaQuerySpan> intervals; // re-voted, longest first,
                                       // <= kRnaExplainMaxIntervals
  int seeds_revoted = 0;
  std::int64_t postings = 0; // postings walked (0 when over budget)
  bool over_budget = false;  // the admissible postings exceeded the budget:
                             // no re-vote at all, read peaks only
  // Copied from RnaRevoteStats.
  std::int64_t emitted_forward = 0;
  std::int64_t emitted_reverse = 0;
  std::int64_t postings_inspected = 0;
  std::int64_t rejected_forward = 0;
  std::int64_t rejected_reverse = 0;
  // E2
  // Residue read peaks (query-disjoint from the realized hull, outside every
  // re-voted interval) ++ the re-vote's peaks. Every unit's peak_indices and
  // every Unit/Revote window's locus index this vector, not the read's.
  std::vector<placement::CoarseDiagonalPeak> unit_peaks;
  int revote_peaks = 0; // how many of unit_peaks came from the re-vote (the
                        // tail)
  std::vector<RnaStaircase> staircases; // over unit_peaks
  // E3
  std::vector<RnaSegmentCandidate> candidates; // [0] = primary; then
                                               // chained; then units
  // E4
  bool solved = false;
  int candidates_offered = 0; // the primary included
  int chained_offered = 0;
  int units_offered = 0; // Unit + Revote
  std::uint64_t dp_cells = 0;
  // E5
  std::vector<RnaExplainWindow> windows; // query order,
                                         // <= kRnaExplainMaxWindows
  int windows_merged = 0;       // blocks that extended a window of another
                                // owner (a continuation)
  int continuation_refused = 0; // windows refused as the primary's own
                                // continuation
  int read_len = 0;             // the read the peaks and spans above belong to
  int max_intron = 0; // the intron ceiling, the continuation rule's gap
};

// A non-elected hypothesis the lifecycle chained: its locus, its bundle, its
// post-repair strand, its coarse vote and its chain span in its own frame.
struct RnaExplainChained {
  int hypothesis_index = -1;
  const placement::CoarseLocus* locus = nullptr;
  const ExactAnchorPath* bundle = nullptr;
  bool reverse = false;
  std::int64_t rank_score = 0;
  int q_begin = -1;
  int q_end = -1;
};

// The re-vote, also run by the backend's second coarse pass on the whole read.

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

// Builds E0..E3 (no solve). `primary` is the realized primary: contig, strand,
// realized reference span [pos, target_end), realized forward-frame hull;
// `primary_rank_score` its locus's coarse vote; `primary_segments` its
// realized segments in the forward frame; `chained` every non-elected chained
// hypothesis. `peaks` is the read's capped peak vector, `forward_seeds` the
// read's extracted forward stream and `views` the fused capture's posting view
// per seed, index-aligned. Leaves `result` empty (no candidates) when the
// primary hull is unmeasurable.
void rna_explain_collect(
    const LongReadSeedContext& seed_ctx,
    const std::vector<QuerySeed>& forward_seeds,
    const std::vector<KmerPostingView>& views,
    const std::vector<placement::CoarseDiagonalPeak>& peaks,
    const placement::CoarseLocusOptions& options,
    const RnaSegmentGeometry& primary, std::int64_t primary_rank_score,
    const std::vector<RnaQuerySpan>& primary_segments,
    const std::vector<RnaExplainChained>& chained, int read_len, int seed_len,
    int min_support, int max_intron, RnaExplainScratch& scratch,
    RnaExplainResult& result);

// Replace a Unit candidate's evidence by its chain after the caller chained it
// (E3'): the hypothesis index, the anchors' tiles minus the primary's, the
// anchor count, and the chain's own geometry. The kind stays Unit / Revote.
void rna_explain_candidate_chained(
    RnaExplainResult& result, std::size_t candidate_index, int hypothesis_index,
    const ExactAnchorPath& bundle, bool reverse, int read_len, int seed_len,
    const ::fa::cpu::voting::QueryTileMask& primary_tiles);

// The unexplained piece of the query `span` overlaps most (unmeasured when
// none does). Valid after rna_explain_collect: the windows' clip intervals
// and every harvest's query window.
RnaQuerySpan rna_explain_piece(const RnaExplainScratch& scratch,
                               RnaQuerySpan span) noexcept;

// Whether unit U competes for tiles with a candidate at another place (mask level 0.5 on
// the shorter popcount); a candidate at the unit's own place is the same segment.
// Returns the first such candidate's index, or -1.
int rna_explain_contender(const RnaExplainResult& result,
                          std::size_t candidate_index);

// A locus built from a unit's staircase peaks in query order (the window
// rule), indexing RnaExplainResult::unit_peaks.
placement::CoarseLocus rna_explain_unit_locus(const RnaExplainResult& result,
                                              int staircase_index);

// E4 + E5: offer, solve once, cut windows.
void rna_explain_solve(const placement::CoarseLocusOptions& options,
                       const RnaSegmentGeometry& primary,
                       const std::vector<RnaQuerySpan>& primary_segments,
                       int read_len, int seed_len, int max_intron,
                       RnaExplainScratch& scratch, RnaExplainResult& result);

// A nominee chains only its own query. harvest_rank1 pads the reference window by a full
// maximum intron on each side, which suits a catalogue rival but not a nominee: a
// residue nominee lies within one intron of the committed envelope, so its padded window
// holds the committed chain's anchors and the anchor path would re-find it. Every anchor
// whose forward-frame interval leaves `forward_window` (widened by one seed each side)
// has its geometry verdict cleared, and the pool is refinalized under the same caps. The
// forward interval is invariant under the frame projection, so the restriction survives
// mirror_candidate_pool_query_frame. `pool_reverse` is the frame the pool was harvested
// in; an unmeasurable `forward_window` refuses the harvest.
void restrict_harvest_to_query_window(Rank1HarvestResult& harvest,
                                      RnaQuerySpan forward_window,
                                      bool pool_reverse, int read_len,
                                      int seed_len);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
