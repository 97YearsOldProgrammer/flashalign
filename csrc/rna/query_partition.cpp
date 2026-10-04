#include "query_partition.h"

#include "anchoring/query_span.h"        // rna_forward_query_span (header-only)
#include "placement/coarse_transition.h" // the selector's transition rule
#include "realization/rival_lifecycle.h" // the chimeric emission geometry
#include "../index/format.h" // packed_ref_contig / packed_ref_local (the re-vote)
#include "../seeding/context.h" // effective_vote_diag_bin_width
#include "../voting/query_tiles.h"
#include "../voting/state.h" // vote_pack_key / vote_floor_div / occ policy

#include <algorithm>
#include <climits>
#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

namespace {

namespace voting = ::fa::cpu::voting;

// Evidence fields are int; saturate rather than wrap, which would reorder the solver's
// tie-breaks.
int clamp_evidence(std::int64_t value) noexcept {
  if (value <= 0)
    return 0;
  return value > static_cast<std::int64_t>(INT_MAX) ? INT_MAX
                                                    : static_cast<int>(value);
}

bool span_measured(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

bool span_contains(RnaQuerySpan interval, int position) noexcept {
  return position >= interval.begin && position < interval.end;
}

bool same_span(RnaQuerySpan a, RnaQuerySpan b) noexcept {
  return a.begin == b.begin && a.end == b.end;
}

RnaQuerySpan clip_to(RnaQuerySpan span, RnaQuerySpan interval) noexcept {
  span.begin = std::max(span.begin, interval.begin);
  span.end = std::min(span.end, interval.end);
  return span;
}

// The block's extent in read bases, forward frame. The tile grid is over seed start
// positions, so the last tile adds the seed's span (query_tile_end), and a block ending
// at tile 128 ends at the read's end.
RnaQuerySpan block_forward_span(int tile_begin, int tile_end, int read_len,
                                int seed_len) noexcept {
  RnaQuerySpan span;
  span.begin = voting::query_tile_begin(tile_begin, read_len, seed_len);
  span.end = voting::query_tile_end(tile_end, read_len, seed_len);
  return span;
}

// Sets every tile a forward-frame interval touches. The mask is over seed start
// positions, so the covered starts are [begin, end - seed_len], and a span shorter than
// one seed still holds the start it begins at. A censored span sets nothing.
void set_span_tiles(voting::QueryTileMask& mask, RnaQuerySpan span,
                    int read_len, int seed_len) noexcept {
  if (span.begin < 0 || span.end <= span.begin)
    return;
  const int last_start = std::max(span.begin, span.end - seed_len);
  const int first_tile =
      voting::query_tile_for_position(span.begin, read_len, seed_len);
  const int last_tile =
      voting::query_tile_for_position(last_start, read_len, seed_len);
  for (int tile = first_tile; tile <= last_tile; ++tile)
    mask.set(tile);
}

// The 128-tile problem for one read. A tile is valid when it holds at least one seed
// start, as in DNA placement.
void reset_partition_problem(voting::QueryPartitionProblem& problem,
                             const voting::QueryPartitionParameters& parameters,
                             int read_len, int seed_len) {
  problem.tile_count = voting::kQueryTileCount;
  problem.valid_tiles = voting::QueryTileMask{};
  for (int tile = 0; tile < problem.tile_count; ++tile) {
    if (voting::query_tile_end(tile + 1, read_len, seed_len) >
        voting::query_tile_begin(tile, read_len, seed_len))
      problem.valid_tiles.set(tile);
  }
  problem.parameters = parameters;
  problem.rival = voting::QueryPartitionRival::WinnerOnly;
  problem.catalogue.candidates.clear();
}

// The valid tiles the primary's mask does not claim: what every other
// candidate's support is confined to.
voting::QueryTileMask
tiles_outside(const voting::QueryTileMask& incumbent,
              const voting::QueryTileMask& valid) noexcept {
  return voting::QueryTileMask{{~incumbent.words[0] & valid.words[0],
                                ~incumbent.words[1] & valid.words[1]}};
}

// One candidate per distinct mask: a later candidate with exactly the mask
// of an earlier one could own a block only by the solver's tie-break, which
// the earlier (better-ranked) one wins. It contributes nothing but tail walks.
bool duplicate_mask(const voting::QueryPartitionProblem& problem,
                    const voting::QueryTileMask& support) noexcept {
  for (const voting::QueryCandidate& held : problem.catalogue.candidates)
    if (held.support == support)
      return true;
  return false;
}

RnaQuerySpan peak_forward_span(const placement::CoarseDiagonalPeak& peak,
                               int read_len) noexcept {
  return rna_forward_query_span(static_cast<int>(peak.oriented_query_begin),
                                static_cast<int>(peak.oriented_query_end),
                                peak.reverse, read_len);
}

// E0: the unexplained pieces of the query: the realized hull's terminal clips, then the
// gaps between realized segments in query order (a running reach, so nested or
// overlapping segments leave no false gap), each at least `min_bp` wide (0 keeps every
// non-empty piece).
void unexplained_pieces(RnaQuerySpan hull,
                        const std::vector<RnaQuerySpan>& segments, int read_len,
                        int min_bp, std::vector<RnaQuerySpan>& sorted,
                        std::vector<RnaQuerySpan>& out) {
  out.clear();
  const int floor = std::max(1, min_bp);
  const auto emit = [&](int begin, int end) {
    if (end - begin >= floor)
      out.push_back(RnaQuerySpan{begin, end});
  };
  emit(0, hull.begin);
  emit(hull.end, read_len);
  sorted.clear();
  for (const RnaQuerySpan& segment : segments)
    if (span_measured(segment))
      sorted.push_back(segment);
  std::stable_sort(sorted.begin(), sorted.end(),
                   [](RnaQuerySpan a, RnaQuerySpan b) {
                     if (a.begin != b.begin)
                       return a.begin < b.begin;
                     return a.end < b.end;
                   });
  int reach = 0;
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    if (i > 0 && sorted[i].begin > reach)
      emit(reach, sorted[i].begin);
    reach = std::max(reach, sorted[i].end);
  }
}

// Two staircases of one contig and strand are colinear when the later one in the
// reference starts past the earlier one's end within one legal intron and the read
// agrees with that order (it begins no more than max_query_overlap before the earlier
// one ends). The post-merge's test.
bool staircases_colinear(
    const RnaStaircase& a, const RnaStaircase& b,
    const placement::CoarseLocusOptions& options) noexcept {
  if (a.reference_id != b.reference_id || a.reverse != b.reverse)
    return false;
  const RnaStaircase& earlier = a.reference_begin <= b.reference_begin ? a : b;
  const RnaStaircase& later = a.reference_begin <= b.reference_begin ? b : a;
  if (&earlier == &later)
    return false;
  if (later.reference_begin <= earlier.reference_end)
    return false;
  if (later.reference_begin - earlier.reference_end >
      static_cast<std::uint64_t>(options.max_intron))
    return false;
  return static_cast<std::int64_t>(later.oriented_query_begin) +
             static_cast<std::int64_t>(options.max_query_overlap) >=
         static_cast<std::int64_t>(earlier.oriented_query_end);
}

// The union of the peaks' forward-frame tiles.
void staircase_tiles(const std::vector<placement::CoarseDiagonalPeak>& peaks,
                     const std::vector<std::uint32_t>& indices,
                     voting::QueryTileMask& mask, int read_len,
                     int seed_len) noexcept {
  for (const std::uint32_t peak_index : indices)
    set_span_tiles(mask, peak_forward_span(peaks[peak_index], read_len),
                   read_len, seed_len);
}

// The selector's own (contig, strand, reference) sort over the peak subset
// `area.order` holds on entry (any order), its colinear recurrence with a
// backtrack, greedy extraction, the one-intron post-merge, forward spans,
// score order. Every peak of the subset ends in exactly one staircase of
// `staircases`, which is cleared first.
void build_staircases(const std::vector<placement::CoarseDiagonalPeak>& peaks,
                      const placement::CoarseLocusOptions& options,
                      int read_len, RnaStaircaseScratch& area,
                      std::vector<RnaStaircase>& staircases) {
  staircases.clear();
  // 1. Group by (contig, strand), reference order within: the selector's own
  // sort, so the recurrence below sees peaks exactly as it sees them there.
  std::vector<std::uint32_t>& order = area.order;
  std::stable_sort(order.begin(), order.end(),
                   [&peaks](std::uint32_t a, std::uint32_t b) {
                     const placement::CoarseDiagonalPeak& pa = peaks[a];
                     const placement::CoarseDiagonalPeak& pb = peaks[b];
                     if (pa.reference_id != pb.reference_id)
                       return pa.reference_id < pb.reference_id;
                     if (pa.reverse != pb.reverse)
                       return pa.reverse < pb.reverse;
                     if (pa.reference_begin != pb.reference_begin)
                       return pa.reference_begin < pb.reference_begin;
                     if (pa.reference_end != pb.reference_end)
                       return pa.reference_end < pb.reference_end;
                     return a < b;
                   });

  // 2. The selector's colinear recurrence (best_colinear_score), with the
  // predecessor recorded: score[i] seeds at the peak's own covered bases and
  // the best of the last max_predecessors valid transitions extends it.
  const std::size_t n = order.size();
  std::vector<std::int64_t>& score = area.score;
  std::vector<std::int32_t>& prev = area.prev;
  score.assign(n, 0);
  prev.assign(n, -1);
  std::size_t group_begin = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const placement::CoarseDiagonalPeak& pi = peaks[order[i]];
    const placement::CoarseDiagonalPeak& first = peaks[order[group_begin]];
    if (pi.reference_id != first.reference_id || pi.reverse != first.reverse)
      group_begin = i;
    score[i] = static_cast<std::int64_t>(pi.query_covered_bases);
    std::uint32_t scanned = 0;
    for (std::size_t j = i;
         j-- > group_begin && scanned < options.max_predecessors;) {
      ++scanned;
      const placement::CoarseDiagonalPeak& pj = peaks[order[j]];
      const placement::detail::Transition t =
          placement::detail::score_transition(pj, pi, options);
      if (t.kind == placement::detail::TransitionKind::kInvalid)
        continue;
      const std::int64_t candidate =
          score[j] + placement::detail::incremental_query_reward(pj, pi) -
          t.penalty;
      if (candidate > score[i]) {
        score[i] = candidate;
        prev[i] = static_cast<std::int32_t>(j);
      }
    }
  }
  // Greedy extraction: the best-scoring unassigned peak heads a staircase and
  // the walk back along `prev` takes every unassigned predecessor, so each
  // peak ends in exactly one staircase, singletons included. The score is
  // what the walked peaks contributed: the head's DP score less the score of
  // the assigned peak the walk stopped at.
  std::vector<std::uint32_t>& extract = area.extract;
  extract.resize(n);
  for (std::size_t i = 0; i < n; ++i)
    extract[i] = static_cast<std::uint32_t>(i);
  std::stable_sort(extract.begin(), extract.end(),
                   [&score](std::uint32_t a, std::uint32_t b) {
                     if (score[a] != score[b])
                       return score[a] > score[b];
                     return a < b;
                   });
  std::vector<std::uint8_t>& assigned = area.assigned;
  assigned.assign(n, 0);
  for (const std::uint32_t head : extract) {
    if (assigned[head])
      continue;
    RnaStaircase staircase;
    std::int32_t cursor = static_cast<std::int32_t>(head);
    while (cursor >= 0 && !assigned[static_cast<std::size_t>(cursor)]) {
      assigned[static_cast<std::size_t>(cursor)] = 1;
      staircase.peak_indices.push_back(order[static_cast<std::size_t>(cursor)]);
      cursor = prev[static_cast<std::size_t>(cursor)];
    }
    staircase.score =
        score[head] -
        (cursor >= 0 ? score[static_cast<std::size_t>(cursor)] : 0);
    // Walked head to root, which is reference-descending; the staircase
    // carries its peaks in reference order.
    std::reverse(staircase.peak_indices.begin(), staircase.peak_indices.end());
    bool first = true;
    for (const std::uint32_t peak_index : staircase.peak_indices) {
      const placement::CoarseDiagonalPeak& peak = peaks[peak_index];
      staircase.seed_support += peak.distinct_seed_support;
      if (first) {
        staircase.reference_id = peak.reference_id;
        staircase.reverse = peak.reverse;
        staircase.reference_begin = peak.reference_begin;
        staircase.reference_end = peak.reference_end;
        staircase.oriented_query_begin = peak.oriented_query_begin;
        staircase.oriented_query_end = peak.oriented_query_end;
        first = false;
        continue;
      }
      staircase.reference_begin =
          std::min(staircase.reference_begin, peak.reference_begin);
      staircase.reference_end =
          std::max(staircase.reference_end, peak.reference_end);
      staircase.oriented_query_begin =
          std::min(staircase.oriented_query_begin, peak.oriented_query_begin);
      staircase.oriented_query_end =
          std::max(staircase.oriented_query_end, peak.oriented_query_end);
    }
    staircases.push_back(std::move(staircase));
  }

  // 3. Post-merge, once per group in reference order: the recurrence's query-gap gate
  // splits one transcript at every unpeaked exon, so two consecutive staircases colinear
  // within one intron become one.
  std::vector<std::uint32_t>& merge_order = area.merge_order;
  merge_order.resize(staircases.size());
  for (std::size_t i = 0; i < staircases.size(); ++i)
    merge_order[i] = static_cast<std::uint32_t>(i);
  std::stable_sort(merge_order.begin(), merge_order.end(),
                   [&staircases](std::uint32_t a, std::uint32_t b) {
                     const RnaStaircase& sa = staircases[a];
                     const RnaStaircase& sb = staircases[b];
                     if (sa.reference_id != sb.reference_id)
                       return sa.reference_id < sb.reference_id;
                     if (sa.reverse != sb.reverse)
                       return sa.reverse < sb.reverse;
                     if (sa.reference_begin != sb.reference_begin)
                       return sa.reference_begin < sb.reference_begin;
                     if (sa.reference_end != sb.reference_end)
                       return sa.reference_end < sb.reference_end;
                     return sa.peak_indices.front() < sb.peak_indices.front();
                   });
  {
    std::int64_t current = -1;
    for (const std::uint32_t index : merge_order) {
      if (current >= 0 &&
          staircases_colinear(staircases[static_cast<std::size_t>(current)],
                              staircases[index], options)) {
        RnaStaircase& target = staircases[static_cast<std::size_t>(current)];
        RnaStaircase& piece = staircases[index];
        // Every peak of the piece begins past the target's reference end, so
        // appending keeps the peaks reference-ordered.
        target.peak_indices.insert(target.peak_indices.end(),
                                   piece.peak_indices.begin(),
                                   piece.peak_indices.end());
        target.reference_end =
            std::max(target.reference_end, piece.reference_end);
        target.oriented_query_begin =
            std::min(target.oriented_query_begin, piece.oriented_query_begin);
        target.oriented_query_end =
            std::max(target.oriented_query_end, piece.oriented_query_end);
        target.score += piece.score;
        target.seed_support += piece.seed_support;
        piece.peak_indices.clear(); // absorbed
        continue;
      }
      current = index;
    }
    staircases.erase(std::remove_if(staircases.begin(), staircases.end(),
                                    [](const RnaStaircase& s) {
                                      return s.peak_indices.empty();
                                    }),
                     staircases.end());
  }
  for (RnaStaircase& staircase : staircases) {
    const RnaQuerySpan forward =
        rna_forward_query_span(static_cast<int>(staircase.oriented_query_begin),
                               static_cast<int>(staircase.oriented_query_end),
                               staircase.reverse, read_len);
    if (forward.begin >= 0 && forward.end > forward.begin) {
      staircase.forward_query_begin = static_cast<std::uint32_t>(forward.begin);
      staircase.forward_query_end = static_cast<std::uint32_t>(forward.end);
    }
  }
  std::stable_sort(staircases.begin(), staircases.end(),
                   [](const RnaStaircase& a, const RnaStaircase& b) {
                     if (a.score != b.score)
                       return a.score > b.score;
                     if (a.reference_id != b.reference_id)
                       return a.reference_id < b.reference_id;
                     if (a.reference_begin != b.reference_begin)
                       return a.reference_begin < b.reference_begin;
                     if (a.reverse != b.reverse)
                       return a.reverse < b.reverse;
                     return a.peak_indices.front() < b.peak_indices.front();
                   });
}

// The tiles a chain's selected anchors touch, forward frame (ignored anchors are not
// part of the chain), plus the chain's extent on the query (own frame) and reference,
// and the anchor count, the solver's chain evidence.
struct ChainExtent {
  std::int64_t anchors = 0;
  int query_begin = -1;
  int query_end = -1;
  std::int64_t reference_begin = -1;
  std::int64_t reference_end = -1;
};

// The mask is the query the chain explains: consecutive anchors (forward order) closer
// than two seed lengths form one run whose tiles are all marked, so a chain does not
// look sparser than the peak-span units. A wider gap (an unpeaked exon, or a sprawled
// chain's bridge) stays unmarked.
ChainExtent anchor_tiles(const ExactAnchorPath& bundle, bool reverse,
                         int read_len, int seed_len,
                         voting::QueryTileMask& mask) {
  ChainExtent extent;
  std::vector<RnaQuerySpan> spans;
  spans.reserve(bundle.selected_raw_indices.size());
  for (const std::uint32_t raw : bundle.selected_raw_indices) {
    if (raw >= bundle.selected_pool_anchors.size())
      continue;
    const SelectedPoolAnchor& anchor = bundle.selected_pool_anchors[raw];
    if (anchor.ignore_flag)
      continue;
    ++extent.anchors;
    const int query_end = anchor.query_begin + anchor.span;
    const RnaQuerySpan forward = rna_forward_query_span(
        anchor.query_begin, query_end, reverse, read_len);
    if (span_measured(forward))
      spans.push_back(forward);
    if (extent.query_begin < 0 || anchor.query_begin < extent.query_begin)
      extent.query_begin = anchor.query_begin;
    extent.query_end = std::max(extent.query_end, query_end);
    const std::int64_t reference_begin = anchor.reference_begin;
    const std::int64_t reference_end = reference_begin + anchor.span;
    if (extent.reference_begin < 0 || reference_begin < extent.reference_begin)
      extent.reference_begin = reference_begin;
    extent.reference_end = std::max(extent.reference_end, reference_end);
  }
  std::stable_sort(spans.begin(), spans.end(),
                   [](RnaQuerySpan a, RnaQuerySpan b) {
                     if (a.begin != b.begin)
                       return a.begin < b.begin;
                     return a.end < b.end;
                   });
  const int fill = 2 * std::max(1, seed_len);
  RnaQuerySpan run;
  for (const RnaQuerySpan& span : spans) {
    if (!span_measured(run)) {
      run = span;
    } else if (span.begin - run.end <= fill) {
      run.end = std::max(run.end, span.end);
    } else {
      set_span_tiles(mask, run, read_len, seed_len);
      run = span;
    }
  }
  if (span_measured(run))
    set_span_tiles(mask, run, read_len, seed_len);
  return extent;
}

// The union of the peaks' forward spans in bases, the meaning rank_score has everywhere
// else. `order` is scratch.
std::int64_t
union_forward_bases(const std::vector<placement::CoarseDiagonalPeak>& peaks,
                    const std::vector<std::uint32_t>& indices, int read_len,
                    std::vector<std::uint32_t>& order) {
  order.assign(indices.begin(), indices.end());
  std::stable_sort(order.begin(), order.end(),
                   [&peaks, read_len](std::uint32_t a, std::uint32_t b) {
                     const RnaQuerySpan sa =
                         peak_forward_span(peaks[a], read_len);
                     const RnaQuerySpan sb =
                         peak_forward_span(peaks[b], read_len);
                     if (sa.begin != sb.begin)
                       return sa.begin < sb.begin;
                     return a < b;
                   });
  std::int64_t union_bases = 0;
  int run_begin = -1;
  int run_end = -1;
  for (const std::uint32_t peak_index : order) {
    const RnaQuerySpan span = peak_forward_span(peaks[peak_index], read_len);
    if (!span_measured(span))
      continue;
    if (run_begin < 0) {
      run_begin = span.begin;
      run_end = span.end;
    } else if (span.begin > run_end) {
      union_bases += run_end - run_begin;
      run_begin = span.begin;
      run_end = span.end;
    } else {
      run_end = std::max(run_end, span.end);
    }
  }
  if (run_begin >= 0)
    union_bases += run_end - run_begin;
  return union_bases;
}

// A locus from a set of peaks: query order (the order the skeleton harvest
// reads peak_indices in), the oriented and reference hulls, and the union of
// the forward spans in bases as rank_score. `indices` is sorted in place;
// `order` is scratch.
void locus_from_peaks(const std::vector<placement::CoarseDiagonalPeak>& peaks,
                      std::vector<std::uint32_t>& indices, int read_len,
                      std::vector<std::uint32_t>& order,
                      placement::CoarseLocus& locus) {
  std::stable_sort(indices.begin(), indices.end(),
                   [&peaks](std::uint32_t a, std::uint32_t b) {
                     const placement::CoarseDiagonalPeak& pa = peaks[a];
                     const placement::CoarseDiagonalPeak& pb = peaks[b];
                     if (pa.oriented_query_begin != pb.oriented_query_begin)
                       return pa.oriented_query_begin < pb.oriented_query_begin;
                     if (pa.reference_begin != pb.reference_begin)
                       return pa.reference_begin < pb.reference_begin;
                     return a < b;
                   });
  locus.peak_indices = indices;
  bool first = true;
  for (const std::uint32_t peak_index : locus.peak_indices) {
    const placement::CoarseDiagonalPeak& peak = peaks[peak_index];
    if (first) {
      locus.reference_id = peak.reference_id;
      locus.reverse = peak.reverse;
      locus.oriented_query_begin = peak.oriented_query_begin;
      locus.oriented_query_end = peak.oriented_query_end;
      locus.reference_begin = peak.reference_begin;
      locus.reference_end = peak.reference_end;
      first = false;
      continue;
    }
    locus.oriented_query_begin =
        std::min(locus.oriented_query_begin, peak.oriented_query_begin);
    locus.oriented_query_end =
        std::max(locus.oriented_query_end, peak.oriented_query_end);
    locus.reference_begin =
        std::min(locus.reference_begin, peak.reference_begin);
    locus.reference_end = std::max(locus.reference_end, peak.reference_end);
  }
  const std::int64_t union_bases =
      union_forward_bases(peaks, locus.peak_indices, read_len, order);
  locus.rank_score = union_bases;
  locus.union_query_coverage_bases = union_bases;
}

} // namespace

std::vector<placement::CoarseDiagonalPeak>
rna_revote_intervals(const LongReadSeedContext& seed_ctx,
                     const std::vector<QuerySeed>& forward_seeds,
                     const std::vector<KmerPostingView>& views,
                     const std::vector<RnaQuerySpan>& intervals, int read_len,
                     int k, int min_support, int n_chr,
                     placement::FusedCaptureScratch& capture,
                     RnaRevoteStats& stats) {
  std::vector<placement::CoarseDiagonalPeak> peaks;
  const auto inside = [&](int forward_position) {
    for (const RnaQuerySpan& interval : intervals)
      if (forward_position >= interval.begin &&
          forward_position + k <= interval.end)
        return true;
    return false;
  };
  const std::size_t seed_count = std::min(forward_seeds.size(), views.size());
  for (std::size_t i = 0; i < seed_count; ++i) {
    if (!inside(forward_seeds[i].read_pos))
      continue;
    const KmerPostingView& view = views[i];
    if (!view.found() || view.count == 0 ||
        !chain_window_seed_occ_allowed(seed_ctx, view))
      continue;
    stats.postings += static_cast<std::int64_t>(view.count);
    ++stats.seeds_revoted;
  }
  if (stats.postings > kRnaExplainMaxPostings) {
    stats.over_budget = true;
    stats.postings = 0;
    stats.seeds_revoted = 0;
    return peaks;
  }
  if (stats.seeds_revoted == 0)
    return peaks;
  capture.fwd_log.clear();
  capture.rc_log.clear();
  capture.fwd_bins.reset(static_cast<std::size_t>(stats.seeds_revoted));
  capture.rc_bins.reset(static_cast<std::size_t>(stats.seeds_revoted));
  const int W = effective_vote_diag_bin_width(seed_ctx, read_len);
  // Strand-compatible re-vote, the capture's rule (fused_capture.cpp): each posting is
  // read once and emits a forward and a reverse record, each only on the lane the
  // posting is compatible with.
  for (std::size_t i = 0; i < seed_count; ++i) {
    const int forward_position = forward_seeds[i].read_pos;
    if (!inside(forward_position))
      continue;
    const KmerPostingView& view = views[i];
    if (!view.found() || view.count == 0 ||
        !chain_window_seed_occ_allowed(seed_ctx, view))
      continue;
    const int reverse_position = read_len - k - forward_position;
    const int seed_index = static_cast<int>(i);
    // A seed's postings ascend, so its postings within one diagonal bin are
    // contiguous: a key change marks the first posting of this seed in that
    // bin (the distinct-seed support signal), as in the fused capture.
    std::int64_t previous_forward = INT64_MIN;
    std::int64_t previous_reverse = INT64_MIN;
    for (std::uint32_t posting = 0; posting < view.count; ++posting) {
      const PackedRefPos packed = view.positions.packed_at(posting);
      const int contig = static_cast<int>(packed_ref_contig(packed));
      if (contig < 0 || contig >= n_chr)
        continue;
      const std::uint32_t local = packed_ref_local(packed);
      if (local > static_cast<std::uint32_t>(INT_MAX))
        continue;
      // One posting, two candidate records; exactly one lane is compatible. The
      // per-lane "first posting of this seed in this bin" signals update only on their
      // own lane's emission, so a declined record cannot suppress the first compatible
      // support.
      const bool emit_forward = packed_ref_orientation_compatible(
          forward_seeds[i].z, packed, /*reverse_lane=*/false);
      const bool emit_reverse = packed_ref_orientation_compatible(
          forward_seeds[i].z, packed, /*reverse_lane=*/true);
      ++stats.postings_inspected;
      if (!emit_forward)
        ++stats.rejected_forward;
      if (!emit_reverse)
        ++stats.rejected_reverse;
      if (emit_forward) {
        const int forward_start = static_cast<int>(local) - forward_position;
        const std::int64_t forward_key =
            vote_pack_key(contig, vote_floor_div(forward_start, W));
        const std::uint32_t forward_bin =
            capture.fwd_bins.add(forward_key, forward_key != previous_forward);
        previous_forward = forward_key;
        capture.fwd_log.push_back(placement::CapturedVote{
            forward_key,
            VoteHit{seed_index, forward_position, forward_start,
                    view.occurrence, view.count},
            forward_bin});
      }
      if (emit_reverse) {
        const int reverse_start = static_cast<int>(local) - reverse_position;
        const std::int64_t reverse_key =
            vote_pack_key(contig, vote_floor_div(reverse_start, W));
        const std::uint32_t reverse_bin =
            capture.rc_bins.add(reverse_key, reverse_key != previous_reverse);
        previous_reverse = reverse_key;
        capture.rc_log.push_back(placement::CapturedVote{
            reverse_key,
            VoteHit{seed_index, reverse_position, reverse_start,
                    view.occurrence, view.count},
            reverse_bin});
      }
    }
  }
  stats.emitted_forward = static_cast<std::int64_t>(capture.fwd_log.size());
  stats.emitted_reverse = static_cast<std::int64_t>(capture.rc_log.size());
  peaks = placement::aggregate_fused_capture_to_peaks(capture, k, min_support);
  return peaks;
}

namespace {

void reset_result(RnaExplainResult& result) {
  result.intervals.clear();
  result.seeds_revoted = 0;
  result.postings = 0;
  result.over_budget = false;
  result.emitted_forward = 0;
  result.emitted_reverse = 0;
  result.postings_inspected = 0;
  result.rejected_forward = 0;
  result.rejected_reverse = 0;
  result.unit_peaks.clear();
  result.revote_peaks = 0;
  result.staircases.clear();
  result.candidates.clear();
  result.solved = false;
  result.candidates_offered = 0;
  result.chained_offered = 0;
  result.units_offered = 0;
  result.dp_cells = 0;
  result.windows.clear();
  result.windows_merged = 0;
  result.continuation_refused = 0;
  result.read_len = 0;
}

} // namespace

bool rna_segment_continuation(const RnaSegmentGeometry& a,
                              const RnaSegmentGeometry& b,
                              int max_intron) noexcept {
  if (a.reference_id < 0 || a.reference_id != b.reference_id ||
      a.reverse != b.reverse)
    return false;
  if (!span_measured(a.forward_span) || !span_measured(b.forward_span))
    return false;
  // Query order and reference order must agree at both ends: the forward
  // strand reads the reference forward, the reverse strand backward. A piece
  // that steps back on the reference is a second place, not the next exon.
  const bool a_begins_first = a.forward_span.begin <= b.forward_span.begin;
  const bool a_ends_first = a.forward_span.end <= b.forward_span.end;
  const bool a_ref_begins_first =
      a.reverse ? a.reference_end >= b.reference_end
                : a.reference_begin <= b.reference_begin;
  const bool a_ref_ends_first =
      a.reverse ? a.reference_begin >= b.reference_begin
                : a.reference_end <= b.reference_end;
  if (a_begins_first != a_ref_begins_first || a_ends_first != a_ref_ends_first)
    return false;
  // Within one intron: reference spans overlapping, or a gap of at most
  // max_intron between them.
  if (a.reference_begin < b.reference_end && b.reference_begin < a.reference_end)
    return true;
  const std::uint64_t intron =
      static_cast<std::uint64_t>(std::max(0, max_intron));
  const std::uint64_t gap = a.reference_end <= b.reference_begin
                                ? b.reference_begin - a.reference_end
                                : a.reference_begin - b.reference_end;
  return gap <= intron;
}

bool rna_same_place(const RnaSegmentGeometry& a, const RnaSegmentGeometry& b,
                    int max_intron) noexcept {
  if (a.reference_id < 0 || a.reference_id != b.reference_id ||
      a.reverse != b.reverse)
    return false;
  if (a.reference_begin < b.reference_end &&
      b.reference_begin < a.reference_end)
    return true;
  return rna_segment_continuation(a, b, max_intron);
}

// E0..E3.
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
    RnaExplainResult& result) {
  reset_result(result);
  const int k = seed_len;
  // A censored primary hull is not evidence that anything is left over; the
  // emission floor declines outright on one, so this does too.
  if (read_len <= 0 || k <= 0 || primary.reference_id < 0 ||
      !span_measured(primary.forward_span) ||
      primary.forward_span.end > read_len)
    return;
  result.read_len = read_len;
  result.max_intron = max_intron;

  // E0. The unexplained pieces of the query (no floor: the windows' clip
  // intervals and the harvest's query windows), then the re-vote's
  // intervals: at least the floor wide, longest first, at most the cap.
  unexplained_pieces(primary.forward_span, primary_segments, read_len, 0,
                     scratch.sorted, scratch.pieces);
  unexplained_pieces(primary.forward_span, primary_segments, read_len,
                     kRnaChimeraMinQueryBases, scratch.sorted,
                     result.intervals);
  std::stable_sort(result.intervals.begin(), result.intervals.end(),
                   [](RnaQuerySpan a, RnaQuerySpan b) {
                     const int wa = a.end - a.begin;
                     const int wb = b.end - b.begin;
                     if (wa != wb)
                       return wa > wb;
                     return a.begin < b.begin;
                   });
  if (result.intervals.size() >
      static_cast<std::size_t>(kRnaExplainMaxIntervals))
    result.intervals.resize(static_cast<std::size_t>(kRnaExplainMaxIntervals));

  // E1. The re-vote over those intervals.
  std::vector<placement::CoarseDiagonalPeak> revote_peaks;
  const int n_chr =
      seed_ctx.chr_names ? static_cast<int>(seed_ctx.chr_names->size()) : 0;
  // A re-vote peak needs kRnaCoprimaryMinAnchors distinct seeds: the re-vote sees every
  // seed of the interval, so two-seed peaks are mostly coincidences.
  {
    RnaRevoteStats revote_stats;
    if (!result.intervals.empty() && seed_ctx.index != nullptr && n_chr > 0 &&
        read_len >= k)
      revote_peaks = rna_revote_intervals(
          seed_ctx, forward_seeds, views, result.intervals, read_len, k,
          std::max(min_support, kRnaCoprimaryMinAnchors), n_chr,
          scratch.revote_capture, revote_stats);
    result.seeds_revoted = revote_stats.seeds_revoted;
    result.postings = revote_stats.postings;
    result.over_budget = revote_stats.over_budget;
    result.emitted_forward = revote_stats.emitted_forward;
    result.emitted_reverse = revote_stats.emitted_reverse;
    result.postings_inspected = revote_stats.postings_inspected;
    result.rejected_forward = revote_stats.rejected_forward;
    result.rejected_reverse = revote_stats.rejected_reverse;
  }
  // Over budget: no re-vote at all, the read's peaks stand alone. A re-vote
  // that ran supersedes every read peak inside its intervals (it saw every
  // seed there, the read's selection included).
  const bool revoted = !result.over_budget && result.seeds_revoted > 0;

  // E2. The unit peaks: every read peak query-disjoint from the realized hull
  // and outside the re-voted intervals, then the re-vote's. Envelope
  // membership is not consulted.
  for (const placement::CoarseDiagonalPeak& peak : peaks) {
    if (peak.reference_id < 0)
      continue;
    const RnaQuerySpan span = peak_forward_span(peak, read_len);
    if (!span_measured(span))
      continue;
    if (!rna_chimera_query_disjoint(span, primary.forward_span))
      continue;
    if (revoted) {
      bool inside = false;
      for (const RnaQuerySpan& interval : result.intervals)
        if (span.begin < interval.end && span.end > interval.begin)
          inside = true;
      if (inside)
        continue;
    }
    result.unit_peaks.push_back(peak);
  }
  const std::size_t first_revote = result.unit_peaks.size();
  result.revote_peaks = static_cast<int>(revote_peaks.size());
  result.unit_peaks.insert(result.unit_peaks.end(), revote_peaks.begin(),
                           revote_peaks.end());
  if (!result.unit_peaks.empty()) {
    scratch.staircase.order.clear();
    for (std::size_t index = 0; index < result.unit_peaks.size(); ++index)
      if (result.unit_peaks[index].reference_id >= 0)
        scratch.staircase.order.push_back(static_cast<std::uint32_t>(index));
    build_staircases(result.unit_peaks, options, read_len, scratch.staircase,
                     result.staircases);
  }

  // E3. The candidates: the primary's realized tiles, pinned; every chained
  // hypothesis on its anchors; every unit on its peaks. Every non-primary
  // mask has the primary's tiles subtracted.
  voting::QueryTileMask primary_tiles;
  for (const RnaQuerySpan& segment : primary_segments)
    set_span_tiles(primary_tiles, segment, read_len, k);
  reset_partition_problem(scratch.problem, voting::QueryPartitionParameters{},
                          read_len, k);
  const voting::QueryTileMask outside =
      tiles_outside(primary_tiles, scratch.problem.valid_tiles);
  result.candidates.reserve(1 + chained.size() + result.staircases.size());
  {
    RnaSegmentCandidate candidate;
    candidate.kind = RnaSegmentKind::Primary;
    candidate.geometry = primary;
    candidate.support = primary_tiles;
    candidate.vote_evidence = primary_rank_score;
    result.candidates.push_back(std::move(candidate));
  }
  for (const RnaExplainChained& entry : chained) {
    if (entry.locus == nullptr || entry.bundle == nullptr)
      continue;
    RnaSegmentCandidate candidate;
    candidate.kind = RnaSegmentKind::Chained;
    candidate.hypothesis_index = entry.hypothesis_index;
    candidate.locus = entry.locus;
    candidate.vote_evidence = entry.rank_score;
    voting::QueryTileMask mask;
    const ChainExtent extent =
        anchor_tiles(*entry.bundle, entry.reverse, read_len, k, mask);
    candidate.support = mask & outside;
    candidate.chain_evidence = extent.anchors;
    candidate.geometry.reference_id = entry.locus->reference_id;
    candidate.geometry.reverse = entry.reverse;
    if (extent.reference_begin >= 0 &&
        extent.reference_end > extent.reference_begin) {
      candidate.geometry.reference_begin =
          static_cast<std::uint64_t>(extent.reference_begin);
      candidate.geometry.reference_end =
          static_cast<std::uint64_t>(extent.reference_end);
    } else {
      candidate.geometry.reference_begin = entry.locus->reference_begin;
      candidate.geometry.reference_end = entry.locus->reference_end;
    }
    candidate.geometry.forward_span = rna_forward_query_span(
        entry.q_begin, entry.q_end, entry.reverse, read_len);
    result.candidates.push_back(std::move(candidate));
  }
  for (std::size_t s = 0; s < result.staircases.size(); ++s) {
    const RnaStaircase& staircase = result.staircases[s];
    RnaSegmentCandidate candidate;
    candidate.staircase_index = static_cast<int>(s);
    bool all_revote = true;
    for (const std::uint32_t peak_index : staircase.peak_indices)
      if (peak_index < first_revote)
        all_revote = false;
    candidate.kind = all_revote ? RnaSegmentKind::Revote : RnaSegmentKind::Unit;
    voting::QueryTileMask mask;
    staircase_tiles(result.unit_peaks, staircase.peak_indices, mask, read_len,
                    k);
    candidate.support = mask & outside;
    candidate.chain_evidence =
        static_cast<std::int64_t>(staircase.seed_support);
    candidate.vote_evidence =
        union_forward_bases(result.unit_peaks, staircase.peak_indices, read_len,
                            scratch.window_peaks);
    candidate.geometry.reference_id = staircase.reference_id;
    candidate.geometry.reverse = staircase.reverse;
    candidate.geometry.reference_begin = staircase.reference_begin;
    candidate.geometry.reference_end = staircase.reference_end;
    if (staircase.forward_query_end > staircase.forward_query_begin) {
      candidate.geometry.forward_span.begin =
          static_cast<int>(staircase.forward_query_begin);
      candidate.geometry.forward_span.end =
          static_cast<int>(staircase.forward_query_end);
    }
    result.candidates.push_back(std::move(candidate));
  }
}

// E3'.
RnaQuerySpan rna_explain_piece(const RnaExplainScratch& scratch,
                               RnaQuerySpan span) noexcept {
  RnaQuerySpan piece;
  int best_overlap = 0;
  for (const RnaQuerySpan& candidate : scratch.pieces) {
    const RnaQuerySpan clipped = clip_to(span, candidate);
    const int overlap = span_measured(clipped) ? clipped.end - clipped.begin : 0;
    if (overlap > best_overlap) {
      best_overlap = overlap;
      piece = candidate;
    }
  }
  return piece;
}

void rna_explain_candidate_chained(RnaExplainResult& result,
                                   std::size_t candidate_index,
                                   int hypothesis_index,
                                   const ExactAnchorPath& bundle, bool reverse,
                                   int read_len, int seed_len,
                                   const voting::QueryTileMask& primary_tiles) {
  if (candidate_index >= result.candidates.size())
    return;
  RnaSegmentCandidate& candidate = result.candidates[candidate_index];
  candidate.hypothesis_index = hypothesis_index;
  voting::QueryTileMask mask;
  const ChainExtent extent =
      anchor_tiles(bundle, reverse, read_len, seed_len, mask);
  candidate.support =
      voting::QueryTileMask{{mask.words[0] & ~primary_tiles.words[0],
                             mask.words[1] & ~primary_tiles.words[1]}};
  candidate.chain_evidence = extent.anchors;
  candidate.geometry.reverse = reverse;
  if (extent.query_begin >= 0 && extent.query_end > extent.query_begin) {
    const RnaQuerySpan forward = rna_forward_query_span(
        extent.query_begin, extent.query_end, reverse, read_len);
    if (span_measured(forward))
      candidate.geometry.forward_span = forward;
  }
  if (extent.reference_begin >= 0 &&
      extent.reference_end > extent.reference_begin) {
    candidate.geometry.reference_begin =
        static_cast<std::uint64_t>(extent.reference_begin);
    candidate.geometry.reference_end =
        static_cast<std::uint64_t>(extent.reference_end);
  }
}

int rna_explain_contender(const RnaExplainResult& result,
                          std::size_t candidate_index) {
  if (candidate_index >= result.candidates.size())
    return -1;
  const RnaSegmentCandidate& unit = result.candidates[candidate_index];
  if (unit.kind != RnaSegmentKind::Unit && unit.kind != RnaSegmentKind::Revote)
    return -1;
  const int unit_tiles = unit.support.count();
  if (unit_tiles < kRnaPartitionMinBlockTiles)
    return -1;
  for (std::size_t c = 1; c < result.candidates.size(); ++c) {
    if (c == candidate_index)
      continue;
    const RnaSegmentCandidate& other = result.candidates[c];
    const int other_tiles = other.support.count();
    if (other_tiles < kRnaPartitionMinBlockTiles)
      continue;
    if (rna_same_place(unit.geometry, other.geometry, result.max_intron))
      continue;
    const int shared = (unit.support & other.support).count();
    if (static_cast<double>(shared) >
        kRnaChainMapqMaskLevel *
            static_cast<double>(std::min(unit_tiles, other_tiles)))
      return static_cast<int>(c);
  }
  return -1;
}

placement::CoarseLocus rna_explain_unit_locus(const RnaExplainResult& result,
                                              int staircase_index) {
  placement::CoarseLocus locus;
  if (staircase_index < 0 ||
      static_cast<std::size_t>(staircase_index) >= result.staircases.size())
    return locus;
  std::vector<std::uint32_t> indices =
      result.staircases[static_cast<std::size_t>(staircase_index)].peak_indices;
  std::vector<std::uint32_t> order;
  locus_from_peaks(result.unit_peaks, indices, result.read_len, order, locus);
  return locus;
}

// E4 + E5.
void rna_explain_solve(const placement::CoarseLocusOptions& options,
                       const RnaSegmentGeometry& primary,
                       const std::vector<RnaQuerySpan>& primary_segments,
                       int read_len, int seed_len, int max_intron,
                       RnaExplainScratch& scratch, RnaExplainResult& result) {
  (void)options;
  result.solved = false;
  result.candidates_offered = 0;
  result.chained_offered = 0;
  result.units_offered = 0;
  result.dp_cells = 0;
  result.windows.clear();
  result.windows_merged = 0;
  result.continuation_refused = 0;
  if (result.candidates.empty() || read_len <= 0 || seed_len <= 0)
    return;
  std::vector<RnaSegmentCandidate>& candidates = result.candidates;
  for (RnaSegmentCandidate& candidate : candidates) {
    candidate.offered = false;
    candidate.owns_block = false;
    candidate.catalogue_slot = -1;
  }

  // 1. The offer: the primary in slot 0, then the chained hypotheses by anchors, then
  // the units by vote evidence. A candidate with fewer than kRnaPartitionMinBlockTiles
  // supported tiles can own no block and is not offered; one candidate per distinct
  // mask; at most the solver's bound.
  voting::QueryPartitionProblem& problem = scratch.problem;
  reset_partition_problem(problem, voting::QueryPartitionParameters{}, read_len,
                          seed_len);
  scratch.offered_index.clear();
  std::vector<std::size_t>& order = scratch.offer_order;
  order.clear();
  for (std::size_t index = 1; index < candidates.size(); ++index)
    order.push_back(index);
  const auto is_chained = [](const RnaSegmentCandidate& c) {
    return c.kind == RnaSegmentKind::Chained;
  };
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t a, std::size_t b) {
                     const RnaSegmentCandidate& ca = candidates[a];
                     const RnaSegmentCandidate& cb = candidates[b];
                     if (is_chained(ca) != is_chained(cb))
                       return is_chained(ca);
                     if (is_chained(ca)) {
                       if (ca.chain_evidence != cb.chain_evidence)
                         return ca.chain_evidence > cb.chain_evidence;
                       return a < b;
                     }
                     if (ca.vote_evidence != cb.vote_evidence)
                       return ca.vote_evidence > cb.vote_evidence;
                     if (ca.chain_evidence != cb.chain_evidence)
                       return ca.chain_evidence > cb.chain_evidence;
                     return a < b;
                   });
  const auto offer = [&](std::size_t index) {
    RnaSegmentCandidate& source = candidates[index];
    voting::QueryCandidate candidate;
    candidate.equivalence_key = static_cast<std::uint64_t>(index) + 1;
    candidate.lane = source.geometry.reverse ? 1 : 0;
    candidate.vote_evidence = clamp_evidence(source.vote_evidence);
    candidate.chain_evidence = clamp_evidence(source.chain_evidence);
    candidate.support = source.support;
    const std::size_t slot = problem.catalogue.candidates.size();
    candidate.id = static_cast<voting::CandidateId>(slot);
    candidate.catalogue_rank = static_cast<int>(slot);
    source.offered = true;
    source.catalogue_slot = static_cast<int>(slot);
    scratch.offered_index.push_back(index);
    problem.catalogue.candidates.push_back(candidate);
  };
  offer(0);
  for (const std::size_t index : order) {
    if (problem.catalogue.candidates.size() >=
        static_cast<std::size_t>(kRnaPartitionCatalogueBound))
      break;
    const RnaSegmentCandidate& candidate = candidates[index];
    if (candidate.support.count() < kRnaPartitionMinBlockTiles)
      continue;
    if (duplicate_mask(problem, candidate.support))
      continue;
    offer(index);
    if (is_chained(candidate))
      ++result.chained_offered;
    else
      ++result.units_offered;
  }
  // Nothing but the primary explains anything: no solve.
  if (problem.catalogue.candidates.size() < 2)
    return;
  result.candidates_offered =
      static_cast<int>(problem.catalogue.candidates.size());

  // 2. The one exact partition.
  const voting::QueryPartitionResult solved =
      voting::solve_query_partition(problem);
  result.solved = true;
  result.dp_cells = solved.dp_cells;

  // 3. Windows: one candidate's explanation of one unexplained piece, the hull of every
  // block the solve gave it inside the piece, clipped to the piece. The hull, because
  // the solve charges nothing for an unsupported tile and a small unit wedged between
  // two blocks of one chain must not cut the chain's window below the floor. A block
  // belongs to the piece it overlaps most, since the tile grid is coarser than the
  // realized boundary.
  std::vector<RnaExplainWindow>& windows = result.windows;
  for (const voting::QueryBlock& block : solved.selected.blocks) {
    const int slot = static_cast<int>(block.candidate);
    if (slot <= 0 ||
        slot >= static_cast<int>(problem.catalogue.candidates.size()))
      continue; // a null block, or the primary's own
    const std::size_t index =
        scratch.offered_index[static_cast<std::size_t>(slot)];
    if (index == 0)
      continue;
    candidates[index].owns_block = true;
    if (block.supporting_tiles < kRnaPartitionMinBlockTiles)
      continue;
    const RnaQuerySpan raw = block_forward_span(
        block.query_tile_begin, block.query_tile_end, read_len, seed_len);
    const RnaQuerySpan piece = rna_explain_piece(scratch, raw);
    if (!span_measured(piece))
      continue;
    std::size_t hull = windows.size();
    for (std::size_t w = 0; w < windows.size(); ++w)
      if (windows[w].candidate_index == static_cast<int>(index) &&
          same_span(windows[w].piece, piece)) {
        hull = w;
        break;
      }
    if (hull == windows.size()) {
      const RnaSegmentCandidate& owner = candidates[index];
      RnaExplainWindow window;
      window.kind = owner.kind;
      window.candidate_index = static_cast<int>(index);
      window.hypothesis_index = owner.hypothesis_index;
      window.piece = piece;
      window.tile_begin = block.query_tile_begin;
      window.tile_end = block.query_tile_end;
      window.owners.assign(1, static_cast<int>(index));
      windows.push_back(std::move(window));
    }
    RnaExplainWindow& window = windows[hull];
    window.tile_begin = std::min(window.tile_begin, block.query_tile_begin);
    window.tile_end = std::max(window.tile_end, block.query_tile_end);
    window.supporting_tiles += block.supporting_tiles;
  }
  const auto by_tile = [](const RnaExplainWindow& a,
                          const RnaExplainWindow& b) {
    if (a.tile_begin != b.tile_begin)
      return a.tile_begin < b.tile_begin;
    return a.tile_end < b.tile_end;
  };
  std::stable_sort(windows.begin(), windows.end(), by_tile);

  // Two hulls in one piece whose owners continue each other (rna_segment_continuation)
  // are one segment the units cut in two: the earlier absorbs the later.
  for (std::size_t w = 0; w + 1 < windows.size();) {
    RnaExplainWindow& a = windows[w];
    const RnaExplainWindow& b = windows[w + 1];
    if (same_span(a.piece, b.piece) &&
        rna_segment_continuation(
            candidates[static_cast<std::size_t>(a.candidate_index)].geometry,
            candidates[static_cast<std::size_t>(b.candidate_index)].geometry,
            max_intron)) {
      a.tile_begin = std::min(a.tile_begin, b.tile_begin);
      a.tile_end = std::max(a.tile_end, b.tile_end);
      a.supporting_tiles += b.supporting_tiles;
      a.owners.insert(a.owners.end(), b.owners.begin(), b.owners.end());
      windows.erase(windows.begin() + static_cast<std::ptrdiff_t>(w) + 1);
      ++result.windows_merged;
      continue;
    }
    ++w;
  }

  // Each hull's window: its span clipped to the piece, at least the floor wide; not a
  // continuation of the primary (an exon the realizer did not reach is not a second
  // place); its locus; its competing vote.
  std::vector<std::uint32_t>& window_peaks = scratch.window_peaks;
  std::size_t kept = 0;
  for (std::size_t w = 0; w < windows.size(); ++w) {
    RnaExplainWindow window = std::move(windows[w]);
    window.span = clip_to(block_forward_span(window.tile_begin, window.tile_end,
                                             read_len, seed_len),
                          window.piece);
    if (!span_measured(window.span) ||
        window.span.end - window.span.begin < kRnaChimeraMinQueryBases)
      continue;
    const RnaSegmentCandidate& owner =
        candidates[static_cast<std::size_t>(window.candidate_index)];
    RnaSegmentGeometry geometry;
    geometry.reference_id = owner.geometry.reference_id;
    geometry.reverse = owner.geometry.reverse;
    geometry.reference_begin = owner.geometry.reference_begin;
    geometry.reference_end = owner.geometry.reference_end;
    for (const int o : window.owners) {
      const RnaSegmentGeometry& g =
          candidates[static_cast<std::size_t>(o)].geometry;
      geometry.reference_begin =
          std::min(geometry.reference_begin, g.reference_begin);
      geometry.reference_end = std::max(geometry.reference_end, g.reference_end);
    }
    geometry.forward_span = window.span;
    if (rna_segment_continuation(primary, geometry, max_intron)) {
      ++result.continuation_refused;
      continue;
    }
    if (owner.kind == RnaSegmentKind::Chained) {
      if (owner.locus != nullptr)
        window.locus = *owner.locus;
      window.unit_peaks = false;
    } else {
      // A window of units holds every unit owner's peaks.
      window_peaks.clear();
      for (const int o : window.owners) {
        const RnaSegmentCandidate& unit =
            candidates[static_cast<std::size_t>(o)];
        if (unit.kind == RnaSegmentKind::Chained || unit.staircase_index < 0)
          continue;
        const std::vector<std::uint32_t>& indices =
            result.staircases[static_cast<std::size_t>(unit.staircase_index)]
                .peak_indices;
        window_peaks.insert(window_peaks.end(), indices.begin(), indices.end());
      }
      locus_from_peaks(result.unit_peaks, window_peaks, read_len,
                       scratch.staircase.order, window.locus);
      window.unit_peaks = true;
    }
    // The demoted vote for this window: the largest vote evidence among the other
    // candidates sharing more than half of the shorter tile mask with it. Tiles rather
    // than spans, because a chain's span sprawls across the breakpoint and its anchors
    // do not. The primary, the window's owners and candidates at the window's own place
    // are excluded.
    voting::QueryTileMask window_tiles;
    set_span_tiles(window_tiles, window.span, read_len, seed_len);
    const int window_tile_count = window_tiles.count();
    window.competing_vote = 0;
    for (std::size_t c = 1; c < candidates.size(); ++c) {
      if (std::find(window.owners.begin(), window.owners.end(),
                    static_cast<int>(c)) != window.owners.end())
        continue;
      const RnaSegmentCandidate& other = candidates[c];
      if (rna_same_place(geometry, other.geometry, max_intron))
        continue;
      const int other_tiles = other.support.count();
      if (other_tiles == 0 || window_tile_count == 0)
        continue;
      const int shared = (other.support & window_tiles).count();
      if (static_cast<double>(shared) >
          kRnaChainMapqMaskLevel *
              static_cast<double>(std::min(other_tiles, window_tile_count)))
        window.competing_vote =
            std::max(window.competing_vote, other.vote_evidence);
    }
    windows[kept++] = std::move(window);
  }
  windows.resize(kept);

  // 4. Query order; at most kRnaExplainMaxWindows by supporting tiles.
  const auto by_begin = [](const RnaExplainWindow& a,
                           const RnaExplainWindow& b) {
    if (a.span.begin != b.span.begin)
      return a.span.begin < b.span.begin;
    return a.span.end < b.span.end;
  };
  if (windows.size() > static_cast<std::size_t>(kRnaExplainMaxWindows)) {
    std::stable_sort(windows.begin(), windows.end(),
                     [&](const RnaExplainWindow& a, const RnaExplainWindow& b) {
                       if (a.supporting_tiles != b.supporting_tiles)
                         return a.supporting_tiles > b.supporting_tiles;
                       return by_begin(a, b);
                     });
    windows.resize(static_cast<std::size_t>(kRnaExplainMaxWindows));
  }
  std::stable_sort(windows.begin(), windows.end(), by_begin);
}

void restrict_harvest_to_query_window(Rank1HarvestResult& harvest,
                                      RnaQuerySpan forward_window,
                                      bool pool_reverse, int read_len,
                                      int seed_len) {
  if (harvest.refused)
    return;
  CandidatePool& pool = harvest.bounded_pool;
  if (pool.refused)
    return;
  if (read_len <= 0 || forward_window.begin < 0 ||
      forward_window.end <= forward_window.begin) {
    // An unmeasurable window must not fall back to the whole padded envelope, so the
    // harvest is refused.
    harvest.refused = true;
    harvest.refusal = AnchoringRefusal::EmptyOrInvalidInput;
    return;
  }
  const int slack = std::max(0, seed_len);
  const int low =
      forward_window.begin > slack ? forward_window.begin - slack : 0;
  const int high = forward_window.end + slack;
  bool narrowed = false;
  for (FineAnchor& anchor : pool.anchors) {
    // Only a live verdict is worth revisiting: an anchor the occurrence cap or
    // the harvest's own geometry already rejected carries its reason already.
    if (!anchor.survived_occurrence || !anchor.survived_geometry)
      continue;
    const RnaQuerySpan span = rna_forward_query_span(
        anchor.query_begin, anchor.query_end, pool_reverse, read_len);
    // A censored projection is not evidence that the anchor is inside the
    // window; censored geometry stays censored, which here means excluded.
    if (span.begin >= low && span.end <= high)
      continue;
    anchor.survived_geometry = false;
    anchor.drop_reason = FineAnchorDropReason::kGeometry;
    narrowed = true;
  }
  // A window that already held the whole pool leaves the pool, ledger and content hash
  // as harvested.
  if (narrowed)
    refinalize_candidate_pool(pool, harvest.regions, harvest.harvest_params,
                              harvest.work.postings_skipped_occurrence);
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
