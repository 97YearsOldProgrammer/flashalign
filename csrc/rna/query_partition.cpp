#include "query_partition.h"

#include "../index/format.h" // packed_ref_contig / packed_ref_local (the re-vote)
#include "../index/seed.h" // query_seed_has_tandem_neighbor (the leftover)
#include "../seeding/context.h" // effective_vote_diag_bin_width
#include "../voting/state.h" // vote_pack_key / vote_floor_div / occ policy

#include <algorithm>
#include <climits>
#include <cstdint>
#include <tuple>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

namespace {

bool span_measured(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

} // namespace

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

void rna_leftover_intervals(RnaQuerySpan primary_hull,
                            const std::vector<RnaQuerySpan>& primary_segments,
                            const std::vector<RnaQuerySpan>& explained,
                            int read_len, int min_bp,
                            std::vector<RnaQuerySpan>& out) {
  RnaQuerySpan hull = primary_hull;
  std::vector<RnaQuerySpan> segments = primary_segments;
  for (const RnaQuerySpan& span : explained) {
    if (!span_measured(span))
      continue;
    segments.push_back(span);
    hull.begin = std::min(hull.begin, span.begin);
    hull.end = std::max(hull.end, span.end);
  }
  std::vector<RnaQuerySpan> sorted;
  unexplained_pieces(hull, segments, read_len, min_bp, sorted, out);
}

bool rna_leftover_pools(const LongReadSeedContext& seed_ctx,
                        const std::vector<QuerySeed>& forward_seeds,
                        const std::vector<KmerPostingView>& views,
                        const std::vector<RnaQuerySpan>& intervals, int read_len,
                        int k, int n_chr, std::vector<RnaLeftoverPool>& pools) {
  pools.clear();
  const auto inside = [&](int forward_position) {
    for (const RnaQuerySpan& interval : intervals)
      if (forward_position >= interval.begin &&
          forward_position + k <= interval.end)
        return true;
    return false;
  };
  const auto admitted = [&](std::size_t i) {
    const KmerPostingView& view = views[i];
    return inside(forward_seeds[i].read_pos) && view.found() && view.count != 0 &&
           chain_window_seed_occ_allowed(seed_ctx, view);
  };
  const std::size_t seed_count = std::min(forward_seeds.size(), views.size());
  std::int64_t postings = 0;
  for (std::size_t i = 0; i < seed_count; ++i)
    if (admitted(i))
      postings += static_cast<std::int64_t>(views[i].count);
  if (postings > kRnaExplainMaxPostings)
    return false;
  struct Hit {
    int contig;
    bool reverse;
    int reference;
    int query;
    std::uint32_t seed;
    std::uint32_t occurrence;
    bool tandem;
  };
  std::vector<Hit> hits;
  const std::uint32_t last_seed =
      static_cast<std::uint32_t>(forward_seeds.size()) - 1u;
  for (std::size_t i = 0; i < seed_count; ++i) {
    if (!admitted(i))
      continue;
    const KmerPostingView& view = views[i];
    const int forward_position = forward_seeds[i].read_pos;
    const bool tandem = query_seed_has_tandem_neighbor(forward_seeds, i);
    for (std::uint32_t posting = 0; posting < view.count; ++posting) {
      const PackedRefPos packed = view.positions.packed_at(posting);
      const int contig = static_cast<int>(packed_ref_contig(packed));
      if (contig < 0 || contig >= n_chr)
        continue;
      const std::uint32_t local = packed_ref_local(packed);
      if (local > static_cast<std::uint32_t>(INT_MAX - k))
        continue;
      for (const bool reverse : {false, true}) {
        if (!packed_ref_orientation_compatible(forward_seeds[i].z, packed, reverse))
          continue;
        const std::uint32_t seed = static_cast<std::uint32_t>(i);
        hits.push_back(Hit{contig, reverse, static_cast<int>(local),
                           reverse ? read_len - k - forward_position
                                   : forward_position,
                           reverse ? last_seed - seed : seed, view.occurrence,
                           tandem});
      }
    }
  }
  const auto key = [](const Hit& h) {
    return std::make_tuple(h.contig, h.reverse, h.reference, h.query);
  };
  std::sort(hits.begin(), hits.end(),
            [&](const Hit& a, const Hit& b) { return key(a) < key(b); });
  hits.erase(std::unique(hits.begin(), hits.end(),
                         [&](const Hit& a, const Hit& b) { return key(a) == key(b); }),
             hits.end());
  for (std::size_t begin = 0; begin < hits.size();) {
    std::size_t end = begin;
    while (end < hits.size() && hits[end].contig == hits[begin].contig &&
           hits[end].reverse == hits[begin].reverse)
      ++end;
    pools.emplace_back();
    RnaLeftoverPool& lane = pools.back();
    lane.reference_id = hits[begin].contig;
    lane.reverse = hits[begin].reverse;
    lane.pool.anchors.reserve(end - begin);
    for (std::size_t j = begin; j < end; ++j) {
      const Hit& hit = hits[j];
      FineAnchor anchor;
      anchor.anchor_id = static_cast<std::uint32_t>(j - begin);
      anchor.query_begin = hit.query;
      anchor.query_end = hit.query + k;
      anchor.reference_begin = hit.reference;
      anchor.reference_end = hit.reference + k;
      anchor.span = k;
      anchor.source_seed_id = hit.seed;
      anchor.source_occurrence = hit.occurrence;
      anchor.tandem_flag = hit.tandem;
      anchor.raw_multiplicity = 1;
      anchor.raw_seen = anchor.survived_occurrence = anchor.survived_geometry =
          anchor.survived_feature_cap = anchor.survived_region_cap =
              anchor.survived_global_cap = anchor.final_pool = true;
      lane.pool.anchors.push_back(std::move(anchor));
    }
    lane.pool.ledger.raw_seen = lane.pool.ledger.final_pool = end - begin;
    begin = end;
  }
  return true;
}

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

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
