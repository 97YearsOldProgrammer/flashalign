// The leftover step: the forward intervals no printed region covers, harvested again from all their syncmers into pools.
#include "rna.h"

#include "../index/format.h" // packed_ref_contig / packed_ref_local (the pools)
#include "../index/seed.h" // query_seed_has_tandem_neighbor (the leftover)
#include "../seed/seed.h" // LongReadSeedContext
#include "../vote/vote.h" // chain_window_seed_occ_allowed

#include <algorithm>
#include <climits>
#include <cstdint>
#include <tuple>

/**** query_partition ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

namespace {

bool span_measured(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

// The unexplained pieces of the query: the hull's terminal clips, then the gaps
// between `segments` in query order (a running reach, so nested or overlapping
// segments leave no false gap), each at least `min_bp` wide (0 keeps every non-empty
// piece). `sorted` is scratch.
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

} // namespace

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
                        int k, int n_chr, const std::uint8_t* readmitted,
                        std::vector<RnaLeftoverPool>& pools) {
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
           (chain_window_seed_occ_allowed(seed_ctx, view) ||
            (readmitted != nullptr &&
             readmitted[forward_seeds[i].read_pos] != 0));
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

namespace {

// Same contig and strand, query order agreeing with reference order at both ends
// (mirrored on the reverse strand), and reference spans overlapping or within
// max_intron. Overlap is allowed because an exon the realizer skipped lies inside the
// primary's span; a piece that steps back on the reference (a duplicated exon, a
// back-splice) is a second place. Query overlap is not consulted; a censored span is
// never a continuation.
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

} // namespace

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

/**** leftover ****/
// Stage 5's leftover pass over the query the committed family and its pieces leave
// unexplained.

namespace fa::cpu::lr::rna {

// The leftover, after the primary's MAPQ: the forward intervals no printed region
// covers, each contig and lane's pool chained and realized as a co-primary, and the
// after-election parent rule rerun with their regions behind the elected family's.
void run_leftover(const ReadFacts& facts, ReadState& rs, Lifecycle& lc,
                  const RnaSpliceRealizationResult& realized,
                  const RnaSpliceRealizationRequest& request) {
  auto& pieces = lc.pieces;
  const auto& read_len = facts.read_len;
  const auto& seed_ctx = facts.seed_ctx;
  const auto& fwd_seeds = facts.fwd_seeds;
  auto& capture = facts.capture;
  const auto& k = facts.k;
  const auto& rctx = facts.rctx;
  const auto& harvest_options = facts.harvest_options;
  const auto& read_key = facts.read_key;
  auto& leftover_rivals = lc.leftover_rivals;
  auto& hypotheses = lc.hypotheses;
  auto& loci = rs.loci;
  auto& lifecycle_elected = lc.lifecycle_elected;
  const auto& max_intron = facts.max_intron;
  auto& pieces_dp_input = lc.pieces_dp_input;
  auto& pieces_left_out = lc.pieces_left_out;
  auto& pieces_chain = lc.pieces_chain;
  std::vector<rna::RnaQuerySpan> intervals;
  if (const std::optional<rna::RnaQuerySpan> primary_hull =
          rna::rna_realized_query_span(realized)) {
    std::vector<rna::RnaQuerySpan> segments;
    if (const rna::RnaSpliceHypothesisResult* winner =
            rna::rna_winning_hypothesis(realized))
      for (const rna::SpliceRegionEvidence& evidence : winner->region_evidence)
        segments.push_back(
            rna::RnaQuerySpan{evidence.query_begin, evidence.query_end});
    std::vector<rna::RnaQuerySpan> explained;
    for (const PieceRegion& piece : pieces)
      explained.push_back(rna::RnaQuerySpan{piece.region.qs, piece.region.qe});
    rna::rna_leftover_intervals(*primary_hull, segments, explained, read_len,
                                rna::kRnaChimeraMinQueryBases, intervals);
  }
  const int n_chr = seed_ctx.chr_names
                        ? static_cast<int>(seed_ctx.chr_names->size())
                        : 0;
  std::vector<rna::RnaLeftoverPool> leftover_pools;
  if (!intervals.empty() && n_chr > 0 &&
      rna::rna_leftover_pools(seed_ctx, *fwd_seeds, capture.views, intervals,
                              read_len, k, n_chr, facts.readmitted,
                              leftover_pools)) {
    auto params =
        rna::rank1_anchor_path_params(rctx.ref.index->k(), harvest_options);
    params.min_cnt = rna::kRnaLeftoverChainMinCount;
    params.min_sc = rna::kRnaLeftoverChainMinScore;
    // Each pool's best chain and the chains that do not compete with it.
    struct LeftoverChain {
      rna::ExactAnchorPath bundle;
      rna::SiblingChainStats siblings;
      int reference_id;
      bool reverse;
    };
    std::vector<LeftoverChain> chains;
    std::vector<rna::DisjointAnchorPath> disjoint;
    for (const rna::RnaLeftoverPool& pool : leftover_pools) {
      rna::SiblingChainStats siblings;
      rna::ExactAnchorPathResult path = rna::select_exact_anchor_path(
          read_key, read_len, pool.reference_id, pool.reverse, pool.pool,
          params, &siblings, &disjoint);
      if (path.refused)
        continue;
      chains.push_back(LeftoverChain{std::move(path.bundle), siblings,
                                     pool.reference_id, pool.reverse});
      for (rna::DisjointAnchorPath& chain : disjoint)
        chains.push_back(LeftoverChain{std::move(chain.bundle), chain.siblings,
                                       pool.reference_id, pool.reverse});
    }
    // No reallocation below: the hypotheses borrow the stored chains.
    leftover_rivals.reserve(chains.size());
    const std::size_t first_leftover = hypotheses.size();
    const int first_index = static_cast<int>(loci.size());
    rna::RnaSegmentGeometry elected_geometry;
    const bool elected_placed =
        shadow_geometry(facts, hypotheses[lifecycle_elected], elected_geometry);
    for (LeftoverChain& chain : chains) {
      ProductionRival entry;
      if (!rna::summarize_exact_anchor_path(chain.bundle, entry.summary, nullptr))
        continue;
      entry.rival.candidate =
          first_index + static_cast<int>(leftover_rivals.size());
      entry.rival.reverse = chain.reverse;
      entry.rival.chained = true;
      entry.rival.chain_score = chain.bundle.anchor_path_score;
      entry.rival.chain_anchors =
          static_cast<int>(chain.bundle.anchor_path_anchor_count);
      entry.rival.q_begin = entry.summary.first_query_begin;
      entry.rival.q_end = entry.summary.last_query_end;
      entry.siblings = chain.siblings;
      entry.reference_id = chain.reference_id;
      entry.chained_path = std::move(chain.bundle);
      entry.owns_path = true;
      leftover_rivals.push_back(std::move(entry));
      const ProductionRival& stored = leftover_rivals.back();
      rna::RnaRealizedHypothesis h;
      h.catalogue_index = stored.rival.candidate;
      h.reference_id = stored.reference_id;
      h.bundle = &stored.chained_path;
      h.summary = &stored.summary;
      h.siblings = stored.siblings;
      h.chained = true;
      h.chain_score = stored.rival.chain_score;
      h.chain_anchors = stored.rival.chain_anchors;
      h.q_begin = stored.rival.q_begin;
      h.q_end = stored.rival.q_end;
      h.reverse = stored.rival.reverse;
      h.klass = rna::RnaRivalClass::CoPrimary;
      h.coverage_ok =
          rna::anchor_union_query_coverage(*h.bundle, h.query_covered_bases);
      // At the elected family's place it is a shadow, left out as the rule leaves
      // them.
      rna::RnaSegmentGeometry geometry;
      h.shadow = elected_placed && shadow_geometry(facts, h, geometry) &&
                 rna::rna_same_place(elected_geometry, geometry, max_intron);
      hypotheses.push_back(std::move(h));
      rna::RnaRealizedHypothesis& candidate = hypotheses.back();
      rna::RnaSpliceRealizationRequest leftover_request = request;
      if (!candidate.coverage_ok ||
          !bind_realization_reference(facts, leftover_request, candidate.reference_id))
        continue;
      candidate.result =
          rna::realize_exact_anchor_path(leftover_request, *candidate.bundle);
      candidate.realized = true;
      price_hypothesis(candidate, candidate.result);
    }
    // Their regions that pass minimap2's region filter, sorted with the regions
    // behind the elected family's, which stay first and unchanged.
    std::vector<rna::RnaParentHit> added;
    for (std::size_t i = first_leftover; i < hypotheses.size(); ++i)
      if (hypotheses[i].realized)
        append_dp_hits(lc, i, hypotheses[i].result, added);
    added.erase(
        std::remove_if(
            added.begin(), added.end(),
            [&](const rna::RnaParentHit& hit) {
              rna::SpliceRegionEvidence evidence =
                  rna::rna_winning_hypothesis(
                      hypotheses[static_cast<std::size_t>(hit.owner)].result)
                      ->region_evidence[static_cast<std::size_t>(hit.segment)];
              evidence.dp_maximum = hit.dp_max;
              return rna::classify_splice_region(evidence, read_len,
                                                 rctx.opts.splice_controller) !=
                     rna::SpliceRegionReject::None;
            }),
        added.end());
    if (!added.empty()) {
      const int elected_owner = static_cast<int>(lifecycle_elected);
      std::vector<rna::RnaParentHit> regions = pieces_dp_input;
      const auto family_end = std::find_if(
          regions.begin(), regions.end(), [&](const rna::RnaParentHit& hit) {
            return hit.owner != elected_owner;
          });
      std::vector<rna::RnaParentHit> rest(family_end, regions.end());
      regions.erase(family_end, regions.end());
      rest.insert(rest.end(), added.begin(), added.end());
      rna::rna_parent_hit_sort(rest);
      regions.insert(regions.end(), rest.begin(), rest.end());
      std::vector<char> left_out = pieces_left_out;
      left_out.resize(hypotheses.size(), 0);
      for (std::size_t i = first_leftover; i < hypotheses.size(); ++i)
        left_out[i] = hypotheses[i].shadow ? 1 : 0;
      decide_pieces(facts, lc, parent_regions(lc, regions, &left_out), pieces_chain);
    }
  }
}

} // namespace fa::cpu::lr::rna
