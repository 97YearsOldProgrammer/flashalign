// The RNA conductor, RnaBackend::map_read: vote capture, the coarse catalogue, harvest and frame election, the window chains, the attempt, the realize loop and the committed record; the second pass is in vote.cpp, the election, the rival lifecycle, the MAPQ evidence, the pieces and the runner-ups in hit.cpp, the leftover in leftover.cpp.
#include "rna.h"

#include "../core/core.h" // reverse_complement_encoded_u8, signed RNA coordinate/count domain
#include "../index/seed.h" // ClosedSyncmerConfig / extract_closed_syncmer_query_seeds_into
#include "../seed/seed.h" // VoteHit / SeedDrop, tie_read_seed (the exact-tie break)

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring> // strcmp (the second pass's own typed refusal)
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

/**** config ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

StrandMode parse_strand_mode(const char* text, StrandMode fallback) {
  if (!text || !*text) return fallback;
  if (std::strcmp(text, "unknown") == 0 || std::strcmp(text, "auto") == 0)
    return StrandMode::Unknown;
  if (std::strcmp(text, "forward") == 0 ||
      std::strcmp(text, "fwd") == 0 || std::strcmp(text, "f") == 0 ||
      std::strcmp(text, "+") == 0)
    return StrandMode::Forward;
  if (std::strcmp(text, "reverse") == 0 ||
      std::strcmp(text, "rev") == 0 || std::strcmp(text, "r") == 0 ||
      std::strcmp(text, "-") == 0)
    return StrandMode::Reverse;
  // minimap2 -u n: score no splice motif at all.
  if (std::strcmp(text, "none") == 0 || std::strcmp(text, "n") == 0)
    return StrandMode::None;
  return fallback;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** output ****/

namespace fa { namespace cpu { namespace lr { namespace rna {

namespace {

void stamp_chain_tags(AlignResult& record, const RnaChainTags& tags) {
    record.chain_anchors = tags.chain_anchors;
    record.chain_score = tags.chain_score;
    record.secondary_chain_score = tags.secondary_chain_score;
}

}  // namespace

// `chain_tags` are the committed family's, stamped on the same records as `chain_mapq`:
// the head and its continuation segments.
static bool commit_realized_alignment(
    const RnaSpliceRealizationResult& realized,
    int read_len,
    const int* chain_mapq,
    const RnaChainTags& chain_tags,
    const std::vector<RnaChimericFamily>& chimeric,
    Result& out) {
    if (realized.refused || realized.segments.empty()) return false;

    AlignResult best = realized.segments.front();
    for (size_t i = 1; i < realized.segments.size(); ++i)
        best.supplementary.push_back(realized.segments[i]);

    static_cast<AlignResult&>(out) = std::move(best);
    out.read_len = read_len;
    if (chain_mapq) {
    // Chain-formula MAPQ, computed upstream in the committed frame.
        out.mapq = *chain_mapq;
    } else {
        // Only the MAPQ input accounting failed: emit a mapped Q0 row, as the formula
        // does for an undecidable placement. The caller counts this.
        out.mapq = 0;
    }
    // Segments of one spliced family share the read-level MAPQ and the chain tags. A
    // chimeric family explains different query bases on its own chain against
    // different rivals, so it carries its own.
    stamp_chain_tags(out, chain_tags);
    for (auto& supplementary : out.supplementary) {
        supplementary.mapq = out.mapq;
        stamp_chain_tags(supplementary, chain_tags);
    }
    // Chimeric families are appended flat (io/io.h requires every
    // supplementary to be a leaf), in the caller's order, each segment stamped with
    // its family's MAPQ.
    for (const RnaChimericFamily& family : chimeric) {
        if (family.result == nullptr) continue;
        if (family.result->refused || family.result->segments.empty()) continue;
        for (const AlignResult& segment : family.result->segments) {
            out.supplementary.push_back(segment);
            out.supplementary.back().mapq = family.mapq;
            stamp_chain_tags(out.supplementary.back(), family.chain_tags);
        }
    }
    if (realized.transcript_strand == 1)
        out.transcript_strand = '+';
    else if (realized.transcript_strand == 2)
        out.transcript_strand = '-';
    return true;
}

void attach_runner_up_family(RnaSpliceRealizationResult& runner_up,
                             const RnaChainTags& chain_tags,
                             AlignResult& out) {
    if (runner_up.refused || runner_up.segments.empty()) return;
    // Flat, as above: segment 0 is the head and its continuations are leaf
    // supplementaries. MAPQ 0 on every record, as DNA does for a retained alternative:
    // the confidence belongs to the committed record. The chain tags are the family's
    // own, as for a chimeric family.
    AlignResult head = std::move(runner_up.segments.front());
    head.mapq = 0;
    stamp_chain_tags(head, chain_tags);
    head.supplementary.reserve(runner_up.segments.size() - 1);
    for (size_t i = 1; i < runner_up.segments.size(); ++i) {
        head.supplementary.push_back(std::move(runner_up.segments[i]));
        head.supplementary.back().mapq = 0;
        stamp_chain_tags(head.supplementary.back(), chain_tags);
    }
    out.secondary.push_back(std::move(head));
}

}}}}  // namespace fa::cpu::lr::rna

/**** placement/output ****/

namespace fa { namespace cpu { namespace lr { namespace rna { namespace placement {

// Fills `out` from the selected locus and its exact anchor path. Returns false when the
// path geometry or the MAPQ is out of range; there is no coarse fallback projection.
static bool project_fine_path_placement(
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

/**** query_span ****/

namespace fa::cpu::lr::rna {

namespace query_span_detail {

// The invalid span (begin < 0 or empty) marks censored geometry.
static bool forward_span_valid(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

// Shared overlap arithmetic; negative when the two spans are disjoint.
static int forward_span_overlap(RnaQuerySpan winner,
                                RnaQuerySpan rival) noexcept {
  return std::min(winner.end, rival.end) - std::max(winner.begin, rival.begin);
}

} // namespace query_span_detail

RnaQuerySpan rna_forward_query_span(int oriented_begin, int oriented_end,
                                    bool reverse, int read_len) noexcept {
  RnaQuerySpan span;
  if (oriented_begin < 0 || oriented_end <= oriented_begin)
    return span;
  if (!reverse) {
    span.begin = oriented_begin;
    span.end = oriented_end;
    return span;
  }
  if (read_len <= 0 || oriented_end > read_len)
    return span;
  span.begin = read_len - oriented_end;
  span.end = read_len - oriented_begin;
  return span;
}

bool rna_query_spans_compete(RnaQuerySpan winner,
                             RnaQuerySpan rival) noexcept {
  using query_span_detail::forward_span_overlap;
  using query_span_detail::forward_span_valid;
  if (!forward_span_valid(winner) || !forward_span_valid(rival))
    return true;
  const int overlap = forward_span_overlap(winner, rival);
  if (overlap <= 0)
    return false;
  const int shorter =
      std::min(winner.end - winner.begin, rival.end - rival.begin);
  return static_cast<double>(overlap) >
         kRnaChainMapqMaskLevel * static_cast<double>(shorter);
}

int rna_uncovered_query_bases(RnaQuerySpan winner,
                              RnaQuerySpan rival) noexcept {
  using query_span_detail::forward_span_overlap;
  using query_span_detail::forward_span_valid;
  if (!forward_span_valid(winner) || !forward_span_valid(rival))
    return 0;
  const int overlap = std::max(0, forward_span_overlap(winner, rival));
  return (rival.end - rival.begin) - overlap;
}

} // namespace fa::cpu::lr::rna

/**** selected_locus_anchoring ****/

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

static placement::CoarseLocus
mirror_locus_query_orientation(const placement::CoarseLocus& locus,
                               std::uint32_t query_length) {
  placement::CoarseLocus mirror = locus;
  mirror.reverse = !locus.reverse;
  mirror.oriented_query_begin = locus.oriented_query_end <= query_length
                                    ? query_length - locus.oriented_query_end
                                    : 0u;
  mirror.oriented_query_end = locus.oriented_query_begin <= query_length
                                  ? query_length - locus.oriented_query_begin
                                  : query_length;
  return mirror;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa

/**** read ****/
// The read's exits and its reverse strand: the single exit and the lazily encoded
// reverse stream.

namespace fa::cpu::lr::rna {

static AlignResult finish(const ReadState& rs,
                          rna::Result result) {
  return static_cast<AlignResult&&>(std::move(result));
}


static const std::vector<uint8_t>& ensure_rc_enc(const ReadFacts& facts) {
  auto& rc_enc = facts.rc_enc;
  const auto& fwd_enc = facts.fwd_enc;
  if (rc_enc.size() != fwd_enc.size())
    rc_enc = ::fa::cpu::lr::reverse_complement_encoded_u8(fwd_enc);
  return rc_enc;
}

} // namespace fa::cpu::lr::rna

/**** catalogue ****/
// The catalogue's rank scores and Stages 1-3: the vote capture and the coarse
// catalogue.

namespace fa::cpu::lr::rna {

static void read_catalogue_rank_scores(ReadState& rs) {
  auto& loci = rs.loci;
  auto& catalogue_rank1_score = rs.catalogue_rank1_score;
  auto& catalogue_rank2_score = rs.catalogue_rank2_score;
  catalogue_rank1_score = loci[0].rank_score;
  catalogue_rank2_score = loci.size() >= 2 ? loci[1].rank_score : 0;
}


// The fine chain is placement evidence in every output mode. A read whose harvest or
// anchor path cannot be formed is left unmapped: a coarse locus is a vote envelope,
// not placement evidence.
static AlignResult finish_fine_failure(const ReadFacts& facts, ReadState& rs,
                                       const AttemptState& as) {
  auto& out = rs.out;
  return finish(rs, std::move(out));
}

// Stage 1: the fused two-strand vote capture over the forward seeds and their reverse
// projection. Sets the projected seeds, the read's posting memo and the capture.
static void capture_votes(ReadFacts& facts, ReadState& rs) {
  const auto& shared_fwd_syncmer_seeds = facts.shared_fwd_syncmer_seeds;
  const auto& read_len = facts.read_len;
  const auto& cfg = facts.cfg;
  const auto& seed_ctx = facts.seed_ctx;
  const auto& rctx = facts.rctx;
  auto& projected_shared_reverse = facts.projected_shared_reverse;
  project_query_seeds_to_rc_coordinates(*shared_fwd_syncmer_seeds, read_len,
                                        cfg.k, projected_shared_reverse);
  // Per-read posting memo for both strands: a read runs several harvests over the same
  // two seed vectors (committed locus, its mirror, rivals), each probing every seed key.
  // The reverse projection keeps canonical keys, so reverse seed j is forward seed
  // n-1-j and one memo serves both. new_read() scopes it to this read in O(1).
  auto& postings = facts.postings;
  postings.new_read();
  postings.bind_forward(*shared_fwd_syncmer_seeds);
  postings.bind_reverse(projected_shared_reverse);
  auto& capture = facts.capture;
  // A context of 0 (every seed) stays 0.
  int vote_seeds_per_strand = seed_ctx.max_query_seeds_per_strand;
  const int vote_spacing = rctx.opts.vote_seed_spacing;
  if (vote_spacing > 0 && vote_seeds_per_strand > 0)
    vote_seeds_per_strand =
        std::max(vote_seeds_per_strand,
                 (read_len + vote_spacing - 1) / vote_spacing);
  // minimap2's seed selection: the capture marks the seeds each over-cap streak
  // re-admits, which every later occurrence gate of the read admits too.
  std::vector<std::uint8_t>* const readmitted =
      rctx.opts.occ_dist > 0 ? &facts.worker_scratch.readmitted : nullptr;
  rna::placement::capture_votes_fused(seed_ctx, *shared_fwd_syncmer_seeds,
                                      projected_shared_reverse, read_len,
                                      vote_seeds_per_strand, capture,
                                      rctx.opts.occ_dist, readmitted);
  facts.readmitted = readmitted != nullptr ? readmitted->data() : nullptr;
}

// Stages 2-3: the captured postings aggregated into exon-shaped peaks and the peaks
// chained into the coarse locus catalogue. Sets the coarse locus options.
static void build_coarse_catalogue(ReadFacts& facts, ReadState& rs) {
  auto& capture = facts.capture;
  const auto& k = facts.k;
  const auto& min_support = facts.min_support;
  const auto& read_len = facts.read_len;
  const auto& min_intron = facts.min_intron;
  const auto& max_intron = facts.max_intron;
  const auto& cfg = facts.cfg;
  const auto& max_loci = facts.max_loci;
  auto& worker_scratch = facts.worker_scratch;
  const auto& tie_seed = facts.tie_seed;
  rs.peaks =
      rna::placement::aggregate_fused_capture_to_peaks(capture, k, min_support);
  auto& peaks = rs.peaks;
  facts.copt =
      rna::placement::make_coarse_locus_options(
          static_cast<uint32_t>(k), static_cast<uint32_t>(min_intron),
          static_cast<uint32_t>(max_intron));
  auto& copt = facts.copt;
  copt.max_predecessors =
      std::max(1u, static_cast<uint32_t>(cfg.max_chain_predecessors));
  copt.max_locus_chains = static_cast<uint32_t>(max_loci);
  // No overflow: nothing downstream reads envelopes past the cut.
  copt.max_overflow_loci = 0;
  // Coalescing needs one piece that is real evidence; the floor scales with the vote
  // threshold.
  copt.coalesce_min_support =
      std::max<uint32_t>(2u, 4u * static_cast<uint32_t>(min_support));
  auto& loci = rs.loci;
  if (!peaks.empty()) {
    // Stage 3: coarse locus envelopes. Locus selection only; the exon path is the
    // spliced anchor chain's job (stage 4).
    // Past the cut, a locus whose seeds cover kRnaChimeraMinQueryBases of query no
    // kept locus's seeds cover is kept too.
    const rna::placement::CoarseLocusQueryCover query_cover{
        &capture.peak_read_pos_offset, &capture.peak_read_pos, read_len, k,
        rna::kRnaChimeraMinQueryBases};
    loci = rna::placement::select_coarse_loci(
        peaks, copt, worker_scratch.coarse, tie_seed, /*coalesced=*/nullptr,
        /*overflow_out=*/nullptr, &query_cover);
  }
}

} // namespace fa::cpu::lr::rna

/**** harvest ****/
// A window's harvest: the mirrored frame's pool and the pool each frame chains.

namespace fa::cpu::lr::rna {

// The mirrored frame's pool, from `harvest` in place; the mirror is always the
// projection, never a second harvest. Its preconditions, equal stream sizes (the
// n-1-j premise) and a single search region, are structural: a window violating one is
// refused, and `harvest`'s nominated pool is left for the escalation to chain.
static void project_mirror_frame(const ReadFacts& facts, ReadState& rs,
                                 rna::Rank1HarvestResult& harvest,
                                 std::size_t harvested_seeds,
                                 std::size_t mirror_seeds) {
  if (mirror_seeds != harvested_seeds || harvest.regions.size() != 1) {
    harvest.refused = true;
    harvest.refusal = rna::AnchoringRefusal::MirrorProjectionDomain;
    return;
  }
  // The opposite lane's pool was collected at harvest time in its own coordinates; a
  // refusal there is carried up.
  if (harvest.opposite_pool.refused) {
    harvest.refused = true;
    harvest.refusal = harvest.opposite_pool.refusal;
  }
}

// Which pool each frame chains.
static const rna::CandidatePool&
pool_for_nominated_chain(const ReadFacts& facts, ReadState& rs,
                         rna::Rank1HarvestResult& harvest) {
  return harvest.nominated_pool;
}

static const rna::CandidatePool&
pool_for_mirror_chain(const ReadFacts& facts, ReadState& rs,
                      rna::Rank1HarvestResult& harvest) {
  return harvest.opposite_pool;
}

} // namespace fa::cpu::lr::rna

/**** frame ****/
// The query frame of a window: the sibling merge of two chained frames and the
// pre-chaining election.

namespace fa::cpu::lr::rna {

// Each fine-stage window may be chained in both query frames, keeping the better
// chain, as minimap2 does: seed keys are canonical, so a window can be nominated on the
// wrong strand, where the true anchors are anti-diagonal but a small spurious chain may
// still pass min_cnt/min_sc.
// With both frames chained, the uncommitted frame's best chain is a runner-up inside
// the same window: chain counts add, and the sibling is the better of the committed
// frame's own runner-up and the other frame's best. A frame whose recurrence never ran
// (-1) contributes nothing.
static void merge_frame_siblings(rna::SiblingChainStats& committed,
                                 const rna::SiblingChainStats& other,
                                 const rna::ExactAnchorPathResult& other_path) {
  if (other.chains_total >= 0)
    committed.chains_total =
        std::max<std::int32_t>(committed.chains_total, 0) +
        other.chains_total;
  if (other_path.refused)
    return;
  const std::int32_t other_best = other_path.bundle.anchor_path_score;
  if (other_best > committed.sib_score) {
    committed.sib_score = other_best;
    committed.sib_count = static_cast<std::int32_t>(
        other_path.bundle.anchor_path_anchor_count);
  }
}

// Pre-chaining query-frame election: where the pool's band margin is decisive only the
// elected frame is chained; an elected frame that refuses escalates to the other, and
// kAmbiguous chains both and lets the chain scores decide.
static rna::FrameElection elect_frame(const ReadFacts& facts, ReadState& rs,
                                      const rna::Rank1HarvestResult& harvest) {
  // Routed: the nominated band comes from the harvested lane's pool and the mirror band
  // from the other lane's (rna.h). A refused lane contributes zeros.
  const bool routed_election = harvest.routed;
  const rna::FrameDiagonalSupport support =
      routed_election
          ? rna::routed_frame_diagonal_support(harvest.nominated_pool,
                                               harvest.opposite_pool,
                                               rna::kFrameElectionTolerance)
          : rna::compute_frame_diagonal_support(harvest.bounded_pool,
                                                rna::kFrameElectionTolerance);
  const rna::FrameElection election = rna::elect_query_frame(support);
  return election;
}


} // namespace fa::cpu::lr::rna

/**** window_chain ****/
// A catalogue rival's window: harvest, frame election, chain and orientation repair, as
// for the committed window, then the summary and the shadow rule; and the scan of the
// catalogue that chains every rival.

namespace fa::cpu::lr::rna {

// One catalogue entry's harvest and chain (harvest, frame election, orientation
// repair, summary). `committed_geometry` is the committed chain's geometry for the
// shadow rule; null while the committed chain is not fixed.
static ProductionRival
chain_one_rival(const ReadFacts& facts, ReadState& rs,
                const rna::placement::CoarseLocus& rival_locus, int catalogue_index,
                const rna::RnaSegmentGeometry* committed_geometry) {
  const auto& reference_count = facts.reference_count;
  const auto& rctx = facts.rctx;
  const auto& fwd_seeds = facts.fwd_seeds;
  auto& peaks = rs.peaks;
  const auto& read_len = facts.read_len;
  const auto& harvest_options = facts.harvest_options;
  auto& postings = facts.postings;
  const auto& read_key = facts.read_key;
  const auto& projected_shared_reverse = facts.projected_shared_reverse;
  const auto& max_intron = facts.max_intron;
  ProductionRival entry;
  entry.rival.candidate = catalogue_index;
  entry.rival.rank_score = rival_locus.rank_score;
  // Frame of q_begin/q_end: the nomination's until a repair adopts the mirror,
  // whose `reverse` is already flipped.
  entry.rival.reverse = rival_locus.reverse;
  {
    const int rival_chr = rival_locus.reference_id;
    entry.reference_id = rival_chr;
    int rival_chr_len = 0;
    if (rival_chr < 0 || rival_chr >= reference_count ||
        !mapping::checked_size_to_int(
            (*rctx.ref.encoded)[static_cast<std::size_t>(rival_chr)]
                .size(),
            rival_chr_len)) {
      entry.refusal = "rival_reference_domain";
    } else {
      const std::vector<QuerySeed>& rival_seeds =
          rival_locus.reverse ? projected_shared_reverse : *fwd_seeds;
      if (rival_seeds.empty()) {
        entry.refusal = "rival_no_strand_seeds";
      } else {
        rna::Rank1HarvestResult rival_harvest = rna::harvest_rank1(
            *rctx.ref.index, rival_seeds, rival_locus, peaks, read_len,
            rival_chr_len, harvest_options, &postings);
        if (rival_harvest.refused) {
          entry.refusal =
              rna::anchoring_refusal_name(rival_harvest.refusal);
        } else {
          // Rivals run the rank-1 window's frame policy, so they compare with the
          // winner they are priced against.
          const rna::FrameElection rival_election =
              elect_frame(facts, rs, rival_harvest);
          const bool rival_mirror_only =
              rival_election == rna::FrameElection::kMirror;
          const bool rival_nominated_only =
              rival_election == rna::FrameElection::kNominated;
          // "Never chained", until something chains it: an accepted
          // empty path is not what a skipped frame means.
          rna::ExactAnchorPathResult rival_path;
          rival_path.refused = true;
          rival_path.refusal = rna::AnchoringRefusal::NoAcceptedChain;
          bool rival_nominated_chained = false;
          if (!rival_mirror_only) {
            rival_path = rna::select_exact_anchor_path(
                read_key, read_len, rival_chr, rival_locus.reverse,
                pool_for_nominated_chain(facts, rs, rival_harvest),
                rival_harvest.anchor_path_params, &entry.siblings);
            rival_nominated_chained = true;
          }
          // Orientation repair, as for the committed window: chain this rival's
          // mirror frame and keep the better chain, unless the election chose the
          // nominated frame and it chained. A refused nominated frame escalates
          // here.
          if (!rival_nominated_only || rival_path.refused) {
            rna::placement::CoarseLocus mirror =
                rna::mirror_locus_query_orientation(
                    rival_locus, static_cast<uint32_t>(read_len));
            const std::vector<QuerySeed>& mirror_seeds =
                mirror.reverse ? projected_shared_reverse : *fwd_seeds;
            if (!mirror_seeds.empty()) {
              entry.retry = 1;
              // The mirrored pool is projected from the one just chained, as for
              // the committed window.
              project_mirror_frame(facts, rs, rival_harvest, rival_seeds.size(),
                                   mirror_seeds.size());
              if (!rival_harvest.refused) {
                rna::SiblingChainStats rival_mirror_siblings;
                rna::ExactAnchorPathResult mirror_path =
                    rna::select_exact_anchor_path(
                        read_key, read_len, rival_chr, mirror.reverse,
                        pool_for_mirror_chain(facts, rs, rival_harvest),
                        rival_harvest.anchor_path_params,
                        &rival_mirror_siblings);
                // The mirror was elected and refused: chain the nominated frame
                // after all and let the comparison decide.
                if (mirror_path.refused && !rival_nominated_chained) {
                  rival_path = rna::select_exact_anchor_path(
                      read_key, read_len, rival_chr, rival_locus.reverse,
                      pool_for_nominated_chain(facts, rs, rival_harvest),
                      rival_harvest.anchor_path_params, &entry.siblings);
                  rival_nominated_chained = true;
                }
                // Strict >: a tie keeps the nominated frame; an unchained
                // nominated frame yields to any mirror path.
                const bool adopt_mirror =
                    !mirror_path.refused &&
                    (!rival_nominated_chained || rival_path.refused ||
                     mirror_path.bundle.anchor_path_score >
                         rival_path.bundle.anchor_path_score);
                if (adopt_mirror) {
                  rna::SiblingChainStats merged = rival_mirror_siblings;
                  // Only a frame that was actually chained is a
                  // runner-up inside this window.
                  if (rival_nominated_chained)
                    merge_frame_siblings(merged, entry.siblings,
                                         rival_path);
                  rival_path = std::move(mirror_path);
                  entry.retry = 2;
                  entry.siblings = merged;
                  entry.locus = std::move(mirror);
                  entry.repaired = true;
                  entry.rival.reverse = entry.locus.reverse;
                } else if (rival_nominated_chained) {
                  merge_frame_siblings(
                      entry.siblings, rival_mirror_siblings, mirror_path);
                }
              }
            }
          }
          // The elected mirror never reached a chain call (empty mirror stream or
          // refused projection), so the pool is still the nominated frame's.
          if (!rival_nominated_chained && rival_path.refused) {
            rival_path = rna::select_exact_anchor_path(
                read_key, read_len, rival_chr, rival_locus.reverse,
                pool_for_nominated_chain(facts, rs, rival_harvest),
                rival_harvest.anchor_path_params, &entry.siblings);
            rival_nominated_chained = true;
          }
          if (rival_path.refused) {
            entry.refusal =
                rna::anchoring_refusal_name(rival_path.refusal);
          } else {
            rna::ExactAnchorPathSummary rival_summary;
            if (!rna::summarize_exact_anchor_path(
                    rival_path.bundle, rival_summary, nullptr)) {
              entry.refusal = "rival_summary_refused";
            } else {
              entry.rival.chained = true;
              entry.rival.chain_score =
                  rival_path.bundle.anchor_path_score;
              entry.rival.chain_anchors = static_cast<int>(
                  rival_path.bundle.anchor_path_anchor_count);
              entry.rival.q_begin = rival_summary.first_query_begin;
              entry.rival.q_end = rival_summary.last_query_end;
              if (!entry.repaired)
                entry.locus = rival_locus;
              entry.summary = rival_summary;
              entry.chained_path = std::move(rival_path.bundle);
              entry.owns_path = true;
              // Shadow rule (RnaChainMapqRival::shadow): a chain at the committed
              // chain's own place, seen from a neighbouring window, is not a rival.
              // Judged on geometry (rna_same_place), never on scores.
              if (committed_geometry != nullptr) {
                rna::RnaSegmentGeometry shadow_geometry;
                shadow_geometry.reference_id = rival_chr;
                shadow_geometry.reverse = entry.rival.reverse;
                shadow_geometry.reference_begin =
                    static_cast<std::uint64_t>(
                        std::max(0, rival_summary.first_reference_begin));
                shadow_geometry.reference_end =
                    static_cast<std::uint64_t>(
                        std::max(0, rival_summary.last_reference_end));
                shadow_geometry.forward_span =
                    rna::rna_forward_query_span(
                        entry.rival.q_begin, entry.rival.q_end,
                        entry.rival.reverse, read_len);
                entry.rival.shadow = rna::rna_same_place(
                    *committed_geometry, shadow_geometry, max_intron);
              }
            }
          }
        }
      }
    }
  }
  return entry;
}

static void chain_production_rivals(
    const ReadFacts& facts, ReadState& rs, rna::RnaChainMapqEvidence& evidence,
    std::vector<ProductionRival>& rivals_out, std::size_t committed_index,
    const rna::RnaSegmentGeometry* committed_geometry = nullptr) {
  auto& loci = rs.loci;
  // The scan covers the whole catalogue, already bounded by max_locus_chains.
  const std::size_t scanned = loci.size();
  for (std::size_t i = 0; i < scanned; ++i) {
    if (i == committed_index)
      continue;
    ProductionRival entry =
        chain_one_rival(facts, rs, loci[i], static_cast<int>(i), committed_geometry);
    evidence.rivals.push_back(entry.rival);
    rivals_out.push_back(std::move(entry));
  }
}

} // namespace fa::cpu::lr::rna

/**** attempt ****/
// A Stage-4 attempt on the catalogue's rank 1: the window checks, the rank-1 window
// chain, the committed path, then the plain lane's record or the CIGAR lane's
// acceptance.

namespace fa::cpu::lr::rna {

namespace {

// How a Stage-4 attempt ended; map_read's attempt loop acts on it: kRefused tries the
// next pass, kEmitted returns the plain lane's record, kAccepted
// leaves the loop for Stage 5.
enum class AttemptOutcome { kRefused, kEmitted, kAccepted };

} // namespace

// The rank-1 window chain: the window's harvest, the frame election, the chain in the
// elected frames with the orientation repair and its escalations. False when the
// harvest refuses or no path survives.
static bool chain_rank1_window(const ReadFacts& facts, ReadState& rs, AttemptState& as,
                               const rna::placement::CoarseLocus& rank1_locus, bool is_rc,
                               const std::vector<QuerySeed>& strand_seeds, int chr_len) {
  const auto& rctx = facts.rctx;
  auto& peaks = rs.peaks;
  const auto& read_len = facts.read_len;
  const auto& harvest_options = facts.harvest_options;
  auto& postings = facts.postings;
  auto& anchor_path_built = as.anchor_path_built;
  const auto& read_key = facts.read_key;
  auto& chr = as.chr;
  auto& rank1_siblings = as.rank1_siblings;
  auto& selected_reverse = as.selected_reverse;
  const auto& projected_shared_reverse = facts.projected_shared_reverse;
  const auto& fwd_seeds = facts.fwd_seeds;
  auto& mirror_siblings = as.mirror_siblings;
  auto& repaired_locus = as.repaired_locus;
  auto& selected_locus = as.selected_locus;
  auto& repaired_orientation = as.repaired_orientation;
  rna::Rank1HarvestResult rank1_harvest =
      rna::harvest_rank1(*rctx.ref.index, strand_seeds, rank1_locus, peaks,
                         read_len, chr_len, harvest_options, &postings);
  if (rank1_harvest.refused) {
    return false;
  }
  // Which frames this window will chain, decided from the pool before
  // either chain runs.
  const rna::FrameElection rank1_election =
      elect_frame(facts, rs, rank1_harvest);
  const bool rank1_mirror_only =
      rank1_election == rna::FrameElection::kMirror;
  const bool rank1_nominated_only =
      rank1_election == rna::FrameElection::kNominated;
  // "Never chained" until something chains it; downstream reads `refused`.
  anchor_path_built = rna::ExactAnchorPathResult{};
  anchor_path_built.refused = true;
  anchor_path_built.refusal = rna::AnchoringRefusal::NoAcceptedChain;
  bool nominated_chained = false;
  if (!rank1_mirror_only) {
    anchor_path_built = rna::select_exact_anchor_path(
        read_key, read_len, chr, is_rc,
        pool_for_nominated_chain(facts, rs, rank1_harvest),
        rank1_harvest.anchor_path_params, &rank1_siblings);
    nominated_chained = true;
  }
  // Orientation repair: chain the rank-1 window in the mirrored query frame too and
  // commit the better chain. Skipped when the election chose the nominated frame and
  // it chained; the only chain under a mirror election; unconditional under kAmbiguous.
  selected_reverse = is_rc;
  if (!rank1_nominated_only || anchor_path_built.refused) {
    const uint32_t query_length = static_cast<uint32_t>(read_len);
    rna::placement::CoarseLocus mirror =
        rna::mirror_locus_query_orientation(rank1_locus, query_length);
    const std::vector<QuerySeed>& mirror_seeds =
        mirror.reverse ? projected_shared_reverse : *fwd_seeds;
    if (!mirror_seeds.empty()) {
      // The mirrored frame's pool is projected in place rather than harvested again;
      // the committed chain above already copied what it keeps.
      project_mirror_frame(facts, rs, rank1_harvest, strand_seeds.size(),
                           mirror_seeds.size());
      if (!rank1_harvest.refused) {
        rna::ExactAnchorPathResult mirror_path =
            rna::select_exact_anchor_path(
                read_key, read_len, chr, mirror.reverse,
                pool_for_mirror_chain(facts, rs, rank1_harvest),
                rank1_harvest.anchor_path_params, &mirror_siblings);
        // The mirror was elected and refused: chain the nominated frame after all and
        // let the comparison decide.
        if (mirror_path.refused && !nominated_chained) {
          anchor_path_built = rna::select_exact_anchor_path(
              read_key, read_len, chr, is_rc,
              pool_for_nominated_chain(facts, rs, rank1_harvest),
              rank1_harvest.anchor_path_params, &rank1_siblings);
          nominated_chained = true;
        }
        // Strict >: a tie keeps the catalogue's frame; an unchained nominated frame
        // yields to any mirror path.
        const bool adopt_mirror =
            !mirror_path.refused &&
            (!nominated_chained || anchor_path_built.refused ||
             mirror_path.bundle.anchor_path_score >
                 anchor_path_built.bundle.anchor_path_score);
        if (adopt_mirror) {
          // Merged before the move: the demoted nominated frame, if chained, is the
          // committed window's other-frame sibling.
          rna::SiblingChainStats merged = mirror_siblings;
          if (nominated_chained)
            merge_frame_siblings(merged, rank1_siblings, anchor_path_built);
          anchor_path_built = std::move(mirror_path);
          repaired_locus = std::move(mirror);
          selected_locus = &repaired_locus;
          selected_reverse = repaired_locus.reverse;
          repaired_orientation = true;
          rank1_siblings = merged;
        } else if (nominated_chained) {
          merge_frame_siblings(rank1_siblings, mirror_siblings, mirror_path);
        }
      }
    }
    // The elected mirror never reached a chain call (empty mirror stream or refused
    // projection), so the pool is still the nominated frame's.
    if (!nominated_chained && anchor_path_built.refused) {
      anchor_path_built = rna::select_exact_anchor_path(
          read_key, read_len, chr, is_rc,
          pool_for_nominated_chain(facts, rs, rank1_harvest),
          rank1_harvest.anchor_path_params, &rank1_siblings);
      nominated_chained = true;
    }
  }
  if (anchor_path_built.refused) {
    return false;
  }
  return true;
}

// The committed path: the path summary, the committed chain's geometry for the shadow
// rule, the rescue floor and the anchor-union query coverage. False when the summary or
// the rescue floor refuses.
static bool commit_attempt_path(const ReadFacts& facts, ReadState& rs, AttemptState& as,
                                const rna::ExactAnchorPath& anchor_path) {
  auto& selected_path_summary = as.selected_path_summary;
  auto& committed_geometry = as.committed_geometry;
  auto& chr = as.chr;
  auto& selected_reverse = as.selected_reverse;
  const auto& read_len = facts.read_len;
  auto& committed_geometry_ptr = as.committed_geometry_ptr;
  auto& coarse_pass = rs.coarse_pass;
  auto& path_coverage_ok = as.path_coverage_ok;
  auto& selected_query_covered_bases = as.selected_query_covered_bases;

  rna::AnchoringRefusal summary_refusal = rna::AnchoringRefusal::None;
  if (!rna::summarize_exact_anchor_path(anchor_path, selected_path_summary,
                                        &summary_refusal)) {
    return false;
  }
  // The committed chain's geometry for the shadow rule, read by every rival's chain
  // call below.
  committed_geometry.reference_id = chr;
  committed_geometry.reverse = selected_reverse;
  committed_geometry.reference_begin = static_cast<std::uint64_t>(
      std::max(0, selected_path_summary.first_reference_begin));
  committed_geometry.reference_end = static_cast<std::uint64_t>(
      std::max(0, selected_path_summary.last_reference_end));
  committed_geometry.forward_span = rna::rna_forward_query_span(
      selected_path_summary.first_query_begin,
      selected_path_summary.last_query_end, selected_reverse, read_len);
  committed_geometry_ptr = &committed_geometry;

  // The rescue floor: a second-pass placement is a rescue, and Stage 4's chain floor
  // accepts very short chains, so it must explain at least kRnaChimeraMinQueryBases
  // query bases. First-pass placements are exempt.
  if (coarse_pass == 2 && selected_path_summary.last_query_end -
                                  selected_path_summary.first_query_begin <
                              rna::kRnaChimeraMinQueryBases) {
    return false;
  }

  // The committed chain's anchor-union query coverage: the map-only PAF matches field
  // and the brake's numerator.
  path_coverage_ok = rna::anchor_union_query_coverage(
      anchor_path, selected_query_covered_bases);
  return true;
}

// The plain lane's record (map-only PAF): the coverage brake, the chain-formula MAPQ
// evidence over the rival chains and the plain election, the projection and its tags.
// kRefused when the coverage or the projection refuses; kEmitted with the record in
// rs.out.
static AttemptOutcome emit_plain_lane(const ReadFacts& facts, ReadState& rs,
                                      AttemptState& as,
                                      const rna::ExactAnchorPath& anchor_path) {
  auto& path_coverage_ok = as.path_coverage_ok;
  auto& out = rs.out;
  auto& selected_locus = as.selected_locus;
  auto& selected_path_summary = as.selected_path_summary;
  auto& rank1_siblings = as.rank1_siblings;
  auto& catalogue_rank2_score = rs.catalogue_rank2_score;
  const auto& read_len = facts.read_len;
  auto& selected_reverse = as.selected_reverse;
  auto& selected_query_covered_bases = as.selected_query_covered_bases;
  const auto& rctx = facts.rctx;
  auto& committed_geometry_ptr = as.committed_geometry_ptr;
  const auto& reference_count = facts.reference_count;
  if (!path_coverage_ok) {
    return AttemptOutcome::kRefused;
  }
  // Chain-formula MAPQ, branch A (no realized DP).
  rna::RnaChainMapqEvidence chain_evidence;
  chain_evidence.f1 = anchor_path.anchor_path_score;
  chain_evidence.cnt =
      static_cast<int>(anchor_path.anchor_path_anchor_count);
  chain_evidence.sib_f2 = rank1_siblings.sib_score;
  chain_evidence.winner_q_begin = selected_path_summary.first_query_begin;
  chain_evidence.winner_q_end = selected_path_summary.last_query_end;
  chain_evidence.rank1_score = selected_locus->rank_score;
  // The disjoint-vote floor's numerator, from the catalogue: rank 2's vote, or
  // rank 1's when the election below adopts a rival.
  chain_evidence.rank2_score = catalogue_rank2_score;
  chain_evidence.read_len = read_len;
  // `selected_reverse` is already the post-repair strand; never XOR'd again.
  chain_evidence.winner_reverse = selected_reverse;
  chain_evidence.selected_query_covered_bases =
      selected_query_covered_bases;
  // The preset MAPQ calibration values ride the evidence (see RnaChainMapqEvidence).
  chain_evidence.qcov_tau = rctx.opts.mapq_qcov_tau;
  chain_evidence.disjoint_vote_damp = rctx.opts.mapq_disjoint_vote_damp;
  std::vector<ProductionRival> production_rivals;
  chain_production_rivals(facts, rs, chain_evidence, production_rivals,
                          /*committed_index=*/0, committed_geometry_ptr);
  const rna::placement::CoarseLocus* committed_locus = selected_locus;
  const rna::ExactAnchorPathSummary* committed_summary =
      &selected_path_summary;
  std::uint64_t committed_qcov = selected_query_covered_bases;
  rna_plain_election(facts, rs, as, anchor_path, production_rivals, chain_evidence,
                     committed_locus, committed_summary, committed_qcov);
  // The breakdown is a pure out-parameter; its f2 is the s2:i tag.
  rna::RnaChainMapqBreakdown plain_breakdown;
  const int chain_mapq =
      rna::rna_chain_mapq(chain_evidence, &plain_breakdown);
  if (!rna::placement::project_fine_path_placement(
          *committed_locus, committed_summary->first_reference_begin,
          committed_summary->last_reference_end,
          committed_summary->first_query_begin,
          committed_summary->last_query_end, committed_qcov, chain_mapq,
          read_len, reference_count, *rctx.ref.names, out)) {
    return AttemptOutcome::kRefused;
  }
  // cm:i and s1:i: the committed chain; s2:i: the best competing chain score the
  // formula weighed, seeded with the committed window's sibling (0 means neither).
  if (chain_evidence.f1 > 0) {
    out.chain_anchors = chain_evidence.cnt;
    out.chain_score = chain_evidence.f1;
    out.secondary_chain_score = plain_breakdown.f2;
  }
  return AttemptOutcome::kEmitted;
}

// One Stage-4 attempt on the catalogue's rank 1: the per-attempt reset, the window
// checks, the window chain and the committed path, then the plain lane's record or,
// with a requested CIGAR, the placement Stage 5 realizes (kAccepted).
static AttemptOutcome run_attempt(const ReadFacts& facts, ReadState& rs,
                                  AttemptState& as) {
  auto& loci = rs.loci;
  auto& chr = as.chr;
  const auto& reference_count = facts.reference_count;
  const auto& fwd_seeds = facts.fwd_seeds;
  const auto& projected_shared_reverse = facts.projected_shared_reverse;
  const auto& rctx = facts.rctx;
  auto& anchor_path_built = as.anchor_path_built;
  auto& fine_accepted = as.fine_accepted;
  const rna::placement::CoarseLocus& rank1_locus = loci.front();
  // Per-attempt reset.
  as.reset_for_attempt(rank1_locus);

  chr = rank1_locus.reference_id;
  if (chr < 0 || chr >= reference_count) {
    return AttemptOutcome::kRefused;
  }
  const bool is_rc = rank1_locus.reverse;
  const std::vector<QuerySeed>& strand_seeds =
      is_rc ? projected_shared_reverse : *fwd_seeds;
  if (strand_seeds.empty()) {
    return AttemptOutcome::kRefused;
  }
  int chr_len = 0;
  if (!mapping::checked_size_to_int(
          (*rctx.ref.encoded)[static_cast<size_t>(chr)].size(), chr_len)) {
    return AttemptOutcome::kRefused;
  }
  const int ref_begin = static_cast<int>(std::min<uint64_t>(
      rank1_locus.reference_begin, static_cast<uint64_t>(chr_len)));
  const int ref_end = static_cast<int>(std::min<uint64_t>(
      rank1_locus.reference_end, static_cast<uint64_t>(chr_len)));
  if (ref_end <= ref_begin) {
    return AttemptOutcome::kRefused;
  }

  if (!chain_rank1_window(facts, rs, as, rank1_locus, is_rc, strand_seeds, chr_len))
    return AttemptOutcome::kRefused;
  const rna::ExactAnchorPath& anchor_path = anchor_path_built.bundle;
  if (!commit_attempt_path(facts, rs, as, anchor_path))
    return AttemptOutcome::kRefused;
  if (!rctx.opts.enable_full_read_cigar) {
    return emit_plain_lane(facts, rs, as, anchor_path);
  }
  // CIGAR output: this placement is the primary and Stage 5 runs on it.
  fine_accepted = true;
  return AttemptOutcome::kAccepted;
}

} // namespace fa::cpu::lr::rna
/**** lifecycle ****/
// Stage 5's rival lifecycle: the conductor of the stages in hit.cpp, keeping the
// rival chains, the realize loop and the adoption.

namespace fa::cpu::lr::rna {

// The rival lifecycle on the CIGAR lane's accepted attempt: (1) the rival chains in
// the rank-1 frame; the hypotheses and the realization queue (hit.cpp); the queue's
// realizations; the rank recalibration, the election and the pieces (hit.cpp); then
// the adoption, which repoints the committed frame and swaps `realized` for the
// adopted rival's family.
static void run_lifecycle(const ReadFacts& facts, ReadState& rs,
                          const AttemptState& as, Lifecycle& lc, Committed& cf,
                          const rna::ExactAnchorPath& anchor_path,
                          const rna::RnaRivalLifecycleConfig& lifecycle_cfg,
                          const rna::RnaSpliceRealizationRequest& request,
                          rna::RnaSpliceRealizationResult& realized) {
  const auto& rank1_siblings = as.rank1_siblings;
  const auto& selected_path_summary = as.selected_path_summary;
  const auto& selected_locus = as.selected_locus;
  auto& catalogue_rank2_score = rs.catalogue_rank2_score;
  const auto& read_len = facts.read_len;
  const auto& selected_reverse = as.selected_reverse;
  const auto& selected_query_covered_bases = as.selected_query_covered_bases;
  const auto& rctx = facts.rctx;
  auto& production_rivals = lc.production_rivals;
  const auto& committed_geometry_ptr = as.committed_geometry_ptr;
  auto& hypotheses = lc.hypotheses;
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& demoted_incumbent = lc.demoted_incumbent;
  auto& committed_locus = cf.committed_locus;
  auto& committed_path = cf.committed_path;
  auto& committed_path_summary = cf.committed_path_summary;
  auto& committed_reverse = cf.committed_reverse;
  auto& committed_chr = cf.committed_chr;
  auto& committed_rank = cf.committed_rank;
  auto& commit_coverage_ok = cf.commit_coverage_ok;
  auto& commit_qcov = cf.commit_qcov;
  auto& committed_siblings = cf.committed_siblings;
  // (1) The rival chains, taken once in the rank-1 frame; the MAPQ below reuses them.
  rna::RnaChainMapqEvidence lifecycle_frame;
  lifecycle_frame.f1 = anchor_path.anchor_path_score;
  lifecycle_frame.cnt =
      static_cast<int>(anchor_path.anchor_path_anchor_count);
  lifecycle_frame.sib_f2 = rank1_siblings.sib_score;
  lifecycle_frame.winner_q_begin = selected_path_summary.first_query_begin;
  lifecycle_frame.winner_q_end = selected_path_summary.last_query_end;
  lifecycle_frame.rank1_score = selected_locus->rank_score;
  // Rank 1 is still the committed frame; the election below may repoint it.
  lifecycle_frame.rank2_score = catalogue_rank2_score;
  lifecycle_frame.read_len = read_len;
  lifecycle_frame.winner_reverse = selected_reverse;
  lifecycle_frame.selected_query_covered_bases = selected_query_covered_bases;
  // Kept equal to the MAPQ evidence, though chaining and classification do not read
  // it.
  lifecycle_frame.qcov_tau = rctx.opts.mapq_qcov_tau;
  lifecycle_frame.disjoint_vote_damp = rctx.opts.mapq_disjoint_vote_damp;
  chain_production_rivals(facts, rs, lifecycle_frame, production_rivals,
                          /*committed_index=*/0, committed_geometry_ptr);

  build_hypotheses(facts, rs, as, lc, anchor_path, lifecycle_cfg, realized);

  std::vector<std::size_t> queue = queue_realizations(facts, rs, lc, lifecycle_cfg);

  for (const std::size_t index : queue) {
    rna::RnaRealizedHypothesis& candidate = hypotheses[index];
    rna::RnaSpliceRealizationRequest rival_request = request;
    // A rival contig that will not bind is skipped; the incumbent is already realized.
    if (!bind_realization_reference(facts, rival_request, candidate.reference_id))
      continue;
    candidate.result = rna::realize_exact_anchor_path(
        rival_request, *candidate.bundle);
    candidate.realized = true;
    price_hypothesis(candidate, candidate.result);
  }

  rank_recalibrate(facts, lc, realized);

  elect_and_decide_pieces(facts, lc, realized);
  if (lifecycle_elected != 0) {
    rna::RnaRealizedHypothesis& winner = hypotheses[lifecycle_elected];
    // Keep the incumbent's realization before it is overwritten: on an adoption it is
    // the strongest loser, and no hypothesis holds its family.
    demoted_incumbent = std::move(realized);
    realized = std::move(winner.result);
    committed_locus = winner.locus;
    committed_path = winner.bundle;
    committed_path_summary = winner.summary;
    committed_reverse = winner.reverse;
    committed_chr = winner.reference_id;
    committed_rank = static_cast<std::size_t>(winner.catalogue_index) + 1;
    commit_coverage_ok = winner.coverage_ok;
    commit_qcov = winner.query_covered_bases;
    committed_siblings = &winner.siblings;
  }
}

} // namespace fa::cpu::lr::rna

/**** backend ****/
// RnaBackend::map_read: vote capture, exon peaks, the coarse locus catalogue and the
// rank-1 exact anchor path, then either the map-only election and projection with the
// chain MAPQ or splice realization with the rival lifecycle and the supplementary
// pieces of minimap2's parent rule.

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

AlignResult
RnaBackend::map_read(const Context& rctx, const LongReadSeedContext& seed_ctx,
                     WorkerScratch& worker_scratch, const std::string& read,
                     const std::vector<uint8_t>& fwd_enc,
                     std::vector<uint8_t>& rc_enc,
                     const std::vector<QuerySeed>* shared_fwd_syncmer_seeds,
                     std::uint32_t read_name_hash) {
  rna::Result out{};
  ReadFacts facts(rctx, seed_ctx, worker_scratch, read, fwd_enc, rc_enc,
                  shared_fwd_syncmer_seeds);
  ReadState rs(out);
  auto& read_len = facts.read_len;
  auto& reference_count = facts.reference_count;
  if (!mapping::checked_size_to_int(read.size(), read_len) ||
      !rctx.ref.valid() ||
      rctx.coordinate_domain_refusal !=
          mapping::CoordinateDomainRefusal::None ||
      rctx.ref.names->size() != rctx.ref.encoded->size() ||
      !mapping::checked_size_to_int(rctx.ref.encoded->size(),
                                    reference_count)) {
    return finish(rs, std::move(out));
  }
  out.read_len = read_len;
  // minimap2's tie seed, from the name hash and the read length. The catalogue order
  // uses it for windows equal on every vote key; every later tie between loci goes to the
  // lower catalogue index, which is that order.
  facts.tie_seed = tie_read_seed(read_name_hash, read_len);

  const auto& cfg = facts.cfg;

  if (shared_fwd_syncmer_seeds == nullptr) {
    return finish(rs, std::move(out));
  }
  // Every consumer below, the vote, the second pass, the harvests and the
  // leftover step, reads the stream without the read's repeated keys, capped
  // at the vote's N.
  {
    std::size_t dropped = 0;
    shared_fwd_syncmer_seeds = drop_query_repeats(
        *shared_fwd_syncmer_seeds,
        effective_long_primary_occ_cap(seed_ctx.occ_policy),
        rctx.opts.query_occ_frac, worker_scratch.query_occ_keys,
        worker_scratch.query_occ_seeds, dropped);
    facts.shared_fwd_syncmer_seeds = shared_fwd_syncmer_seeds;
  }
  // The read's key: the immutable exact-anchor-path bundle identity.
  facts.read_key = rna::hash_rna_read_key(read);

  if (read_len < std::max(2 * cfg.k, cfg.k)) {
    return finish(rs, std::move(out));
  }

  facts.k = cfg.k;
  facts.min_support = cfg.min_support;
  facts.min_intron = cfg.min_intron;
  facts.max_intron = cfg.max_intron;
  facts.max_loci = cfg.max_locus_chains;
  const auto& max_intron = facts.max_intron;

  // Stage 1: fused two-strand vote capture (rna.h). The
  // reverse-projected seed stream is also what the fine harvests below read.
  capture_votes(facts, rs);

  // Stages 2-3: the captured postings aggregate into exon-shaped peaks and the peaks
  // into the coarse locus catalogue.
  build_coarse_catalogue(facts, rs);
  auto& loci = rs.loci;

  // The second pass, when the first gives no catalogue; a read with no catalogue from
  // either pass is refused.
  auto& coarse_pass = rs.coarse_pass;
  if (loci.empty() && !build_second_pass_catalogue(facts, rs)) {
    return finish(rs, std::move(out));
  }

  // The catalogue's rank_score at ranks 1 and 2, the vote evidence the MAPQ reads;
  // re-read when the second pass replaces the catalogue.
  read_catalogue_rank_scores(rs);
  AttemptState as;
  // Repointed at rank 1 at the top of every Stage-4 attempt, and at
  // `repaired_locus` when the orientation repair adopts the mirror.
  as.selected_locus = &loci.front();
  auto& selected_locus = as.selected_locus;
  // Sibling-chain evidence for the committed window, always requested: the recurrence
  // has already built the chains, and the MAPQ reads them.
  auto& rank1_siblings = as.rank1_siblings;

  // Strand seed streams for the fine harvest.
  facts.fwd_seeds = shared_fwd_syncmer_seeds;

  auto& harvest_options = facts.harvest_options;
  harvest_options.global_occurrence_cap = rctx.opts.cigar_local_global_occ;
  harvest_options.readmitted = facts.readmitted;
  harvest_options.chain_query_gap =
      std::max(1, rctx.opts.fine_chain_query_gap);
  harvest_options.chain_band = rctx.opts.fine_chain_band;
  harvest_options.chain_min_count = 2;
  harvest_options.chain_min_score = 20;
  harvest_options.max_intron = max_intron;


  // Stage 4: rank 1 of the catalogue, one attempt per coarse pass: the budgeted
  // catalogue's rank 1 and, if that refuses, the second pass's. Everything downstream
  // (rivals, MAPQ, realization, lifecycle, pieces) runs on whichever catalogue placed
  // the read. All per-attempt state is reset at the top of the loop.
  auto& chr = as.chr;
  auto& fine_accepted = as.fine_accepted;
  auto& selected_reverse = as.selected_reverse;
  auto& path_coverage_ok = as.path_coverage_ok;
  auto& selected_query_covered_bases = as.selected_query_covered_bases;
  auto& anchor_path_built = as.anchor_path_built;
  auto& selected_path_summary = as.selected_path_summary;
  for (int attempt = 1; attempt <= 2; ++attempt) {
    if (attempt == 2) {
      // Attempt 1 placed nothing: the second pass gets one attempt, unless attempt 1
      // already ran on it (the budgeted vote produced no catalogue).
      if (coarse_pass == 2 || !build_second_pass_catalogue(facts, rs))
        break;
      read_catalogue_rank_scores(rs);
    }
    const AttemptOutcome outcome = run_attempt(facts, rs, as);
    if (outcome == AttemptOutcome::kRefused) {
      continue;
    }
    if (outcome == AttemptOutcome::kEmitted)
      return finish(rs, std::move(out));
    break;
  }
  if (!fine_accepted) {
    return finish_fine_failure(facts, rs, as);
  }
  const rna::ExactAnchorPath& anchor_path = anchor_path_built.bundle;

  // Stage 5: splice realization and rival arbitration.
  rna::RnaSpliceRealizationRequest request;
  request.query_forward = &fwd_enc;
  request.query_reverse = &ensure_rc_enc(facts);
  if (!bind_realization_reference(facts, request, chr)) {
    return finish(rs, std::move(out));
  }
  request.index_k = rctx.ref.index->k();
  request.options = rctx.opts.splice_controller;
  // The rival, chimeric and runner-up requests copy this one, so every emitted segment
  // carries the cs/MD strings.
  request.cigar_replay_request = rctx.opts.cigar_replay_request;
  request.scratch = &worker_scratch.realization;
  if (cfg.strand_mode == rna::StrandMode::Forward)
    request.forced_transcript_orientation =
        static_cast<int>(rna::TranscriptOrientation::Forward);
  else if (cfg.strand_mode == rna::StrandMode::Reverse)
    request.forced_transcript_orientation =
        static_cast<int>(rna::TranscriptOrientation::Reverse);
  else if (cfg.strand_mode == rna::StrandMode::None)
    // minimap2 -un: one motif-free pass, like a forced orientation.
    request.forced_transcript_orientation =
        static_cast<int>(rna::TranscriptOrientation::None);

  rna::RnaSpliceRealizationResult realized = rna::realize_exact_anchor_path(
      request, anchor_path);

  // The CIGAR lane's arbitration over the bounded catalogue's rivals. The incumbent
  // family is kept unless the election adopts a rival; ties and failures keep the
  // incumbent.
  Committed cf;
  cf.committed_locus = selected_locus;
  cf.committed_path = &anchor_path;
  auto& committed_path = cf.committed_path;
  cf.committed_path_summary = &selected_path_summary;
  cf.committed_reverse = selected_reverse;
  auto& committed_reverse = cf.committed_reverse;
  cf.committed_chr = chr;
  auto& committed_chr = cf.committed_chr;
  auto& committed_rank = cf.committed_rank;
  // The committed record's own MAPQ inputs, re-pointed on adoption.
  cf.commit_coverage_ok = path_coverage_ok;
  auto& commit_coverage_ok = cf.commit_coverage_ok;
  cf.commit_qcov = selected_query_covered_bases;
  // The committed window's sibling evidence, repointed by the election.
  cf.committed_siblings = &rank1_siblings;
  const std::optional<int> rank1_realized_dp =
      rna::rna_segment0_dp_maximum(realized);

  const rna::RnaRivalLifecycleConfig& lifecycle_cfg = rctx.opts.rival_lifecycle;
  Lifecycle lc;
  auto& hypotheses = lc.hypotheses;
  // Per hypothesis: whether a printed piece already put that locus on the wire.
  auto& locus_on_wire = lc.locus_on_wire;
  auto& lifecycle_elected = lc.lifecycle_elected;
  auto& lifecycle_dp2_owner = lc.lifecycle_dp2_owner;
  auto& pieces = lc.pieces;
  // minimap2's sub_diff (2a + b) for the parent rule.
  lc.parent_sub_diff = 2 * rctx.opts.splice_controller.match +
                       rctx.opts.splice_controller.mismatch;

  // The rival realization lifecycle, in minimap2's order: chain the rivals, classify
  // each against the incumbent by forward-frame query geometry, admit by the retention
  // band (competitors) or the co-primary floor, realize a bounded number of competitors
  // and the chain-stage co-primaries, elect the primary on the realized DP maximum, and
  // hand the demoted maxima to the MAPQ as dp2.
  if (path_coverage_ok) {
    run_lifecycle(facts, rs, as, lc, cf, anchor_path, lifecycle_cfg, request, realized);
  }
  // Set exactly when the lifecycle ran, which it declines only for a rank-1 chain
  // without query coverage.
  const bool lifecycle_ran = !hypotheses.empty();

  // The committed frame's realized DP maxima: dp1 the elected hypothesis's, dp2 the
  // lifecycle's dp2 owner's. Without the lifecycle only the incumbent's maximum exists.
  std::optional<int> committed_dp_max;
  std::optional<int> rival_dp_max;
  if (lifecycle_ran) {
    committed_dp_max = hypotheses[lifecycle_elected].dp_maximum;
    if (lifecycle_dp2_owner.has_value()) {
      const rna::RnaRealizedHypothesis& owner =
          hypotheses[*lifecycle_dp2_owner];
      rival_dp_max = owner.dp_maximum;
    }
  } else {
    committed_dp_max = rank1_realized_dp;
  }


  // The chain MAPQ in the committed frame. On adoption the adopted window owns f1,
  // sib_f2, rs1 and dp1, and the demoted incumbent becomes a chained rival owning dp2.
  int cigar_chain_mapq = 0;
  const int* chain_mapq_ptr = nullptr;
  // The pieces' store, outside the coverage guard because `chimeric_families` borrows
  // into it and both must outlive the commit. Filled only inside the guard: a read that
  // cannot be priced gets no piece.
  std::vector<rna::RnaSpliceRealizationResult> chimeric_results;
  std::vector<rna::RnaChimericFamily> chimeric_families;
  // Meaningful exactly when chain_mapq_ptr is set.
  rna::RnaChainMapqEvidence chain_evidence;
  // The breakdown of the same call, for the committed record's s2:i; stamped after the
  // commit, which overwrites `out`.
  rna::RnaChainMapqBreakdown primary_breakdown;
  if (commit_coverage_ok) {
    committed_chain_evidence(facts, rs, lc, cf, committed_dp_max, rival_dp_max,
                             realized, chain_evidence);
    cigar_chain_mapq = rna::rna_chain_mapq(chain_evidence, &primary_breakdown);
    chain_mapq_ptr = &cigar_chain_mapq;

    // The leftover, after the primary's MAPQ: every forward interval of at least
    // kRnaChimeraMinQueryBases that neither the primary's segments nor a piece covers is
    // harvested from all its syncmers (rna_leftover_pools); each contig and lane's pool
    // is chained by the lane's fine chain at minimap2's chain floor, each chain realized
    // as a co-primary, and the after-election parent rule runs again with their regions
    // behind the elected family's, in mm_hit_sort's order, so a primary among them is a
    // piece. One pass.
    if (lifecycle_ran) {
      run_leftover(facts, rs, lc, realized, request);
    }

    locus_on_wire.assign(hypotheses.size(), 0);
    // The pieces, each one flat supplementary record: its region's segment, priced by
    // the chain MAPQ on its own region (dp_max, dp_max2, n_sub, identity),
    // its hypothesis's chain and vote, its children as its rivals and its region's
    // query span as the coverage brake's denominator.
    if (lifecycle_ran && !pieces.empty()) {
      price_pieces(facts, lc, realized, chain_evidence, chimeric_results,
                   chimeric_families);
    }
  }

  // The commit: the realized family with its MAPQ and chain tags, then the pieces,
  // become the read's record.
  // cm:i, s1:i and s2:i for the committed family, stamped on the same records as its
  // MAPQ. Without committed-frame evidence (f1 == 0) the record carries none.
  rna::RnaChainTags primary_chain_tags;
  if (chain_evidence.f1 > 0) {
    primary_chain_tags.chain_anchors = chain_evidence.cnt;
    primary_chain_tags.chain_score = chain_evidence.f1;
    primary_chain_tags.secondary_chain_score = primary_breakdown.f2;
  }
  const bool commit_ok = rna::commit_realized_alignment(
      realized, read_len, chain_mapq_ptr, primary_chain_tags, chimeric_families,
      out);
  if (!commit_ok) {
    return finish(rs, std::move(out));
  }

  // The runner-ups beside the committed record: one complete realized family per
  // surviving competitor, attached as alternative hypotheses. The realizations are
  // moved out; nothing reads `hypotheses[i].result` past the commit.
  attach_runner_ups(facts, rs, lc);

  return finish(rs, std::move(out));
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
