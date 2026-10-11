// RnaBackend::map_read: vote capture, exon peaks, the coarse locus catalogue and the
// rank-1 exact anchor path, then either the map-only election and projection with the
// chain MAPQ or splice realization with the rival lifecycle and the supplementary
// pieces of minimap2's parent rule.
#include "backend.h"

#include "result.h"
#include "output.h"
#include "chain_mapq.h"
#include "rival_pricing.h" // segment-0 and overlap prices of a realization
#include "plain_elect.h"   // the map-only lane's catalogue election
#include "../core/sequence.h"                   // reverse_complement_encoded_u8
#include "../index/seed.h" // ClosedSyncmerConfig / extract_closed_syncmer_query_seeds_into
#include "../core/checked_range.h" // signed RNA coordinate/count domain
#include "../seeding/types.h"      // VoteHit
#include "../seeding/tie_hash.h"   // tie_read_seed (the exact-tie break)
#include "placement/fused_capture.h" // fused single-pass two-strand capture
#include "placement/output.h"        // project_fine_path_placement
#include "anchoring/exact_anchor_path.h" // immutable raw pool + one exact-anchor-path traceback
#include "anchoring/select_exact_anchor_path.h" // one generic chain + immutable bundle build
#include "anchoring/frame_diagonal_support.h" // pool diagonal coherence + frame election
#include "anchoring/selected_locus_anchoring.h" // checked rank-1 skeleton + one harvest
#include "realization/splice_realizer.h" // native fixed exact-anchor-path splice realization
#include "realization/rank2_arbitration.h" // realized_primary_dp_maximum
#include "realization/rival_lifecycle.h" // rival realization lifecycle
#include "realization/rna_parent.h"     // minimap2's parent rule
#include "placement/coarse_chain.h"      // coarse locus-envelope selection
#include "query_partition.h" // the re-vote, the leftover and the segment geometry

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring> // strcmp (the second pass's own typed refusal)
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

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
  auto finish = [](rna::Result result) {
    return static_cast<AlignResult&&>(std::move(result));
  };
  int read_len = 0;
  int reference_count = 0;
  if (!mapping::checked_size_to_int(read.size(), read_len) ||
      !rctx.ref.valid() ||
      rctx.coordinate_domain_refusal !=
          mapping::CoordinateDomainRefusal::None ||
      rctx.ref.names->size() != rctx.ref.encoded->size() ||
      !mapping::checked_size_to_int(rctx.ref.encoded->size(),
                                    reference_count)) {
    return finish(std::move(out));
  }
  out.read_len = read_len;
  // minimap2's tie seed, from the name hash and the read length. The catalogue order
  // uses it for windows equal on every vote key; every later tie between loci goes to the
  // lower catalogue index, which is that order.
  const std::uint32_t tie_seed = tie_read_seed(read_name_hash, read_len);

  const rna::RnaConfig& cfg = rctx.opts.cfg;

  if (shared_fwd_syncmer_seeds == nullptr) {
    return finish(std::move(out));
  }
  // The read's key: the immutable exact-anchor-path bundle identity.
  const uint64_t read_key = rna::hash_rna_read_key(read);

  if (read_len < std::max(2 * cfg.k, cfg.k)) {
    return finish(std::move(out));
  }

  const int k = cfg.k;
  const int min_support = cfg.min_support;
  const int min_intron = cfg.min_intron;
  const int max_intron = cfg.max_intron;
  const int max_loci = cfg.max_locus_chains;
  auto ensure_rc_enc = [&]() -> const std::vector<uint8_t>& {
    if (rc_enc.size() != fwd_enc.size())
      rc_enc = ::fa::cpu::lr::reverse_complement_encoded_u8(fwd_enc);
    return rc_enc;
  };

  // Stage 1: fused two-strand vote capture (placement/fused_capture.h). The
  // reverse-projected seed stream is also what the fine harvests below read.
  std::vector<QuerySeed> projected_shared_reverse;
  project_query_seeds_to_rc_coordinates(*shared_fwd_syncmer_seeds, read_len,
                                        cfg.k, projected_shared_reverse);
  // Per-read posting memo for both strands: a read runs several harvests over the same
  // two seed vectors (committed locus, its mirror, rivals), each probing every seed key.
  // The reverse projection keeps canonical keys, so reverse seed j is forward seed
  // n-1-j and one memo serves both. new_read() scopes it to this read in O(1).
  rna::SeedPostingMemo& postings = worker_scratch.postings;
  postings.new_read();
  postings.bind_forward(*shared_fwd_syncmer_seeds);
  postings.bind_reverse(projected_shared_reverse);
  rna::placement::FusedCaptureScratch& capture = worker_scratch.capture;
  rna::placement::capture_votes_fused(seed_ctx, *shared_fwd_syncmer_seeds,
                                      projected_shared_reverse, read_len,
                                      capture);

  // Stage 2: aggregate captured postings into exon-shaped peaks.
  std::vector<rna::placement::CoarseDiagonalPeak> peaks =
      rna::placement::aggregate_fused_capture_to_peaks(capture, k, min_support);
  rna::placement::CoarseLocusOptions copt =
      rna::placement::make_coarse_locus_options(
          static_cast<uint32_t>(k), static_cast<uint32_t>(min_intron),
          static_cast<uint32_t>(max_intron));
  copt.max_predecessors =
      std::max(1u, static_cast<uint32_t>(cfg.max_chain_predecessors));
  copt.max_locus_chains = static_cast<uint32_t>(max_loci);
  // No overflow: nothing downstream reads envelopes past the cut.
  copt.max_overflow_loci = 0;
  // Coalescing needs one piece that is real evidence; the floor scales with the vote
  // threshold.
  copt.coalesce_min_support =
      std::max<uint32_t>(2u, 4u * static_cast<uint32_t>(min_support));
  std::vector<rna::placement::CoarseLocus> loci;
  if (!peaks.empty()) {
    // Stage 3: coarse locus envelopes. Locus selection only; the exon path is the
    // spliced anchor chain's job (stage 4).
    loci = rna::placement::select_coarse_loci(
        peaks, copt, worker_scratch.coarse, tie_seed, /*coalesced=*/nullptr,
        /*overflow_out=*/nullptr);
  }

  // The second coarse pass. The budgeted vote needs two distinct seeds in one diagonal
  // bin, so a spliced read with one exact hit per exon may have no peak at its true
  // locus, leaving an empty catalogue or only a decoy that Stage 4 refuses. Such a read
  // gets one more pass: every extracted seed of the whole read re-voted from the
  // capture's posting views (within kRnaExplainMaxPostings), drained at
  // kRnaCoprimaryMinAnchors (a re-vote sees every seed, so two-seed agreement is a
  // coincidence), and the catalogue rebuilt by select_coarse_loci with the same options.
  int coarse_pass = 1;
  const auto build_second_pass_catalogue = [&]() -> bool {
    // Every extracted seed of the whole read, both strand frames: -u is a
    // realization-time hypothesis and constrains no vote.
    const std::vector<rna::RnaQuerySpan> whole_read{
        rna::RnaQuerySpan{0, read_len}};
    rna::RnaRevoteStats revote_stats;
    std::vector<rna::placement::CoarseDiagonalPeak> pass_peaks =
        rna::rna_revote_intervals(
            seed_ctx, *shared_fwd_syncmer_seeds, capture.views, whole_read,
            read_len, k, rna::kRnaCoprimaryMinAnchors, reference_count,
            worker_scratch.explain.revote_capture, revote_stats);
    if (revote_stats.over_budget)
      return false;
    std::vector<rna::placement::CoarseLocus> pass_loci =
        rna::placement::select_coarse_loci(
            pass_peaks, copt, worker_scratch.coarse, tie_seed,
            /*coalesced=*/nullptr, /*overflow_out=*/nullptr);
    if (pass_loci.empty())
      return false;
    peaks = std::move(pass_peaks);
    loci = std::move(pass_loci);
    coarse_pass = 2;
    return true;
  };

  if (loci.empty() && !build_second_pass_catalogue()) {
    return finish(std::move(out));
  }

  // The catalogue's rank_score at ranks 1 and 2, the vote evidence the MAPQ reads;
  // re-read when the second pass replaces the catalogue.
  int64_t catalogue_rank1_score = 0;
  int64_t catalogue_rank2_score = 0;
  const auto read_catalogue_rank_scores = [&]() {
    catalogue_rank1_score = loci[0].rank_score;
    catalogue_rank2_score = loci.size() >= 2 ? loci[1].rank_score : 0;
  };
  read_catalogue_rank_scores();
  // Orientation repair may replace the nomination with its mirror locus, which needs
  // storage the catalogue does not own; both output modes project the same one.
  rna::placement::CoarseLocus repaired_locus;
  // Repointed at rank 1 at the top of every Stage-4 attempt, and at
  // `repaired_locus` when the orientation repair adopts the mirror.
  const rna::placement::CoarseLocus* selected_locus = &loci.front();
  // Sibling-chain evidence for the committed window, always requested: the recurrence
  // has already built the chains, and the MAPQ reads them.
  rna::SiblingChainStats rank1_siblings;
  rna::SiblingChainStats mirror_siblings;
  // The fine chain is placement evidence in every output mode. A read whose harvest or
  // anchor path cannot be formed is left unmapped: a coarse locus is a vote envelope,
  // not placement evidence.
  auto finish_fine_failure = [&]() -> AlignResult {
    return finish(std::move(out));
  };

  // Strand seed streams for the fine harvest.
  const std::vector<QuerySeed>* fwd_seeds = shared_fwd_syncmer_seeds;

  rna::Rank1HarvestOptions harvest_options;
  harvest_options.global_occurrence_cap = rctx.opts.cigar_local_global_occ;
  harvest_options.chain_query_gap =
      std::max(1, rctx.opts.cigar_local_interval_anchor_chain_max_gap);
  harvest_options.chain_band = rctx.opts.fine_chain_band;
  harvest_options.chain_min_count = 2;
  harvest_options.chain_min_score = 20;
  harvest_options.max_intron = max_intron;

  // Each fine-stage window may be chained in both query frames, keeping the better
  // chain, as minimap2 does: seed keys are canonical, so a window can be nominated on the
  // wrong strand, where the true anchors are anti-diagonal but a small spurious chain may
  // still pass min_cnt/min_sc.
  // With both frames chained, the uncommitted frame's best chain is a runner-up inside
  // the same window: chain counts add, and the sibling is the better of the committed
  // frame's own runner-up and the other frame's best. A frame whose recurrence never ran
  // (-1) contributes nothing.
  const auto merge_frame_siblings =
      [](rna::SiblingChainStats& committed, const rna::SiblingChainStats& other,
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
      };

  // The mirrored frame's pool, from `harvest` in place; the mirror is always the
  // projection, never a second harvest. Its preconditions, equal stream sizes (the
  // n-1-j premise) and a single search region, are structural: a window violating one is
  // refused, and `harvest`'s nominated pool is left for the escalation to chain.
  const auto project_mirror_frame = [&](rna::Rank1HarvestResult& harvest,
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
  };
  // Which pool each frame chains.
  const auto pool_for_nominated_chain =
      [](rna::Rank1HarvestResult& harvest) -> const rna::CandidatePool& {
    return harvest.nominated_pool;
  };
  const auto pool_for_mirror_chain =
      [](rna::Rank1HarvestResult& harvest) -> const rna::CandidatePool& {
    return harvest.opposite_pool;
  };
  // Pre-chaining query-frame election: where the pool's band margin is decisive only the
  // elected frame is chained; an elected frame that refuses escalates to the other, and
  // kAmbiguous chains both and lets the chain scores decide.
  const auto elect_frame = [&](const rna::Rank1HarvestResult& harvest) {
    // Routed: the nominated band comes from the harvested lane's pool and the mirror band
    // from the other lane's (frame_diagonal_support.h). A refused lane contributes zeros.
    const bool routed_election = harvest.routed;
    const rna::FrameDiagonalSupport support =
        routed_election
            ? rna::routed_frame_diagonal_support(harvest.nominated_pool,
                                                 harvest.opposite_pool,
                                                 rna::kFrameElectionTolerance)
            : rna::compute_frame_diagonal_support(harvest.bounded_pool,
                                                  rna::kFrameElectionTolerance);
    return rna::elect_query_frame(support);
  };

  // A catalogue rival's chain evidence for the MAPQ formula. A rival that could not be
  // chained still competes through its vote ratio (censored evidence is not zero
  // evidence) and never counts toward n_sub.
  struct ProductionRival {
    rna::RnaChainMapqRival rival;
    const char* refusal = "none";
    int retry = 0;
    // The rival's chain, kept so the rival can be realized without
    // re-harvesting its window; `owns_path` says whether `chained_path` holds one.
    rna::ExactAnchorPath chained_path;
    bool owns_path = false;
    // The locus the chain was taken on (the catalogue entry, or its mirror when
    // `repaired`), with that window's summary and sibling evidence; `reference_id` is the
    // harvested contig.
    rna::placement::CoarseLocus locus;
    rna::ExactAnchorPathSummary summary;
    rna::SiblingChainStats siblings;
    int reference_id = -1;
    bool repaired = false;
  };
  // One catalogue entry's harvest and chain (harvest, frame election, orientation
  // repair, summary). `committed_geometry` is the committed chain's geometry for the
  // shadow rule; null while the committed chain is not fixed.
  const auto chain_one_rival =
      [&](const rna::placement::CoarseLocus& rival_locus, int catalogue_index,
          const rna::RnaSegmentGeometry* committed_geometry) {
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
                    elect_frame(rival_harvest);
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
                      pool_for_nominated_chain(rival_harvest),
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
                    project_mirror_frame(rival_harvest, rival_seeds.size(),
                                         mirror_seeds.size());
                    if (!rival_harvest.refused) {
                      rna::SiblingChainStats rival_mirror_siblings;
                      rna::ExactAnchorPathResult mirror_path =
                          rna::select_exact_anchor_path(
                              read_key, read_len, rival_chr, mirror.reverse,
                              pool_for_mirror_chain(rival_harvest),
                              rival_harvest.anchor_path_params,
                              &rival_mirror_siblings);
                      // The mirror was elected and refused: chain the nominated frame
                      // after all and let the comparison decide.
                      if (mirror_path.refused && !rival_nominated_chained) {
                        rival_path = rna::select_exact_anchor_path(
                            read_key, read_len, rival_chr, rival_locus.reverse,
                            pool_for_nominated_chain(rival_harvest),
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
                      pool_for_nominated_chain(rival_harvest),
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
      };
  const auto chain_production_rivals =
      [&](rna::RnaChainMapqEvidence& evidence,
          std::vector<ProductionRival>& rivals_out,
          std::size_t committed_index,
          const rna::RnaSegmentGeometry* committed_geometry = nullptr) {
        // The scan covers the whole catalogue, already bounded by max_locus_chains.
        const std::size_t scanned = loci.size();
        for (std::size_t i = 0; i < scanned; ++i) {
          if (i == committed_index)
            continue;
          ProductionRival entry =
              chain_one_rival(loci[i], static_cast<int>(i), committed_geometry);
          evidence.rivals.push_back(entry.rival);
          rivals_out.push_back(std::move(entry));
        }
      };

  // Stage 4: rank 1 of the catalogue, one attempt per coarse pass: the budgeted
  // catalogue's rank 1 and, if that refuses, the second pass's. Everything downstream
  // (rivals, MAPQ, realization, lifecycle, pieces) runs on whichever catalogue placed
  // the read. All per-attempt state is reset at the top of the loop.
  int chr = -1;
  bool fine_accepted = false;
  bool selected_reverse = false;
  bool repaired_orientation = false;
  bool path_coverage_ok = false;
  std::uint64_t selected_query_covered_bases = 0;
  rna::ExactAnchorPathResult anchor_path_built;
  rna::ExactAnchorPathSummary selected_path_summary;
  // The committed chain's geometry for the shadow rule; null until this attempt's chain
  // and summary are fixed.
  rna::RnaSegmentGeometry committed_geometry;
  const rna::RnaSegmentGeometry* committed_geometry_ptr = nullptr;
  for (int attempt = 1; attempt <= 2; ++attempt) {
    if (attempt == 2) {
      // Attempt 1 placed nothing: the second pass gets one attempt, unless attempt 1
      // already ran on it (the budgeted vote produced no catalogue).
      if (coarse_pass == 2 || !build_second_pass_catalogue())
        break;
      read_catalogue_rank_scores();
    }
    const rna::placement::CoarseLocus& rank1_locus = loci.front();
    // Per-attempt reset.
    selected_locus = &rank1_locus;
    selected_reverse = false;
    repaired_orientation = false;
    path_coverage_ok = false;
    selected_query_covered_bases = 0;
    rank1_siblings = rna::SiblingChainStats{};
    mirror_siblings = rna::SiblingChainStats{};

    chr = rank1_locus.reference_id;
    if (chr < 0 || chr >= reference_count) {
      continue;
    }
    const bool is_rc = rank1_locus.reverse;
    const std::vector<QuerySeed>& strand_seeds =
        is_rc ? projected_shared_reverse : *fwd_seeds;
    if (strand_seeds.empty()) {
      continue;
    }
    int chr_len = 0;
    if (!mapping::checked_size_to_int(
            (*rctx.ref.encoded)[static_cast<size_t>(chr)].size(), chr_len)) {
      continue;
    }
    const int ref_begin = static_cast<int>(std::min<uint64_t>(
        rank1_locus.reference_begin, static_cast<uint64_t>(chr_len)));
    const int ref_end = static_cast<int>(std::min<uint64_t>(
        rank1_locus.reference_end, static_cast<uint64_t>(chr_len)));
    if (ref_end <= ref_begin) {
      continue;
    }

    rna::Rank1HarvestResult rank1_harvest =
        rna::harvest_rank1(*rctx.ref.index, strand_seeds, rank1_locus, peaks,
                           read_len, chr_len, harvest_options, &postings);
    if (rank1_harvest.refused) {
      continue;
    }
    // Which frames this window will chain, decided from the pool before
    // either chain runs.
    const rna::FrameElection rank1_election = elect_frame(rank1_harvest);
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
          pool_for_nominated_chain(rank1_harvest),
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
        project_mirror_frame(rank1_harvest, strand_seeds.size(),
                             mirror_seeds.size());
        if (!rank1_harvest.refused) {
          rna::ExactAnchorPathResult mirror_path =
              rna::select_exact_anchor_path(
                  read_key, read_len, chr, mirror.reverse,
                  pool_for_mirror_chain(rank1_harvest),
                  rank1_harvest.anchor_path_params, &mirror_siblings);
          // The mirror was elected and refused: chain the nominated frame after all and
          // let the comparison decide.
          if (mirror_path.refused && !nominated_chained) {
            anchor_path_built = rna::select_exact_anchor_path(
                read_key, read_len, chr, is_rc,
                pool_for_nominated_chain(rank1_harvest),
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
            pool_for_nominated_chain(rank1_harvest),
            rank1_harvest.anchor_path_params, &rank1_siblings);
        nominated_chained = true;
      }
    }
    if (anchor_path_built.refused) {
      continue;
    }
    const rna::ExactAnchorPath& anchor_path = anchor_path_built.bundle;

    rna::AnchoringRefusal summary_refusal = rna::AnchoringRefusal::None;
    if (!rna::summarize_exact_anchor_path(anchor_path, selected_path_summary,
                                          &summary_refusal)) {
      continue;
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
      continue;
    }

    // The committed chain's anchor-union query coverage: the map-only PAF matches field
    // and the brake's numerator.
    path_coverage_ok = rna::anchor_union_query_coverage(
        anchor_path, selected_query_covered_bases);
    if (!rctx.opts.enable_full_read_cigar) {
      if (!path_coverage_ok) {
        continue;
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
      chain_production_rivals(chain_evidence, production_rivals,
                              /*committed_index=*/0, committed_geometry_ptr);
      // The election (plain_elect.h): the incumbent and every chained rival are scored,
      // and the best eligible rival takes the primary when it scores strictly higher; a
      // tie between rivals goes to the lower catalogue index. An adoption swaps the roles
      // as the CIGAR lane's lifecycle does: the elected chain owns the MAPQ evidence and
      // the record, the demoted incumbent joins the rivals in catalogue order, rs2 is
      // rank 1's vote, and every shadow verdict is re-read against the elected chain.
      const rna::placement::CoarseLocus* committed_locus = selected_locus;
      const rna::ExactAnchorPathSummary* committed_summary =
          &selected_path_summary;
      std::uint64_t committed_qcov = selected_query_covered_bases;
      {
        const std::int64_t incumbent_score =
            rna::rna_plain_elect_score(anchor_path, min_intron);
        const std::size_t rivals = production_rivals.size();
        std::vector<std::int64_t> rival_scores(rivals, 0);
        std::vector<std::uint64_t> rival_qcov(rivals, 0);
        std::optional<std::size_t> elected;
        for (std::size_t i = 0; i < rivals; ++i) {
          const ProductionRival& entry = production_rivals[i];
          if (!entry.rival.chained)
            continue;
          const int qspan = entry.rival.q_end - entry.rival.q_begin;
          rival_scores[i] =
              rna::rna_plain_elect_score(entry.chained_path, min_intron);
          // Eligible on the incumbent's terms: an anchor union the record can report,
          // and on a second-pass catalogue the rescue floor.
          const bool eligible =
              rna::anchor_union_query_coverage(entry.chained_path,
                                               rival_qcov[i]) &&
              (coarse_pass != 2 || qspan >= rna::kRnaChimeraMinQueryBases);
          if (eligible &&
              (!elected || rival_scores[i] > rival_scores[*elected]))
            elected = i;
        }
        if (elected && rival_scores[*elected] > incumbent_score) {
          const ProductionRival& winner = production_rivals[*elected];
          committed_locus = &winner.locus;
          committed_summary = &winner.summary;
          committed_qcov = rival_qcov[*elected];
          chain_evidence.f1 = winner.rival.chain_score;
          chain_evidence.cnt = winner.rival.chain_anchors;
          chain_evidence.sib_f2 = winner.siblings.sib_score;
          chain_evidence.winner_q_begin = winner.rival.q_begin;
          chain_evidence.winner_q_end = winner.rival.q_end;
          chain_evidence.rank1_score = winner.locus.rank_score;
          chain_evidence.rank2_score = catalogue_rank1_score;
          chain_evidence.winner_reverse = winner.rival.reverse;
          chain_evidence.selected_query_covered_bases = committed_qcov;
          // The shadow rule's geometry, as chain_one_rival builds it.
          const auto chain_geometry = [&](const ProductionRival& entry) {
            rna::RnaSegmentGeometry geometry;
            geometry.reference_id = entry.reference_id;
            geometry.reverse = entry.rival.reverse;
            geometry.reference_begin = static_cast<std::uint64_t>(
                std::max(0, entry.summary.first_reference_begin));
            geometry.reference_end = static_cast<std::uint64_t>(
                std::max(0, entry.summary.last_reference_end));
            geometry.forward_span = rna::rna_forward_query_span(
                entry.rival.q_begin, entry.rival.q_end, entry.rival.reverse,
                read_len);
            return geometry;
          };
          const rna::RnaSegmentGeometry winner_geometry =
              chain_geometry(winner);
          chain_evidence.rivals.clear();
          rna::RnaChainMapqRival demoted;
          demoted.candidate = 0;
          demoted.chained = true;
          demoted.chain_score = anchor_path.anchor_path_score;
          demoted.chain_anchors =
              static_cast<int>(anchor_path.anchor_path_anchor_count);
          demoted.q_begin = selected_path_summary.first_query_begin;
          demoted.q_end = selected_path_summary.last_query_end;
          demoted.rank_score = selected_locus->rank_score;
          demoted.reverse = selected_reverse;
          // committed_geometry is still the incumbent's chain.
          demoted.shadow = rna::rna_same_place(winner_geometry,
                                               committed_geometry, max_intron);
          chain_evidence.rivals.push_back(demoted);
          for (std::size_t i = 0; i < rivals; ++i) {
            if (i == *elected)
              continue;
            rna::RnaChainMapqRival rival = production_rivals[i].rival;
            rival.shadow =
                rival.chained &&
                rna::rna_same_place(winner_geometry,
                                    chain_geometry(production_rivals[i]),
                                    max_intron);
            chain_evidence.rivals.push_back(rival);
          }
        }
      }
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
        continue;
      }
      // cm:i and s1:i: the committed chain; s2:i: the best competing chain score the
      // formula weighed, seeded with the committed window's sibling (0 means neither).
      if (chain_evidence.f1 > 0) {
        out.chain_anchors = chain_evidence.cnt;
        out.chain_score = chain_evidence.f1;
        out.secondary_chain_score = plain_breakdown.f2;
      }
      return finish(std::move(out));
    }
    // CIGAR output: this placement is the primary and Stage 5 runs on it.
    fine_accepted = true;
    break;
  }
  if (!fine_accepted) {
    return finish_fine_failure();
  }
  const rna::ExactAnchorPath& anchor_path = anchor_path_built.bundle;

  // Stage 5: splice realization and rival arbitration.
  rna::RnaSpliceRealizationRequest request;
  request.query_forward = &fwd_enc;
  request.query_reverse = &ensure_rc_enc();
  const auto bind_realization_reference =
      [&](rna::RnaSpliceRealizationRequest& target, int reference_id) {
        return rna::bind_splice_realization_reference(
            target, *rctx.ref.encoded, *rctx.ref.names,
            rctx.opts.known_junctions.get(), reference_id);
      };
  if (!bind_realization_reference(request, chr)) {
    return finish(std::move(out));
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

  rna::RnaSpliceRealizationResult realized =
      rna::realize_exact_anchor_path(request, anchor_path);

  // The CIGAR lane's arbitration over the bounded catalogue's rivals. The incumbent
  // family is kept unless the election adopts a rival; ties and failures keep the
  // incumbent.
  const rna::placement::CoarseLocus* committed_locus = selected_locus;
  const rna::ExactAnchorPath* committed_path = &anchor_path;
  const rna::ExactAnchorPathSummary* committed_path_summary =
      &selected_path_summary;
  bool committed_reverse = selected_reverse;
  int committed_chr = chr;
  size_t committed_rank = 1;
  // The committed record's own MAPQ inputs, re-pointed on adoption.
  bool commit_coverage_ok = path_coverage_ok;
  std::uint64_t commit_qcov = selected_query_covered_bases;
  // The committed window's sibling evidence, repointed by the election.
  const rna::SiblingChainStats* committed_siblings = &rank1_siblings;
  const std::optional<int> rank1_realized_dp =
      rna::realized_primary_dp_maximum(realized);

  const rna::RnaRivalLifecycleConfig& lifecycle_cfg = rctx.opts.rival_lifecycle;
  // Prices one realization for the lifecycle: the election price is segment 0's
  // dp_maximum, and the realized hull feeds the recalibration gate and the emission
  // floor.
  const auto price_hypothesis =
      [](rna::RnaRealizedHypothesis& h,
         const rna::RnaSpliceRealizationResult& result) {
        h.dp_segment0 = rna::rna_segment0_dp_maximum(result);
        if (const auto hull = rna::rna_realized_query_span(result)) {
          h.aq_begin = hull->begin;
          h.aq_end = hull->end;
        }
        h.dp_maximum = h.dp_segment0;
      };
  // The rival set, chained once for the lifecycle and reused by the MAPQ.
  std::vector<ProductionRival> production_rivals;
  // The leftover's chains, which their hypotheses borrow.
  std::vector<ProductionRival> leftover_rivals;
  std::vector<rna::RnaRealizedHypothesis> hypotheses;
  // The demoted incumbent's realization, parked here on an adoption. Element 0's
  // `result` stays empty by the borrowed-storage contract (rival_lifecycle.h). Read by
  // the pieces and the runner-up attach.
  rna::RnaSpliceRealizationResult demoted_incumbent;
  // Per hypothesis: whether a printed piece already put that locus on the wire.
  std::vector<char> locus_on_wire;
  std::size_t lifecycle_elected = 0;
  std::optional<std::size_t> lifecycle_dp2_owner;
  // Whether the rank recalibration repriced this read; selects the election rule.
  bool rank_recal_fired = false;
  // The parent rule's chain-stage hits (owner = hypothesis index), in the queue's order.
  std::vector<rna::RnaParentHit> parent_chain_hits;
  // The supplementary pieces: an after-DP primary outside the elected family, its count
  // of children and dp_max2 without its shadows, the hypotheses owning its children at
  // either stage, and the first child at dp_max2.
  struct PieceRegion {
    rna::RnaParentHit region;
    std::vector<std::size_t> children;
    std::optional<std::size_t> dp2_owner;
    // Its hypothesis and every hypothesis at its place, kept off the runner-ups.
    std::vector<std::size_t> place;
  };
  std::vector<PieceRegion> pieces;
  // minimap2's mask_level and sub_diff (2a + b) for the parent rule.
  constexpr float kParentMaskLevel = 0.5f;
  const int parent_sub_diff = 2 * rctx.opts.splice_controller.match +
                              rctx.opts.splice_controller.mismatch;
  // mm_set_parent over hits already in their stage's order, less the hits whose owner
  // `left_out` marks.
  const auto parent_regions = [&](const std::vector<rna::RnaParentHit>& ordered,
                                  const std::vector<char>* left_out) {
    std::vector<rna::RnaParentHit> regions;
    regions.reserve(ordered.size());
    for (const rna::RnaParentHit& hit : ordered) {
      const std::size_t owner = static_cast<std::size_t>(hit.owner);
      if (left_out != nullptr && hit.owner >= 0 && owner < left_out->size() &&
          (*left_out)[owner] != 0)
        continue;
      regions.push_back(hit);
    }
    rna::rna_parent_set(regions, kParentMaskLevel, parent_sub_diff);
    return regions;
  };
  // A chained hypothesis's place for the shadow rule; false without a chain.
  const auto shadow_geometry = [&](const rna::RnaRealizedHypothesis& h,
                                   rna::RnaSegmentGeometry& g) {
    if (!h.chained || h.summary == nullptr ||
        h.summary->last_reference_end <= h.summary->first_reference_begin)
      return false;
    g.reference_id = h.reference_id;
    g.reverse = h.reverse;
    g.reference_begin = static_cast<std::uint64_t>(
        std::max(0, h.summary->first_reference_begin));
    g.reference_end = static_cast<std::uint64_t>(
        std::max(0, h.summary->last_reference_end));
    g.forward_span = rna::rna_forward_query_span(h.q_begin, h.q_end,
                                                 h.reverse, read_len);
    return true;
  };
  const auto same_place = [&](std::size_t a, std::size_t b) {
    rna::RnaSegmentGeometry ga;
    rna::RnaSegmentGeometry gb;
    return shadow_geometry(hypotheses[a], ga) &&
           shadow_geometry(hypotheses[b], gb) &&
           rna::rna_same_place(ga, gb, max_intron);
  };
  // Hypothesis i's realized regions, one per segment of `result`, as the parent rule
  // reads them.
  const auto append_dp_hits = [&](std::size_t i,
                                  const rna::RnaSpliceRealizationResult& result,
                                  std::vector<rna::RnaParentHit>& hits) {
    const rna::RnaRealizedHypothesis& h = hypotheses[i];
    const rna::RnaSpliceHypothesisResult* winner =
        rna::rna_winning_hypothesis(result);
    if (winner == nullptr)
      return;
    for (std::size_t s = 0; s < winner->segments.size(); ++s) {
      const rna::SpliceRegionEvidence& evidence = winner->region_evidence[s];
      const AlignResult& segment = winner->segments[s];
      const std::vector<std::uint32_t>& cigar = winner->segment_cigars[s];
      rna::RnaParentHit hit;
      hit.qs = evidence.query_begin;
      hit.qe = evidence.query_end;
      hit.score = evidence.chain_score;
      hit.cnt = evidence.region_anchor_count;
      hit.rid = h.reference_id;
      hit.rs = segment.pos;
      hit.re = segment.target_end;
      hit.has_dp = true;
      hit.dp_max = evidence.dp_maximum;
      hit.accounting.matches = segment.matches;
      hit.accounting.mismatches = segment.mismatches;
      hit.accounting.ambiguities = segment.ambiguities;
      hit.accounting.cigar = cigar.data();
      hit.accounting.cigar_len = cigar.size();
      hit.accounting_ok = segment.alignment_accounting_valid && !cigar.empty();
      hit.owner = static_cast<int>(i);
      hit.catalogue = h.catalogue_index;
      hit.segment = static_cast<int>(s);
      hits.push_back(hit);
    }
  };
  // The after-election parent rule's input, kept for the leftover's rerun: its regions
  // in order (the elected family's first), the hypotheses it leaves out, and the chain
  // stage's regions.
  std::vector<rna::RnaParentHit> pieces_dp_input;
  std::vector<char> pieces_left_out;
  std::vector<rna::RnaParentHit> pieces_chain;
  // The pieces from the after-election regions `dp` and the chain stage's `chain`,
  // rebuilt from scratch: every after-DP primary but the elected family's own regions, in
  // hit order, unless its hypothesis is at the place of a printed record. Its children
  // less its shadows are counted as mm_set_parent counts them.
  const auto decide_pieces = [&](const std::vector<rna::RnaParentHit>& dp,
                                 const std::vector<rna::RnaParentHit>& chain) {
    pieces.clear();
    for (std::size_t j = 1; j < dp.size(); ++j) {
      const rna::RnaParentHit& region = dp[j];
      const std::size_t owner = static_cast<std::size_t>(region.owner);
      if (region.parent != static_cast<int>(j) || owner == lifecycle_elected ||
          same_place(owner, lifecycle_elected))
        continue;
      if (std::any_of(pieces.begin(), pieces.end(),
                      [&](const PieceRegion& printed) {
                        const std::size_t other =
                            static_cast<std::size_t>(printed.region.owner);
                        return other != owner && same_place(other, owner);
                      }))
        continue;
      PieceRegion piece;
      piece.region = region;
      piece.region.n_sub = 0;
      piece.region.dp_max2 = 0;
      const auto add_child = [&](std::size_t child_owner) {
        if (child_owner != owner &&
            std::find(piece.children.begin(), piece.children.end(),
                      child_owner) == piece.children.end())
          piece.children.push_back(child_owner);
      };
      for (std::size_t c = j + 1; c < dp.size(); ++c) {
        const rna::RnaParentHit& child = dp[c];
        const std::size_t child_owner = static_cast<std::size_t>(child.owner);
        if (child.parent != static_cast<int>(j) ||
            (child_owner != owner && same_place(owner, child_owner)))
          continue;
        const int si = child.qs;
        const int ei = child.qe;
        const int sj = region.qs;
        const int ej = region.qe;
        const int min = std::min(ej - sj, ei - si);
        const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                               : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
        bool cnt_sub = child.cnt >= region.cnt;
        if (region.has_dp && child.has_dp &&
            (region.rid != child.rid || region.rs != child.rs ||
             region.re != child.re || ol != min)) {
          piece.region.dp_max2 = std::max(piece.region.dp_max2, child.dp_max);
          if (region.dp_max - child.dp_max <= parent_sub_diff)
            cnt_sub = true;
        }
        if (cnt_sub)
          ++piece.region.n_sub;
        add_child(child_owner);
      }
      for (std::size_t c = j + 1; c < dp.size() && !piece.dp2_owner; ++c) {
        const rna::RnaParentHit& child = dp[c];
        const std::size_t child_owner = static_cast<std::size_t>(child.owner);
        if (child.parent == static_cast<int>(j) && piece.region.dp_max2 > 0 &&
            child.dp_max == piece.region.dp_max2 && child_owner != owner &&
            !same_place(owner, child_owner))
          piece.dp2_owner = child_owner;
      }
      // Its chain's chain-stage children, when that chain was a primary there.
      for (std::size_t e = 0; e < chain.size(); ++e) {
        if (static_cast<std::size_t>(chain[e].owner) != owner)
          continue;
        if (chain[e].parent == static_cast<int>(e))
          for (std::size_t c = e + 1; c < chain.size(); ++c) {
            const std::size_t child_owner =
                static_cast<std::size_t>(chain[c].owner);
            if (chain[c].parent == static_cast<int>(e) &&
                !same_place(owner, child_owner))
              add_child(child_owner);
          }
        break;
      }
      piece.place.push_back(owner);
      for (std::size_t i = 0; i < hypotheses.size(); ++i)
        if (i != owner && same_place(owner, i))
          piece.place.push_back(i);
      pieces.push_back(std::move(piece));
    }
  };

  // The rival realization lifecycle, in minimap2's order: chain the rivals, classify
  // each against the incumbent by forward-frame query geometry, admit by the retention
  // band (competitors) or the co-primary floor, realize a bounded number of competitors
  // and the chain-stage co-primaries, elect the primary on the realized DP maximum, and
  // hand the demoted maxima to the MAPQ as dp2.
  if (path_coverage_ok) {
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
    chain_production_rivals(lifecycle_frame, production_rivals,
                            /*committed_index=*/0, committed_geometry_ptr);

    // (2) The hypothesis vector, incumbent first. Every pointer borrows storage that
    // outlives it (the rank-1 chain, or a rival chain in `production_rivals`, which no
    // longer grows), so reordering `hypotheses` moves no bundle.
    hypotheses.reserve(production_rivals.size() + 1);
    {
      rna::RnaRealizedHypothesis incumbent;
      incumbent.catalogue_index = 0;
      incumbent.reference_id = chr;
      incumbent.locus = selected_locus;
      incumbent.bundle = &anchor_path;
      incumbent.summary = &selected_path_summary;
      incumbent.siblings = rank1_siblings;
      incumbent.rank_score = selected_locus->rank_score;
      incumbent.chained = true;
      incumbent.chain_score = anchor_path.anchor_path_score;
      incumbent.chain_anchors =
          static_cast<int>(anchor_path.anchor_path_anchor_count);
      incumbent.q_begin = selected_path_summary.first_query_begin;
      incumbent.q_end = selected_path_summary.last_query_end;
      incumbent.reverse = selected_reverse;
      incumbent.repaired = repaired_orientation;
      incumbent.query_covered_bases = selected_query_covered_bases;
      incumbent.coverage_ok = path_coverage_ok;
      incumbent.klass = rna::RnaRivalClass::Incumbent;
      incumbent.admitted = true;
      // Rank 1 was realized above into `realized`; the hypothesis carries only its
      // price.
      incumbent.realized = true;
      price_hypothesis(incumbent, realized);
      hypotheses.push_back(std::move(incumbent));
    }

    // (3) Classification and admission, in the read's forward frame: chains on opposite
    // strands can carry identical oriented spans over opposite ends of the read.
    // `reverse` is the post-repair strand on both sides.
    const rna::RnaQuerySpan incumbent_span = rna::rna_forward_query_span(
        selected_path_summary.first_query_begin,
        selected_path_summary.last_query_end, selected_reverse, read_len);
    for (const ProductionRival& entry : production_rivals) {
      const std::size_t cat =
          static_cast<std::size_t>(std::max(0, entry.rival.candidate));
      rna::RnaRealizedHypothesis h;
      h.catalogue_index = entry.rival.candidate;
      h.reference_id = entry.reference_id;
      // A chained rival owns its window (the catalogue entry or its adopted mirror); an
      // unchained one has only the catalogue's nomination.
      h.locus = entry.rival.chained || entry.repaired
                    ? &entry.locus
                    : (cat < loci.size() ? &loci[cat] : nullptr);
      h.rank_score = entry.rival.rank_score;
      h.bundle = entry.owns_path ? &entry.chained_path : nullptr;
      h.summary = &entry.summary;
      h.siblings = entry.siblings;
      h.chained = entry.rival.chained;
      h.chain_score = entry.rival.chain_score;
      h.chain_anchors = entry.rival.chain_anchors;
      h.q_begin = entry.rival.q_begin;
      h.q_end = entry.rival.q_end;
      h.reverse = entry.rival.reverse;
      h.repaired = entry.repaired;
      h.shadow = entry.rival.shadow;
      const rna::RnaQuerySpan rival_span =
          rna::rna_forward_query_span(entry.rival.q_begin, entry.rival.q_end,
                                      entry.rival.reverse, read_len);
      h.klass = rna::rna_classify_rival(
          rna::rna_query_spans_compete(incumbent_span, rival_span));
      if (h.chained && h.bundle != nullptr) {
        h.coverage_ok =
            rna::anchor_union_query_coverage(*h.bundle, h.query_covered_bases);
        if (h.klass == rna::RnaRivalClass::CoPrimary) {
          // A co-primary explains different query bases, so it is neither an election
          // candidate nor dp2/n_sub evidence; its admission is read by the traces only.
          h.admitted = rna::rna_coprimary_admitted(
              h.chain_score, h.chain_anchors,
              rna::rna_uncovered_query_bases(incumbent_span, rival_span));
        } else {
          h.admitted =
              h.coverage_ok &&
              rna::rna_competitor_admitted(
                  h.chain_score, anchor_path.anchor_path_score, lifecycle_cfg);
        }
      }
      hypotheses.push_back(std::move(h));
    }

    // (4) The realization budget, spent on admitted competitors in a total order (chain
    // score, anchor count, catalogue index) so the choice does not depend on the sort.
    std::vector<std::size_t> queue;
    for (std::size_t i = 1; i < hypotheses.size(); ++i) {
      const rna::RnaRealizedHypothesis& candidate = hypotheses[i];
      if (candidate.admitted &&
          candidate.klass == rna::RnaRivalClass::Competitor &&
          candidate.bundle != nullptr)
        queue.push_back(i);
    }
    std::sort(queue.begin(), queue.end(),
              [&hypotheses](std::size_t a, std::size_t b) {
                return rna::rna_rival_queue_before(hypotheses[a],
                                                   hypotheses[b]);
              });
    const std::size_t realize_budget =
        static_cast<std::size_t>(std::max(0, lifecycle_cfg.realize_max));
    if (queue.size() > realize_budget)
      queue.resize(realize_budget);
    // The chain-stage co-primaries: mm_set_parent over the incumbent and every chained
    // rival with a reportable anchor union, shadows left out, in the queue's order. Each
    // primary the queue does not hold is realized after it, outside its budget, and
    // classed a co-primary: it takes no part in the election, the recalibration, the
    // dp2 owner, the MAPQ's prices or the runner-ups.
    {
      std::vector<std::size_t> order;
      for (std::size_t i = 0; i < hypotheses.size(); ++i) {
        const rna::RnaRealizedHypothesis& h = hypotheses[i];
        if (i == 0 || (h.chained && h.bundle != nullptr && h.coverage_ok &&
                       !h.shadow))
          order.push_back(i);
      }
      std::stable_sort(order.begin(), order.end(),
                       [&hypotheses](std::size_t a, std::size_t b) {
                         return rna::rna_rival_queue_before(hypotheses[a],
                                                            hypotheses[b]);
                       });
      parent_chain_hits.reserve(order.size());
      for (const std::size_t i : order) {
        const rna::RnaRealizedHypothesis& h = hypotheses[i];
        rna::RnaParentHit hit;
        const rna::RnaQuerySpan span =
            rna::rna_forward_query_span(h.q_begin, h.q_end, h.reverse, read_len);
        hit.qs = span.begin;
        hit.qe = span.end;
        hit.score = h.chain_score;
        hit.cnt = h.chain_anchors;
        hit.rid = h.reference_id;
        hit.rs = h.summary->first_reference_begin;
        hit.re = h.summary->last_reference_end;
        hit.owner = static_cast<int>(i);
        hit.catalogue = h.catalogue_index;
        parent_chain_hits.push_back(hit);
      }
      const std::vector<rna::RnaParentHit> regions =
          parent_regions(parent_chain_hits, nullptr);
      for (std::size_t j = 0; j < regions.size(); ++j) {
        const std::size_t i = static_cast<std::size_t>(regions[j].owner);
        if (regions[j].parent != static_cast<int>(j) || i == 0 ||
            std::find(queue.begin(), queue.end(), i) != queue.end())
          continue;
        hypotheses[i].klass = rna::RnaRivalClass::CoPrimary;
        queue.push_back(i);
      }
    }
    if (!queue.empty()) {
    }

    for (const std::size_t index : queue) {
      rna::RnaRealizedHypothesis& candidate = hypotheses[index];
      rna::RnaSpliceRealizationRequest rival_request = request;
      // A rival contig that will not bind is skipped; the incumbent is already realized.
      if (!bind_realization_reference(rival_request, candidate.reference_id))
        continue;
      candidate.result =
          rna::realize_exact_anchor_path(rival_request, *candidate.bundle);
      candidate.realized = true;
      price_hypothesis(candidate, candidate.result);
    }

    // Rank recalibration, as minimap2's mm_update_dp_max. An intronless retrocopy of a
    // spliced gene can win the raw-DP election, because the spliced hit pays a long-gap
    // open at every junction. minimap2 reprices every aligned hit on an intron-free,
    // divergence-amplified scale and elects on that alone; here the election combines
    // the raw and recalibrated prices (rna_elect_primary_rank_recal), since the
    // recalibrated scale alone is blind to the aligned-mass differences the raw one
    // carries. dp_maximum stays the raw segment-0 price, the currency the MAPQ is
    // calibrated on, so when a flip leaves the demoted incumbent's raw dp2 above the
    // winner's dp1, branch C sees rho > 1 and lowers the MAPQ. Everything here fails
    // open: a gate that does not fire, a missing divergence estimate or missing
    // accounting leaves the read on raw prices.
    {
      // One hypothesis's primary-segment accounting. The incumbent's realization is in
      // `realized` (element 0's `result` is empty by the borrowed-storage contract). The
      // numeric CIGAR is the winning transcript hypothesis's segment_cigars entry, found
      // by orientation, rather than a re-parse of the clip-padded public string.
      const auto recal_accounting =
          [&](std::size_t index,
              rna::RnaRankRecalAccounting& accounting) -> bool {
        const rna::RnaSpliceRealizationResult& result =
            index == 0 ? realized : hypotheses[index].result;
        if (result.refused || result.segments.empty())
          return false;
        const rna::RnaSpliceHypothesisResult* winner = nullptr;
        for (const rna::RnaSpliceHypothesisResult& candidate :
             result.hypotheses)
          if (candidate.orientation == result.winning_hypothesis) {
            winner = &candidate;
            break;
          }
        if (winner == nullptr ||
            winner->segment_cigars.size() != winner->segments.size() ||
            winner->segment_cigars.empty())
          return false;
        const AlignResult& primary = result.segments.front();
        if (!primary.alignment_accounting_valid)
          return false;
        accounting.matches = primary.matches;
        accounting.mismatches = primary.mismatches;
        accounting.ambiguities = primary.ambiguities;
        // Borrowed, never copied: the vector lives in the realization result,
        // which outlives this block.
        accounting.cigar = winner->segment_cigars.front().data();
        accounting.cigar_len = winner->segment_cigars.front().size();
        return accounting.cigar_len != 0;
      };

      // The dp leader and runner-up, as the election orders them: argmax on the realized
      // maximum, ties to the lower catalogue index. Only the realized competitors count,
      // not the co-primaries, so this is minimap2's aligned-hit set.
      const auto off_scale = [](const rna::RnaRealizedHypothesis& h) {
        return !h.realized || !h.dp_maximum.has_value() || *h.dp_maximum <= 0 ||
               h.klass == rna::RnaRivalClass::CoPrimary;
      };
      std::size_t leader = hypotheses.size();
      for (std::size_t i = 0; i < hypotheses.size(); ++i) {
        const rna::RnaRealizedHypothesis& h = hypotheses[i];
        if (off_scale(h))
          continue;
        if (leader == hypotheses.size() ||
            *h.dp_maximum > *hypotheses[leader].dp_maximum ||
            (*h.dp_maximum == *hypotheses[leader].dp_maximum &&
             h.catalogue_index < hypotheses[leader].catalogue_index))
          leader = i;
      }
      int second_dp = 0;
      if (leader != hypotheses.size())
        for (std::size_t i = 0; i < hypotheses.size(); ++i) {
          if (i == leader)
            continue;
          const rna::RnaRealizedHypothesis& h = hypotheses[i];
          if (off_scale(h))
            continue;
          second_dp = std::max(second_dp, *h.dp_maximum);
        }
      // Fewer than two aligned hits is minimap2's max2 == 0: nothing to recalibrate
      // against.
      if (leader != hypotheses.size() && second_dp > 0) {
        const rna::RnaRealizedHypothesis& top = hypotheses[leader];
        rna::RnaRankRecalAccounting top_accounting;
        double identity = -1.0;
        if (rna::rna_rank_recal_gate(read_len, top.aq_begin, top.aq_end,
                                     *top.dp_maximum, second_dp) &&
            recal_accounting(leader, top_accounting))
          identity = rna::rna_event_identity(top_accounting);
        // A divergence the top hit cannot supply is no scale to reprice on.
        if (identity >= 0.0) {
          const int match_sc = rctx.opts.splice_controller.match;
          const double b2 = rna::rna_rank_recal_b2(
              identity, match_sc, rctx.opts.splice_controller.mismatch);
          // Repricing is atomic per read: a hypothesis left on its raw price beside
          // repriced rivals would be compared across currencies.
          std::vector<std::pair<std::size_t, int>> repriced;
          repriced.reserve(hypotheses.size());
          bool recal_complete = true;
          for (std::size_t i = 0; i < hypotheses.size(); ++i) {
            const rna::RnaRealizedHypothesis& h = hypotheses[i];
            if (off_scale(h))
              continue;
            rna::RnaRankRecalAccounting accounting;
            if (!recal_accounting(i, accounting)) {
              recal_complete = false; // fail open for the whole read
              break;
            }
            // A negative price floors to 0, as in minimap2.
            const int recal = rna::rna_recal_max_dp(accounting, b2, match_sc);
            repriced.emplace_back(i, recal < 0 ? 0 : recal);
          }
          if (recal_complete) {
            // dp_maximum keeps the raw price; the recalibrated one lands beside it.
            for (const auto& [index, price] : repriced)
              hypotheses[index].dp_recal = price;
            rank_recal_fired = true;
          }
        }
      }
    }

    // (5) Election and adoption. The raw election fails closed (a tie or a missing
    // maximum keeps the incumbent); a repriced read runs the combined rule with a raw
    // veto band (rival_lifecycle.h).
    lifecycle_elected = rank_recal_fired
                            ? rna::rna_elect_primary_rank_recal(hypotheses)
                            : rna::rna_elect_primary_hypothesis(hypotheses);
    // The shadow rule again, against the elected hypothesis: the winner's shadows may
    // include the demoted incumbent or a tie rival. The election is not touched.
    rna::RnaSegmentGeometry winner_geometry;
    if (shadow_geometry(hypotheses[lifecycle_elected], winner_geometry)) {
      for (std::size_t i = 0; i < hypotheses.size(); ++i) {
        if (i == lifecycle_elected)
          continue;
        rna::RnaSegmentGeometry other_geometry;
        hypotheses[i].shadow =
            shadow_geometry(hypotheses[i], other_geometry) &&
            rna::rna_same_place(winner_geometry, other_geometry, max_intron);
      }
    }
    lifecycle_dp2_owner =
        rna::rna_select_dp2_owner(hypotheses, lifecycle_elected);
    // The pieces, after the election: every realized region, repriced by
    // mm_update_dp_max and in mm_hit_sort's order with the elected family's regions
    // first (they own the query they cover), through mm_set_parent less the elected
    // family's shadows. A region outside the elected family stays only if it passes
    // minimap2's region filter, before and after the repricing (mm_filter_regs in
    // mm_align_skeleton). The election and the MAPQ do not read it.
    {
      const int elected_owner = static_cast<int>(lifecycle_elected);
      // Every segment of every realized hypothesis; a refused or empty realization has
      // no region.
      std::vector<rna::RnaParentHit> dp_hits;
      for (std::size_t i = 0; i < hypotheses.size(); ++i)
        if (i == 0 || hypotheses[i].realized)
          append_dp_hits(i, i == 0 ? realized : hypotheses[i].result, dp_hits);
      const auto filter_regions = [&]() {
        const auto passes = [&](const rna::RnaParentHit& hit) {
          if (hit.owner == elected_owner)
            return true;
          const std::size_t i = static_cast<std::size_t>(hit.owner);
          const rna::RnaSpliceHypothesisResult* winner =
              rna::rna_winning_hypothesis(i == 0 ? realized : hypotheses[i].result);
          rna::SpliceRegionEvidence evidence =
              winner->region_evidence[static_cast<std::size_t>(hit.segment)];
          evidence.dp_maximum = hit.dp_max;
          return rna::classify_splice_region(evidence, read_len,
                                             rctx.opts.splice_controller) ==
                 rna::SpliceRegionReject::None;
        };
        dp_hits.erase(
            std::stable_partition(dp_hits.begin(), dp_hits.end(), passes),
            dp_hits.end());
      };
      filter_regions();
      rna::rna_parent_update_dp_max(
          read_len, rna::kRnaRankRecalMinReadLen, dp_hits,
          static_cast<float>(rna::kRnaRankRecalFrac),
          rctx.opts.splice_controller.match,
          rctx.opts.splice_controller.mismatch);
      filter_regions();
      rna::rna_parent_hit_sort(dp_hits);
      std::stable_partition(dp_hits.begin(), dp_hits.end(),
                            [elected_owner](const rna::RnaParentHit& hit) {
                              return hit.owner == elected_owner;
                            });
      std::vector<char> left_out(hypotheses.size(), 0);
      for (std::size_t i = 0; i < hypotheses.size(); ++i)
        left_out[i] = i != lifecycle_elected && hypotheses[i].shadow ? 1 : 0;
      pieces_chain = parent_regions(parent_chain_hits, &left_out);
      decide_pieces(parent_regions(dp_hits, &left_out), pieces_chain);
      pieces_dp_input = std::move(dp_hits);
      pieces_left_out = std::move(left_out);
    }
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
    chain_evidence.f1 = committed_path->anchor_path_score;
    chain_evidence.cnt =
        static_cast<int>(committed_path->anchor_path_anchor_count);
    chain_evidence.sib_f2 = committed_siblings->sib_score;
    chain_evidence.winner_q_begin = committed_path_summary->first_query_begin;
    chain_evidence.winner_q_end = committed_path_summary->last_query_end;
    chain_evidence.rank1_score = committed_locus->rank_score;
    // The largest vote at any locus but the committed one: rank 2's without adoption,
    // rank 1's on any adoption (the catalogue is ordered).
    chain_evidence.rank2_score =
        committed_rank >= 2 ? catalogue_rank1_score : catalogue_rank2_score;
    chain_evidence.read_len = read_len;
    // As for map-only, in the committed frame: an adopted rival carries its own
    // post-repair strand.
    chain_evidence.winner_reverse = committed_reverse;
    chain_evidence.selected_query_covered_bases = commit_qcov;
    // Every catalogue rival is chained here, so a vote-only rival contributes x strength
    // 0 (see unchained_zero).
    chain_evidence.unchained_zero = true;
    if (committed_dp_max && *committed_dp_max > 0)
      chain_evidence.dp1 = static_cast<double>(*committed_dp_max);
    if (chain_evidence.dp1 > 0.0 && rival_dp_max && *rival_dp_max > 0)
      chain_evidence.dp2 = static_cast<double>(*rival_dp_max);
    if (!realized.segments.empty()) {
      const AlignResult& primary = realized.segments.front();
      chain_evidence.identity =
          primary.alignment_accounting_valid && primary.block_len > 0
              ? static_cast<double>(primary.matches) /
                    static_cast<double>(primary.block_len)
              : 1.0;
    }
    chain_evidence.match_sc = rctx.opts.splice_controller.match;
    // minimap2 map.c: sub_diff = 2a + b, per preset (4 splice / 6 splice:hq),
    // following any -A/-B override.
    chain_evidence.sub_diff = 2 * rctx.opts.splice_controller.match +
                              rctx.opts.splice_controller.mismatch;
    // The preset MAPQ calibration values ride the evidence.
    chain_evidence.qcov_tau = rctx.opts.mapq_qcov_tau;
    chain_evidence.disjoint_vote_damp = rctx.opts.mapq_disjoint_vote_damp;
    // The rivals are assembled from the hypothesis vector rather than chained again, in
    // catalogue order, with a demoted incumbent among them. owns_dp2 lands on the
    // lifecycle's dp2 owner, never on a co-primary.
    for (std::size_t i = 0; i < hypotheses.size(); ++i) {
      if (i == lifecycle_elected)
        continue;
      const rna::RnaRealizedHypothesis& other = hypotheses[i];
      rna::RnaChainMapqRival rival;
      rival.candidate = other.catalogue_index;
      rival.chained = other.chained;
      rival.chain_score = other.chain_score;
      rival.chain_anchors = other.chain_anchors;
      rival.q_begin = other.q_begin;
      rival.q_end = other.q_end;
      rival.rank_score = other.rank_score;
      rival.reverse = other.reverse;
      // The shadow verdict, re-read against the elected winner.
      rival.shadow = other.shadow;
      rival.owns_dp2 = chain_evidence.dp2 > 0.0 &&
                       lifecycle_dp2_owner.has_value() &&
                       *lifecycle_dp2_owner == i;
      // A realized rival is already bounded by rho = dp2/dp1 and stays out of x_floor;
      // a chained-only or censored rival enters it, and so does a co-primary, realized
      // for the pieces rather than as a price.
      rival.realized_dp = other.realized && other.dp_maximum.has_value() &&
                          *other.dp_maximum > 0 &&
                          other.klass != rna::RnaRivalClass::CoPrimary;
      chain_evidence.rivals.push_back(rival);
    }
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
                                  read_len, k, n_chr, leftover_pools)) {
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
            shadow_geometry(hypotheses[lifecycle_elected], elected_geometry);
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
          h.shadow = elected_placed && shadow_geometry(h, geometry) &&
                     rna::rna_same_place(elected_geometry, geometry, max_intron);
          hypotheses.push_back(std::move(h));
          rna::RnaRealizedHypothesis& candidate = hypotheses.back();
          rna::RnaSpliceRealizationRequest leftover_request = request;
          if (!candidate.coverage_ok ||
              !bind_realization_reference(leftover_request, candidate.reference_id))
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
            append_dp_hits(i, hypotheses[i].result, added);
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
          decide_pieces(parent_regions(regions, &left_out), pieces_chain);
        }
      }
    }

    locus_on_wire.assign(hypotheses.size(), 0);
    // The pieces, each one flat supplementary record: its region's segment, priced by
    // the chain MAPQ on its own region (dp_max, dp_max2, n_sub, identity),
    // its hypothesis's chain and vote, its children as its rivals and its region's
    // query span as the coverage brake's denominator.
    if (lifecycle_ran && !pieces.empty()) {
      const auto realization_of =
          [&](std::size_t i) -> const rna::RnaSpliceRealizationResult& {
        if (i == lifecycle_elected)
          return realized;
        if (i == 0)
          return demoted_incumbent;
        return hypotheses[i].result;
      };
      // Reserved once: `chimeric_families` borrows into `chimeric_results`.
      chimeric_results.reserve(pieces.size());
      chimeric_families.reserve(pieces.size());
      for (const PieceRegion& piece : pieces) {
        const rna::RnaParentHit& region = piece.region;
        const std::size_t owner = static_cast<std::size_t>(region.owner);
        const rna::RnaSpliceRealizationResult& result = realization_of(owner);
        const std::size_t s = static_cast<std::size_t>(region.segment);
        if (result.refused || s >= result.segments.size())
          continue;
        const AlignResult& segment = result.segments[s];
        const rna::RnaRealizedHypothesis& candidate = hypotheses[owner];
        rna::RnaChainMapqEvidence piece_evidence;
        piece_evidence.f1 = candidate.chain_score;
        piece_evidence.cnt = candidate.chain_anchors;
        piece_evidence.sib_f2 = candidate.siblings.sib_score;
        piece_evidence.winner_q_begin = candidate.q_begin;
        piece_evidence.winner_q_end = candidate.q_end;
        piece_evidence.rank1_score = candidate.rank_score;
        // The demoted vote: the largest vote among its children.
        for (const std::size_t c : piece.children)
          piece_evidence.rank2_score =
              std::max(piece_evidence.rank2_score, hypotheses[c].rank_score);
        piece_evidence.read_len = read_len;
        piece_evidence.winner_reverse = candidate.reverse;
        piece_evidence.selected_query_covered_bases =
            candidate.query_covered_bases;
        piece_evidence.qcov_bases = region.qe - region.qs;
        piece_evidence.unchained_zero = true;
        if (region.dp_max > 0)
          piece_evidence.dp1 = static_cast<double>(region.dp_max);
        if (piece_evidence.dp1 > 0.0 && region.dp_max2 > 0)
          piece_evidence.dp2 = static_cast<double>(region.dp_max2);
        piece_evidence.identity =
            segment.alignment_accounting_valid && segment.block_len > 0
                ? static_cast<double>(segment.matches) /
                      static_cast<double>(segment.block_len)
                : 1.0;
        piece_evidence.match_sc = chain_evidence.match_sc;
        piece_evidence.sub_diff = chain_evidence.sub_diff;
        piece_evidence.qcov_tau = chain_evidence.qcov_tau;
        piece_evidence.disjoint_vote_damp = chain_evidence.disjoint_vote_damp;
        piece_evidence.rivals.reserve(piece.children.size());
        for (const std::size_t c : piece.children) {
          const rna::RnaRealizedHypothesis& other = hypotheses[c];
          rna::RnaChainMapqRival rival;
          rival.candidate = other.catalogue_index;
          rival.chained = other.chained;
          rival.chain_score = other.chain_score;
          rival.chain_anchors = other.chain_anchors;
          rival.q_begin = other.q_begin;
          rival.q_end = other.q_end;
          rival.rank_score = other.rank_score;
          rival.reverse = other.reverse;
          rival.owns_dp2 = piece_evidence.dp2 > 0.0 && piece.dp2_owner == c;
          rival.realized_dp = other.realized && other.dp_maximum.has_value() &&
                              *other.dp_maximum > 0;
          rival.parented = true;
          piece_evidence.rivals.push_back(rival);
        }
        piece_evidence.parent_rule = true;
        piece_evidence.parent_n_sub = region.n_sub;
        rna::RnaChainMapqBreakdown piece_breakdown;
        const int piece_mapq =
            rna::rna_chain_mapq(piece_evidence, &piece_breakdown);
        rna::RnaSpliceRealizationResult piece_result;
        piece_result.segments.push_back(segment);
        chimeric_results.push_back(std::move(piece_result));
        rna::RnaChimericFamily family;
        family.result = &chimeric_results.back();
        family.mapq = piece_mapq;
        if (piece_evidence.f1 > 0) {
          family.chain_tags.chain_anchors = piece_evidence.cnt;
          family.chain_tags.chain_score = piece_evidence.f1;
          family.chain_tags.secondary_chain_score = piece_breakdown.f2;
        }
        chimeric_families.push_back(family);
        for (const std::size_t i : piece.place)
          locus_on_wire[i] = 1;
      }
    }
  }

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
    return finish(std::move(out));
  }

  // The runner-ups beside the committed record: one complete realized family per
  // surviving competitor, attached as alternative hypotheses. The realizations are
  // moved out; nothing reads `hypotheses[i].result` past the commit.
  {
    // One survivor: what to move, where it sits, and the keys the printed
    // order uses.
    struct RunnerUp {
      rna::RnaSpliceRealizationResult* result;
      rna::RnaSegmentGeometry place;
      int dp_maximum;
      int catalogue_index;
      // cm:i and s1:i from the hypothesis's own chain, as a chimeric family's.
      rna::RnaChainTags chain_tags;
    };
    std::vector<RunnerUp> runner_ups;
    const auto on_wire = [&](std::size_t j) {
      return j < locus_on_wire.size() && locus_on_wire[j] != 0;
    };
    // The realized place, in the geometry every duplicate test here
    // reads: the segments' reference hull on their own strand, with the
    // realized forward-frame query hull the lifecycle already computed.
    const auto collect = [&](rna::RnaSpliceRealizationResult& result,
                             const rna::RnaRealizedHypothesis& h) {
      if (result.refused || result.segments.empty() ||
          !h.dp_maximum.has_value() || *h.dp_maximum <= 0)
        return;
      rna::RnaSegmentGeometry place;
      place.reference_id = h.reference_id;
      place.reverse = result.segments.front().is_reverse;
      place.forward_span = rna::RnaQuerySpan{h.aq_begin, h.aq_end};
      for (std::size_t s = 0; s < result.segments.size(); ++s) {
        const AlignResult& segment = result.segments[s];
        const std::uint64_t begin = static_cast<std::uint64_t>(
            std::max<std::int64_t>(0, static_cast<std::int64_t>(segment.pos)));
        const std::uint64_t end =
            static_cast<std::uint64_t>(std::max<std::int64_t>(
                0, static_cast<std::int64_t>(segment.target_end)));
        place.reference_begin =
            s == 0 ? begin : std::min(place.reference_begin, begin);
        place.reference_end = s == 0 ? end : std::max(place.reference_end, end);
      }
      rna::RnaChainTags chain_tags;
      if (h.chain_score > 0) {
        chain_tags.chain_anchors = h.chain_anchors;
        chain_tags.chain_score = h.chain_score;
      }
      runner_ups.push_back(
          {&result, place, *h.dp_maximum, h.catalogue_index, chain_tags});
    };
    // Whether hypothesis i is an alternative placement of the bases the committed
    // record explains: not a shadow, and competing with the committed chain's span at
    // mask_level in the forward frame. `klass` was read against the incumbent, so after
    // an adoption it cannot answer this; rna_chain_mapq_admissible compares raw oriented
    // spans and would drop an opposite-strand alternative. With no lifecycle the span
    // stays censored and the vector is empty.
    rna::RnaQuerySpan committed_span;
    if (lifecycle_elected < hypotheses.size()) {
      const rna::RnaRealizedHypothesis& committed =
          hypotheses[lifecycle_elected];
      committed_span = rna::rna_forward_query_span(
          committed.q_begin, committed.q_end, committed.reverse, read_len);
    }
    const auto weighed = [&](std::size_t i) {
      const rna::RnaRealizedHypothesis& h = hypotheses[i];
      return !h.shadow &&
             rna::rna_query_spans_compete(
                 committed_span, rna::rna_forward_query_span(
                                     h.q_begin, h.q_end, h.reverse, read_len));
    };
    for (std::size_t i = 1; i < hypotheses.size(); ++i) {
      rna::RnaRealizedHypothesis& h = hypotheses[i];
      // Rivals of the committed record only: not the elected hypothesis (it is the
      // record), not a co-primary, not a shadow.
      if (i == lifecycle_elected || h.klass != rna::RnaRivalClass::Competitor ||
          !h.realized || on_wire(i) || !weighed(i))
        continue;
      collect(h.result, h);
    }
    // The demoted incumbent, from its own local; it was among the MAPQ's rivals, so the
    // same test applies.
    if (lifecycle_elected != 0 && !on_wire(0) && weighed(0))
      collect(demoted_incumbent, hypotheses[0]);
    if (!runner_ups.empty()) {
      // Richest realized maximum first, ties to the lower catalogue index
      // (rna_select_dp2_owner's rule), so the first alternative is the dp2 the MAPQ was
      // priced against.
      std::stable_sort(runner_ups.begin(), runner_ups.end(),
                       [](const RunnerUp& a, const RunnerUp& b) {
                         if (a.dp_maximum != b.dp_maximum)
                           return a.dp_maximum > b.dp_maximum;
                         return a.catalogue_index < b.catalogue_index;
                       });
      // Two catalogue windows can re-find one locus; the sort makes the kept copy the
      // better priced one.
      out.secondary.reserve(runner_ups.size());
      std::vector<std::size_t> kept;
      kept.reserve(runner_ups.size());
      for (std::size_t i = 0; i < runner_ups.size(); ++i) {
        bool duplicate = false;
        for (const std::size_t j : kept) {
          if (rna::rna_same_place(runner_ups[j].place, runner_ups[i].place,
                                  max_intron)) {
            duplicate = true;
            break;
          }
        }
        if (duplicate)
          continue;
        kept.push_back(i);
        rna::attach_runner_up_family(*runner_ups[i].result,
                                     runner_ups[i].chain_tags, out);
      }
    }
  }


  return finish(std::move(out));
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
