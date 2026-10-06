// RnaBackend::map_read: vote capture, exon peaks, the coarse locus catalogue and the
// rank-1 exact anchor path, then either the map-only election and projection with the
// chain MAPQ or splice realization with the rival lifecycle and the query partition's
// second families.
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
#include "placement/coarse_chain.h"      // coarse locus-envelope selection
#include "query_partition.h" // the query partition below the MAPQ line

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
  // A query-partition nominee's window restriction, applied to each routed pool in its
  // own lane by the same routine as a single pool (restrict_harvest_to_query_window).
  const auto restrict_routed_pools = [&](rna::Rank1HarvestResult& harvest,
                                         const rna::RnaQuerySpan& window,
                                         bool nominated_reverse) {
    if (!harvest.routed)
      return;
    const auto narrow = [&](rna::CandidatePool& pool, bool pool_reverse) {
      if (pool.refused)
        return;
      rna::Rank1HarvestResult shim;
      shim.bounded_pool = std::move(pool);
      shim.regions = harvest.regions;
      shim.harvest_params = harvest.harvest_params;
      shim.work.postings_skipped_occurrence =
          harvest.work.postings_skipped_occurrence;
      rna::restrict_harvest_to_query_window(shim, window, pool_reverse,
                                            read_len, k);
      pool = std::move(shim.bounded_pool);
      if (shim.refused) {
        pool.refused = true;
        pool.refusal = shim.refusal;
      }
    };
    narrow(harvest.nominated_pool, nominated_reverse);
    narrow(harvest.opposite_pool, !nominated_reverse);
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
  // One catalogue entry's harvest and chain, shared by the rival scan and the query
  // partition's nominees so every window goes through the same routine (harvest, frame
  // election, orientation repair, summary). `query_window`, when non-null, restricts the
  // pool to the anchors explaining that forward-frame stretch of the read (nominees
  // only). `skeleton_peaks` is the peak vector `peak_indices` names (the read's own by
  // default). `committed_geometry` is the committed chain's geometry for the shadow
  // rule; null while the committed chain is not fixed.
  const auto chain_one_rival =
      [&](const rna::placement::CoarseLocus& rival_locus, int catalogue_index,
          const rna::RnaQuerySpan* query_window = nullptr,
          const std::vector<rna::placement::CoarseDiagonalPeak>*
              skeleton_peaks = nullptr,
          const rna::RnaSegmentGeometry* committed_geometry = nullptr) {
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
                  *rctx.ref.index, rival_seeds, rival_locus,
                  skeleton_peaks != nullptr ? *skeleton_peaks : peaks, read_len,
                  rival_chr_len, harvest_options, &postings);
              // The chain reads only the nominee's window.
              if (query_window != nullptr) {
                rna::restrict_harvest_to_query_window(
                    rival_harvest, *query_window, rival_locus.reverse, read_len,
                    k);
                restrict_routed_pools(rival_harvest, *query_window,
                                      rival_locus.reverse);
              }
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
              chain_one_rival(loci[i], static_cast<int>(i), nullptr, nullptr,
                              committed_geometry);
          evidence.rivals.push_back(entry.rival);
          rivals_out.push_back(std::move(entry));
        }
      };

  // Stage 4: rank 1 of the catalogue, one attempt per coarse pass: the budgeted
  // catalogue's rank 1 and, if that refuses, the second pass's. Everything downstream
  // (rivals, MAPQ, realization, lifecycle, partition) runs on whichever catalogue placed
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
    // query bases, the floor every partition family clears. First-pass placements are
    // exempt.
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
  // The query partition's chain storage (query_partition.h). A separate vector because
  // hypotheses borrow raw pointers into `production_rivals`, which must not grow;
  // reserved once before anything points into it.
  std::vector<ProductionRival> partition_rivals;
  // Parallel to `partition_rivals`: the best vote competing with each harvested window
  // (RnaExplainWindow::competing_vote; 0 for a contender), read by that family's MAPQ.
  std::vector<std::int64_t> partition_competing;
  // Where the partition's hypotheses begin in `hypotheses`. They are appended last and
  // contiguously, so this indexes `partition_competing`.
  std::size_t partition_first_hypothesis = 0;
  std::vector<rna::RnaRealizedHypothesis> hypotheses;
  // The demoted incumbent's realization, parked here on an adoption. Element 0's
  // `result` stays empty by the borrowed-storage contract (rival_lifecycle.h), which the
  // second-family pricing relies on. Read only by the runner-up attach.
  rna::RnaSpliceRealizationResult demoted_incumbent;
  // Per hypothesis: whether an emitted second family already put that locus on the
  // wire. Empty when no family was emitted.
  std::vector<char> locus_on_wire;
  std::size_t lifecycle_elected = 0;
  std::optional<std::size_t> lifecycle_dp2_owner;
  // Whether the rank recalibration repriced this read; selects the election rule.
  bool rank_recal_fired = false;

  // The rival realization lifecycle, in minimap2's order: chain the rivals, classify
  // each against the incumbent by forward-frame query geometry, admit by the retention
  // band (competitors) or the co-primary floor, realize a bounded number of competitors,
  // elect the primary on the realized DP maximum, and hand the demoted maxima to the
  // MAPQ as dp2.
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
          // candidate nor dp2/n_sub evidence. Admission makes it eligible for the
          // second-family block, realized after the primary's MAPQ is fixed.
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
      // maximum, ties to the lower catalogue index. Only admitted competitors are
      // realized, so this is minimap2's aligned-hit set.
      std::size_t leader = hypotheses.size();
      for (std::size_t i = 0; i < hypotheses.size(); ++i) {
        const rna::RnaRealizedHypothesis& h = hypotheses[i];
        if (!h.realized || !h.dp_maximum.has_value() || *h.dp_maximum <= 0)
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
          if (!h.realized || !h.dp_maximum.has_value() || *h.dp_maximum <= 0)
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
            if (!h.realized || !h.dp_maximum.has_value() || *h.dp_maximum <= 0)
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
  // The second-family store, outside the coverage guard because `chimeric_families`
  // borrows into it and both must outlive the commit. Filled only inside the guard: a
  // read that cannot be priced gets no second family.
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
      // a chained-only or censored rival enters it.
      rival.realized_dp = other.realized && other.dp_maximum.has_value() &&
                          *other.dp_maximum > 0;
      chain_evidence.rivals.push_back(rival);
    }
    cigar_chain_mapq = rna::rna_chain_mapq(chain_evidence, &primary_breakdown);
    chain_mapq_ptr = &cigar_chain_mapq;

    // The query partition, run only after the primary's MAPQ is fixed: above this line
    // every decision reads `hypotheses`, and a segment the catalogue never returned must
    // not enter the election, the dp2 selection or the primary's rivals. The primary is
    // its realized segments; the lifecycle's other chains, colinear staircases of the
    // residue peaks and a re-vote over the unexplained query are candidates, and one
    // exact partition of the read's tiles decides which segment explains which stretch
    // (query_partition.h). Each window it returns is harvested and chained by
    // chain_one_rival on that stretch, then realized, priced and emitted like a
    // co-primary.
    rna::RnaExplainResult explain;
    // Every hypothesis the partition nominated (a window's harvest, or a unit chained as
    // a contender), emitted strongest chain first.
    struct ExplainNominee {
      std::size_t hypothesis;
      rna::RnaSegmentKind kind;
      bool contender;
    };
    std::vector<ExplainNominee> explain_nominees;
    // A contender's locus, built from its unit's peaks and borrowed by its hypothesis;
    // reserved once before the first pointer is taken.
    std::vector<rna::placement::CoarseLocus> contender_loci;
    // The realized primary's geometry: the partition's pinned candidate, and
    // the emission's other side of the colinearity rule.
    rna::RnaSegmentGeometry primary_geometry;
    if (lifecycle_ran) {
      // The realized primary hull, wider than the chain's span since the DP extends past
      // its end anchors; the emission floor measures against it. Without one no second
      // family is emitted, so nothing is harvested.
      const std::optional<rna::RnaQuerySpan> realized_primary_hull =
          rna::rna_realized_query_span(realized);
      if (realized_primary_hull.has_value()) {
        const AlignResult& realized_primary = realized.segments.front();
        primary_geometry.reference_id = committed_chr;
        primary_geometry.reverse = realized_primary.is_reverse;
        primary_geometry.reference_begin =
            static_cast<std::uint64_t>(std::max<std::int64_t>(
                0, static_cast<std::int64_t>(realized_primary.pos)));
        primary_geometry.reference_end =
            static_cast<std::uint64_t>(std::max<std::int64_t>(
                0, static_cast<std::int64_t>(realized_primary.target_end)));
        primary_geometry.forward_span = *realized_primary_hull;
        // One forward-frame query interval per realized primary segment, so the mask is
        // what the primary aligned, not the gaps between.
        std::vector<rna::RnaQuerySpan>& primary_segments =
            worker_scratch.explain.segments;
        primary_segments.clear();
        if (const rna::RnaSpliceHypothesisResult* winner =
                rna::rna_winning_hypothesis(realized)) {
          for (const rna::SpliceRegionEvidence& evidence :
               winner->region_evidence)
            primary_segments.push_back(
                rna::RnaQuerySpan{evidence.query_begin, evidence.query_end});
        }
        // Every non-elected hypothesis the lifecycle chained, the demoted incumbent
        // included, is a candidate on its anchors.
        std::vector<rna::RnaExplainChained> chained;
        chained.reserve(hypotheses.size());
        for (std::size_t i = 0; i < hypotheses.size(); ++i) {
          if (i == lifecycle_elected)
            continue;
          const rna::RnaRealizedHypothesis& h = hypotheses[i];
          if (!h.chained || h.bundle == nullptr || h.locus == nullptr)
            continue;
          rna::RnaExplainChained entry;
          entry.hypothesis_index = static_cast<int>(i);
          entry.locus = h.locus;
          entry.bundle = h.bundle;
          entry.reverse = h.reverse;
          entry.rank_score = h.rank_score;
          entry.q_begin = h.q_begin;
          entry.q_end = h.q_end;
          chained.push_back(entry);
        }
        rna::rna_explain_collect(
            seed_ctx, *fwd_seeds, capture.views, peaks, copt, primary_geometry,
            committed_locus->rank_score, primary_segments, chained, read_len, k,
            min_support, max_intron, worker_scratch.explain, explain);
        // Storage for every harvest this stage makes (at most one contender per unit,
        // then the windows), reserved once: `hypotheses` holds raw pointers into these
        // vectors.
        std::size_t units = 0;
        for (const rna::RnaSegmentCandidate& candidate : explain.candidates)
          if (candidate.kind == rna::RnaSegmentKind::Unit ||
              candidate.kind == rna::RnaSegmentKind::Revote)
            ++units;
        const std::size_t harvest_capacity =
            units + static_cast<std::size_t>(rna::kRnaExplainMaxWindows);
        partition_rivals.reserve(harvest_capacity);
        partition_competing.reserve(harvest_capacity);
        contender_loci.reserve(units);
        partition_first_hypothesis = hypotheses.size();
        // One running catalogue index past the kept range, so the queue's total order
        // still resolves ties by index.
        int nominee_counter = 0;
        // One nominee window, harvested and appended through chain_one_rival like a
        // catalogue rival. Returns whether the co-primary floor admitted it.
        // `skeleton_peaks`: the vector the locus's peak_indices name when not the
        // read's own (the partition's unit peaks).
        const auto harvest_nominee =
            [&](rna::placement::CoarseLocus& nominated, int catalogue_index,
                std::int64_t competing, rna::RnaQuerySpan window,
                std::size_t twin_hypothesis,
                const std::vector<rna::placement::CoarseDiagonalPeak>*
                    skeleton_peaks) {
              const std::vector<rna::placement::CoarseDiagonalPeak>&
                  harvest_peaks =
                      skeleton_peaks != nullptr ? *skeleton_peaks : peaks;
              // The skeleton reads peak_indices as query-ordered anchor priority and
              // the coarse sweep leaves them in reference order, so sort them as the
              // selector does for kept loci. The colinear score and intron estimate
              // stay unset.
              std::stable_sort(
                  nominated.peak_indices.begin(), nominated.peak_indices.end(),
                  [&harvest_peaks](uint32_t a, uint32_t b) {
                    const auto& pa = harvest_peaks[a];
                    const auto& pb = harvest_peaks[b];
                    if (pa.oriented_query_begin != pb.oriented_query_begin)
                      return pa.oriented_query_begin < pb.oriented_query_begin;
                    if (pa.reference_begin != pb.reference_begin)
                      return pa.reference_begin < pb.reference_begin;
                    return a < b;
                  });
              // Restricted to the nominee's window: otherwise the harvest's
              // +/- max_intron pad hands the chain the committed locus's anchors and
              // the path re-finds the primary.
              partition_rivals.push_back(
                  chain_one_rival(nominated, catalogue_index, &window,
                                  skeleton_peaks, committed_geometry_ptr));
              partition_competing.push_back(competing);
              const ProductionRival& entry = partition_rivals.back();
              rna::RnaRealizedHypothesis h;
              h.catalogue_index = catalogue_index;
              h.reference_id = entry.reference_id;
              h.locus = entry.rival.chained || entry.repaired ? &entry.locus
                                                              : &nominated;
              // The vote evidence the partition computed for this window.
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
              h.partition_nominated = true;
              h.twin_hypothesis = twin_hypothesis;
              h.klass = rna::RnaRivalClass::CoPrimary;
              if (h.chained && h.bundle != nullptr) {
                h.coverage_ok = rna::anchor_union_query_coverage(
                    *h.bundle, h.query_covered_bases);
                // minimap2's chain floor only; the query the family adds is judged on
                // its realized hull by the emission floor.
                h.admitted = h.chain_score >= rna::kRnaCoprimaryMinChainScore &&
                             h.chain_anchors >= rna::kRnaCoprimaryMinAnchors;
              }
              const bool admitted = h.admitted;
              hypotheses.push_back(std::move(h));
              return admitted;
            };
        int contenders_seen = 0;
        // Contenders: a unit whose tiles compete with a chained hypothesis's is chained
        // before the solve, on its own span, so a contested interval is decided chain
        // against chain. At most kRnaExplainMaxContenders units, strongest vote first.
        std::vector<std::size_t>& contender_order =
            worker_scratch.explain.offer_order;
        contender_order.clear();
        for (std::size_t u = 0; u < explain.candidates.size(); ++u) {
          const rna::RnaSegmentCandidate& unit = explain.candidates[u];
          if (unit.kind != rna::RnaSegmentKind::Unit &&
              unit.kind != rna::RnaSegmentKind::Revote)
            continue;
          const rna::RnaQuerySpan unit_span = unit.geometry.forward_span;
          if (unit_span.begin < 0 || unit_span.end <= unit_span.begin)
            continue;
          if (rna::rna_explain_contender(explain, u) < 0)
            continue;
          contender_order.push_back(u);
        }
        std::stable_sort(contender_order.begin(), contender_order.end(),
                         [&explain](std::size_t a, std::size_t b) {
                           const auto& ca = explain.candidates[a];
                           const auto& cb = explain.candidates[b];
                           if (ca.vote_evidence != cb.vote_evidence)
                             return ca.vote_evidence > cb.vote_evidence;
                           if (ca.chain_evidence != cb.chain_evidence)
                             return ca.chain_evidence > cb.chain_evidence;
                           return a < b;
                         });
        // One contender per place: a unit at a chained hypothesis's place is that
        // hypothesis's own evidence, and a second unit at a chosen contender's place is
        // the same segment again.
        {
          std::size_t chosen = 0;
          for (std::size_t i = 0; i < contender_order.size(); ++i) {
            if (chosen >=
                static_cast<std::size_t>(rna::kRnaExplainMaxContenders))
              break;
            const rna::RnaSegmentGeometry& geometry =
                explain.candidates[contender_order[i]].geometry;
            bool same_place = false;
            for (std::size_t c = 1;
                 c < explain.candidates.size() && !same_place; ++c)
              if (explain.candidates[c].kind == rna::RnaSegmentKind::Chained)
                same_place = rna::rna_same_place(
                    geometry, explain.candidates[c].geometry, max_intron);
            for (std::size_t j = 0; j < chosen && !same_place; ++j)
              same_place = rna::rna_same_place(
                  geometry, explain.candidates[contender_order[j]].geometry,
                  max_intron);
            if (same_place)
              continue;
            contender_order[chosen++] = contender_order[i];
          }
          contender_order.resize(chosen);
        }
        for (const std::size_t u : contender_order) {
          const rna::RnaSegmentCandidate& unit = explain.candidates[u];
          const rna::RnaQuerySpan unit_span = unit.geometry.forward_span;
          contender_loci.push_back(
              rna::rna_explain_unit_locus(explain, unit.staircase_index));
          const int catalogue_index =
              static_cast<int>(loci.size()) + nominee_counter++;
          ++contenders_seen;
          // Chained over the whole unexplained piece the unit lies in, like
          // a window's harvest: the chain, not the unit's peaks, says how
          // far the segment reaches.
          const rna::RnaQuerySpan contender_piece =
              rna::rna_explain_piece(worker_scratch.explain, unit_span);
          harvest_nominee(
              contender_loci.back(), catalogue_index, 0,
              contender_piece.begin >= 0 ? contender_piece : unit_span,
              static_cast<std::size_t>(-1), &explain.unit_peaks);
          const std::size_t chained_index = hypotheses.size() - 1;
          const rna::RnaRealizedHypothesis& h = hypotheses[chained_index];
          if (h.chained && h.bundle != nullptr) {
            rna::rna_explain_candidate_chained(
                explain, u, static_cast<int>(chained_index), *h.bundle,
                h.reverse, read_len, k, explain.candidates.front().support);
            explain_nominees.push_back({chained_index, unit.kind, true});
          }
        }
        rna::rna_explain_solve(copt, primary_geometry, primary_segments,
                               read_len, k, max_intron, worker_scratch.explain,
                               explain);
        // The windows, in query order. A window owned by a contender reuses the
        // contender's chain; one owned by a chained hypothesis is re-chained inside the
        // window with that hypothesis as its twin (the same placement, never its rival).
        for (std::size_t w = 0; w < explain.windows.size(); ++w) {
          rna::RnaExplainWindow& window = explain.windows[w];
          const int catalogue_index =
              static_cast<int>(loci.size()) + nominee_counter++;
          const bool contender = window.kind != rna::RnaSegmentKind::Chained &&
                                 window.hypothesis_index >= 0;
          if (contender) {
            // The owner was chained on its own span before the solve; that chain is
            // the family, and its competing vote is the window's.
            const std::size_t owner =
                static_cast<std::size_t>(window.hypothesis_index);
            if (owner >= partition_first_hypothesis &&
                owner - partition_first_hypothesis < partition_competing.size())
              partition_competing[owner - partition_first_hypothesis] =
                  window.competing_vote;
            --nominee_counter;
            continue;
          }
          const std::size_t twin =
              window.hypothesis_index >= 0
                  ? static_cast<std::size_t>(window.hypothesis_index)
                  : static_cast<std::size_t>(-1);
          harvest_nominee(window.locus, catalogue_index, window.competing_vote,
                          window.piece, twin,
                          window.unit_peaks ? &explain.unit_peaks : nullptr);
          const std::size_t harvested = hypotheses.size() - 1;
          // A contender chained over this window's stretch also answers to the window's
          // competing vote; the emission takes the strongest chain either way.
          for (const ExplainNominee& nominee : explain_nominees) {
            if (!nominee.contender)
              continue;
            const rna::RnaRealizedHypothesis& ch =
                hypotheses[nominee.hypothesis];
            const rna::RnaQuerySpan contender_span =
                rna::rna_forward_query_span(ch.q_begin, ch.q_end, ch.reverse,
                                            read_len);
            if (!rna::rna_query_spans_compete(window.span, contender_span))
              continue;
            if (nominee.hypothesis >= partition_first_hypothesis &&
                nominee.hypothesis - partition_first_hypothesis <
                    partition_competing.size()) {
              std::int64_t& competing =
                  partition_competing[nominee.hypothesis -
                                      partition_first_hypothesis];
              competing = std::max(competing, window.competing_vote);
            }
          }
          explain_nominees.push_back({harvested, window.kind, false});
        }
      }
    }
    // The partition has finished appending, so the hypothesis vector has its final
    // size.
    locus_on_wire.assign(hypotheses.size(), 0);

    // The second families: the partition's nominees, realized and emitted. Everything
    // here is downstream of the primary's price and writes nothing back onto
    // `hypotheses`. A family is emitted when its realized hull explains at least
    // kRnaChimeraMinQueryBases query bases no accepted segment explains, at the
    // alignment quality the primary must meet. When no window survives, a contender
    // (chained before the solve) is emitted only if it won a partition block.
    if (lifecycle_ran && !explain_nominees.empty()) {
      const std::optional<rna::RnaQuerySpan> primary_hull =
          rna::rna_realized_query_span(realized);
      if (primary_hull.has_value()) {
        constexpr std::size_t kMaxFamilies =
            static_cast<std::size_t>(rna::kRnaChimeraMaxFamilies);
        // The accepted segments, primary first: their hulls, segment-0 prices and
        // owning hypotheses. An accepted family is a co-family member of every later
        // one, never its rival.
        std::vector<rna::RnaQuerySpan> accepted_hulls;
        std::vector<int> accepted_dp0;
        std::vector<std::size_t> accepted_indices;
        accepted_hulls.reserve(kMaxFamilies + 1);
        accepted_dp0.reserve(kMaxFamilies);
        accepted_indices.reserve(kMaxFamilies);
        accepted_hulls.push_back(*primary_hull);
        // Reserved once: `chimeric_families` holds raw pointers into this vector.
        chimeric_results.reserve(kMaxFamilies);
        chimeric_families.reserve(kMaxFamilies);
        // Strongest chain first, minimap2's primary order. A chain competing at
        // mask_level with an accepted hull is refused before it is realized.
        std::vector<std::size_t> emission_order;
        emission_order.reserve(explain_nominees.size());
        for (std::size_t q = 0; q < explain_nominees.size(); ++q)
          emission_order.push_back(q);
        std::stable_sort(emission_order.begin(), emission_order.end(),
                         [&](std::size_t a, std::size_t b) {
                           const rna::RnaRealizedHypothesis& ha =
                               hypotheses[explain_nominees[a].hypothesis];
                           const rna::RnaRealizedHypothesis& hb =
                               hypotheses[explain_nominees[b].hypothesis];
                           if (ha.chain_score != hb.chain_score)
                             return ha.chain_score > hb.chain_score;
                           if (ha.chain_anchors != hb.chain_anchors)
                             return ha.chain_anchors > hb.chain_anchors;
                           // A tie between a window's own harvest and a
                           // contender goes to the harvest: the solve nominated
                           // it on the tiles it explains.
                           if (explain_nominees[a].contender !=
                               explain_nominees[b].contender)
                             return !explain_nominees[a].contender;
                           return a < b;
                         });
        for (const std::size_t q : emission_order) {
          if (chimeric_families.size() >= kMaxFamilies)
            break;
          const std::size_t index = explain_nominees[q].hypothesis;
          const rna::RnaRealizedHypothesis& candidate = hypotheses[index];
          if (!candidate.admitted || candidate.bundle == nullptr)
            continue;
          // With no window left, a contender the solve gave no block was not chosen by
          // the partition: it stays a hypothesis, a rival and runner-up of any other
          // family. One that won a block stays eligible though the window filters
          // dropped its window.
          if (explain_nominees[q].contender && explain.windows.empty() &&
              std::none_of(explain.candidates.begin(), explain.candidates.end(),
                           [index](const rna::RnaSegmentCandidate& c) {
                             return c.owns_block &&
                                    c.hypothesis_index == static_cast<int>(index);
                           }))
            continue;
          const rna::RnaQuerySpan candidate_span = rna::rna_forward_query_span(
              candidate.q_begin, candidate.q_end, candidate.reverse, read_len);
          {
            bool competes = false;
            for (const rna::RnaQuerySpan& accepted : accepted_hulls)
              if (rna::rna_query_spans_compete(candidate_span, accepted))
                competes = true;
            if (competes) {
              continue;
            }
          }
          rna::RnaSpliceRealizationRequest chimeric_request = request;
          // A contig that will not bind is a family skipped; the primary is unaffected.
          if (!bind_realization_reference(chimeric_request,
                                          candidate.reference_id)) {
            continue;
          }
          rna::RnaSpliceRealizationResult chimeric_result;
          if (candidate.realized) {
            chimeric_result = candidate.result; // realized above, as a price
          } else {
            chimeric_result = rna::realize_exact_anchor_path(
                chimeric_request, *candidate.bundle);
          }
          if (chimeric_result.refused || chimeric_result.segments.empty()) {
            continue;
          }
          const std::optional<int> chimeric_dp0 =
              rna::rna_segment0_dp_maximum(chimeric_result);
          const std::optional<rna::RnaQuerySpan> chimeric_hull =
              rna::rna_realized_query_span(chimeric_result);
          const rna::RnaQuerySpan hull =
              chimeric_hull.value_or(rna::RnaQuerySpan{});
          if (!rna::rna_chimera_family_admitted(
                  chimeric_dp0, hull,
                  rctx.opts.splice_controller.minimum_dp_maximum)) {
            continue;
          }
          // The floor: the query the family explains that no accepted segment does. No
          // accepted segment may compete with it at mask_level either: then it is a
          // second description of those bases, a rival, not a second primary.
          std::int64_t exclusive = hull.end - hull.begin;
          bool competes = false;
          for (const rna::RnaQuerySpan& accepted : accepted_hulls) {
            exclusive -= std::max(0, std::min(hull.end, accepted.end) -
                                         std::max(hull.begin, accepted.begin));
            if (rna::rna_query_spans_compete(hull, accepted))
              competes = true;
          }
          if (competes || exclusive < rna::kRnaChimeraMinQueryBases) {
            continue;
          }

          // This family's MAPQ, from the same formula with the family as the committed
          // hypothesis: its chain supplies f1, cnt, sib_f2 and the winner span, its
          // realization dp1 and identity, its hull the brake's denominator, and every
          // other competing hypothesis becomes a rival. The primary's evidence is only
          // read. Another hypothesis competes on its realized hull when it has one (the
          // elected primary's is the primary hull), else on its chain span. Excluded
          // throughout: the family itself, its twin and every accepted family.
          const auto competition_span =
              [&](std::size_t j) -> rna::RnaQuerySpan {
            if (j == lifecycle_elected)
              return *primary_hull;
            const rna::RnaRealizedHypothesis& other = hypotheses[j];
            if (other.aq_begin >= 0 && other.aq_end > other.aq_begin)
              return rna::RnaQuerySpan{other.aq_begin, other.aq_end};
            return rna::rna_forward_query_span(other.q_begin, other.q_end,
                                               other.reverse, read_len);
          };
          // The family's place: its realized segments' reference hull on its strand. A
          // hypothesis at the same place (rna_same_place) is a second description of
          // this segment, not a rival.
          const auto reference_hull =
              [](const std::vector<AlignResult>& segments,
                 rna::RnaSegmentGeometry& geometry) {
                bool first_segment = true;
                for (const AlignResult& segment : segments) {
                  const std::uint64_t begin =
                      static_cast<std::uint64_t>(std::max<std::int64_t>(
                          0, static_cast<std::int64_t>(segment.pos)));
                  const std::uint64_t end =
                      static_cast<std::uint64_t>(std::max<std::int64_t>(
                          0, static_cast<std::int64_t>(segment.target_end)));
                  if (first_segment) {
                    geometry.reference_begin = begin;
                    geometry.reference_end = end;
                    first_segment = false;
                  } else {
                    geometry.reference_begin =
                        std::min(geometry.reference_begin, begin);
                    geometry.reference_end =
                        std::max(geometry.reference_end, end);
                  }
                }
              };
          rna::RnaSegmentGeometry family_geometry;
          family_geometry.reference_id = candidate.reference_id;
          family_geometry.reverse = chimeric_result.segments.front().is_reverse;
          family_geometry.forward_span = hull;
          reference_hull(chimeric_result.segments, family_geometry);
          // A family that continues the primary (same contig and strand, query order
          // agreeing with reference order at both ends, reference spans overlapping or
          // within one legal intron) is an exon the realizer skipped or did not reach,
          // not a second place. A family that steps back on the reference (a duplicated
          // exon, a back-splice) stays.
          if (rna::rna_segment_continuation(primary_geometry, family_geometry,
                                            max_intron)) {
            continue;
          }
          const auto hypothesis_geometry = [&](std::size_t j) {
            const rna::RnaRealizedHypothesis& other = hypotheses[j];
            rna::RnaSegmentGeometry geometry;
            geometry.reference_id = other.reference_id;
            geometry.reverse = other.reverse;
            if (other.realized && !other.result.segments.empty()) {
              geometry.reverse = other.result.segments.front().is_reverse;
              reference_hull(other.result.segments, geometry);
            } else if (other.chained && other.summary != nullptr &&
                       other.summary->last_reference_end >
                           other.summary->first_reference_begin) {
              geometry.reference_begin = static_cast<std::uint64_t>(
                  std::max(0, other.summary->first_reference_begin));
              geometry.reference_end = static_cast<std::uint64_t>(
                  std::max(0, other.summary->last_reference_end));
            } else if (other.locus != nullptr) {
              geometry.reference_begin = other.locus->reference_begin;
              geometry.reference_end = other.locus->reference_end;
            }
            geometry.forward_span = competition_span(j);
            return geometry;
          };
          const auto excluded = [&](std::size_t j) {
            if (j == index || j == candidate.twin_hypothesis)
              return true;
            for (const std::size_t accepted : accepted_indices)
              if (accepted == j)
                return true;
            return j != lifecycle_elected &&
                   rna::rna_same_place(family_geometry, hypothesis_geometry(j),
                                       max_intron);
          };
          // The family's runner-up, realized for the price only: the strongest admitted
          // nominee whose chain competes with this hull at another place, as the
          // lifecycle does for the primary (dp2 is realized evidence, never a chain
          // score). One per family.
          for (const std::size_t r : emission_order) {
            const std::size_t j = explain_nominees[r].hypothesis;
            if (excluded(j))
              continue;
            rna::RnaRealizedHypothesis& other = hypotheses[j];
            if (!other.admitted || other.bundle == nullptr || other.realized)
              continue;
            if (!rna::rna_query_spans_compete(
                    rna::rna_forward_query_span(other.q_begin, other.q_end,
                                                other.reverse, read_len),
                    hull))
              continue;
            rna::RnaSpliceRealizationRequest runner_request = request;
            if (bind_realization_reference(runner_request,
                                           other.reference_id)) {
              other.result = rna::realize_exact_anchor_path(
                  runner_request, *other.bundle);
              other.realized = true;
              price_hypothesis(other, other.result);
            }
            break;
          }
          // dp2: the largest positive realized DP maximum among the competing
          // hypotheses and competing accepted families (priced by segment 0). Ascending
          // order with a strict `>` gives ties to the lower index, as
          // rna_select_dp2_owner does. 0 means nothing competes (branch B).
          int family_dp2 = 0;
          std::optional<std::size_t> family_dp2_owner;
          for (std::size_t j = 0; j < hypotheses.size(); ++j) {
            if (excluded(j))
              continue;
            const rna::RnaRealizedHypothesis& other = hypotheses[j];
            if (!other.realized || !other.dp_maximum.has_value() ||
                *other.dp_maximum <= 0)
              continue;
            // A shadow of the primary would count the primary's own maximum again.
            if (other.shadow)
              continue;
            if (!rna::rna_query_spans_compete(candidate_span,
                                              competition_span(j)))
              continue;
            if (*other.dp_maximum > family_dp2) {
              family_dp2 = *other.dp_maximum;
              family_dp2_owner = j;
            }
          }
          for (std::size_t a = 0; a < accepted_dp0.size(); ++a) {
            if (!rna::rna_query_spans_compete(candidate_span,
                                              accepted_hulls[a + 1]))
              continue;
            if (accepted_dp0[a] > family_dp2) {
              family_dp2 = accepted_dp0[a];
              // Owned by no hypothesis, so no rival carries owns_dp2; the ratio
              // still damps x.
              family_dp2_owner.reset();
            }
          }
          rna::RnaChainMapqEvidence family_evidence;
          family_evidence.f1 = candidate.chain_score;
          family_evidence.cnt = candidate.chain_anchors;
          family_evidence.sib_f2 = candidate.siblings.sib_score;
          family_evidence.winner_q_begin = candidate.q_begin;
          family_evidence.winner_q_end = candidate.q_end;
          family_evidence.rank1_score = candidate.rank_score;
          // The demoted vote for this family: the largest vote competing with its
          // window, among the partition's own candidates (competing_vote) and the
          // competing hypotheses. 0 when nothing competes.
          std::int64_t family_rank2 = 0;
          if (index >= partition_first_hypothesis &&
              index - partition_first_hypothesis < partition_competing.size())
            family_rank2 =
                partition_competing[index - partition_first_hypothesis];
          for (std::size_t j = 0; j < hypotheses.size(); ++j) {
            if (excluded(j))
              continue;
            if (rna::rna_query_spans_compete(candidate_span,
                                             competition_span(j)))
              family_rank2 = std::max(family_rank2, hypotheses[j].rank_score);
          }
          family_evidence.rank2_score = family_rank2;
          family_evidence.read_len = read_len;
          family_evidence.winner_reverse = candidate.reverse;
          family_evidence.selected_query_covered_bases =
              candidate.query_covered_bases;
          // The brake divides by the family's hull, not the read length: a family
          // explains a piece of the read by construction.
          family_evidence.qcov_bases = hull.end - hull.begin;
          // Every catalogue rival is chained, as for the primary.
          family_evidence.unchained_zero = true;
          family_evidence.dp1 = static_cast<double>(*chimeric_dp0);
          if (family_dp2 > 0)
            family_evidence.dp2 = static_cast<double>(family_dp2);
          {
            const AlignResult& family_primary =
                chimeric_result.segments.front();
            family_evidence.identity =
                family_primary.alignment_accounting_valid &&
                        family_primary.block_len > 0
                    ? static_cast<double>(family_primary.matches) /
                          static_cast<double>(family_primary.block_len)
                    : 1.0;
          }
          // Scoring and calibration constants are the read's, taken from the
          // primary's evidence.
          family_evidence.match_sc = chain_evidence.match_sc;
          family_evidence.sub_diff = chain_evidence.sub_diff;
          family_evidence.qcov_tau = chain_evidence.qcov_tau;
          family_evidence.disjoint_vote_damp =
              chain_evidence.disjoint_vote_damp;
          // The rivals: other hypotheses competing with the family for its query bases
          // in the forward frame, the test that chose its dp2 and demoted vote. The
          // formula's raw-span admission still applies.
          family_evidence.rivals.reserve(hypotheses.size());
          for (std::size_t j = 0; j < hypotheses.size(); ++j) {
            if (excluded(j))
              continue;
            const rna::RnaRealizedHypothesis& other = hypotheses[j];
            if (!rna::rna_query_spans_compete(candidate_span,
                                              competition_span(j)))
              continue;
            rna::RnaChainMapqRival rival;
            rival.candidate = other.catalogue_index;
            rival.chained = other.chained;
            rival.chain_score = other.chain_score;
            rival.chain_anchors = other.chain_anchors;
            rival.q_begin = other.q_begin;
            rival.q_end = other.q_end;
            rival.rank_score = other.rank_score;
            rival.reverse = other.reverse;
            rival.shadow = other.shadow;
            rival.owns_dp2 = family_evidence.dp2 > 0.0 &&
                             family_dp2_owner.has_value() &&
                             *family_dp2_owner == j;
            rival.realized_dp = other.realized &&
                                other.dp_maximum.has_value() &&
                                *other.dp_maximum > 0;
            family_evidence.rivals.push_back(rival);
          }
          // The same pure out-parameter as the primary's.
          rna::RnaChainMapqBreakdown family_breakdown;
          const int family_mapq =
              rna::rna_chain_mapq(family_evidence, &family_breakdown);

          chimeric_results.push_back(std::move(chimeric_result));
          rna::RnaChimericFamily family;
          family.result = &chimeric_results.back();
          family.mapq = family_mapq;
          // cm:i and s1:i from this family's own chain, s2:i from its own rival set,
          // never the primary's.
          if (family_evidence.f1 > 0) {
            family.chain_tags.chain_anchors = family_evidence.cnt;
            family.chain_tags.chain_score = family_evidence.f1;
            family.chain_tags.secondary_chain_score = family_breakdown.f2;
          }
          chimeric_families.push_back(family);
          accepted_hulls.push_back(hull);
          accepted_dp0.push_back(*chimeric_dp0);
          accepted_indices.push_back(index);
          // Record which loci this family put on the wire (its hypothesis, its twin,
          // accepted families and everything at its place), so the runner-up attach
          // does not print one again. For the incumbent the geometry is its chain span,
          // not its realized hull.
          for (std::size_t j = 0; j < locus_on_wire.size(); ++j)
            if (excluded(j))
              locus_on_wire[j] = 1;
        }
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
      // record), not a partition nominee (a co-primary realized for a price), not a
      // shadow.
      if (i == lifecycle_elected || h.klass != rna::RnaRivalClass::Competitor ||
          h.partition_nominated || !h.realized || on_wire(i) || !weighed(i))
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
