// Placement chaining. The whole-read vote catalogue is the only source of
// candidates; this stage harvests and chains anchors for the candidates the
// partition selects, maps each exact chain back to the 128 forward-query
// tiles and re-solves the partition. Both map-only projection and CIGAR
// realization consume the resulting family and anchor paths. No DP runs here.
#pragma once

#include "context.h"
#include "retained_seed_density.h"
#include "placement_family_adapter.h"
#include "alternative_hypothesis.h"
#include "../chaining/anchor.h"
#include "../chaining/colinear_params.h"
#include "../chaining/dense_chain.h"
#include "../chaining/partition.h"
#include "../index/seed.h"

#include <cstdint>
#include <cstddef>
#include <vector>

namespace fa::cpu::lr {

enum class DnaPlacementChainStatus : std::uint8_t {
  NotSelected,
  Sparse,
  Accepted,
  InvalidReference,
  MissingSeeds,
  NoIntervalHits,
  NoChain,
};

struct DnaPlacementCandidateChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  ::fa::cpu::voting::QueryTileMask sparse_support;
  ::fa::cpu::voting::QueryTileMask dense_support;
  std::vector<chaining::Anchor> primary;
  // Sibling paths from the same exact DP. Realization and the map-only
  // projection may use one when the primary path has no anchors in a selected
  // block (dna_selected_sibling_path).
  std::vector<std::vector<chaining::Anchor>> sibling_paths;
  // Parallel to sibling_paths: each sibling's chain score. rival_sibling is
  // the index of the sibling rival_chain_score came from
  // (ChainPartition::f2_index), or -1 when there is none or it was not
  // materialized. Whole-query passes only; read by dna_sibling_rival_chain,
  // and by the map-only projection for a record projected from a sibling.
  std::vector<int> sibling_scores;
  int rival_sibling = -1;
  DnaPlacementChainStatus status = DnaPlacementChainStatus::NotSelected;
  std::int64_t interval_hits = 0;
  int filtered_hits = 0;
  int chain_score = 0;
  int rival_chain_score = 0;
  int overlapping_rivals = 0;
  int chain_anchors = 0;
  // The screening pass's result, kept because the exact restore later
  // overwrites the chain fields of selected candidates.
  int screening_chain_score = 0;
  int screening_chain_anchors = 0;
  int screening_forward_query_begin = -1;
  int screening_forward_query_end = -1;
  int oriented_query_begin = -1;
  int oriented_query_end = -1;
  int forward_query_begin = -1;
  int forward_query_end = -1;
  int reference_begin = -1;
  int reference_end = -1;
  int initial_selected_tiles = 0;
  int final_selected_tiles = 0;
  bool exact = false;
  // Harvest counters. sparse_anchors: anchors built inside the harvest loop.
  // deferred_anchors: postings of the slices deferred as metadata, within
  // this candidate's window. rescued_anchors: anchors built by the deferred
  // restore. The pool is sparse + rescued. The whole-query pass defers every
  // slice, so there sparse_anchors is 0 and rescued_anchors is the pool; the
  // screening pass defers nothing.
  std::uint64_t sparse_anchors = 0;
  std::uint64_t deferred_anchors = 0;
  std::uint64_t rescued_anchors = 0;
  // Runs produced by the whole-query pass's collapse; 0 on the screening pass.
  std::int64_t dense_runs = 0;
  // The stage 1 trim (see kDnaSkipPoolBudget), run on the deferred slices'
  // metadata before any is restored. pool_trimmed: some tile was over its
  // share. trimmed_tiles: how many. trimmed_postings: postings never
  // appended. retained_slices / retained_anchors: the admissible keepers in
  // those tiles and their anchors, a subset of rescued_anchors. All 0 on the
  // screening pass.
  bool pool_trimmed = false;
  int trimmed_tiles = 0;
  std::uint64_t trimmed_postings = 0;
  int retained_slices = 0;
  std::uint64_t retained_anchors = 0;
};

// Residue recovery: query tiles the stable partition owns without support are
// searched again, from cached postings only. With no candidate or diagonal to
// start from, anchors carry their contig and are clustered by a coarse
// (contig, diagonal) histogram before the colinear DP scores them.

// Diagonal histogram bin width, in reference base pairs.
inline constexpr int kDnaResidueDiagonalWidth = 256;
// A histogram cluster is chained only with at least this many anchors, and the
// same floor gates the cached fine-seed supply of a candidate interval.
inline constexpr int kDnaResidueMinClusterAnchors = 3;
// Primary-path anchor floor for admission, an absolute bar.
inline constexpr int kDnaResidueMinChainAnchors = 10;
// The chain score floor in exact seed matches, independent of the anchor
// floor. N exact k-mer matches score at most N * min(seed_length, 255).
inline constexpr int kDnaResidueScoreFloorMatches = 3;
// Catalogue rivals chained over the whole query for the MAPQ, per ranking.
// mm_set_mapq2 reads one best rival and a count; two let the best rival by
// vote lose to the runner-up on chain score. Under
// kDnaChainMapqStudyBlockRivals each block also ranks its own.
inline constexpr int kDnaMapqRivalChains = 2;
// Only rivals with at least 1/kDnaMapqRivalVoteDenominator of the strongest
// owner's vote are chained; the rest still enter the MAPQ through their vote.
inline constexpr int kDnaMapqRivalVoteDenominator = 4;

// Stage 1 of the pool admission: a posting budget per whole-query pass, spent
// on the deferred slices' metadata before any anchor is built.
// kDnaSkipPoolBudget postings are shared among the chaining::kDenseAdmitTileBp
// oriented-query tiles that hold a slice, each getting max(kDnaSkipTileFloor,
// budget / occupied tiles). A tile within its share is restored whole. A tile
// over it keeps only slices whose key is under cigar_local_global_occ, rarest
// first by (global_count, count, read_pos, slot), while their total stays
// within the share (the rarest admissible slice always); a tile with no
// admissible slice keeps nothing, and the DP crosses it as a gap. The floor
// and the rarest admissible slice can take a long read's pool above the
// budget, so the bound grows with read length. The budget spans the whole
// query because a per-tile bar does not bound a sum of many ordinary tiles,
// and rarity is genome-wide because a key rare in the genome localizes the
// read while one merely rare in the window may not.
inline constexpr std::uint64_t kDnaSkipPoolBudget = 65536;
inline constexpr std::uint64_t kDnaSkipTileFloor = 128;

// Top histogram clusters considered per strand.
inline constexpr int kDnaResidueMaxClustersPerStrand = 2;
// Hard bounds. Exceeding the posting budget refuses the whole interval.
inline constexpr int kDnaResidueMaxIntervalPostings = 65536;
inline constexpr int kDnaResidueMaxIntervalsPerRead = 2;
inline constexpr int kDnaResidueMaxAdmissionsPerRead = 2;

bool dna_residue_posting_budget_allows(std::int64_t used,
                                       std::uint32_t next) noexcept;

// One posting expanded into an exact anchor, tagged with its contig.
struct DnaResidueAnchor {
  int contig = -1;
  chaining::Anchor anchor;
};

struct DnaResidueCluster {
  int contig = -1;
  // floor((r - q) / kDnaResidueDiagonalWidth) of the modal bin.
  std::int64_t diagonal_bin = 0;
  // Median (r - q) over the collected anchors: the cluster's peak diagonal.
  std::int64_t peak_diagonal = 0;
  std::vector<chaining::Anchor> anchors;
};

struct DnaResidueChainOutcome {
  bool admitted = false;
  int chain_score = 0;
  int chain_anchors = 0;
  // Query bp spanned by the primary path. The anchor-density bar divides by it.
  int chain_query_span = 0;
  std::vector<chaining::Anchor> primary;
};

// The admission bar of one residue chain. The default is the absolute bar; a
// post-commit pass may scale the floors to the query interval it opened.
struct DnaResidueAdmissionBar {
  int min_chain_anchors = kDnaResidueMinChainAnchors;
  // Negative uses context.opts.residue_min_anchor_density_per_100bp; zero
  // disables the density test.
  int min_anchor_density_per_100bp = -1;
};

// Bins anchors by (contig, floor((r - q) / kDnaResidueDiagonalWidth)), takes
// the `max_clusters` densest bins (never two adjacent on one contig) and
// collects each cluster from its bin and both neighbours, so a split at a bin
// edge loses nothing. Clusters under `min_anchors` are dropped. Ties on count
// break by the read-seeded hash of (contig, reverse, bin), as minimap2, then
// by (contig, bin).
std::vector<DnaResidueCluster>
dna_residue_diagonal_clusters(const std::vector<DnaResidueAnchor>& anchors,
                              int max_clusters, int min_anchors,
                              std::uint32_t tie_seed, bool reverse);

// Sorts and deduplicates the anchors, runs the dense chain with
// dna_candidate_chain_params and applies the admission bar to the expanded
// chain: at least bar.min_chain_anchors anchors, a score of at least
// kDnaResidueScoreFloorMatches * min(seed_length, 255), and the anchor
// density per 100 query bases.
DnaResidueChainOutcome dna_residue_chain_cluster(
    const DnaContext& context, std::vector<chaining::Anchor> anchors,
    int seed_length, int read_length,
    const DnaResidueAdmissionBar& bar = DnaResidueAdmissionBar{});

struct DnaResidueAdmission {
  DnaResidueCluster cluster;
  bool reverse = false;
  DnaResidueChainOutcome outcome;
};

// Helpers for post-commit terminal recovery. They read only the whole-query
// fine seeds and cached posting views.
int dna_residue_cached_supply(const DnaPlacementFamily& family,
                              const std::vector<RetainedSeedRef>& fine_forward,
                              const std::vector<RetainedSeedRef>& fine_reverse,
                              const ChainSeedLookupCache& lookup_cache,
                              std::uint32_t occurrence_cap, int query_begin_bp,
                              int query_end_bp);

bool dna_residue_collect_anchors(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<RetainedSeedRef>& fine_forward,
    const std::vector<RetainedSeedRef>& fine_reverse,
    const ChainSeedLookupCache& lookup_cache, std::uint32_t occurrence_cap,
    int query_begin_bp, int query_end_bp,
    std::vector<DnaResidueAnchor>& forward_anchors,
    std::vector<DnaResidueAnchor>& reverse_anchors);

// dna_residue_collect_anchors for a higher occurrence ceiling. Seeds are
// expanded rarest first (ties by key, read position, strand), so a frequent
// seed cannot use up the budget before a rare one at the true locus, and the
// walk stops at the posting budget instead of refusing the interval. An
// uncached key contributes nothing.
void dna_residue_collect_anchors_rarest_first(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<RetainedSeedRef>& fine_forward,
    const std::vector<RetainedSeedRef>& fine_reverse,
    const ChainSeedLookupCache& lookup_cache, std::uint32_t occurrence_cap,
    int query_begin_bp, int query_end_bp,
    std::vector<DnaResidueAnchor>& forward_anchors,
    std::vector<DnaResidueAnchor>& reverse_anchors);

bool dna_residue_best_admission(
    const DnaContext& context, const DnaPlacementFamily& family,
    const std::vector<DnaResidueAnchor>& forward,
    const std::vector<DnaResidueAnchor>& reverse, DnaResidueAdmission& best,
    const DnaResidueAdmissionBar& bar = DnaResidueAdmissionBar{},
    int max_clusters_per_strand = kDnaResidueMaxClustersPerStrand);

// The committed chain's anchor density: the primary-path anchors and
// forward-query span of the exact whole-query chain of the candidate that
// dominates the stable partition. A length-aware bar scales it by the interval
// length, so each read is held to its own anchor density. Zero when there is
// no such chain.
struct DnaResidueObservedDensity {
  int anchors = 0;
  int query_span = 0;
};
DnaResidueObservedDensity
dna_residue_observed_density(const DnaPlacementChainingResult& placement);

struct DnaResidueDetachedChain {
  DnaPlacementCandidate candidate;
  DnaPlacementCandidateChain chain;
};

DnaResidueDetachedChain
dna_residue_detached_chain(const DnaPlacementFamily& family,
                           const DnaResidueCluster& cluster, bool reverse,
                           const DnaResidueChainOutcome& outcome,
                           ::fa::cpu::voting::CandidateId id);

// A catalogue rival chained over the whole query for the MAPQ; a candidate's
// own chain partition only sees rivals inside its window. Evidence only: kept
// out of `candidates`, which is parallel to family.candidates and drives
// selection.
struct DnaRivalExactChain {
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  DnaPlacementCandidateChain chain;
  // Taken from a record that already had a whole-query chain instead of
  // chained again.
  bool reused = false;
};

struct DnaPlacementChainingResult {
  DnaPlacementFamily family;
  std::vector<DnaPlacementCandidateChain> candidates;
  DnaAlternativeSelection alternative;
  DnaPlacementCandidateChain alternative_exact;
  // The MAPQ's cross-locus rivals; empty when none is worth chaining.
  std::vector<DnaRivalExactChain> rival_exact;
  // Whole-query chains built for the MAPQ rivals, reused ones excluded.
  int mapq_rival_chains = 0;
  std::vector<RetainedSeedRef> residue_fine_forward;
  std::vector<RetainedSeedRef> residue_fine_reverse;
  std::int64_t initial_score = 0;
  int initial_blocks = 0;
  bool selection_changed = false;
  bool accepted = false;
  // Some stabilization under DnaTileOwnership::Span gave an accepted owner a
  // tile beyond its anchor tiles. Only then can a retry under AnchorTiles
  // produce a different family.
  bool span_widened = false;

  const DnaPlacementCandidateChain*
  find(::fa::cpu::voting::CandidateId candidate) const noexcept;

  // A candidate's whole-query exact chain: its MAPQ rival chain (any status),
  // else the retained alternative's exact restore, else its own record if it
  // was exact-restored. Null when there is none. Callers that need an
  // accepted chain check the status.
  const DnaPlacementCandidateChain*
  whole_query_chain(::fa::cpu::voting::CandidateId candidate) const noexcept;
};

// The winner's rival sibling as a chain record of its own, for
// build_dna_rival_placement: the sibling path as `primary`, exact and
// Accepted, with its own score, anchor count, dense_support and spans; other
// fields are left empty. Empty (NotSelected) when the winner has no rival
// sibling.
DnaPlacementCandidateChain
dna_sibling_rival_chain(const DnaPlacementCandidateChain& winner, bool reverse,
                        int read_length, int seed_length);

// The whole-query chain of the committed hypothesis. A promoted alternative
// uses its restricted rerun, not the stable run's record for that candidate.
const DnaPlacementCandidateChain* dna_committed_winner_chain(
    const DnaPlacementChainingResult* stable,
    const DnaPlacementChainingResult* promoted_alternative,
    ::fa::cpu::voting::CandidateId primary_candidate,
    bool primary_is_alternative) noexcept;

// Chain parameters for the candidate paths.
chaining::ColinearChainParams
dna_candidate_chain_params(const DnaContext& context, int seed_length,
                           int read_length);

// The same parameters for the dense run chain: the scoring is unchanged, and
// max_iter and max_skip have no counterpart. The diagonal-keyed search runs
// when diag_min_runs >= 0 && runs >= diag_min_runs, the linear scan
// otherwise; both are exact.
chaining::DenseChainParams
dna_dense_chain_params(const chaining::ColinearChainParams& params,
                       int seed_length, int diag_min_runs);

// Runs the restore/stabilize pipeline over {retained alternative, null}. The
// returned family uses solver id 0 and records the catalogue id in
// family.original_candidate_id.
DnaPlacementChainingResult
build_dna_alternative_placement(const DnaContext& context,
                                const DnaPlacementFamily& stable_family,
                                const DnaPlacementChainingResult& stable,
                                const std::vector<std::uint8_t>& forward_query,
                                const std::vector<std::uint8_t>& reverse_query);

// The same over {one catalogue candidate, null}, from a given exact
// whole-query chain.
DnaPlacementChainingResult
build_dna_rival_placement(const DnaContext& context,
                          const DnaPlacementFamily& stable_family,
                          ::fa::cpu::voting::CandidateId original,
                          const DnaPlacementCandidateChain& exact_chain,
                          const std::vector<std::uint8_t>& forward_query,
                          const std::vector<std::uint8_t>& reverse_query);

// Which query tiles an accepted whole-query chain owns when stabilization
// re-solves the partition (see stabilize_selected_family).
enum class DnaTileOwnership : std::uint8_t {
  // Its anchor tiles and every tile between its first and last anchor,
  // except those a rival at the same locus keeps. Used first.
  Span,
  // Its anchor tiles alone; used for the retry when a spanned family fails
  // to realize.
  AnchorTiles,
};

DnaPlacementChainingResult build_dna_placement_chains(
    const DnaContext& context, DnaPlacementFamily family,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query,
    const std::vector<ChainWindowRetainedSeed>* forward_seeds,
    const std::vector<ChainWindowRetainedSeed>* reverse_seeds,
    const std::vector<QuerySeed>* fine_forward_seeds,
    const std::vector<QuerySeed>* fine_reverse_seeds,
    ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots,
    DnaTileOwnership ownership);

} // namespace fa::cpu::lr
