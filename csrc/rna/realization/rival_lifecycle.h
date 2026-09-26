// The RNA rival realization lifecycle, in minimap2's shape:
//   1. Classify every retained chain against the primary by query geometry (the
//      mask_level test of mm_set_parent): an overlapping chain competes for the same
//      bases, a disjoint one is a co-primary, a different part of the read.
//   2. Admit competitors by chain score (pri_ratio / min_diff), co-primaries by their
//      own chain-quality floor.
//   3. Realize the admitted competitors, bounded.
//   4. Elect the primary on the realized DP maximum (mm_hit_sort on dp_max).
//   5. Feed the demoted maxima to the MAPQ as dp2 / n_sub evidence.
// This header is the pure policy; the backend owns storage, harvests, realizations and
// the adoption.
//
// The election fails closed: a competitor takes the primary slot only when its own and
// the incumbent's realized DP maxima are both present and positive and its own is
// strictly greater. Exact ties (queue keys, both elections, the dp2 owner, and in the
// backend the recalibration leader and runner-up order) go to the lower catalogue index,
// the incumbent's (0) first. The catalogue is ordered by vote, then by minimap2's
// read-seeded hash, so loci tied on votes too are decided by the hash.
//
// A co-primary never takes part in the election and never supplies dp2 or n_sub. An
// admitted co-primary may be realized after the primary's MAPQ is fixed and emitted as a
// second family with its own chain MAPQ; that writes nothing back onto the hypotheses.
#pragma once

#include "splice_realizer.h"
#include "rival_lifecycle_config.h"                // RnaRivalLifecycleConfig
#include "../anchoring/exact_anchor_path.h"        // ExactAnchorPath(+Summary)
#include "../anchoring/query_span.h"               // RnaQuerySpan
#include "../anchoring/select_exact_anchor_path.h" // SiblingChainStats

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Borrowed by pointer only; declared here to avoid a dependency cycle with the coarse
// catalogue.
namespace placement {
struct CoarseLocus;
}

// What a hypothesis IS relative to the rank-1 incumbent, decided once before
// any realization on forward-frame query geometry.
enum class RnaRivalClass : std::uint8_t {
  // The rank-1 chain, already realized by the caller. Exactly one hypothesis carries
  // this, and it is always element 0 of the hypothesis vector.
  Incumbent = 0,
  // Competes for the same query bases (mask_level overlap, or censored geometry): it may
  // win the primary slot, and if it loses it is what dp2 / n_sub measure.
  Competitor = 1,
  // Explains a different part of the read: never elected, never dp2 or n_sub.
  CoPrimary = 2,
};

// One hypothesis in the lifecycle. Its bundle, summary and locus are borrowed from the
// backend's per-read storage, which outlives the vector, so reordering the vector moves
// only pointers. `result` is the one owned member, the realization run for a
// non-incumbent; it stays empty for element 0, whose realization the caller holds.
struct RnaRealizedHypothesis {
  int catalogue_index = -1;
  int reference_id = -1;
  const placement::CoarseLocus* locus = nullptr;
  const ExactAnchorPath* bundle = nullptr;
  const ExactAnchorPathSummary* summary = nullptr;
  SiblingChainStats siblings;
  // The catalogue's order key, so the committed-frame MAPQ rivals can be rebuilt from
  // this vector.
  std::int64_t rank_score = 0;
  // The chain's own values, in its own mapping frame; `reverse` is the post-repair
  // strand, read straight.
  bool chained = false;
  int chain_score = 0;
  int chain_anchors = 0;
  int q_begin = -1;
  int q_end = -1;
  bool reverse = false;
  bool repaired = false;
  // Anchor-union query coverage of `bundle`, the commit's qcov input.
  std::uint64_t query_covered_bases = 0;
  bool coverage_ok = false;
  // True only for a hypothesis the query partition appended after the primary's MAPQ
  // was fixed; false for the incumbent and every catalogue rival.
  bool partition_nominated = false;
  // For a nominee re-chained from a hypothesis already in the vector: that hypothesis's
  // index. The twin is the same placement, not a rival, and the family's pricing skips
  // it. SIZE_MAX when there is none.
  std::size_t twin_hypothesis = static_cast<std::size_t>(-1);
  // Lifecycle verdicts.
  RnaRivalClass klass = RnaRivalClass::Competitor;
  bool admitted = false;
  bool realized = false;
  // This chain is the committed placement re-found from a neighbouring catalogue window
  // (RnaChainMapqRival::shadow). It stays in the lifecycle but never supplies dp2 and is
  // inadmissible to the MAPQ.
  bool shadow = false;
  // The election price: segment 0's realized DP maximum, minimap2's dp_max currency.
  // The election and the dp2 selection read only this field.
  std::optional<int> dp_maximum;
  // Segment 0 alone, and the realized segments' forward-frame query hull
  // [aq_begin, aq_end), -1/-1 when unmeasured (rival_pricing.h).
  std::optional<int> dp_segment0;
  int aq_begin = -1;
  int aq_end = -1;
  // The rank-recalibration price, set only when the recalibration fired for this read
  // and this hypothesis's accounting could be formed. dp_maximum keeps the raw price.
  std::optional<int> dp_recal;
  RnaSpliceRealizationResult result;
};

// Step 1: minimap2's mask_level verdict on the forward-frame spans
// (rna_query_spans_compete, evaluated by the caller). Censored geometry competes: an
// unmeasurable span is not evidence of disjointness.
RnaRivalClass rna_classify_rival(bool spans_compete) noexcept;

// Step 2a: minimap2's secondary retention band, as an OR: within pri_ratio of the
// primary chain score, or within min_diff of it. Both inclusive; an unchained rival is
// never admitted.
bool rna_competitor_admitted(int chain_score, int primary_chain_score,
                             const RnaRivalLifecycleConfig& config) noexcept;

// Step 2b: a co-primary's floor, on chain score, anchor count and the query bases it adds
// over the incumbent (rna_uncovered_query_bases). minimap2's numbers, fixed.
inline constexpr int kRnaCoprimaryMinChainScore = 40;
inline constexpr int kRnaCoprimaryMinAnchors = 3;
inline constexpr int kRnaCoprimaryMinQueryBases = 50;

bool rna_coprimary_admitted(int chain_score, int chain_anchors,
                            int uncovered_query_bases) noexcept;

// Step 2c: the chimeric emission limits, stricter than the co-primary floor and applied
// to realized geometry, since an emitted family must have been aligned, not merely
// chained: at most kRnaChimeraMaxFamilies second families per read, each with a realized
// hull of at least kRnaChimeraMinQueryBases. Fixed, not tunable.
inline constexpr int kRnaChimeraMaxFamilies = 2;   // extra families per read
                                                   // beyond the primary
inline constexpr int kRnaChimeraOverlapNumerator = 1;
inline constexpr int kRnaChimeraOverlapDenominator = 10;
inline constexpr int kRnaChimeraMinQueryBases = 100;

// Both spans in the forward frame. Disjoint when the overlap is at most
// floor(shorter / 10); a censored span on either side is not disjoint.
bool rna_chimera_query_disjoint(RnaQuerySpan a, RnaQuerySpan b) noexcept;

// A realized second family's quality bar: a segment-0 DP maximum of at least
// `minimum_dp_maximum`, the floor the primary realization is held to, and a hull
// covering at least kRnaChimeraMinQueryBases.
bool rna_chimera_family_admitted(std::optional<int> dp_segment0,
                                 RnaQuerySpan hull,
                                 int minimum_dp_maximum) noexcept;

// Step 3: the realization queue's total order: chain score, then anchor count, then
// catalogue index. Total, so the bounded budget never depends on the sort.
bool rna_rival_queue_before(const RnaRealizedHypothesis& a,
                            const RnaRealizedHypothesis& b) noexcept;

// Step 3b: minimap2's rank recalibration (mm_update_dp_max / mm_recal_max_dp), between
// realization and the election. On raw DP a spliced alignment pays the long-gap open at
// every junction while an intronless retrocopy pays nothing, so the retrocopy can
// outscore the gene. When two aligned hits both cover most of the read, minimap2 ranks
// every hit by an intron-free, divergence-amplified score instead:
//
//     match_sc * (mlen - b2*n_mis - gap_cost)
//
// Introns contribute nothing, a mismatch costs b2 = 0.5/divergence match units, and a
// gap costs b2 + log2(1+len) per operation with no per-base cost. A hypothesis whose
// accounting cannot be formed keeps its raw price.

// minimap2's rank_min_len, rank_frac and divergence floor.
inline constexpr int kRnaRankRecalMinReadLen = 500;
inline constexpr double kRnaRankRecalFrac = 0.9;
inline constexpr double kRnaRankRecalMinDivergence = 0.02;

// One realized hypothesis's primary-segment accounting and numeric CIGAR, borrowed from
// RnaSpliceHypothesisResult::segment_cigars. The counters follow minimap2's
// mm_update_extra: `matches` is mlen (unambiguous equal residues), `ambiguities` is
// n_ambi (excluded from the block length), and `mismatches` is blen - mlen - n_gap. The
// CIGAR uses the realizer's encoding, op = value & 0xf over "MIDNSHP=XB".
struct RnaRankRecalAccounting {
  int matches = 0;
  int mismatches = 0;
  int ambiguities = 0;
  const std::uint32_t* cigar = nullptr;
  std::size_t cigar_len = 0;
};

// minimap2's mm_event_identity, mlen / (blen + n_ambi - n_gap + n_gapo): identity over
// events, each gap counting once however long it is. blen is matches + mismatches +
// n_gap; introns are neither gap nor block. Returns -1 on an empty or malformed CIGAR,
// negative counters or a non-positive denominator, and the caller then skips the
// recalibration.
double rna_event_identity(const RnaRankRecalAccounting& accounting) noexcept;

// minimap2's divergence-to-mismatch weight:
//
//     div = 1 - identity;  if (div < 0.02) div = 0.02;
//     b2 = .5 / div;  if (b2 * a < b) b2 = (double)a / b;
//
// The second line is transcribed as minimap2 writes it, reciprocal included, so the same
// hypothesis is elected. `mismatch_sc` is minimap2's positive `b`; a non-positive a or b
// skips that line.
double rna_rank_recal_b2(double top_event_identity, int match_sc,
                         int mismatch_sc) noexcept;

// minimap2's mm_recal_max_dp. The caller floors the result at 0, as mm_update_dp_max
// does.
int rna_recal_max_dp(const RnaRankRecalAccounting& accounting, double b2,
                     int match_sc) noexcept;

// minimap2's gate: the read is long enough to rank, the top hit is aligned across most of
// it (its realized forward-frame hull [aq_begin, aq_end)), and the runner-up is within
// `frac` of the top. Both DP prices must be positive.
bool rna_rank_recal_gate(int read_len, int top_aq_begin, int top_aq_end,
                         int top_dp, int second_dp) noexcept;

// Raw-price veto band for the combined election, in match units. Not a minimap2 number:
// electing on the recalibrated price alone flips reads the raw price placed correctly,
// so a flip may overturn at most this much raw margin.
inline constexpr int kRnaRankRecalRawVetoMargin = 20;

// Step 4: the elected primary's index; element 0 must be the incumbent. Only realized
// incumbent and competitor entries with a positive DP maximum take part. Returns 0 unless
// a competitor strictly beat a present incumbent maximum.
std::size_t rna_elect_primary_hypothesis(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept;

// Step 4 on a recalibrated read: argmax over raw segment-0 plus recalibrated price with
// the same tie discipline. A winner that differs from the raw election's takes the slot
// only when the raw margin it overturns is within kRnaRankRecalRawVetoMargin
// (inclusive). Missing evidence falls back to the raw verdict.
std::size_t rna_elect_primary_rank_recal(
    const std::vector<RnaRealizedHypothesis>& hypotheses) noexcept;

// Step 5: the hypothesis supplying the MAPQ's dp2: the largest realized DP maximum among
// the non-elected competitors (the demoted incumbent included), ties to the lower
// catalogue index. Co-primaries and shadows are excluded.
std::optional<std::size_t>
rna_select_dp2_owner(const std::vector<RnaRealizedHypothesis>& hypotheses,
                     std::size_t elected) noexcept;

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
