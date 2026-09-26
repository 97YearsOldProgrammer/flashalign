// Post-DP scoring for the HiFi presets: the dp1 / dp2 the MAPQ reads. As in
// minibwa, each record's CIGAR is re-walked under a log-gap cost
// (mm_update_extra), and when the best hit covers the read and a runner-up is
// comparable, every record is rescored with one b2 derived from the best hit's
// identity (mb_recal_max_dp: a gap costs b2 * log2(1 + len) and the clipped
// query is charged as mismatches).
//
// Runs once per read with CIGAR output, after the incumbent and alternative
// families are realized. It does not decide the promotion, which uses the raw
// decision score; the backend swaps dp1 / dp2 on a promotion. One minibwa hit
// corresponds to one emitted record: the trigger compares the two primaries,
// and rescoring touches every record of both families.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include the DNA post-DP scoring stage"
#endif

#include "context.h"
#include "dp_runner.h" // DpScoringParams
#include "family_realization.h"

#include <cstdint>
#include <utility>
#include <vector>

namespace fa::cpu::lr {

struct DnaPostDpOutcome {
  bool ran = false;
  // Records priced per family: the primary, then the supplementaries.
  int incumbent_records = 0;
  int alternative_records = 0; // 0 when no alternative was realized
  // In the incumbent's frame: dp1 prices the incumbent records that compete
  // with the alternative's query span, dp2 the alternative primary (0 without
  // one). The backend swaps them on a promotion.
  int mapq_dp1 = 0;
  int mapq_dp2 = 0;
  // Family totals: `*_dp_raw` is the decision score the commit uses, the
  // other two are sums over the records.
  int incumbent_dp_raw = 0;
  int incumbent_dp_sweep = 0;
  int incumbent_dp_rescored = 0;
  int alternative_dp_raw = 0;
  int alternative_dp_sweep = 0;
  int alternative_dp_rescored = 0;
  // Rescoring inputs.
  bool triggered = false;
  double b2 = 0.0;
  double best_identity = 0.0;
  int best_span = 0;
  int runner_span = 0;
  // Query bases outside the best hit: max(0, read_len - best_span).
  int clip_bp = 0;
};

// Prices one read. `alternative` is null when no alternative was realized.
// `fwd` / `rc` are the encoded strands the CIGARs were committed against.
DnaPostDpOutcome
run_dna_postdp_scoring(const DnaContext& context,
                       const DnaFamilyRealizationOutcome& incumbent,
                       const DnaFamilyRealizationOutcome* alternative,
                       const std::vector<std::uint8_t>& fwd,
                       const std::vector<std::uint8_t>& rc, int read_len);

// The ms:i value of one record: the maximum-scoring segment of its CIGAR under
// the log-gap cost, as minimap2 computes ms. Independent of postdp_rescoring.
// Returns -1 (no tag) for an unmapped record, an empty CIGAR or an unknown
// contig.
int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc);

// One record's dual-affine score split by a forward-query interval [lo, hi):
// what it earns on the bases inside and outside. With [lo, hi) the overlap of
// two hypotheses' primary spans, the inside margin is the part of
// dp1_raw - dp2_raw the two actually contest.
//
// Scoring is ksw2's: +match, -mismatch, -ambi for an ambiguous base, and a gap
// of length L costs min(o1 + e1 * L, o2 + e2 * L). There is no end bonus or
// Z-drop, so inside + outside need not equal the record's score. An aligned
// base counts on the side of its forward coordinate; a gap counts once, an
// insertion at the midpoint of its run and a deletion at the query base that
// follows it (clamped to the last base). An empty interval puts everything
// outside.
//
// `ops` is the CIGAR (parse_cigar_ops form, clips included); `query` is the
// oriented encoded query from base 0 and `reference` the encoded contig from
// the record's pos.
struct AffineIntervalScore {
  long long inside = 0;
  long long outside = 0;
  int aligned_inside = 0;  // M/=/X query bases attributed inside
  int aligned_outside = 0;
};

AffineIntervalScore affine_score_over_query_interval(
    const std::vector<std::pair<int, char>>& ops, const std::uint8_t* query,
    const std::uint8_t* reference, int read_len, bool is_reverse, int lo,
    int hi, const DpScoringParams& dp) noexcept;

} // namespace fa::cpu::lr
