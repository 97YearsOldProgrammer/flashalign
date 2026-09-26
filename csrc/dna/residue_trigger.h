#pragma once

#include "geometry_work.h"
#include "context.h"
#include "placement_chaining.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

struct DnaQuerySpan {
  int begin = 0;
  int end = 0;
  int length() const noexcept { return end - begin; }
};

struct DnaCommittedCoverage {
  int read_length = 0;
  std::vector<DnaQuerySpan> spans;
};

// A terminal-clip chain must span at least this share of the clip and carry
// at least this many anchors per 100 bases of it.
inline constexpr int kDnaTerminalClipMinSpanPercent = 85;
inline constexpr int kDnaTerminalClipMinAnchorsPer100Bp = 9;

DnaCommittedCoverage dna_committed_query_coverage(const dna::Result& primary);
std::vector<DnaQuerySpan> dna_terminal_clip_intervals(
    const DnaCommittedCoverage& coverage, int minimum_bp);
bool dna_residue_interval_non_interposing(const DnaCommittedCoverage& coverage,
                                          int query_begin,
                                          int query_end) noexcept;
// True when [query_begin, query_end) lies inside the committed query
// coverage, so a chain there would only restate an existing record.
bool dna_residue_span_committed(const DnaCommittedCoverage& coverage,
                                int query_begin, int query_end) noexcept;
bool dna_terminal_clip_interval_bars(int interval_length, int chain_query_span,
                                     int chain_anchors) noexcept;

// --dna-clip-nominate: a terminal clip of at least
// kDnaClipNominateMinIntervalBp gets a second attempt after the regular one
// declined it.
inline constexpr int kDnaClipNominateMinIntervalBp = 500;
// The admitted chain must span at least this share of the clip, so a repeat
// fragment cannot nominate a whole arm.
inline constexpr int kDnaClipNominateMinSpanPercent = 50;
// Anchor floor: 3/10 of the anchors the read's committed chain would place in
// an interval this long, and at least kDnaResidueMinChainAnchors.
inline constexpr int kDnaClipNominateAnchorNumerator = 3;
inline constexpr int kDnaClipNominateAnchorDenominator = 10;
// Occurrence ceiling of the rarest-first sweep over the clip's cached views.
// Under the whole-query pool gate the ceiling is capped at the gate, so the
// sweep cannot re-admit keys the gate kept out of the chain.
inline constexpr std::uint32_t kDnaClipNominateOccurrenceCeiling = 4095;
// Diagonal clusters tried per strand. A clipped arm in a satellite has many
// decoy diagonals, so its true locus often misses the regular top two.
inline constexpr int kDnaClipNominateMaxClustersPerStrand = 8;

int dna_clip_nominate_min_anchors(const DnaResidueObservedDensity& observed,
                                  int interval_length) noexcept;
bool dna_clip_nominate_interval_bars(int interval_length,
                                     int chain_query_span) noexcept;

struct DnaResidueTriggerOutcome {
  int emitted = 0;
  int ksw2_attempts = 0;
  std::int64_t estimated_cells = 0;
  DnaGeometryWork geometry;
};

// Recovers the terminal clips of the committed family as extra supplementary
// records of `primary`. The family, catalogue and partition are unchanged.
DnaResidueTriggerOutcome dna_run_terminal_clip_recovery(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    const ChainSeedLookupCache& lookup_cache,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query, bool cigar_lane,
    dna::Result& primary);

}  // namespace fa::cpu::lr
