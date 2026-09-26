#include "residue_trigger.h"

#include "residue_emission.h"

#include <algorithm>
#include <limits>

namespace fa::cpu::lr {
namespace {

void accumulate(const AlignResult& record, DnaCommittedCoverage& coverage) {
  if (!record.mapped()) return;
  const int begin = std::max(0, record.query_start);
  const int end = std::min(coverage.read_length, record.query_end);
  if (end > begin) coverage.spans.push_back({begin, end});
  for (const AlignResult& child : record.supplementary)
    accumulate(child, coverage);
}

}  // namespace

DnaCommittedCoverage dna_committed_query_coverage(const dna::Result& primary) {
  DnaCommittedCoverage coverage;
  coverage.read_length = primary.read_len;
  if (coverage.read_length <= 0 || !primary.mapped()) return coverage;
  accumulate(primary, coverage);
  std::sort(coverage.spans.begin(), coverage.spans.end(),
            [](const DnaQuerySpan& left, const DnaQuerySpan& right) {
              if (left.begin != right.begin) return left.begin < right.begin;
              return left.end < right.end;
            });
  std::vector<DnaQuerySpan> merged;
  for (const DnaQuerySpan& span : coverage.spans) {
    if (!merged.empty() && span.begin <= merged.back().end) {
      merged.back().end = std::max(merged.back().end, span.end);
    } else {
      merged.push_back(span);
    }
  }
  coverage.spans = std::move(merged);
  return coverage;
}

std::vector<DnaQuerySpan> dna_terminal_clip_intervals(
    const DnaCommittedCoverage& coverage, int minimum_bp) {
  std::vector<DnaQuerySpan> intervals;
  if (coverage.spans.empty() || coverage.read_length <= 0 || minimum_bp <= 0)
    return intervals;
  const DnaQuerySpan head{0, coverage.spans.front().begin};
  if (head.length() >= minimum_bp) intervals.push_back(head);
  const DnaQuerySpan tail{coverage.spans.back().end, coverage.read_length};
  if (tail.length() >= minimum_bp) intervals.push_back(tail);
  return intervals;
}

bool dna_residue_interval_non_interposing(const DnaCommittedCoverage& coverage,
                                          int query_begin,
                                          int query_end) noexcept {
  return !coverage.spans.empty() && query_end > query_begin &&
         (query_end <= coverage.spans.front().begin ||
          query_begin >= coverage.spans.back().end);
}

bool dna_residue_span_committed(const DnaCommittedCoverage& coverage,
                                int query_begin, int query_end) noexcept {
  if (query_end <= query_begin) return true;
  for (const DnaQuerySpan& span : coverage.spans)
    if (span.begin <= query_begin && query_end <= span.end) return true;
  return false;
}

int dna_clip_nominate_min_anchors(const DnaResidueObservedDensity& observed,
                                  int interval_length) noexcept {
  if (observed.anchors <= 0 || observed.query_span <= 0 ||
      interval_length <= 0)
    return kDnaResidueMinChainAnchors;
  const std::int64_t scaled =
      static_cast<std::int64_t>(kDnaClipNominateAnchorNumerator) *
      static_cast<std::int64_t>(observed.anchors) *
      static_cast<std::int64_t>(interval_length) /
      (static_cast<std::int64_t>(kDnaClipNominateAnchorDenominator) *
       static_cast<std::int64_t>(observed.query_span));
  const std::int64_t floored =
      std::max<std::int64_t>(kDnaResidueMinChainAnchors, scaled);
  return static_cast<int>(
      std::min<std::int64_t>(floored, std::numeric_limits<int>::max()));
}

bool dna_clip_nominate_interval_bars(int interval_length,
                                     int chain_query_span) noexcept {
  if (interval_length <= 0) return false;
  return static_cast<std::int64_t>(std::max(0, chain_query_span)) * 100 >=
         static_cast<std::int64_t>(kDnaClipNominateMinSpanPercent) *
             static_cast<std::int64_t>(interval_length);
}

bool dna_terminal_clip_interval_bars(int interval_length, int chain_query_span,
                                     int chain_anchors) noexcept {
  if (interval_length <= 0) return false;
  const std::int64_t length = interval_length;
  return static_cast<std::int64_t>(std::max(0, chain_query_span)) * 100 >=
             static_cast<std::int64_t>(kDnaTerminalClipMinSpanPercent) *
                 length &&
         static_cast<std::int64_t>(std::max(0, chain_anchors)) * 100 >=
             static_cast<std::int64_t>(
                 kDnaTerminalClipMinAnchorsPer100Bp) * length;
}

DnaResidueTriggerOutcome dna_run_terminal_clip_recovery(
    const DnaContext& context, const DnaPlacementChainingResult& placement,
    const ChainSeedLookupCache& lookup_cache,
    const std::vector<std::uint8_t>& forward_query,
    const std::vector<std::uint8_t>& reverse_query, bool cigar_lane,
    dna::Result& primary) {
  DnaResidueTriggerOutcome result;
  // Reference bases are not required here: map-only projection reads none,
  // and the realizer checks for them itself.
  if (!placement.accepted || !placement.family.valid || !primary.mapped() ||
      context.ref.names == nullptr || context.ref.index == nullptr)
    return result;

  const DnaPlacementFamily& family = placement.family;
  const DnaCommittedCoverage coverage = dna_committed_query_coverage(primary);
  const int minimum_bp = std::max(1, context.opts.residue_min_interval_bp);
  struct TriggerInterval {
    DnaQuerySpan span;
    // A --dna-clip-nominate attempt; these come after all regular ones.
    bool nominate = false;
  };
  std::vector<TriggerInterval> intervals;
  for (const DnaQuerySpan& span :
       dna_terminal_clip_intervals(coverage, minimum_bp))
    intervals.push_back({span, false});

  // --dna-clip-nominate reopens the same clips after the regular attempts. A
  // clip that already produced a record is refused by the committed-span
  // guard, so nomination only adds records.
  const bool clip_nominate = context.opts.clip_nominate;
  if (clip_nominate) {
    for (const DnaQuerySpan& span :
         dna_terminal_clip_intervals(coverage, kDnaClipNominateMinIntervalBp))
      intervals.push_back({span, true});
  }
  int budget = kDnaResidueMaxAdmissionsPerRead;
  const std::uint32_t occurrence_cap = static_cast<std::uint32_t>(
      std::max(0, context.opts.cigar_local_global_occ));
  // Nomination never exceeds the whole-query pool's occurrence gate, so a key
  // the gate kept out of the chain cannot return as a clip record.
  const std::uint32_t nominate_ceiling =
      context.opts.dna_pool_gate_occ > 0
          ? std::min(kDnaClipNominateOccurrenceCeiling,
                     static_cast<std::uint32_t>(
                         context.opts.dna_pool_gate_occ))
          : kDnaClipNominateOccurrenceCeiling;
  // Intervals come from the snapshot above; `committed` is refreshed after
  // each emitted record, so no later interval restates it.
  DnaCommittedCoverage committed = coverage;
  // The committed chain's anchor density, which scales nomination's floors.
  const DnaResidueObservedDensity observed =
      dna_residue_observed_density(placement);

  for (const TriggerInterval& opened : intervals) {
    const DnaQuerySpan& interval = opened.span;
    if (budget <= 0) break;
    const bool nominate = opened.nominate;
    DnaResidueAdmissionBar bar;
    if (nominate) {
      // An interval-scaled anchor floor, the same score floor and no density
      // floor: a second breakend arm is sparse by nature.
      bar.min_chain_anchors =
          dna_clip_nominate_min_anchors(observed, interval.length());
      bar.min_anchor_density_per_100bp = 0;
    }
    const int clusters_per_strand =
        nominate ? kDnaClipNominateMaxClustersPerStrand
                 : kDnaResidueMaxClustersPerStrand;
    const std::uint32_t interval_cap =
        nominate ? nominate_ceiling : occurrence_cap;
    if (!dna_residue_interval_non_interposing(coverage, interval.begin,
                                              interval.end))
      continue;
    if (dna_residue_cached_supply(
            family, placement.residue_fine_forward,
            placement.residue_fine_reverse, lookup_cache, interval_cap,
            interval.begin, interval.end) < kDnaResidueMinClusterAnchors)
      continue;
    std::vector<DnaResidueAnchor> forward_anchors;
    std::vector<DnaResidueAnchor> reverse_anchors;
    if (!dna_residue_collect_anchors(
            context, family, placement.residue_fine_forward,
            placement.residue_fine_reverse, lookup_cache, occurrence_cap,
            interval.begin, interval.end, forward_anchors, reverse_anchors)) {
      if (!nominate) continue;
      // The collector refuses an interval whose seeds overrun the posting
      // budget; nomination goes on to the rarest-first sweep instead.
      forward_anchors.clear();
      reverse_anchors.clear();
    }
    DnaResidueAdmission admission;
    // Picks the best (chain_score, chain_anchors) chain regardless of
    // geometry; the geometry test comes later.
    const auto pick = [&](const std::vector<DnaResidueAnchor>& forward,
                          const std::vector<DnaResidueAnchor>& reverse) {
      return dna_residue_best_admission(context, family, forward, reverse,
                                        admission, bar, clusters_per_strand);
    };
    bool admitted = pick(forward_anchors, reverse_anchors);
    if (!admitted && nominate) {
      // Nothing admissible under the regular cap: sweep this interval's
      // cached views rarest-first up to the higher ceiling. No index lookup.
      forward_anchors.clear();
      reverse_anchors.clear();
      dna_residue_collect_anchors_rarest_first(
          context, family, placement.residue_fine_forward,
          placement.residue_fine_reverse, lookup_cache, nominate_ceiling,
          interval.begin, interval.end, forward_anchors, reverse_anchors);
      admitted = pick(forward_anchors, reverse_anchors);
    }
    if (!admitted)
      continue;
    if (admission.cluster.contig < 0 ||
        admission.cluster.contig >= static_cast<int>(context.ref.names->size()))
      continue;

    DnaPlacementFamily window = family;
    const auto id = static_cast<::fa::cpu::voting::CandidateId>(
        window.candidates.size());
    DnaResidueDetachedChain detached = dna_residue_detached_chain(
        family, admission.cluster, admission.reverse, admission.outcome, id);
    window.candidates.push_back(detached.candidate);
    // A chain inside committed query only restates an existing record.
    if (dna_residue_span_committed(committed,
                                   detached.chain.forward_query_begin,
                                   detached.chain.forward_query_end))
      continue;
    DnaPlacementChainingResult scratch;
    scratch.accepted = true;
    scratch.candidates.push_back(detached.chain);

    const bool same_contig =
        (*context.ref.names)[static_cast<std::size_t>(admission.cluster.contig)] ==
        primary.chromosome;
    const std::int64_t gap =
        same_contig
            ? dna_residue_reference_gap(detached.chain.reference_begin,
                                        detached.chain.reference_end,
                                        primary.pos, primary.target_end)
            : -1;
    // A regular clip must lie on the primary's contig within
    // kDnaResidueEmissionLocalBp of it, on either strand. Nomination skips the
    // test: a translocation's other breakend is on another contig.
    const bool geometry =
        nominate || (same_contig && gap >= 0 &&
                     gap < kDnaResidueEmissionLocalBp);
    if (!geometry)
      continue;
    const int chain_span =
        std::max(0, detached.chain.forward_query_end -
                        detached.chain.forward_query_begin);
    const bool interval_bars =
        nominate ? dna_clip_nominate_interval_bars(interval.length(),
                                                   chain_span)
                 : dna_terminal_clip_interval_bars(
                       interval.length(), chain_span,
                       admission.outcome.chain_anchors);
    if (!interval_bars)
      continue;
    ::fa::cpu::voting::QueryBlock block;
    if (!dna_residue_chain_block(window, detached.candidate, detached.chain,
                                 block,
                                 nominate ? bar.min_chain_anchors
                                          : kDnaResidueEmissionMinAnchors))
      continue;
    if (dna_emit_residue_record(
            context, window, scratch, detached.chain, block, forward_query,
            reverse_query, cigar_lane, primary, result.ksw2_attempts,
            result.estimated_cells, result.geometry)) {
      ++result.emitted;
      --budget;
      committed = dna_committed_query_coverage(primary);
    }
  }
  return result;
}

}  // namespace fa::cpu::lr
