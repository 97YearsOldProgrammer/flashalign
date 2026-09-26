// Prices a realized RNA hypothesis for the rival lifecycle. A realized hypothesis is a
// family (a primary plus its continuation segments), while minimap2's dp_max is per hit.
// Every segment's dp_maximum is already in region_evidence, in the read's forward frame,
// so pricing runs no DP.
#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include RNA rival pricing"
#endif

#include "chain_mapq.h" // RnaQuerySpan, rna_query_spans_compete
#include "realization/splice_realizer.h" // RnaSpliceRealizationResult

#include <optional>

namespace fa::cpu::lr::rna {

// The winning transcript-orientation hypothesis, or null when the realization was refused
// or produced no segment (the predicate realized_primary_dp_maximum uses).
const RnaSpliceHypothesisResult*
rna_winning_hypothesis(const RnaSpliceRealizationResult& result) noexcept;

// Segment 0's dp_maximum; equals realized_primary_dp_maximum.
std::optional<int>
rna_segment0_dp_maximum(const RnaSpliceRealizationResult& result) noexcept;

// The family's realized query footprint in the read's forward frame: the
// half-open hull [min begin, max end) over its segments. Null when there is
// no winning hypothesis or no segment carries a valid span.
std::optional<RnaQuerySpan>
rna_realized_query_span(const RnaSpliceRealizationResult& result) noexcept;

// Sum of dp_maximum over the segments whose forward-frame span competes with `rival_span`
// under mask_level; segment 0's value when none competes or the rival span is invalid.
std::optional<int>
rna_overlap_dp_maximum(const RnaSpliceRealizationResult& result,
                       RnaQuerySpan rival_span) noexcept;

} // namespace fa::cpu::lr::rna
