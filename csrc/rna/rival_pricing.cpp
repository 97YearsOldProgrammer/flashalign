#include "rival_pricing.h"

#include <algorithm>

namespace fa::cpu::lr::rna {

namespace {

RnaQuerySpan segment_span(const SpliceRegionEvidence& evidence) noexcept {
  RnaQuerySpan span;
  span.begin = evidence.query_begin;
  span.end = evidence.query_end;
  return span;
}

bool span_valid(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

} // namespace

const RnaSpliceHypothesisResult*
rna_winning_hypothesis(const RnaSpliceRealizationResult& result) noexcept {
  if (result.refused || result.segments.empty())
    return nullptr;
  for (const RnaSpliceHypothesisResult& hypothesis : result.hypotheses) {
    if (hypothesis.orientation == result.winning_hypothesis &&
        !hypothesis.refused && !hypothesis.segments.empty())
      return &hypothesis;
  }
  return nullptr;
}

std::optional<int>
rna_segment0_dp_maximum(const RnaSpliceRealizationResult& result) noexcept {
  const RnaSpliceHypothesisResult* winner = rna_winning_hypothesis(result);
  if (winner == nullptr)
    return std::nullopt;
  return winner->primary_dp_maximum;
}

std::optional<RnaQuerySpan>
rna_realized_query_span(const RnaSpliceRealizationResult& result) noexcept {
  const RnaSpliceHypothesisResult* winner = rna_winning_hypothesis(result);
  if (winner == nullptr)
    return std::nullopt;
  std::optional<RnaQuerySpan> hull;
  for (const SpliceRegionEvidence& evidence : winner->region_evidence) {
    const RnaQuerySpan span = segment_span(evidence);
    if (!span_valid(span))
      continue;
    if (!hull) {
      hull = span;
    } else {
      hull->begin = std::min(hull->begin, span.begin);
      hull->end = std::max(hull->end, span.end);
    }
  }
  return hull;
}

std::optional<int>
rna_overlap_dp_maximum(const RnaSpliceRealizationResult& result,
                       RnaQuerySpan rival_span) noexcept {
  const RnaSpliceHypothesisResult* winner = rna_winning_hypothesis(result);
  if (winner == nullptr)
    return std::nullopt;
  if (winner->region_evidence.empty() || !span_valid(rival_span))
    return winner->primary_dp_maximum;
  int total = 0;
  bool any = false;
  for (const SpliceRegionEvidence& evidence : winner->region_evidence) {
    const RnaQuerySpan span = segment_span(evidence);
    if (!span_valid(span))
      continue;
    if (!rna_query_spans_compete(span, rival_span))
      continue;
    total += std::max(0, evidence.dp_maximum);
    any = true;
  }
  if (!any)
    return winner->primary_dp_maximum;
  return total;
}

} // namespace fa::cpu::lr::rna
