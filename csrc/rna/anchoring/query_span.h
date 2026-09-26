// Forward-frame query spans and minimap2's mask_level competition test (mm_set_parent
// with uncov_len = 0), shared by the fine chain's sibling scan and the chain MAPQ.
#pragma once

#include <algorithm>

namespace fa::cpu::lr::rna {

// minimap2's mask_level: a rival overlapping the committed hypothesis over
// more than half the shorter query span competes with it.
inline constexpr double kRnaChainMapqMaskLevel = 0.5;

// A query span in its chain's own mapping frame: the read for a forward chain, its
// reverse complement for a reverse one. Spans from opposite strands are comparable only
// after rna_forward_query_span.
struct RnaQuerySpan {
  int begin = -1;
  int end = -1;
};

namespace query_span_detail {

// The invalid span (begin < 0 or empty) marks censored geometry.
inline bool forward_span_valid(RnaQuerySpan span) noexcept {
  return span.begin >= 0 && span.end > span.begin;
}

// Shared overlap arithmetic; negative when the two spans are disjoint.
inline int forward_span_overlap(RnaQuerySpan winner,
                                RnaQuerySpan rival) noexcept {
  return std::min(winner.end, rival.end) - std::max(winner.begin, rival.begin);
}

} // namespace query_span_detail

// Projects a mapping-oriented span into the read's forward frame; a reverse span mirrors
// to [read_len - end, read_len - begin). An invalid span, or a reverse span leaving the
// read, projects to the invalid span.
inline RnaQuerySpan rna_forward_query_span(int oriented_begin,
                                           int oriented_end, bool reverse,
                                           int read_len) noexcept {
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

// minimap2's mask_level test on forward-frame spans: an invalid span on either side
// competes (censored evidence is not disjoint evidence), a non-positive overlap does not,
// otherwise the overlap must exceed half the shorter span.
inline bool rna_query_spans_compete(RnaQuerySpan winner,
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

// Rival query bases the winner's span does not cover; 0 when either span is invalid.
inline int rna_uncovered_query_bases(RnaQuerySpan winner,
                                     RnaQuerySpan rival) noexcept {
  using query_span_detail::forward_span_overlap;
  using query_span_detail::forward_span_valid;
  if (!forward_span_valid(winner) || !forward_span_valid(rival))
    return 0;
  const int overlap = std::max(0, forward_span_overlap(winner, rival));
  return (rival.end - rival.begin) - overlap;
}

} // namespace fa::cpu::lr::rna
