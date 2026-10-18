// The anchor path of one block and the pre-DP control minimap2's mm_align1 does before the kernels run.
#include "dna.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

/**** ordered_anchor_path ****/

namespace fa {
namespace cpu {
namespace lr {
namespace ordered_anchor {
namespace {

constexpr std::uint32_t kKnownAnchorFlags =
    ::fa::cpu::chaining::ANCHOR_IGNORE | ::fa::cpu::chaining::ANCHOR_TANDEM |
    ::fa::cpu::chaining::ANCHOR_LONG_JOIN;

bool anchor_in_bounds(const ::fa::cpu::chaining::Anchor& anchor,
                      int query_length, int target_length) {
  if (anchor.span <= 0 || anchor.q < 0 || anchor.r < 0)
    return false;
  const std::int64_t query_end =
      static_cast<std::int64_t>(anchor.q) + anchor.span;
  const std::int64_t target_end =
      static_cast<std::int64_t>(anchor.r) + anchor.span;
  return query_end <= query_length && target_end <= target_length;
}

bool anchor_in_path_bounds(const ::fa::cpu::chaining::Anchor& anchor,
                           const OrderedAnchorPath& path) {
  return anchor.q >= path.query_bound_begin &&
         anchor.r >= path.target_bound_begin &&
         static_cast<std::int64_t>(anchor.q) + anchor.span <=
             path.query_bound_end &&
         static_cast<std::int64_t>(anchor.r) + anchor.span <=
             path.target_bound_end;
}

ValidationResult validate_anchor_vector(const AnchorWindow& anchors,
                                        int query_length, int target_length) {
  for (std::size_t i = 0; i < anchors.size(); ++i) {
    const auto& anchor = anchors[i];
    if (anchor.span <= 0)
      return {ValidationCode::InvalidAnchorSpan, i};
    if (!anchor_in_bounds(anchor, query_length, target_length))
      return {ValidationCode::AnchorOutOfBounds, i};
    if ((anchor.flags & ~kKnownAnchorFlags) != 0)
      return {ValidationCode::InvalidAnchorFlags, i};
    if (i == 0)
      continue;
    const auto& previous = anchors[i - 1];
    if (anchor.q <= previous.q || anchor.r <= previous.r)
      return {ValidationCode::NonMonotoneSelectedPath, i};
  }
  return {};
}

int score_expanded_reference_span(int query_span,
                                  const TerminalWindowControl& control) {
  std::int64_t span = query_span;
  const std::int64_t reward = span * static_cast<std::int64_t>(control.match);
  if (reward > control.gap_open) {
    span += (reward - control.gap_open) / control.gap_extend;
  }
  span = std::min<std::int64_t>(span, control.max_gap);
  return static_cast<int>(
      std::min<std::int64_t>(span, std::numeric_limits<int>::max()));
}

int fuzzy_match_length(const OrderedAnchorPath& path, std::size_t begin,
                       std::size_t end) {
  if (begin >= end)
    return 0;
  int length = path.selected[begin].span;
  for (std::size_t i = begin + 1; i < end; ++i) {
    const auto& anchor = path.selected[i];
    const auto& previous = path.selected[i - 1];
    const int target_length = anchor.r - previous.r;
    const int query_length = anchor.q - previous.q;
    length += target_length > anchor.span && query_length > anchor.span
                  ? anchor.span
                  : std::min(target_length, query_length);
  }
  return length;
}

std::vector<std::size_t> collect_long_gaps(const OrderedAnchorPath& path,
                                           std::size_t begin, std::size_t end,
                                           int minimum_gap) {
  std::vector<std::size_t> positions;
  for (std::size_t offset = 1; begin + offset < end; ++offset) {
    const auto& anchor = path.selected[begin + offset];
    const auto& previous = path.selected[begin + offset - 1];
    const int gap = (anchor.q - previous.q) - (anchor.r - previous.r);
    if (gap < -minimum_gap || gap > minimum_gap)
      positions.push_back(offset);
  }
  // As minimap2's collect_long_gaps: fewer than two long gaps count as none.
  if (positions.size() <= 1)
    positions.clear();
  return positions;
}

void filter_bad_seeds(OrderedAnchorPath& path, int minimum_gap,
                      int difference_threshold, int maximum_extension_length,
                      int maximum_extension_count) {
  const std::size_t begin = path.realization_begin;
  const std::vector<std::size_t> gaps =
      collect_long_gaps(path, begin, path.realization_end, minimum_gap);
  if (gaps.empty())
    return;

  int maximum = 0;
  int maximum_start = -1;
  int maximum_end = -1;
  for (int k = 0;; ++k) {
    int gap = 0;
    int insertions = 0;
    int deletions = 0;
    int maximum_difference = 0;
    int maximum_difference_end = -1;

    if (k == static_cast<int>(gaps.size()) || k >= maximum_end) {
      if (maximum_end > 0) {
        for (std::size_t i = gaps[static_cast<std::size_t>(maximum_start)];
             i < gaps[static_cast<std::size_t>(maximum_end)]; ++i) {
          path.selected[begin + i].flags |= ::fa::cpu::chaining::ANCHOR_IGNORE;
        }
      }
      maximum = 0;
      maximum_start = maximum_end = -1;
      if (k == static_cast<int>(gaps.size()))
        break;
    }

    const std::size_t i = gaps[static_cast<std::size_t>(k)];
    const auto& anchor = path.selected[begin + i];
    const auto& previous = path.selected[begin + i - 1];
    gap = (anchor.q - previous.q) - (anchor.r - previous.r);
    if (gap > 0)
      insertions += gap;
    else
      deletions += -gap;
    const int query_start = previous.q;
    const int target_start = previous.r;
    for (int l = k + 1;
         l < static_cast<int>(gaps.size()) && l <= k + maximum_extension_count;
         ++l) {
      const std::size_t j = gaps[static_cast<std::size_t>(l)];
      const auto& next = path.selected[begin + j];
      if (next.q - query_start > maximum_extension_length ||
          next.r - target_start > maximum_extension_length)
        break;
      const auto& next_previous = path.selected[begin + j - 1];
      gap = (next.q - next_previous.q) - (next.r - next_previous.r);
      if (gap > 0)
        insertions += gap;
      else
        deletions += -gap;
      const int difference =
          insertions + deletions - std::abs(insertions - deletions);
      if (maximum_difference < difference) {
        maximum_difference = difference;
        maximum_difference_end = l;
      }
    }
    if (maximum_difference > difference_threshold &&
        maximum_difference > maximum) {
      maximum = maximum_difference;
      maximum_start = k;
      maximum_end = maximum_difference_end;
    }
  }
}

void filter_bad_seeds_alt(OrderedAnchorPath& path, int minimum_gap,
                          int maximum_extension) {
  const std::size_t begin = path.realization_begin;
  const std::vector<std::size_t> gaps =
      collect_long_gaps(path, begin, path.realization_end, minimum_gap);
  if (gaps.empty())
    return;

  for (std::size_t k = 0; k < gaps.size();) {
    const std::size_t i = gaps[k];
    const auto& anchor = path.selected[begin + i];
    const auto& previous = path.selected[begin + i - 1];
    int gap1 = std::abs((anchor.q - previous.q) - (anchor.r - previous.r));
    int target_end1 = anchor.r;
    int query_end1 = anchor.q;
    std::size_t l = k + 1;
    for (; l < gaps.size(); ++l) {
      const std::size_t j = gaps[l];
      const auto& next = path.selected[begin + j];
      if (next.q - query_end1 > maximum_extension ||
          next.r - target_end1 > maximum_extension)
        break;
      const auto& next_previous = path.selected[begin + j - 1];
      const int gap2 =
          std::abs((next.q - next_previous.q) - (next.r - next_previous.r));
      const int target_start2 = next_previous.r + next_previous.span;
      const int query_start2 = next_previous.q + next_previous.span;
      const int matching_span =
          std::min(target_start2 - target_end1, query_start2 - query_end1);
      if (matching_span > gap1 + gap2)
        break;
      target_end1 = next.r;
      query_end1 = next.q;
      gap1 = gap2;
    }
    if (l > k + 1) {
      const std::size_t end = gaps[l - 1];
      for (std::size_t j = gaps[k]; j < end; ++j) {
        path.selected[begin + j].flags |= ::fa::cpu::chaining::ANCHOR_IGNORE;
      }
      path.selected[begin + end].flags |= ::fa::cpu::chaining::ANCHOR_LONG_JOIN;
    }
    k = l;
  }
}

} // namespace

ValidationResult validate(const OrderedAnchorPath& path) {
  if (path.selected.empty())
    return {ValidationCode::EmptySelectedPath, 0};
  if (path.reference_id < 0)
    return {ValidationCode::InvalidReference, 0};
  if (path.query_length <= 0 || path.target_length <= 0 ||
      path.minimizer_k <= 0)
    return {ValidationCode::InvalidBounds, 0};
  if (path.query_bound_begin < 0 ||
      path.query_bound_begin >= path.query_bound_end ||
      path.query_bound_end > path.query_length || path.target_bound_begin < 0 ||
      path.target_bound_begin >= path.target_bound_end ||
      path.target_bound_end > path.target_length)
    return {ValidationCode::InvalidBounds, 0};
  if (path.realization_begin >= path.realization_end ||
      path.realization_end > path.selected.size())
    return {ValidationCode::InvalidSelectedSubrange, 0};
  if (path.coordinate_convention != CoordinateConvention::ExactKmerStart)
    return {ValidationCode::UnsupportedCoordinateConvention, 0};
  if (path.homopolymer_compressed)
    return {ValidationCode::UnsupportedHomopolymerCompression, 0};

  const ValidationResult result = validate_anchor_vector(
      path.selected, path.query_length, path.target_length);
  if (!result)
    return result;
  for (std::size_t i = 0; i < path.selected.size(); ++i) {
    if (path.selected[i].span != path.minimizer_k)
      return {ValidationCode::InconsistentMinimizerSpan, i};
    if (!anchor_in_path_bounds(path.selected[i], path))
      return {ValidationCode::AnchorOutOfBounds, i};
  }
  return {};
}

AdjustedEndpoint adjusted_endpoint(const OrderedAnchorPath& path,
                                   const ::fa::cpu::chaining::Anchor& anchor) {
  // minimap2 stores the inclusive k-mer end and subtracts floor(k/2);
  // validate() guarantees span == k.
  const int shift = path.minimizer_k - 1 - (path.minimizer_k >> 1);
  return {anchor.q + shift, anchor.r + shift};
}

ValidationResult normalize_for_realization(OrderedAnchorPath& path,
                                           const NormalizationControl& control) {
  if (control.bandwidth <= 0 || control.min_chain_score <= 0 ||
      control.max_gap <= 0)
    return {ValidationCode::InvalidControl, 0};

  const std::size_t original_begin = path.realization_begin;
  const std::size_t original_end = path.realization_end;
  if (original_end - original_begin >= 3) {
    const int fuzzy_length =
        fuzzy_match_length(path, original_begin, original_end);
    std::size_t begin = original_begin;
    int length = path.selected[original_begin].span;
    int matched = length;
    for (std::size_t i = original_begin + 1; i + 1 < original_end; ++i) {
      const auto& anchor = path.selected[i];
      const auto& previous = path.selected[i - 1];
      if ((anchor.flags & ::fa::cpu::chaining::ANCHOR_LONG_JOIN) != 0)
        break;
      const int target_length = anchor.r - previous.r;
      const int query_length = anchor.q - previous.q;
      const int minimum = std::min(target_length, query_length);
      const int maximum = std::max(target_length, query_length);
      if (maximum - minimum > (length >> 1))
        begin = i;
      length += minimum;
      matched += std::min(minimum, anchor.span);
      if (length >= (control.bandwidth << 1) ||
          (matched >= control.min_chain_score * 2 &&
           matched >= control.bandwidth) ||
          matched >= (fuzzy_length >> 1))
        break;
    }

    std::size_t end = original_end;
    length = path.selected[original_end - 1].span;
    matched = length;
    for (std::size_t i = original_end - 2; i > begin; --i) {
      const auto& anchor = path.selected[i];
      const auto& next = path.selected[i + 1];
      if ((next.flags & ::fa::cpu::chaining::ANCHOR_LONG_JOIN) != 0)
        break;
      const int target_length = next.r - anchor.r;
      const int query_length = next.q - anchor.q;
      const int minimum = std::min(target_length, query_length);
      const int maximum = std::max(target_length, query_length);
      if (maximum - minimum > (length >> 1))
        end = i + 1;
      length += minimum;
      matched += std::min(minimum, next.span);
      if (length >= (control.bandwidth << 1) ||
          (matched >= control.min_chain_score * 2 &&
           matched >= control.bandwidth) ||
          matched >= (fuzzy_length >> 1))
        break;
    }
    if (control.trim_left)
      path.realization_begin = begin;
    if (control.trim_right)
      path.realization_end = end;
  }

  // minimap2's constants for non-splice long reads.
  filter_bad_seeds(path, /*minimum_gap=*/10,
                   /*difference_threshold=*/40,
                   /*maximum_extension_length=*/control.max_gap >> 1,
                   /*maximum_extension_count=*/10);
  filter_bad_seeds_alt(path, /*minimum_gap=*/30,
                       /*maximum_extension=*/control.max_gap >> 1);
  return {};
}

namespace {

// adjusted_endpoint() shifts both axes equally, so raw and adjusted
// coordinates lie on the same diagonal.
int anchor_diagonal(const ::fa::cpu::chaining::Anchor& anchor) {
  return anchor.r - anchor.q;
}

// The ungapped certificate over one gap of `length` bases (GapCertificate).
bool certified_gap(const std::uint8_t* query, const std::uint8_t* target,
                   int length, const GapCertificate& certificate) {
  if (length < 2)
    return true;
  std::int64_t ungapped = 0;
  for (int i = 0; i < length; ++i) {
    const std::size_t q = query[i] > 4 ? 4 : query[i];
    const std::size_t t = target[i] > 4 ? 4 : target[i];
    ungapped += certificate.matrix[q * 5 + t];
  }
  const std::int64_t max_gapped =
      static_cast<std::int64_t>(length - 1) * certificate.match -
      static_cast<std::int64_t>(2) * certificate.one_base_gap;
  return ungapped >= max_gapped;
}

} // namespace

GeometryPlanner::GeometryPlanner(const OrderedAnchorPath& path,
                                 const GeometryControl& control,
                                 const GapCertificate& certificate,
                                 GeometryPlan& plan)
    : path_(path), control_(control), certificate_(certificate), plan_(plan),
      begin_(path.realization_begin), last_(path.realization_end - 1),
      next_(begin_ + 1), open_{begin_, begin_, 1, false, 0}, corner_(begin_) {
  plan_ = GeometryPlan{};
  const AdjustedEndpoint opening =
      adjusted_endpoint(path_, path_.selected[begin_]);
  const AdjustedEndpoint closing =
      adjusted_endpoint(path_, path_.selected[last_]);
  plan_.initial_query_cursor = opening.query;
  plan_.initial_target_cursor = opening.target;
  plan_.final_query_cursor = closing.query;
  plan_.final_target_cursor = closing.target;
}

bool GeometryPlanner::reach(std::size_t s) {
  while (plan_.steps.size() <= s)
    if (!advance())
      return false;
  return true;
}

// Skipped anchors, as in minimap2; never the first or the last.
bool GeometryPlanner::skipped(std::size_t i) const {
  return i != begin_ && i != last_ &&
         (path_.selected[i].flags & (::fa::cpu::chaining::ANCHOR_IGNORE |
                                     ::fa::cpu::chaining::ANCHOR_TANDEM)) != 0;
}

// One fill step [center(left), center(right)): a seam or a piece.
GeometryStep GeometryPlanner::fill(GeometryStepKind kind, std::size_t left,
                                   std::size_t right) const {
  const ::fa::cpu::chaining::Anchor& prev = path_.selected[left];
  const ::fa::cpu::chaining::Anchor& next = path_.selected[right];
  const AdjustedEndpoint lo = adjusted_endpoint(path_, prev);
  const AdjustedEndpoint hi = adjusted_endpoint(path_, next);
  GeometryStep step;
  step.kind = kind;
  step.selected_index = right;
  step.query_begin = lo.query;
  step.query_end = hi.query;
  step.target_begin = lo.target;
  step.target_end = hi.target;
  step.long_join = (next.flags & ::fa::cpu::chaining::ANCHOR_LONG_JOIN) != 0;
  const int query_span = hi.query - lo.query;
  const int target_span = hi.target - lo.target;
  step.selected_band = step.long_join ? std::max(query_span, target_span)
                                      : control_.normal_long_band;
  step.anchor_free = std::max(next.q - prev.q_end(), next.r - prev.r_end());
  return step;
}

// A stretch from corner `start` to corner `end`, cut as minimap2's loop
// does: a piece ends at the stretch end, at a long join, or once both spans
// since the previous corner reach min_ksw_len.
void GeometryPlanner::emit_stretch(std::size_t start, std::size_t end,
                                   GeometryStretch stretch) {
  stretch.first_step = plan_.steps.size();
  stretch.query_bases = adjusted_endpoint(path_, path_.selected[end]).query -
                        adjusted_endpoint(path_, path_.selected[start]).query;
  const int index = static_cast<int>(plan_.stretches.size());
  std::size_t corner = start;
  AdjustedEndpoint corner_center =
      adjusted_endpoint(path_, path_.selected[start]);
  int crossed = 0;
  for (std::size_t i = start + 1; i <= end; ++i) {
    if (skipped(i))
      continue;
    const ::fa::cpu::chaining::Anchor& anchor = path_.selected[i];
    const AdjustedEndpoint center = adjusted_endpoint(path_, anchor);
    const bool long_join =
        (anchor.flags & ::fa::cpu::chaining::ANCHOR_LONG_JOIN) != 0;
    if (i == end || long_join ||
        (center.query - corner_center.query >= control_.min_ksw_len &&
         center.target - corner_center.target >= control_.min_ksw_len)) {
      GeometryStep step = fill(GeometryStepKind::Piece, corner, i);
      step.stretch = index;
      step.piece = stretch.pieces++;
      step.anchor_count = crossed;
      plan_.steps.push_back(step);
      corner = i;
      corner_center = center;
      crossed = 0;
      continue;
    }
    ++crossed;
  }
  plan_.stretches.push_back(stretch);
}

bool GeometryPlanner::advance() {
  using ::fa::cpu::chaining::Anchor;
  if (done_)
    return false;

  // 1. The open region grows over the retained anchors on its diagonal; one
  // off it, or a long join, opens the next region.
  GeometryRegion region = open_;
  bool more = false;
  for (; next_ <= last_; ++next_) {
    if (skipped(next_))
      continue;
    const Anchor& anchor = path_.selected[next_];
    const bool long_join =
        (anchor.flags & ::fa::cpu::chaining::ANCHOR_LONG_JOIN) != 0;
    if (!long_join && anchor_diagonal(anchor) ==
                          anchor_diagonal(path_.selected[region.last])) {
      region.last = next_;
      ++region.anchors;
      continue;
    }
    open_ = {next_, next_, 1, false, 0};
    ++next_;
    more = true;
    break;
  }

  // 2. Verification, gap by gap. Offsets are relative to the region's first
  // center, where its Verified step begins; a failed region keeps no gaps.
  // Bases under a skipped anchor are tested too: that anchor may lie on
  // another diagonal.
  const std::size_t gap_begin = plan_.gaps.size();
  if (region.anchors >= 2) {
    const int origin =
        adjusted_endpoint(path_, path_.selected[region.first]).query;
    std::size_t p = region.first;
    for (std::size_t n = region.first + 1; n <= region.last; ++n) {
      if (skipped(n))
        continue;
      const Anchor& left = path_.selected[p];
      const Anchor& right = path_.selected[n];
      p = n;
      const int length = right.q - left.q_end();
      if (length <= 0)
        continue;
      plan_.gaps.push_back({left.q_end() - origin, right.q - origin});
      if (!certified_gap(certificate_.query + left.q_end(),
                         certificate_.target + left.r_end(), length,
                         certificate_))
        ++region.failing_gaps;
    }
    region.verified = region.failing_gaps == 0;
    if (!region.verified)
      plan_.gaps.resize(gap_begin);
  }
  const std::size_t gap_count = plan_.gaps.size() - gap_begin;

  // 3. The layout.
  if (!region.verified) {
    previous_verified_ = false;
    if (region.anchors < 2) {
      ++pending_.lone_anchors;
    } else {
      ++pending_.failed_regions;
      pending_.failed_anchors += region.anchors;
      pending_.failed_bases +=
          adjusted_endpoint(path_, path_.selected[region.last]).query -
          adjusted_endpoint(path_, path_.selected[region.first]).query;
      pending_.failing_gaps += region.failing_gaps;
    }
  } else {
    if (region.first != corner_) {
      // A neighbour has no retained anchor between the two regions.
      if (previous_verified_)
        plan_.steps.push_back(
            fill(GeometryStepKind::Seam, corner_, region.first));
      else
        emit_stretch(corner_, region.first, pending_);
      pending_ = GeometryStretch{};
    }
    const AdjustedEndpoint lo =
        adjusted_endpoint(path_, path_.selected[region.first]);
    const AdjustedEndpoint hi =
        adjusted_endpoint(path_, path_.selected[region.last]);
    GeometryStep step;
    step.kind = GeometryStepKind::Verified;
    step.selected_index = region.last;
    step.query_begin = lo.query;
    step.query_end = hi.query;
    step.target_begin = lo.target;
    step.target_end = hi.target;
    step.selected_band = control_.normal_long_band;
    step.gap_begin = gap_begin;
    step.gap_count = gap_count;
    step.anchor_count = region.anchors;
    plan_.steps.push_back(step);
    corner_ = region.last;
    previous_verified_ = true;
  }
  if (!more) {
    if (corner_ != last_)
      emit_stretch(corner_, last_, pending_);
    done_ = true;
  }
  return true;
}

TerminalWindowPlan plan_terminal_windows(const OrderedAnchorPath& path,
                                         const TerminalWindowControl& control) {
  TerminalWindowPlan plan;
  if (control.max_gap <= 0 || control.match <= 0 || control.gap_open < 0 ||
      control.gap_extend <= 0) {
    plan.validation = {ValidationCode::InvalidControl, 0};
    return plan;
  }

  const auto& first_anchor = path.selected[path.realization_begin];
  const auto& last_anchor = path.selected[path.realization_end - 1];
  plan.first_adjusted = adjusted_endpoint(path, first_anchor);
  plan.last_adjusted = adjusted_endpoint(path, last_anchor);

  // Left window, as minimap2's mm_align1. The raw bound is the first seed's
  // exact start; the adjusted endpoint is its center.
  int target_start0 = path.selected.front().r;
  int query_start0 = path.selected.front().q;
  int target_start1 = path.target_bound_begin;
  int query_start1 = path.query_bound_begin;
  if (control.left_cap) {
    target_start1 = std::max(target_start1, control.left_cap_target);
    query_start1 = std::max(query_start1, control.left_cap_query);
  }
  const int adjusted_query_start = plan.first_adjusted.query;
  const int adjusted_target_start = plan.first_adjusted.target;
  if (adjusted_query_start > path.query_bound_begin &&
      adjusted_target_start > path.target_bound_begin) {
    int query_span = std::min(adjusted_query_start - path.query_bound_begin,
                              control.max_gap);
    query_start1 = std::max(query_start1, adjusted_query_start - query_span);
    query_start0 = std::min(query_start0, query_start1);
    int target_span = score_expanded_reference_span(query_span, control);
    target_span =
        std::min(target_span, adjusted_target_start - path.target_bound_begin);
    target_start1 =
        std::max(target_start1, adjusted_target_start - target_span);
    target_start0 = std::min(target_start0, target_start1);
    target_start0 = std::min(target_start0, adjusted_target_start);
  } else {
    target_start0 = adjusted_target_start;
    query_start0 = adjusted_query_start;
  }

  // Right window; the raw bound is the last seed's exact k-mer end.
  int target_end0 = path.selected.back().r + path.selected.back().span;
  int query_end0 = path.selected.back().q + path.selected.back().span;
  int target_end1 = path.target_bound_end;
  int query_end1 = path.query_bound_end;
  if (control.right_cap) {
    target_end1 = std::min(target_end1, control.right_cap_target);
    query_end1 = std::min(query_end1, control.right_cap_query);
  }
  const int adjusted_query_end = plan.last_adjusted.query;
  const int adjusted_target_end = plan.last_adjusted.target;
  if (adjusted_query_end < path.query_bound_end &&
      adjusted_target_end < path.target_bound_end) {
    int query_span =
        std::min(path.query_bound_end - adjusted_query_end, control.max_gap);
    query_end1 = std::min(query_end1, adjusted_query_end + query_span);
    query_end0 = std::max(query_end0, query_end1);
    int target_span = score_expanded_reference_span(query_span, control);
    target_span =
        std::min(target_span, path.target_bound_end - adjusted_target_end);
    target_end1 = std::min(target_end1, adjusted_target_end + target_span);
    target_end0 = std::max(target_end0, target_end1);
  } else {
    target_end0 = adjusted_target_end;
    query_end0 = adjusted_query_end;
  }

  plan.query_start = query_start0;
  plan.target_start = target_start0;
  plan.query_end = query_end0;
  plan.target_end = target_end0;
  return plan;
}

} // namespace ordered_anchor
} // namespace lr
} // namespace cpu
} // namespace fa
