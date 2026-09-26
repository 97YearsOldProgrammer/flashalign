// Query intervals and strands for split mapping. An interval passed between stages is in
// forward-read coordinates; the strand frame is used only inside realization.
#pragma once

#include <cstdint>

namespace fa {
namespace cpu {
namespace split {

// Strand of a placement: Forward is the read as given, Reverse its reverse complement.
enum class Strand : uint8_t { Forward = 0, Reverse = 1 };

// Conversions to and from the `bool is_reverse` form.
constexpr Strand strand_from_is_reverse(bool is_reverse) noexcept {
  return is_reverse ? Strand::Reverse : Strand::Forward;
}
constexpr bool is_reverse_strand(Strand s) noexcept {
  return s == Strand::Reverse;
}

// Half-open query interval [begin, end).
struct QueryInterval {
  int begin = 0;
  int end = 0;

  // 0 for an empty or inverted interval.
  constexpr int length() const noexcept { return end > begin ? end - begin : 0; }
  constexpr bool valid() const noexcept { return end > begin; }

  friend constexpr bool operator==(const QueryInterval &a,
                                   const QueryInterval &b) noexcept {
    return a.begin == b.begin && a.end == b.end;
  }
  friend constexpr bool operator!=(const QueryInterval &a,
                                   const QueryInterval &b) noexcept {
    return !(a == b);
  }
};

// [b, e) on a read of length L maps to [L - e, L - b). Its own inverse.
constexpr QueryInterval reverse_complement_interval(QueryInterval qi,
                                                    int read_len) noexcept {
  return QueryInterval{read_len - qi.end, read_len - qi.begin};
}

// Forward-read interval to the strand frame used for chaining and realization.
constexpr QueryInterval to_strand_frame(QueryInterval forward, int read_len,
                                        Strand s) noexcept {
  return s == Strand::Reverse ? reverse_complement_interval(forward, read_len)
                              : forward;
}

// Strand-frame interval back to forward-read coordinates.
constexpr QueryInterval to_forward(QueryInterval strand_frame, int read_len,
                                   Strand s) noexcept {
  return s == Strand::Reverse
             ? reverse_complement_interval(strand_frame, read_len)
             : strand_frame;
}

} // namespace split
} // namespace cpu
} // namespace fa
