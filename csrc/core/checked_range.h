// Checked conversions for `int` sequence coordinates and counts. The supported domain is
// [0, INT_MAX], and a half-open end may equal INT_MAX. Check values while still 64-bit;
// narrowing first and validating afterwards is wrong.
#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

namespace fa {
namespace cpu {
namespace mapping {

constexpr std::int64_t kSignedCoordinateMinimum = 0;
constexpr std::int64_t kSignedCoordinateMaximum =
    static_cast<std::int64_t>(std::numeric_limits<int>::max());

enum class CoordinateDomainRefusal : std::uint8_t {
  None = 0,
  UnrepresentableCount,
  NegativeCoordinate,
  InvalidRange,
  EndpointOverflow,
  ProductOverflow,
  ReferenceShapeMismatch,
};

inline bool checked_size_to_int(std::size_t value, int& converted) noexcept {
  if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return false;
  converted = static_cast<int>(value);
  return true;
}

inline bool checked_nonnegative_to_int(std::int64_t value,
                                       int& converted) noexcept {
  if (value < kSignedCoordinateMinimum ||
      value > kSignedCoordinateMaximum)
    return false;
  converted = static_cast<int>(value);
  return true;
}

inline bool checked_half_open_range(std::int64_t begin, std::int64_t end,
                                    std::uint64_t length) noexcept {
  return begin >= 0 && end >= begin &&
      static_cast<std::uint64_t>(end) <= length;
}

inline bool checked_count_range(std::int64_t begin, std::uint64_t count,
                                std::uint64_t size) noexcept {
  if (begin < 0) return false;
  const std::uint64_t unsigned_begin = static_cast<std::uint64_t>(begin);
  return unsigned_begin <= size && count <= size - unsigned_begin;
}

inline bool checked_signed_count_range(std::int64_t begin,
                                       std::uint64_t count,
                                       std::uint64_t size) noexcept {
  if (size > static_cast<std::uint64_t>(kSignedCoordinateMaximum) ||
      count > static_cast<std::uint64_t>(kSignedCoordinateMaximum))
    return false;
  return checked_count_range(begin, count, size);
}

inline bool checked_add(std::int64_t left, std::int64_t right,
                        std::int64_t& sum) noexcept {
  const std::int64_t maximum = std::numeric_limits<std::int64_t>::max();
  const std::int64_t minimum = std::numeric_limits<std::int64_t>::min();
  if ((right > 0 && left > maximum - right) ||
      (right < 0 && left < minimum - right))
    return false;
  sum = left + right;
  return true;
}

inline bool checked_signed_endpoint(std::int64_t begin, std::int64_t span,
                                    std::int64_t length,
                                    int& endpoint) noexcept {
  if (begin < kSignedCoordinateMinimum || span < 0 ||
      length < kSignedCoordinateMinimum ||
      begin > kSignedCoordinateMaximum ||
      length > kSignedCoordinateMaximum)
    return false;
  if (span > kSignedCoordinateMaximum - begin) return false;
  const std::int64_t widened_endpoint = begin + span;
  if (widened_endpoint > length) return false;
  endpoint = static_cast<int>(widened_endpoint);
  return true;
}

inline bool checked_u64_product(std::uint64_t left, std::uint64_t right,
                                std::uint64_t& product) noexcept {
  const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
  if (left != 0 && right > maximum / left) return false;
  product = left * right;
  return true;
}

}  // namespace mapping
}  // namespace cpu
}  // namespace fa
