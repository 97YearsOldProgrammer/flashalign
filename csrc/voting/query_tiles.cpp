#include "query_tiles.h"

#include <algorithm>

namespace fa::cpu::voting {

namespace {
int word_count(std::uint64_t word) noexcept {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(word);
#else
  int count = 0;
  for (; word != 0; word &= word - 1)
    ++count;
  return count;
#endif
}
} // namespace

void QueryTileMask::set(int tile) {
  if (tile < 0 || tile >= kMaxQueryTiles)
    return;
  const std::size_t index = static_cast<std::size_t>(tile / 64);
  if (index < words.size()) {
    words[index] |= std::uint64_t{1} << (tile % 64);
  } else {
    if (extra.size() <= index - words.size())
      extra.resize(index - words.size() + 1);
    extra[index - words.size()] |= std::uint64_t{1} << (tile % 64);
  }
}

bool QueryTileMask::test(int tile) const noexcept {
  if (tile < 0 || tile >= kMaxQueryTiles)
    return false;
  const std::size_t index = static_cast<std::size_t>(tile / 64);
  const std::uint64_t word =
      index < words.size()
          ? words[index]
          : (index - words.size() < extra.size() ? extra[index - words.size()]
                                                 : 0);
  return (word & (std::uint64_t{1} << (tile % 64))) != 0;
}

int QueryTileMask::count() const noexcept {
  int count = word_count(words[0]) + word_count(words[1]);
  for (const auto word : extra)
    count += word_count(word);
  return count;
}

bool QueryTileMask::empty() const noexcept {
  if (words[0] != 0 || words[1] != 0)
    return false;
  for (const auto word : extra)
    if (word != 0)
      return false;
  return true;
}

bool operator==(const QueryTileMask& left,
                const QueryTileMask& right) noexcept {
  if (left.words != right.words)
    return false;
  const std::size_t common = std::min(left.extra.size(), right.extra.size());
  for (std::size_t i = 0; i < common; ++i)
    if (left.extra[i] != right.extra[i])
      return false;
  for (std::size_t i = common; i < left.extra.size(); ++i)
    if (left.extra[i] != 0)
      return false;
  for (std::size_t i = common; i < right.extra.size(); ++i)
    if (right.extra[i] != 0)
      return false;
  return true;
}

QueryTileMask operator&(const QueryTileMask& left, const QueryTileMask& right) {
  QueryTileMask result{
      {left.words[0] & right.words[0], left.words[1] & right.words[1]}, {}};
  result.extra.resize(std::min(left.extra.size(), right.extra.size()));
  for (std::size_t i = 0; i < result.extra.size(); ++i)
    result.extra[i] = left.extra[i] & right.extra[i];
  return result;
}

QueryTileMask operator|(const QueryTileMask& left, const QueryTileMask& right) {
  QueryTileMask result{
      {left.words[0] | right.words[0], left.words[1] | right.words[1]}, {}};
  result.extra.resize(std::max(left.extra.size(), right.extra.size()));
  for (std::size_t i = 0; i < result.extra.size(); ++i)
    result.extra[i] = (i < left.extra.size() ? left.extra[i] : 0) |
                      (i < right.extra.size() ? right.extra[i] : 0);
  return result;
}

QueryTileMask valid_query_tiles(int tile_count) {
  QueryTileMask result;
  for (int tile = 0; tile < std::clamp(tile_count, 0, kMaxQueryTiles); ++tile)
    result.set(tile);
  return result;
}

int query_tile_for_position(int query_position, int query_length,
                            int seed_length, int tile_count) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int position = std::clamp(query_position, 0, span - 1);
  return std::min(
      tile_count - 1,
      static_cast<int>((static_cast<std::int64_t>(position) * tile_count) /
                       span));
}

int query_tile_begin(int tile, int query_length, int seed_length,
                     int tile_count) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int bounded = std::clamp(tile, 0, tile_count);
  return static_cast<int>((static_cast<std::int64_t>(bounded) * span) /
                          tile_count);
}

int query_tile_end(int tile_end, int query_length, int seed_length,
                   int tile_count) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int bounded = std::clamp(tile_end, 0, tile_count);
  const int exclusive = static_cast<int>(
      (static_cast<std::int64_t>(bounded) * span + tile_count - 1) /
      tile_count);
  return std::min(query_length,
                  exclusive + std::max(0, seed_length - 1));
}

}  // namespace fa::cpu::voting
