#include "query_tiles.h"

#include <algorithm>

namespace fa::cpu::voting {

void QueryTileMask::set(int tile) noexcept {
  if (tile < 0 || tile >= kFixedQueryTiles) return;
  words[static_cast<std::size_t>(tile / 64)] |=
      std::uint64_t{1} << (tile % 64);
}

bool QueryTileMask::test(int tile) const noexcept {
  return tile >= 0 && tile < kFixedQueryTiles &&
         (words[static_cast<std::size_t>(tile / 64)] &
          (std::uint64_t{1} << (tile % 64))) != 0;
}

int QueryTileMask::count() const noexcept {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(words[0]) + __builtin_popcountll(words[1]);
#else
  int result = 0;
  for (std::uint64_t word : words) {
    while (word != 0) {
      word &= word - 1;
      ++result;
    }
  }
  return result;
#endif
}

QueryTileMask valid_query_tiles(int tile_count) noexcept {
  QueryTileMask result;
  for (int tile = 0; tile < std::clamp(tile_count, 0, kFixedQueryTiles);
       ++tile) {
    result.set(tile);
  }
  return result;
}

int query_tile_for_position(int query_position, int query_length,
                            int seed_length) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int position = std::clamp(query_position, 0, span - 1);
  return std::min(kFixedQueryTiles - 1,
                  static_cast<int>((static_cast<std::int64_t>(position) *
                                    kFixedQueryTiles) /
                                   span));
}

int query_tile_begin(int tile, int query_length, int seed_length) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int bounded = std::clamp(tile, 0, kFixedQueryTiles);
  return static_cast<int>(
      (static_cast<std::int64_t>(bounded) * span) / kFixedQueryTiles);
}

int query_tile_end(int tile_end, int query_length, int seed_length) noexcept {
  const int span = std::max(1, query_length - std::max(1, seed_length) + 1);
  const int bounded = std::clamp(tile_end, 0, kFixedQueryTiles);
  const int exclusive = static_cast<int>(
      (static_cast<std::int64_t>(bounded) * span +
       kFixedQueryTiles - 1) /
      kFixedQueryTiles);
  return std::min(query_length,
                  exclusive + std::max(0, seed_length - 1));
}

}  // namespace fa::cpu::voting
