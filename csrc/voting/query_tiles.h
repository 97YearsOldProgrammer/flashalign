#pragma once

#include <array>
#include <cstdint>

namespace fa::cpu::voting {

inline constexpr int kFixedQueryTiles = 128;

struct QueryTileMask {
  std::array<std::uint64_t, 2> words{};

  void set(int tile) noexcept;
  bool test(int tile) const noexcept;
  int count() const noexcept;
  bool empty() const noexcept { return words[0] == 0 && words[1] == 0; }

  friend bool operator==(const QueryTileMask& left,
                         const QueryTileMask& right) noexcept {
    return left.words == right.words;
  }
  friend QueryTileMask operator&(const QueryTileMask& left,
                                 const QueryTileMask& right) noexcept {
    return {{left.words[0] & right.words[0],
             left.words[1] & right.words[1]}};
  }
  friend QueryTileMask operator|(const QueryTileMask& left,
                                 const QueryTileMask& right) noexcept {
    return {{left.words[0] | right.words[0],
             left.words[1] | right.words[1]}};
  }
};

QueryTileMask valid_query_tiles(int tile_count) noexcept;
int query_tile_for_position(int query_position, int query_length,
                            int seed_length) noexcept;
int query_tile_begin(int tile, int query_length, int seed_length) noexcept;
int query_tile_end(int tile_end, int query_length, int seed_length) noexcept;

}  // namespace fa::cpu::voting
