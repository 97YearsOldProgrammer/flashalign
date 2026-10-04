#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace fa::cpu::voting {

inline constexpr int kFixedQueryTiles = 128;
// The fewest tiles that can place a read: an owned block needs two
// supported tiles (QueryPartitionParameters).
inline constexpr int kMinQueryTiles = 2;
// The largest tile count a mask holds (--tiles).
inline constexpr int kMaxQueryTiles = 4096;
// The finest grid on which DNA ratio admission compares masks; a read
// partitioned into more tiles is admitted on this many.
inline constexpr int kMaxAdmissionQueryTiles = 2048;

struct QueryTileMask {
  std::array<std::uint64_t, 2> words{};
  // Tiles from 128 up; a mask of at most 128 tiles allocates nothing.
  std::vector<std::uint64_t> extra;

  void set(int tile);
  bool test(int tile) const noexcept;
  int count() const noexcept;
  bool empty() const noexcept;

  friend bool operator==(const QueryTileMask& left,
                         const QueryTileMask& right) noexcept;
  friend QueryTileMask operator&(const QueryTileMask& left,
                                 const QueryTileMask& right);
  friend QueryTileMask operator|(const QueryTileMask& left,
                                 const QueryTileMask& right);
};

QueryTileMask valid_query_tiles(int tile_count);
int query_tile_for_position(int query_position, int query_length,
                            int seed_length,
                            int tile_count = kFixedQueryTiles) noexcept;
int query_tile_begin(int tile, int query_length, int seed_length,
                     int tile_count = kFixedQueryTiles) noexcept;
int query_tile_end(int tile_end, int query_length, int seed_length,
                   int tile_count = kFixedQueryTiles) noexcept;

}  // namespace fa::cpu::voting
