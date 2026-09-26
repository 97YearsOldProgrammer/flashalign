#pragma once

#include <cstdint>

namespace fa {
namespace cpu {
namespace chaining {

// minimap2 single-segment DNA chaining defaults (map-ont).
struct ColinearChainParams {
  int32_t max_dist_x = 5000;
  int32_t max_dist_y = 5000;
  int32_t bw = 500;
  int32_t max_skip = 25;
  int32_t max_iter = 5000;
  int32_t min_cnt = 3;
  int32_t min_sc = 40;
  float chn_pen_gap = 0.168f;
  float chn_pen_skip = 0.0f;
};

// Spliced (cDNA) chaining, as minimap2's: forward reference gaps pay the cheaper of the
// linear and log penalties, spans are uncapped, max_dist_y is not floored to bw, and the
// chain Z-drop is off.
struct SplicedChainParams {
  int32_t max_dist_x = 5000;
  int32_t max_dist_y = 5000;
  int32_t bw = 500;
  int32_t max_skip = 25;
  int32_t max_iter = 5000;
  int32_t min_cnt = 3;
  int32_t min_sc = 40;
  float chn_pen_gap = 0.168f;
  float chn_pen_skip = 0.0f;
};

}  // namespace chaining
}  // namespace cpu
}  // namespace fa
