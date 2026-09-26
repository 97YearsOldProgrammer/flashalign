// A gap inside a verified region: the bases between one anchor's k-mer and
// the next one's, which passed the ungapped certificate. Offsets are relative
// to the region start and apply to both axes, since a verified region lies on
// one diagonal.
#pragma once

namespace fa::cpu::lr {

struct RegionGap {
  int begin = 0;
  int end = 0;
};

} // namespace fa::cpu::lr
