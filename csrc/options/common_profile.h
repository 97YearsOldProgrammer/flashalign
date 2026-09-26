// Options shared by DNA and RNA mapping.
#pragma once

namespace fa {
namespace cpu {
namespace options {

// Seeding parameters. from_index means they were read from a prebuilt index;
// otherwise they come from the preset.
struct IndexIdentity {
  bool from_index = false;
  int k = 15;
  int syncmer_s = 9;
  int syncmer_downsample = 1;
};

// cs:Z form (minimap2 --cs[=short|long]); None emits no cs:Z.
enum class CsMode { None, Short, Long };

struct CommonOptions {
  int num_threads = 1;
  int min_support = 3;
  int max_query_seeds_per_strand = 128;
  bool enable_full_read_cigar = true;
  // minimap2 --cs / --MD; both imply a realized CIGAR.
  CsMode cs = CsMode::None;
  bool emit_md = false;
  // =/X CIGAR operators instead of M (minimap2 --eqx).
  bool emit_eqx = false;
};

} // namespace options
} // namespace cpu
} // namespace fa
