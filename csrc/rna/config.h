// Runtime configuration of the RNA pipeline, built once per configure by
// LongReadEngine::compose_rna_runtime and read by RnaBackend::map_read.
#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include RNA configuration"
#endif

#include <cstdint>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Transcript-strand policy, read-relative like minimap2's -u: Unknown runs both
// hypotheses and elects one (-ub); Forward treats every read as the sense strand (-uf),
// Reverse as antisense (-ur). It never fixes the mapping strand: seeding, voting and
// chaining stay two-stranded and only the second realization pass is skipped.
// None is minimap2's -un: no splice motif is scored, one hypothesis runs and no ts:A tag
// is emitted.
enum class StrandMode : uint8_t {
  Unknown = 0,
  Forward = 1,
  Reverse = 2,
  None = 3,
};

StrandMode parse_strand_mode(const char* text, StrandMode fallback);

struct RnaConfig {
  // Seeding and intron geometry, mirrored from the shared long-read config.
  int k = 15;
  int min_support = 2;
  int min_intron = 20;
  int max_intron = 200000;
  int max_locus_chains = 6;
  int max_chain_predecessors = 64;

  // Unknown runs minimap2's two transcript hypotheses; a forward/reverse hint runs one.
  StrandMode strand_mode = StrandMode::Unknown;
};

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
