// RNA (spliced) mapping options: the DNA options as `base`, plus intron,
// strand, rival and splice-scoring settings.
#pragma once

#include "dna_profile.h"
#include <string>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct RnaLongOptions {
  // Seeding, occurrence, chaining and DP scoring shared with DNA.
  DnaLongOptions base;
  // --min-intron / -G
  int min_intron = 20;
  int max_intron = 200000;
  int max_locus_chains = 6;
  int max_chain_predecessors = 64;
  // Rival realization (rna/realization/rival_lifecycle.h). A rival is
  // retained when its chain score is within rival_pri_ratio of the primary's,
  // or within rival_min_diff of it absolutely (default 30 = 2k).
  double rival_pri_ratio = 0.8;
  int rival_min_diff = 30;
  // At most this many rival realizations per read; -N sets it to the loci
  // kept minus one.
  int rival_realize_max = 5;
  // Splice-DP terms with no DNA analogue, as in minimap2's splice presets.
  int splice_transition = 0;
  int splice_junction_bonus = 9;
  // Read only by the scored-site mode, which no option selects.
  int splice_junction_penalty = 5;
  std::string junction_bed;
  int splice_inversion_zdrop = 100;
  // Transcript strand as rna::StrandMode: -1 unset, 0 auto, 1 forward,
  // 2 reverse, 3 none (no splice motif scored).
  int strand_mode = -1;
  // RNA MAPQ parameters (rna/chain_mapq.h), set by the preset only.
  //
  // Coverage threshold, in (0,1]: the penalty's coverage term is
  // clip((sel_qcov/read_len) / tau, 0, 1).
  double rna_mapq_qcov_tau = 0.7;
  // In [0,1]: a rival that fails the mask_level test adds
  // min(1, damp * rs_i/rs1) to the ambiguity term; 0 disables this.
  double rna_mapq_disjoint_vote_damp = 0.0;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
