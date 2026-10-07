// The -c lane's controller in minimap2's form, from the final ownership
// selection to records. Every owner and the attached
// chains mm_select_sub keeps are regions; each is realized from its whole
// anchor vector as mm_align1 realizes one (both ends trimmed, terminal windows
// to the read ends capped at nearby seeds, a Z-drop splitting off the rest as
// a region of its own), the records are filtered (mm_filter_regs,
// mm_update_dp_max) and roles are assigned again on the aligned spans
// (mm_set_parent, mm_select_sub). Records may overlap.
#pragma once

#include "family_realization.h"

#include <vector>

namespace fa::cpu::lr {

// What the MAPQ reads for one primary record: the SAM primary first, then each
// supplementary in order.
struct DnaRegionPrice {
  // The catalogue candidate of its chain; -1 for an inversion middle.
  int candidate = -1;
  bool inversion = false;
  // minimap2's score and cnt (the region's share of its chain after Z-drop
  // splits, and its anchors), and score0 with the whole chain's anchors.
  int score = 0;
  int cnt = 0;
  int score0 = 0;
  int anchors0 = 0;
  // The whole chain's forward-query and reference spans.
  int chain_q_begin = -1;
  int chain_q_end = -1;
  int chain_ref_begin = -1;
  int chain_ref_end = -1;
  int dp_max = 0;
  int dp_max2 = 0;
  // The catalogue candidate of the record dp_max2 came from, -1 for none.
  int dp_max2_candidate = -1;
  // dp_max and dp_max2 before the rank rescale: the sweep scores (ms:i) of
  // the record and of the record dp_max2 came from.
  int dp_max0 = 0;
  int dp_max2_0 = 0;
  // mm_set_parent's subsc and n_sub, and the same two over the secondaries
  // from its own candidate's chains.
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct DnaRegionOutcome {
  // `output` holds the SAM primary, the other primaries as supplementaries
  // (supplementary_candidates parallel, no parts) and each secondary as a
  // hypothesis of its own. When no record survives it is unmapped, with
  // `failure` and `failed_block` naming the first region refused, or `floored`
  // set when records were realized and the filter dropped them all.
  // block_count counts the regions, segment_count the segments realized and
  // unit_count the records before the filter.
  DnaFamilyRealizationOutcome family;
  std::vector<DnaRegionPrice> prices;
  bool floored = false;
};

DnaRegionOutcome realize_dna_regions(const DnaContext& context,
                                     const DnaFamilyRealizationRequest& request);

}  // namespace fa::cpu::lr
