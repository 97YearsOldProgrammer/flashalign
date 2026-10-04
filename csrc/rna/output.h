// Commits the selected realized RNA family to the output record.
#pragma once

#include "realization/splice_realizer.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa { namespace cpu { namespace lr { namespace rna {

// A family's chain tags, each -1 when none: cm:i and s1:i, the anchors and score of its
// own chain, and s2:i, the best competing chain score its MAPQ was weighed against.
struct RnaChainTags {
  int chain_anchors = -1;
  int chain_score = -1;
  int secondary_chain_score = -1;
};

// A second family emitted beside the committed one: a realized co-primary that cleared
// the chimeric emission floor (realization/rival_lifecycle.h), with its own chain MAPQ.
// `result` is borrowed from the caller's realization store and must outlive the commit.
struct RnaChimericFamily {
  const RnaSpliceRealizationResult* result = nullptr;
  int mapq = 0;
  // This family's own chain tags.
  RnaChainTags chain_tags;
};

// `chain_tags` are the committed family's, stamped on the same records as `chain_mapq`:
// the head and its continuation segments.
bool commit_realized_alignment(
    const RnaSpliceRealizationResult& realized,
    int read_len,
    const int* chain_mapq,
    const RnaChainTags& chain_tags,
    const std::vector<RnaChimericFamily>& chimeric,
    Result& out);

// Attaches one runner-up family to `out`'s alternative hypotheses, moving its segments
// out. Call only after commit_realized_alignment, which overwrites `out`. `chain_tags`,
// the family's own chain (s2:i left -1), go on every record of it.
void attach_runner_up_family(RnaSpliceRealizationResult& runner_up,
                             const RnaChainTags& chain_tags,
                             AlignResult& out);

}}}}  // namespace fa::cpu::lr::rna
