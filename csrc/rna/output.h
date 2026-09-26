// Commits the selected realized RNA family to the output record.
#pragma once

#include "realization/splice_realizer.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa { namespace cpu { namespace lr { namespace rna {

// A second family emitted beside the committed one: a realized co-primary that cleared
// the chimeric emission floor (realization/rival_lifecycle.h), with its own chain MAPQ.
// `result` is borrowed from the caller's realization store and must outlive the commit.
struct RnaChimericFamily {
  const RnaSpliceRealizationResult* result = nullptr;
  int mapq = 0;
  // This family's own s2:i (-1 when none): the best competing chain score its MAPQ was
  // weighed against.
  int secondary_chain_score = -1;
};

// `secondary_chain_score` is the committed family's s2:i (-1 when none), stamped on the
// same records as `chain_mapq`: the head and its continuation segments.
bool commit_realized_alignment(
    const RnaSpliceRealizationResult& realized,
    int read_len,
    const int* chain_mapq,
    int secondary_chain_score,
    const std::vector<RnaChimericFamily>& chimeric,
    Result& out);

// Attaches one runner-up family to `out`'s alternative hypotheses, moving its segments
// out. Call only after commit_realized_alignment, which overwrites `out`.
void attach_runner_up_family(RnaSpliceRealizationResult& runner_up,
                             AlignResult& out);

}}}}  // namespace fa::cpu::lr::rna
