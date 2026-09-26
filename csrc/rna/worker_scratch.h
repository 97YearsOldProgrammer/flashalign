// Caller-owned per-worker RNA mapping scratch, one per worker thread (minimap2's
// mm_tbuf_t). Buffers are reused across reads.
#pragma once

#include "../seeding/scratch.h"         // ChainWindowPeakScratch
#include "anchoring/skeleton_harvest.h" // SeedPostingMemo (per-read posting memo)
#include "placement/coarse_chain.h"  // CoarseLocusScratch (envelope selection)
#include "query_partition.h"         // RnaExplainScratch (the query partition)
#include "placement/fused_capture.h" // FusedCaptureScratch (fused vote capture)
#include "realization/splice_scratch.h" // SpliceRealizationScratch (DP arena)

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

struct WorkerScratch {
  ChainWindowPeakScratch fwd;
  ChainWindowPeakScratch rc;
  placement::FusedCaptureScratch capture;
  placement::CoarseLocusScratch coarse;
  // The query partition's buffers: the solver's problem, the re-vote capture, the
  // staircase DP arrays.
  RnaExplainScratch explain;
  // Per-read posting memo shared by both strands (the reverse projection keeps canonical
  // keys and indexes it through n-1-j). Scoped to a read by new_read()'s generation bump:
  // the stream identity check cannot tell two reads apart when the allocator recycles a
  // block.
  SeedPostingMemo postings;
  // Realization's DP arena and packet buffers, shared by rank 1, every rival and both
  // transcript hypotheses, which run one after another on this thread.
  SpliceRealizationScratch realization;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
