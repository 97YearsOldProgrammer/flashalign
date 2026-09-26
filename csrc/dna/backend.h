// The DNA mapping backend. Internal header, not public API.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include the DNA backend"
#endif

#include "context.h"            // DnaContext
#include "worker_scratch.h"     // DnaWorkerScratch
#include "result.h"
#include "../index/seed.h"      // QuerySeed
#include "../seeding/context.h" // LongReadSeedContext

#include <cstdint>
#include <string>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

namespace dna {

// Maps one read: whole-read vote, bounded-candidate chaining, projection
// (map-only) or realization (CIGAR), then MAPQ and output assembly.
AlignResult map_read(const DnaContext &dctx, const LongReadSeedContext &seed_ctx,
                  DnaWorkerScratch &worker_scratch, const std::string &read,
                  const std::vector<uint8_t> &fwd_enc,
                  std::vector<uint8_t> &rc_enc,
                  const std::vector<QuerySeed> *shared_fwd_syncmer_seeds);

// Stateless backend tag for the templated engine; the RNA counterpart is
// rna::RnaBackend.
struct DnaBackend {
  using WorkerScratch = DnaWorkerScratch;
  using Context = DnaContext;

  static AlignResult map_read(const DnaContext &dctx,
                           const LongReadSeedContext &seed_ctx,
                           DnaWorkerScratch &worker_scratch,
                           const std::string &read,
                           const std::vector<uint8_t> &fwd_enc,
                           std::vector<uint8_t> &rc_enc,
                           const std::vector<QuerySeed> *shared_fwd_syncmer_seeds) {
    return ::fa::cpu::lr::dna::map_read(dctx, seed_ctx, worker_scratch, read,
                                        fwd_enc, rc_enc,
                                        shared_fwd_syncmer_seeds);
  }
};

} // namespace dna
} // namespace lr
} // namespace cpu
} // namespace fa
