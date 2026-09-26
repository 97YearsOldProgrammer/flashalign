// The RNA spliced backend. The shared engine dispatches each read as
// Backend::map_read(...); this header must not include engine/aligner.h.
#pragma once

#ifdef FLASHALIGN_BUILDING_DNA
#error "flashalign_dna may not include the RNA backend"
#endif

#include "result.h"
#include "../index/seed.h"      // QuerySeed
#include "../seeding/context.h" // LongReadSeedContext
#include "context.h"            // rna::Context
#include "worker_scratch.h"     // rna::WorkerScratch

#include <cstdint>
#include <string>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Maps one read: vote capture, exon peaks, coarse locus catalogue, exact anchor path,
// then map-only projection or splice realization with the rival lifecycle. rctx and
// seed_ctx are composed by the engine (compose_rna_runtime); worker_scratch is the
// caller's per-worker scratch. `read_name_hash` is tie_name_hash of the read name (0
// without one, as minimap2 without a qname); with the read length it seeds the tie-break
// between two equally scored loci.
struct RnaBackend {
  using WorkerScratch = ::fa::cpu::lr::rna::WorkerScratch;

  static AlignResult map_read(const Context &rctx,
                           const LongReadSeedContext &seed_ctx,
                           WorkerScratch &worker_scratch,
                           const std::string &read,
                           const std::vector<uint8_t> &fwd_enc,
                           std::vector<uint8_t> &rc_enc,
                           const std::vector<QuerySeed> *shared_fwd_syncmer_seeds,
                           std::uint32_t read_name_hash);
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
