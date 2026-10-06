// Dependency bundle handed to RnaBackend::map_read (the RNA counterpart of dna/context.h):
// reference views plus an options snapshot composed once per configure.
#pragma once

#include "annotation/known_junctions.h"
#include "../core/checked_range.h"           // CoordinateDomainRefusal
#include "../core/cigar.h"                   // output::CigarReplayRequest (cs/MD)
#include "../index/reference_context.h"      // engine::ReferenceContext
#include "chain_mapq.h"                      // kRnaChainMapqQcovTau (compat)
#include "config.h"                          // rna::RnaConfig
#include "realization/splice_controller.h"   // SpliceControllerOptions
#include "realization/rival_lifecycle_config.h" // RnaRivalLifecycleConfig
#include "../voting/query_partition.h"          // QueryPartitionParameters

#include <memory>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

// Execution snapshot composed once from the engine's resolved options.
struct ResolvedOptions {
  // The RNA config, with the entry-point clamps applied.
  rna::RnaConfig cfg;
  // The resolved splice-controller scoring row (defaults: minimap2's `splice`), used
  // unchanged by every hypothesis and continuation.
  SpliceControllerOptions splice_controller =
      ordinary_long_read_splice_options();

  // Query-tile objective of the RNA query partition (query_partition.h), its only RNA
  // reader. RNA keeps the struct's defaults.
  ::fa::cpu::voting::QueryPartitionParameters query_partition;

  // Rank-1 harvest and exact-anchor-path chain options: the fine chain's query gap, which
  // -g sets, and its band, -r's first value.
  int cigar_local_global_occ = 200;
  int cigar_local_interval_anchor_chain_max_gap = 5000;
  int fine_chain_band = 200000;

  // Rival realization lifecycle options; the lifecycle is the requested-CIGAR path's
  // only arbitration.
  RnaRivalLifecycleConfig rival_lifecycle;
  // The MAPQ formula's two preset calibration constants (chain_mapq.h), copied onto
  // every RnaChainMapqEvidence by backend.cpp. The defaults reproduce the formula when
  // nothing was composed.
  double mapq_qcov_tau = kRnaChainMapqQcovTau;
  double mapq_disjoint_vote_damp = 0.0;
  // Gates the placement-only return versus splice realization in backend.cpp.
  bool enable_full_read_cigar = true;
  // cs:Z / MD:Z request (minimap2 --cs / --MD), copied onto every realization request.
  ::fa::cpu::output::CigarReplayRequest cigar_replay_request;
  std::shared_ptr<const KnownJunctionStore> known_junctions;
};

// Non-owning views plus the options snapshot, composed once per config; the engine owns
// the storage.
struct Context {
  engine::ReferenceContext ref;
  ResolvedOptions opts;
  // Computed once when the RNA runtime is composed. Per-read code refuses a reference
  // whose counts or contig lengths do not fit the signed coordinate domain.
  mapping::CoordinateDomainRefusal coordinate_domain_refusal =
      mapping::CoordinateDomainRefusal::None;
};

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
