// The long-read engine: owns the index and the reference, holds the resolved options and
// maps reads through the DNA or RNA backend. WindowedAlignSession runs it over a stream of
// batches.
#pragma once

#include "../core/cigar.h"
#include "../core/sequence.h"
#include "../core/types.h"
#include "../dp/kernel.h"
#include "../index/format.h"
#include "../index/index.h"
#include "../index/occ_stats.h"
#include "../index/seed.h"
#include "../core/checked_range.h"
#include "../rna/backend.h"
#include "../rna/context.h"
#include "../rna/worker_scratch.h"
#include "../threading/batch_window.h"
#include "../threading/parallel_for.h"
#include "../voting/vote.h"
#include "../options/resolved_options.h"
#include "../options/presets.h"
#include "../dna/dp_runner.h"
#include "../dna/worker_scratch.h"
#include "../dna/backend.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

// Fills `order` with read indices, longest read first (LPT scheduling with one-item claims).
// Small batches keep the input order.
inline void lpt_read_order(const std::vector<std::string> &reads, int n_threads,
                           std::vector<int> &order) {
  order.resize(reads.size());
  std::iota(order.begin(), order.end(), 0);
  if (reads.size() >= static_cast<size_t>(std::max(8, n_threads * 4))) {
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
      return reads[static_cast<size_t>(a)].size() >
             reads[static_cast<size_t>(b)].size();
    });
  }
}

class LongReadEngine {
public:
  LongReadEngine() = default;

  const ::fa::cpu::options::ResolvedOptions& config() const { return cfg_; }
  void reconfigure(const ::fa::cpu::options::ResolvedOptions& config) {
    install_config(config);
  }
  LongReadEngine(std::unordered_map<std::string, std::string> genome,
                 int k = 15, int min_support = 3, int chain_syncmer_s = 9,
                 int chain_syncmer_downsample = 1, int build_threads = 0) {
    cfg_.index.k = k;
    cfg_.common.min_support = min_support;
    apply_chain_seed_constructor_config(chain_syncmer_s,
                                        chain_syncmer_downsample);
    // Only k and s are needed to build the index; the caller installs the resolved map-time
    // options with reconfigure() afterwards.
    ingest_genome(std::move(genome));
    FaixBuildConfig cfg;
    cfg.k = cfg_.index.k;
    cfg.syncmer_s = effective_chain_syncmer_s();
    // build_threads <= 0 uses every available core.
    cfg.build_threads = build_threads > 0
                            ? build_threads
                            : ::fa::cpu::threading::default_thread_count();
    ::fa::cpu::FaixBuildStatus build_status;
    index_ = std::make_shared<SeedIndex>(
        build_seed_index(chr_encs_, cfg, &build_status));
    if (index_->empty()) {
      throw std::runtime_error(
          "failed to build the index: " +
          ::fa::cpu::faix_build_error_message(build_status));
    }
    cfg_.common.num_threads =
        ::fa::cpu::threading::default_thread_count();
    // Provisional; reconfigure() redoes both with the final options.
    resolve_index_occ_cap_into(cfg_);
    compose_rna_runtime();
  }

  // Loads an engine over a prebuilt .faix. The seed length k and syncmer s come from the
  // index, as in minimap2.
  //
  // `n_threads` bounds the threads used to read the file and unpack the reference; <= 0
  // uses the default. When `reference_bases_needed` is false (DNA map-only PAF) the
  // reference stays packed and only names and lengths are read from the index.
  // `image_offset` / `image_bytes` select one part of a multi-part index; both 0 means the
  // whole file.
  static std::unique_ptr<LongReadEngine>
  from_index(std::unordered_map<std::string, std::string> genome,
             const std::string &index_path, int min_support = 3,
             int chain_syncmer_downsample = 1, int n_threads = 0,
             bool reference_bases_needed = true, uint64_t image_offset = 0,
             uint64_t image_bytes = 0) {
    auto obj = std::make_unique<LongReadEngine>();
    obj->cfg_.common.min_support = min_support;
    // s is set from the index below.
    obj->apply_chain_seed_constructor_config(0, chain_syncmer_downsample);
    FaixLoadStatus load_status;
    obj->index_ = std::make_shared<SeedIndex>(load_seed_index(
        index_path, &load_status, n_threads, image_offset, image_bytes));
    if (obj->index_->empty()) {
      throw std::runtime_error("cannot load index " + index_path + ": " +
                               faix_load_error_message(load_status));
    }
    // Names and lengths come from the index. The bases are unpacked from the index only when
    // needed; an index without them falls back to encoding `genome`.
    bool reference_bound = false;
    if (reference_bases_needed) {
      if (obj->index_->has_reference_payload() &&
          obj->index_->reference_chromosomes_u8(obj->chr_encs_, n_threads)) {
        obj->bind_index_reference_metadata();
        reference_bound = true;
      }
    } else {
      obj->bind_index_reference_metadata();
      reference_bound = true;
    }
    if (!reference_bound) {
      obj->ingest_genome(std::move(genome));
    }
    obj->cfg_.index.k = obj->index_->k();
    obj->cfg_.index.syncmer_s = obj->index_->build_syncmer_s();
    obj->cfg_.index.from_index = true;
    obj->cfg_.common.num_threads =
        ::fa::cpu::threading::default_thread_count();
    obj->resolve_index_occ_cap_into(obj->cfg_);
    obj->compose_rna_runtime();
    return obj;
  }

  // Builds an engine over an already-loaded index, which may be shared by several engines.
  // The index must embed its reference. Options come only from `config`.
  static std::unique_ptr<LongReadEngine>
  from_shared_index(std::shared_ptr<SeedIndex> index,
                    const ::fa::cpu::options::ResolvedOptions& config) {
    if (!index || index->empty()) {
      throw std::runtime_error("the index is empty");
    }
    auto obj = std::make_unique<LongReadEngine>();
    obj->index_ = std::move(index);
    if (!(obj->index_->has_reference_payload() &&
          obj->index_->reference_chromosomes_u8(
              obj->chr_encs_, config.common.num_threads))) {
      throw std::runtime_error(
          "the index does not embed its reference sequence (built with "
          "--idx-no-seq); an Aligner needs it");
    }
    obj->bind_index_reference_metadata();
    // Configure after binding the reference names: loading a junction BED resolves its
    // contigs against them.
    obj->reconfigure(config);
    if (obj->index_->k() != obj->cfg_.index.k) {
      throw std::invalid_argument(
          "k=" + std::to_string(obj->cfg_.index.k) +
          " does not match the index (k=" + std::to_string(obj->index_->k()) +
          ")");
    }
    const int stored_s = obj->index_->build_syncmer_s();
    const int want_s = obj->effective_chain_syncmer_s();
    if (stored_s != want_s) {
      throw std::invalid_argument(
          "s=" + std::to_string(want_s) + " does not match the index (s=" +
          std::to_string(stored_s) + ")");
    }
    return obj;
  }

  // Maps one read with the given Backend (dna::DnaBackend or rna::RnaBackend). This overload
  // uses a fresh scratch; the next one takes a per-worker scratch reused across reads, like
  // minimap2's mm_tbuf_t.
  //
  // `read_name_hash` is tie_name_hash() of the read name and seeds the tie-break between two
  // equally good loci; 0 when there is no name. `self_contig`, where taken, is contig_id()
  // of the read name, read only where the preset leaves a read out of its own vote
  // (DnaLongOptions::skip_self); -1 when it is not known. `name_rank`, where taken, is
  // contig_name_rank() of the read name, read only where the preset prints a pair once
  // (DnaLongOptions::dual); 0 when it is not known.
  template <class Backend>
  AlignResult align(const std::string &read,
                    std::uint32_t read_name_hash = 0) const {
    typename Backend::WorkerScratch worker_scratch;
    return align<Backend>(read, worker_scratch, read_name_hash);
  }

  template <class Backend>
  AlignResult align(const std::string &read,
                    typename Backend::WorkerScratch &worker_scratch,
                    std::uint32_t read_name_hash = 0,
                    int self_contig = -1, int name_rank = 0) const {
    AlignResult out{};
    int read_len = 0;
    if (!::fa::cpu::mapping::checked_size_to_int(read.size(), read_len)) {
      return out;
    }
    out.read_len = read_len;
    if (index_->empty() || read_len < cfg_.index.k) {
      return out;
    }

    std::vector<uint8_t> fwd_enc;
    std::vector<QuerySeed> shared_fwd_syncmer_seeds;
    const std::vector<QuerySeed> *shared_fwd_syncmer_seed_ptr = nullptr;
    ClosedSyncmerConfig stream_cfg;
    stream_cfg.k = cfg_.index.k;
    stream_cfg.s = cfg_.index.syncmer_s;
    stream_cfg.downsample = 1;
    encode_and_extract_closed_syncmer_query_seeds_into(
        read.data(), read_len, stream_cfg, fwd_enc, shared_fwd_syncmer_seeds);
    shared_fwd_syncmer_seed_ptr = &shared_fwd_syncmer_seeds;
    return align_encoded<Backend>(read, fwd_enc, worker_scratch,
                                  shared_fwd_syncmer_seed_ptr, read_name_hash,
                                  self_contig, name_rank);
  }

  template <class Backend>
  AlignResult align_encoded(
      const std::string &read, const std::vector<uint8_t> &fwd_enc,
      typename Backend::WorkerScratch &worker_scratch,
      const std::vector<QuerySeed> *precomputed_shared_fwd_syncmer_seeds =
          nullptr,
      std::uint32_t read_name_hash = 0, int self_contig = -1,
      int name_rank = 0) const {
    AlignResult out{};
    int read_len = 0;
    if (!::fa::cpu::mapping::checked_size_to_int(read.size(), read_len)) {
      return out;
    }
    out.read_len = read_len;
    if (index_->empty() || read_len < cfg_.index.k)
      return out;

    std::vector<QuerySeed> shared_fwd_syncmer_seeds;
    const std::vector<QuerySeed> *shared_fwd_syncmer_seed_ptr = nullptr;
    if (precomputed_shared_fwd_syncmer_seeds) {
      shared_fwd_syncmer_seed_ptr = precomputed_shared_fwd_syncmer_seeds;
    } else {
      ClosedSyncmerConfig stream_cfg;
      stream_cfg.k = cfg_.index.k;
      stream_cfg.s = cfg_.index.syncmer_s;
      stream_cfg.downsample = 1;
      extract_closed_syncmer_query_seeds_into(
          fwd_enc.data(), read_len, stream_cfg, shared_fwd_syncmer_seeds);
      shared_fwd_syncmer_seed_ptr = &shared_fwd_syncmer_seeds;
    }

    // Filled by the backend only when a read needs it.
    std::vector<uint8_t> rc_enc;
    return align_anchor_chain<Backend>(read, fwd_enc, rc_enc, worker_scratch,
                                       shared_fwd_syncmer_seed_ptr,
                                       read_name_hash, self_contig, name_rank);
  }

  // Maps a batch in parallel. `read_name_hashes` is null or holds one hash per read,
  // `self_contigs` is null or holds one contig_id() per read, and `name_ranks` is null or
  // holds one contig_name_rank() per read.
  template <class Backend>
  std::vector<AlignResult>
  align_batch(const std::vector<std::string> &reads,
              const std::vector<std::uint32_t> *read_name_hashes = nullptr,
              const std::vector<int> *self_contigs = nullptr,
              const std::vector<int> *name_ranks = nullptr) const {
    const int n_threads = cfg_.common.num_threads;
    std::vector<AlignResult> out(reads.size());
    std::vector<int> order;
    lpt_read_order(reads, n_threads, order);
    // One scratch per worker, indexed by tid.
    std::vector<typename Backend::WorkerScratch> worker_scratch(
        static_cast<size_t>(std::max(1, n_threads)));
    ::fa::cpu::threading::parallel_for(
        n_threads, static_cast<int64_t>(order.size()),
        [&](int64_t oi, int tid) {
          const int i = order[static_cast<size_t>(oi)];
          out[static_cast<size_t>(i)] =
              align<Backend>(
                  reads[static_cast<size_t>(i)],
                  worker_scratch[static_cast<size_t>(tid)],
                  read_name_hashes != nullptr
                      ? (*read_name_hashes)[static_cast<size_t>(i)]
                      : std::uint32_t{0},
                  self_contigs != nullptr
                      ? (*self_contigs)[static_cast<size_t>(i)]
                      : -1,
                  name_ranks != nullptr
                      ? (*name_ranks)[static_cast<size_t>(i)]
                      : 0);
        });
    return out;
  }

  std::vector<std::string> chromosome_names() const { return chr_names_; }
  // The contig named exactly `name`, or -1. Contigs are numbered in name order.
  int contig_id(std::string_view name) const {
    const auto it = std::lower_bound(
        chr_names_.begin(), chr_names_.end(), name,
        [](const std::string &contig, std::string_view wanted) {
          return contig < wanted;
        });
    return it != chr_names_.end() && *it == name
               ? static_cast<int>(it - chr_names_.begin())
               : -1;
  }
  // The number of contigs whose name sorts before `name`, which are the contigs
  // numbered below it.
  int contig_name_rank(std::string_view name) const {
    return static_cast<int>(
        std::lower_bound(chr_names_.begin(), chr_names_.end(), name,
                         [](const std::string &contig, std::string_view wanted) {
                           return contig < wanted;
                         }) -
        chr_names_.begin());
  }
  // Lengths in chromosome_names() order, available even when the bases were not unpacked.
  std::vector<int64_t> chromosome_lengths() const { return chr_lengths_; }
  int get_index_syncmer_s() const { return index_->build_syncmer_s(); }

private:
  // Sets the DNA vote occurrence cap from the attached index under LongOccPolicy::Platform:
  // max(preset floor, the index's kLongOccQuantileF occurrence quantile as minimap2's
  // mm_idx_cal_max_occ), then at most the preset ceiling. Records how it was derived in
  // long_occ_index_resolution. Does nothing for RNA, for a fixed cap (--max-vote-occ) or
  // without an index. Repeated calls start again from the recorded floor.
  void resolve_index_occ_cap_into(
      ::fa::cpu::options::ResolvedOptions& config) {
    if (config.is_rna()) return;
    auto& mapping = config.long_read();
    if (mapping.long_occ_policy != LongOccPolicy::Platform) return;
    if (!index_ || index_->empty()) return;
    if (!occ_dist_ || occ_dist_source_ != index_.get()) {
      occ_dist_ = std::make_shared<const ::fa::cpu::FaixOccDistribution>(
          ::fa::cpu::faix_occ_distribution(*index_));
      occ_dist_source_ = index_.get();
    }
    const auto& previous = mapping.long_occ_index_resolution;
    LongOccIndexResolution record;
    record.quantile_f = ::fa::cpu::lr::kLongOccQuantileF;
    record.distinct_keys = occ_dist_->n_distinct;
    record.total_postings = occ_dist_->total_postings;
    record.floor_cap = std::max(
        0, previous.resolved ? previous.floor_cap : mapping.long_occ_cap);
    // The ceiling applies after the floor, as minimap2's mm_mapopt_update clamps mid_occ, and
    // not at all when below the floor. A ceiling equal to the floor still applies, so a
    // resolved config fed back through reconfigure() resolves to the same cap.
    record.ceiling_cap = mapping.long_occ_ceiling >= record.floor_cap
                             ? std::max(0, mapping.long_occ_ceiling)
                             : 0;
    if (occ_dist_->n_distinct == 0) {
      record.resolved_cap = record.floor_cap;
      mapping.long_occ_cap = record.resolved_cap;
      mapping.long_occ_index_resolution = record;
      return;
    }
    record.raw_cap =
        ::fa::cpu::faix_cal_max_occ(*occ_dist_, record.quantile_f);
    std::int64_t cap = std::max<std::int64_t>(
        record.raw_cap, static_cast<std::int64_t>(record.floor_cap));
    if (record.ceiling_cap > 0)
      cap = std::min<std::int64_t>(cap, record.ceiling_cap);
    record.resolved_cap = static_cast<int>(cap);
    record.resolved = true;
    mapping.long_occ_cap = record.resolved_cap;
    mapping.long_occ_index_resolution = record;
  }

  static LongOccPolicyConfig
  long_occ_policy_config(const ::fa::cpu::lr::DnaLongOptions& mapping) {
    LongOccPolicyConfig cfg;
    cfg.policy = mapping.long_occ_policy;
    cfg.primary_occ_cap = mapping.long_primary_occ_cap;
    cfg.platform_occ_cap = mapping.long_occ_cap;
    return cfg;
  }
  int effective_chain_syncmer_s() const {
    return effective_closed_syncmer_s(
        std::max(1, cfg_.index.k), cfg_.index.syncmer_s);
  }

  // The seeding parameters the shared seeding code reads.
  LongReadSeedContext
  make_seed_context(
      const ::fa::cpu::options::ResolvedOptions& config) const {
    const auto& mapping = config.long_read();
    LongReadSeedContext ctx;
    ctx.k = config.index.k;
    ctx.index = index_.get();
    ctx.syncmer_s = config.index.syncmer_s;
    ctx.syncmer_downsample = config.index.syncmer_downsample;
    ctx.max_query_seeds_per_strand =
        config.common.max_query_seeds_per_strand;
    ctx.nested_vote_seeds = !config.is_rna();
    ctx.syncmer_occ_aware_enabled = true; // DNA: per-strip single-cap selector
    ctx.occ_policy = long_occ_policy_config(mapping);
    ctx.vote_diag_bin_width = mapping.vote_diag_bin_width;
    ctx.vote_diag_bin_width_adaptive =
        mapping.vote_diag_bin_width_adaptive;
    ctx.vote_diag_slope_den = mapping.vote_diag_slope_den;
    ctx.vote_diag_width_max = mapping.vote_diag_width_max;
    ctx.min_support = config.common.min_support;
    ctx.chain_max_candidates_per_window =
        config.is_rna() ? mapping.chain_max_candidates_per_window
                        : ::fa::cpu::lr::dna_chain_max_candidates(mapping);
    ctx.vote_batched_refine =
        !config.is_rna() &&
        (mapping.vote_admission_ratio > 0.0 || mapping.all_chains);
    ctx.tile_rescue_occ = mapping.dna_tile_rescue_occ;
    ctx.chr_names = &chr_names_;
    return ctx;
  }
  LongReadSeedContext make_seed_context() const {
    return make_seed_context(cfg_);
  }

  // The --cs / --MD / --eqx request passed to both realizers.
  static ::fa::cpu::output::CigarReplayRequest cigar_replay_request_for(
      const ::fa::cpu::options::CommonOptions& common) {
    ::fa::cpu::output::CigarReplayRequest request;
    switch (common.cs) {
    case ::fa::cpu::options::CsMode::Short:
      request.cs = ::fa::cpu::output::CigarReplayRequest::Cs::Short;
      break;
    case ::fa::cpu::options::CsMode::Long:
      request.cs = ::fa::cpu::output::CigarReplayRequest::Cs::Long;
      break;
    case ::fa::cpu::options::CsMode::None:
      request.cs = ::fa::cpu::output::CigarReplayRequest::Cs::None;
      break;
    }
    request.md = common.emit_md;
    request.eqx = common.emit_eqx;
    return request;
  }

  // Builds the RNA backend's context and seed context from `config` and the reference.
  // Unused in DNA mode.
  void compose_rna_runtime(
                           const ::fa::cpu::options::ResolvedOptions& config,
                           rna::Context &next_rna_ctx,
                           LongReadSeedContext &next_rna_seed_ctx) const {
    const auto* rna_mapping = config.rna();
    const auto& mapping = config.long_read();
    rna::RnaConfig rc;
    rc.k = config.index.k;
    rc.min_support = std::max(0, config.common.min_support);
    rc.min_intron =
        std::max(1, rna_mapping ? rna_mapping->min_intron : 20);
    rc.max_intron = std::max(
        rc.min_intron, rna_mapping ? rna_mapping->max_intron : 200000);
    rc.max_locus_chains =
        std::max(1, rna_mapping ? rna_mapping->max_locus_chains : 6);
    rc.max_chain_predecessors =
        std::max(
            1,
            rna_mapping ? rna_mapping->max_chain_predecessors : 64);
    const int strand_mode = rna_mapping ? rna_mapping->strand_mode : -1;
    rc.strand_mode = strand_mode >= 0
                         ? static_cast<rna::StrandMode>(strand_mode)
                         : rna::StrandMode::Unknown;

    rna::SpliceControllerOptions controller =
        rna::ordinary_long_read_splice_options();
    if (rna_mapping) {
      controller.k = config.index.k;
      controller.match = mapping.cigar_dp_match;
      controller.mismatch = mapping.cigar_dp_mismatch;
      controller.gap_open = mapping.cigar_dp_gap_open1;
      controller.gap_extend = mapping.cigar_dp_gap_extend1;
      controller.long_gap_open = mapping.cigar_dp_gap_open2;
      controller.long_gap_extend = mapping.cigar_dp_gap_extend2;
      controller.ambiguous = mapping.cigar_dp_ambi;
      controller.zdrop = mapping.cigar_dp_tail_zdrop;
      controller.end_bonus = mapping.cigar_dp_tail_end_bonus;
      controller.maximum_gap = mapping.cigar_dp_max_gap;
      controller.maximum_reference_gap = rc.max_intron;
      controller.minimum_dp_maximum = mapping.cigar_dp_min_dp_max;
      controller.transition = rna_mapping->splice_transition;
      controller.junction_bonus = rna_mapping->splice_junction_bonus;
      controller.junction_penalty = rna_mapping->splice_junction_penalty;
      controller.inversion_zdrop =
          rna_mapping->splice_inversion_zdrop;
      controller.minimum_intron = rc.min_intron;
      if (!rna::splice_controller_options_supported(controller)) {
        throw std::invalid_argument(
            "these alignment scores are not supported by the splice presets");
      }
    }

    rna::ResolvedOptions opts;
    opts.cfg = rc;
    opts.splice_controller = controller;
    // RNA uses the preset value; DNA takes the resolved vote cap (make_dna_context).
    opts.cigar_local_global_occ = mapping.cigar_local_global_occ;
    opts.cigar_local_interval_anchor_chain_max_gap =
        mapping.cigar_local_interval_anchor_chain_max_gap;
    opts.fine_chain_band = mapping.cigar_dp_bw;
    opts.query_partition = mapping.query_partition;
    if (rna_mapping) {
      opts.rival_lifecycle.pri_ratio = rna_mapping->rival_pri_ratio;
      opts.rival_lifecycle.min_diff = rna_mapping->rival_min_diff;
      opts.rival_lifecycle.realize_max = rna_mapping->rival_realize_max;
      opts.mapq_qcov_tau = rna_mapping->rna_mapq_qcov_tau;
      opts.mapq_disjoint_vote_damp = rna_mapping->rna_mapq_disjoint_vote_damp;
    }
    opts.enable_full_read_cigar = config.common.enable_full_read_cigar;
    opts.cigar_replay_request = cigar_replay_request_for(config.common);
    if (rna_mapping && !rna_mapping->junction_bed.empty()) {
      opts.known_junctions = std::make_shared<const rna::KnownJunctionStore>(
          rna::KnownJunctionStore::load_bed(
              rna_mapping->junction_bed, chr_names_, chr_encs_));
    }

    next_rna_ctx.ref = ::fa::cpu::engine::ReferenceContext{
        index_.get(), &chr_names_, &chr_lengths_, &chr_encs_};
    next_rna_ctx.opts = std::move(opts);
    next_rna_ctx.coordinate_domain_refusal =
        ::fa::cpu::mapping::CoordinateDomainRefusal::None;
    int checked_reference_count = 0;
    if (chr_names_.size() != chr_encs_.size()) {
      next_rna_ctx.coordinate_domain_refusal =
          ::fa::cpu::mapping::CoordinateDomainRefusal::ReferenceShapeMismatch;
    } else if (!::fa::cpu::mapping::checked_size_to_int(
                   chr_encs_.size(), checked_reference_count)) {
      next_rna_ctx.coordinate_domain_refusal =
          ::fa::cpu::mapping::CoordinateDomainRefusal::UnrepresentableCount;
    } else {
      for (const auto& encoded_reference : chr_encs_) {
        int checked_reference_length = 0;
        if (!::fa::cpu::mapping::checked_size_to_int(
                encoded_reference.size(), checked_reference_length)) {
          next_rna_ctx.coordinate_domain_refusal =
              ::fa::cpu::mapping::CoordinateDomainRefusal::UnrepresentableCount;
          break;
        }
      }
    }
    next_rna_seed_ctx = make_seed_context(config);
  }

  void validate_reconfiguration(
      const ::fa::cpu::options::ResolvedOptions& config) const {
    if (config.index.k <= 0) {
      throw std::invalid_argument("k must be at least 1");
    }
    if (!index_->empty() && config.index.k != index_->k()) {
      throw std::invalid_argument(
          "k=" + std::to_string(config.index.k) +
          " does not match the index (k=" + std::to_string(index_->k()) +
          "); rebuild the index to change k");
    }
    if (!index_->empty()) {
      const int requested_s =
          effective_closed_syncmer_s(
              std::max(1, config.index.k), config.index.syncmer_s);
      if (requested_s != index_->build_syncmer_s()) {
        throw std::invalid_argument(
            "s=" + std::to_string(requested_s) +
            " does not match the index (s=" +
            std::to_string(index_->build_syncmer_s()) +
            "); rebuild the index to change s");
      }
    }
  }

  void install_config(
      const ::fa::cpu::options::ResolvedOptions& config) {
    validate_reconfiguration(config);
    // Resolve the occurrence cap before composing the contexts that read it.
    ::fa::cpu::options::ResolvedOptions next_config = config;
    resolve_index_occ_cap_into(next_config);
    rna::Context next_rna_ctx;
    LongReadSeedContext next_rna_seed_ctx;
    compose_rna_runtime(next_config, next_rna_ctx, next_rna_seed_ctx);

    // Commit only once everything above has succeeded.
    cfg_ = std::move(next_config);
    rna_ctx_ = std::move(next_rna_ctx);
    rna_seed_ctx_ = std::move(next_rna_seed_ctx);
  }

  void compose_rna_runtime() {
    rna::Context next_rna_ctx;
    LongReadSeedContext next_rna_seed_ctx;
    compose_rna_runtime(cfg_, next_rna_ctx, next_rna_seed_ctx);
    rna_ctx_ = std::move(next_rna_ctx);
    rna_seed_ctx_ = std::move(next_rna_seed_ctx);
  }

private:
  static std::vector<uint8_t>
  reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc) {
    return ::fa::cpu::lr::reverse_complement_encoded_u8(fwd_enc);
  }

  void apply_chain_seed_constructor_config(int chain_syncmer_s,
                                           int chain_syncmer_downsample) {
    cfg_.index.syncmer_s = std::max(0, chain_syncmer_s);
    cfg_.index.syncmer_downsample =
        std::max(1, chain_syncmer_downsample);
  }

  // Per-read dispatch to the backend. DNA builds its contexts per call from cfg_; RNA reads
  // the ones compose_rna_runtime() built.
  template <class Backend>
  AlignResult align_anchor_chain(
      const std::string &read, const std::vector<uint8_t> &fwd_enc,
      std::vector<uint8_t> &rc_enc,
      typename Backend::WorkerScratch &worker_scratch,
      const std::vector<QuerySeed> *shared_fwd_syncmer_seeds = nullptr,
      std::uint32_t read_name_hash = 0, int self_contig = -1,
      int name_rank = 0) const {
    if constexpr (std::is_same_v<Backend, ::fa::cpu::lr::rna::RnaBackend>) {
      return ::fa::cpu::lr::rna::RnaBackend::map_read(
          rna_ctx_, rna_seed_ctx_, worker_scratch, read, fwd_enc, rc_enc,
          shared_fwd_syncmer_seeds, read_name_hash);
    } else {
      ::fa::cpu::lr::DnaContext dctx = make_dna_context();
      dctx.read_name_hash = read_name_hash;
      if (cfg_.long_read().all_chains && !cfg_.long_read().dual)
        dctx.dual_rank = name_rank;
      LongReadSeedContext seed_ctx = make_seed_context();
      if (cfg_.long_read().skip_self)
        seed_ctx.self_contig = self_contig;
      return Backend::map_read(dctx, seed_ctx, worker_scratch, read, fwd_enc,
                               rc_enc, shared_fwd_syncmer_seeds);
    }
  }

private:
  void ingest_genome(std::unordered_map<std::string, std::string> genome) {
    // Encodes contigs in name order, freeing each ASCII sequence as soon as it is encoded
    // so the reference is never held twice. Empty sequences are skipped.
    std::vector<std::string> names;
    names.reserve(genome.size());
    for (const auto &[name, seq] : genome) {
      if (!seq.empty())
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    chr_names_.reserve(chr_names_.size() + names.size());
    chr_encs_.reserve(chr_encs_.size() + names.size());
    chr_lengths_.reserve(chr_lengths_.size() + names.size());
    for (auto &name : names) {
      auto it = genome.find(name);
      chr_encs_.push_back(encode_sequence_u8(it->second));
      chr_lengths_.push_back(
          static_cast<int64_t>(chr_encs_.back().size()));
      std::string().swap(it->second);
      chr_names_.push_back(std::move(name));
    }
  }

  // Reference names and lengths from the index's name and offset tables.
  void bind_index_reference_metadata() {
    chr_names_ = index_->chromosome_names();
    chr_lengths_.clear();
    const uint64_t *offsets = index_->chrom_offsets_data();
    const std::size_t count = static_cast<std::size_t>(index_->chrom_count());
    if (offsets == nullptr)
      return;
    chr_lengths_.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
      chr_lengths_.push_back(static_cast<int64_t>(offsets[i + 1] - offsets[i]));
  }

  // The DNA backend's per-read context: the options it reads and views of the reference.
  ::fa::cpu::lr::DnaContext make_dna_context() const {
    // dna/context.h repeats these defaults rather than include options/dna_profile.h.
    static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.dna_dense_diag_min_runs ==
                      ::fa::cpu::lr::kDnaDenseDiagMinRuns,
                  "dna/context.h dna_dense_diag_min_runs mirror drifted");
    static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.dna_tandem_window ==
                      ::fa::cpu::lr::kDnaTandemWindow,
                  "dna/context.h dna_tandem_window mirror drifted");
    static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.cigar_dp_min_ksw_len ==
                      ::fa::cpu::lr::kDnaMinKswLen,
                  "dna/context.h cigar_dp_min_ksw_len mirror drifted");
    static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.dna_pool_gate_occ ==
                      ::fa::cpu::lr::kDnaPoolGateOcc,
                  "dna/context.h dna_pool_gate_occ mirror drifted");
    const auto& mapping = cfg_.long_read();
    // The chain's occurrence thresholds follow the resolved vote cap, unless --max-chain-occ
    // sets the pool gate.
    const ::fa::cpu::lr::DnaChainOccThresholds chain_occ =
        ::fa::cpu::lr::dna_chain_occ_thresholds(mapping);
    ::fa::cpu::lr::DnaContext dctx;
    dctx.ref = ::fa::cpu::engine::ReferenceContext{
        index_.get(), &chr_names_, &chr_lengths_, &chr_encs_};
    dctx.opts.cigar_local_interval_anchor_interval_pad =
        mapping.cigar_local_interval_anchor_interval_pad;
    dctx.opts.cigar_local_interval_anchor_chain_max_gap =
        mapping.cigar_local_interval_anchor_chain_max_gap;
    dctx.opts.screen_diag_band = mapping.screen_diag_band;
    dctx.opts.dna_dense_diag_min_runs = mapping.dna_dense_diag_min_runs;
    dctx.opts.dna_tandem_window = mapping.dna_tandem_window;
    dctx.opts.dna_pool_gate_occ = chain_occ.pool_gate_occ;
    dctx.opts.cigar_local_global_occ = chain_occ.global_occ;
    dctx.opts.k = cfg_.index.k;
    dctx.opts.cigar_dp_ambi = mapping.cigar_dp_ambi;
    dctx.opts.cigar_dp_bw = mapping.cigar_dp_bw;
    dctx.opts.cigar_dp_bw_long = mapping.cigar_dp_bw_long;
    dctx.opts.cigar_dp_gap_extend1 = mapping.cigar_dp_gap_extend1;
    dctx.opts.cigar_dp_gap_extend2 = mapping.cigar_dp_gap_extend2;
    dctx.opts.cigar_dp_gap_open1 = mapping.cigar_dp_gap_open1;
    dctx.opts.cigar_dp_gap_open2 = mapping.cigar_dp_gap_open2;
    dctx.opts.cigar_dp_match = mapping.cigar_dp_match;
    dctx.opts.cigar_dp_max_gap = mapping.cigar_dp_max_gap;
    dctx.opts.cigar_dp_min_dp_max = mapping.cigar_dp_min_dp_max;
    dctx.opts.cigar_dp_inversion_zdrop = mapping.cigar_dp_inversion_zdrop;
    dctx.opts.min_chain_score = mapping.min_chain_score;
    dctx.opts.pri_ratio = mapping.pri_ratio;
    dctx.opts.cigar_dp_split_min_anchors =
        mapping.cigar_dp_split_min_anchors;
    dctx.opts.cigar_dp_min_ksw_len = mapping.cigar_dp_min_ksw_len;
    dctx.opts.cigar_dp_mismatch = mapping.cigar_dp_mismatch;
    dctx.opts.cigar_dp_tail_end_bonus = mapping.cigar_dp_tail_end_bonus;
    dctx.opts.cigar_dp_tail_zdrop = mapping.cigar_dp_tail_zdrop;
    dctx.opts.fill_dp_match = mapping.fill_dp_match;
    dctx.opts.fill_dp_mismatch = mapping.fill_dp_mismatch;
    dctx.opts.fill_dp_ambi = mapping.fill_dp_ambi;
    dctx.opts.fill_dp_gap_open1 = mapping.fill_dp_gap_open1;
    dctx.opts.fill_dp_gap_extend1 = mapping.fill_dp_gap_extend1;
    dctx.opts.fill_dp_gap_open2 = mapping.fill_dp_gap_open2;
    dctx.opts.fill_dp_gap_extend2 = mapping.fill_dp_gap_extend2;
    dctx.opts.fill_dp_tail_zdrop = mapping.fill_dp_tail_zdrop;
    dctx.opts.fill_dp_inversion_zdrop = mapping.fill_dp_inversion_zdrop;
    dctx.opts.fill_dp_min_dp_max = mapping.fill_dp_min_dp_max;
    dctx.opts.enable_full_read_cigar =
        cfg_.common.enable_full_read_cigar;
    dctx.opts.cigar_replay_request = cigar_replay_request_for(cfg_.common);
    dctx.opts.residue_min_interval_bp = mapping.residue_min_interval_bp;
    dctx.opts.residue_min_anchor_density_per_100bp =
        mapping.residue_min_anchor_density_per_100bp;
    dctx.opts.postdp_rescoring = mapping.postdp_rescoring;
    dctx.opts.chain_mapq_hifi_margin = mapping.chain_mapq_hifi_margin;
    dctx.opts.inversion_probe_local_gate = mapping.inversion_probe_local_gate;
    dctx.opts.clip_nominate = mapping.dna_clip_nominate;
    dctx.opts.alternative_realize_max = mapping.alternative_realize_max;
    dctx.opts.chain_syncmer_s = cfg_.index.syncmer_s;
    dctx.opts.chain_syncmer_downsample = cfg_.index.syncmer_downsample;
    dctx.opts.vote_admission_ratio = mapping.vote_admission_ratio;
    dctx.opts.query_partition = mapping.query_partition;
    dctx.opts.query_tiles = mapping.query_tiles;
    dctx.opts.tile_owner_anchors = mapping.tile_owner_anchors;
    dctx.opts.all_chains = mapping.all_chains;
    dctx.opts.catalogue_lane_bound =
        ::fa::cpu::lr::dna_chain_max_candidates(mapping);
    return dctx;
  }

  ::fa::cpu::options::ResolvedOptions cfg_;
  std::vector<std::string> chr_names_;
  std::vector<std::vector<uint8_t>> chr_encs_;
  // Always filled; chr_encs_ is empty when the run needs no reference bases.
  std::vector<int64_t> chr_lengths_;
  std::shared_ptr<SeedIndex> index_ = std::make_shared<SeedIndex>();
  // Occurrence distribution of the attached index, computed once per index.
  std::shared_ptr<const ::fa::cpu::FaixOccDistribution> occ_dist_;
  const SeedIndex* occ_dist_source_ = nullptr;
  // Built by compose_rna_runtime(); read-only while mapping.
  rna::Context rna_ctx_;
  LongReadSeedContext rna_seed_ctx_;
};

// Maps a stream of read batches on a BatchWindow, so workers move on to the next batch
// instead of waiting for the slowest read of the current one. Output does not depend on the
// window size: results are stored at each read's original index and batches come back in
// submission order. The calling thread submits and collects; it does not map reads.
template <class Backend>
class WindowedAlignSession {
public:
  // One batch in flight. The session owns the reads until collect() returns them.
  struct Batch {
    std::vector<std::string> reads;
    std::vector<AlignResult> results; // one per read, at its original index
    std::vector<int> order;        // LPT position -> original index
    // Empty, or one tie_name_hash() per read.
    std::vector<std::uint32_t> read_name_hashes;
    // Empty, or one contig_id() per read.
    std::vector<int> self_contigs;
    // Empty, or one contig_name_rank() per read.
    std::vector<int> name_ranks;
  };

  // `engine` must outlive the session and must not be reconfigured while it exists.
  // `n_threads` is the number of scratch slots, i.e. threads that may map reads. The window
  // owns tids [0, n_threads - n_helpers); the other `n_helpers` tids belong to threads of
  // another stage that call try_help(). At least one worker is always the window's own.
  WindowedAlignSession(const LongReadEngine &engine, int n_threads, int window,
                       int n_helpers = 0)
      : engine_(&engine), n_threads_(n_threads < 1 ? 1 : n_threads),
        n_helpers_(clamp_helpers(n_threads_, n_helpers)),
        worker_scratch_(static_cast<size_t>(n_threads_)),
        window_(n_threads_ - n_helpers_,
                static_cast<size_t>(window < 1 ? 1 : window)) {}

  WindowedAlignSession(const WindowedAlignSession &) = delete;
  WindowedAlignSession &operator=(const WindowedAlignSession &) = delete;

  // Queues a batch and returns at once. Precondition: in_flight() < window().
  // `read_name_hashes` is empty or holds one hash per read, `self_contigs` is empty or
  // holds one contig_id() per read, and `name_ranks` is empty or holds one
  // contig_name_rank() per read.
  void submit(std::vector<std::string> reads,
              std::vector<std::uint32_t> read_name_hashes = {},
              std::vector<int> self_contigs = {},
              std::vector<int> name_ranks = {}) {
    std::unique_ptr<Batch> batch(new Batch());
    batch->reads = std::move(reads);
    batch->read_name_hashes = std::move(read_name_hashes);
    batch->self_contigs = std::move(self_contigs);
    batch->name_ranks = std::move(name_ranks);
    batch->results.resize(batch->reads.size());
    lpt_read_order(batch->reads, n_threads_, batch->order);
    // The window holds the Batch until collect(), so the pointer stays valid.
    Batch *state = batch.get();
    const std::int64_t items = static_cast<std::int64_t>(batch->order.size());
    window_.submit(items,
                   [this, state](std::int64_t position, int tid) {
                     map_one(*state, position, tid);
                   },
                   std::move(batch));
  }

  // Waits for the oldest batch and returns it, or rethrows its first exception.
  // Precondition: in_flight() > 0.
  std::unique_ptr<Batch> collect() {
    std::unique_ptr<Batch> batch = window_.collect();
    return batch;
  }

  // Maps at most one read from a thread outside the session (BatchWindow::run_one). `tid`
  // must be one of the helper tids reserved by the constructor, so each scratch slot has
  // one writer.
  bool try_help(int tid) {
    assert(tid >= n_threads_ - n_helpers_ && tid < n_threads_);
    return window_.run_one(tid);
  }

  std::size_t in_flight() const { return window_.in_flight(); }
  std::size_t window() const { return window_.window(); }
  // Includes the helper tids.
  int thread_count() const { return n_threads_; }

private:
  // Leaves the window at least one worker, so batches finish even if helpers stop calling.
  static int clamp_helpers(int n_threads, int n_helpers) {
    if (n_helpers < 0)
      return 0;
    return n_helpers > n_threads - 1 ? n_threads - 1 : n_helpers;
  }

  // Runs on a worker. Each position is claimed once and each tid has its own scratch, so no
  // lock is needed.
  void map_one(Batch &batch, std::int64_t position, int tid) {
    const int i = batch.order[static_cast<size_t>(position)];
    batch.results[static_cast<size_t>(i)] =
        engine_->template align<Backend>(
            batch.reads[static_cast<size_t>(i)],
            worker_scratch_[static_cast<size_t>(tid)],
            batch.read_name_hashes.empty()
                ? std::uint32_t{0}
                : batch.read_name_hashes[static_cast<size_t>(i)],
            batch.self_contigs.empty()
                ? -1
                : batch.self_contigs[static_cast<size_t>(i)],
            batch.name_ranks.empty()
                ? 0
                : batch.name_ranks[static_cast<size_t>(i)]);
  }

  const LongReadEngine *engine_;
  int n_threads_;
  int n_helpers_; // tids reserved for try_help()
  std::vector<typename Backend::WorkerScratch> worker_scratch_;
  // Declared last: its workers use the members above, so it must be destroyed first.
  ::fa::cpu::threading::BatchWindow<std::unique_ptr<Batch>> window_;
};

} // namespace lr
} // namespace cpu
} // namespace fa

namespace fa {
namespace cpu {

using LongReadEngine = ::fa::cpu::lr::LongReadEngine;

} // namespace cpu
} // namespace fa
