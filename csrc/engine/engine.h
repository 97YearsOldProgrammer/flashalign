// The long-read engine: owns the index, the reference and the resolved options; its per-read and windowed align templates over a lane backend.
#pragma once

#include "../core/core.h"
#include "../dna/dna.h"
#include "../index/format.h"
#include "../index/index.h"
#include "../index/seed.h"
#include "../options/options.h"
#include "../rna/rna.h"
#include "../threading/threading.h"


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
#include <variant>
#include <vector>

/**** aligner ****/
// The long-read engine: owns the index and the reference, holds the resolved options and
// maps reads through the DNA or RNA backend. WindowedAlignSession runs it over a stream of
// batches.

namespace fa {
namespace cpu {
namespace lr {

// Fills `order` with read indices, longest read first (LPT scheduling with one-item claims).
// Small batches keep the input order.
void lpt_read_order(const std::vector<std::string> &reads, int n_threads,
                    std::vector<int> &order);

// What the DNA lane derives from the read length L under the current options.
struct DnaReadGeometry {
  int vote_diag_bin_width = 0;
};

// The whole-query (dense) chain pass's band and reach (max_dist_x/y), the
// same for every read.
struct DnaChainGeometry {
  int band = 0;
  int reach = 0;
};

class LongReadEngine {
public:
  LongReadEngine() = default;

  const ::fa::cpu::options::ResolvedOptions& config() const;
  void reconfigure(const ::fa::cpu::options::ResolvedOptions& config);
  // Builds the index over `genome`. With `packed_reference` (the DNA presets) the bases
  // are packed into the index as a .faix holds them and the unpacked copy is freed; the
  // DNA context unpacks each region's window.
  LongReadEngine(std::unordered_map<std::string, std::string> genome,
                 int k = 15, int min_support = 3, int chain_syncmer_s = 9,
                 int chain_syncmer_downsample = 1, int build_threads = 0,
                 bool packed_reference = false);

  // Loads an engine over a prebuilt .faix. The seed length k and syncmer s come from the
  // index, as in minimap2.
  //
  // `n_threads` bounds the threads used to read the file and unpack the reference; <= 0
  // uses the default. When `reference_bases_needed` is false (DNA map-only PAF) the
  // reference stays packed and only names and lengths are read from the index.
  // `image_offset` / `image_bytes` select one part of a multi-part index; both 0 means the
  // whole file. `packed_reference` (the DNA presets) keeps the needed bases packed in the
  // index, for the DNA context to unpack per region; the caller passes it only for an
  // index that embeds its reference.
  static std::unique_ptr<LongReadEngine>
  from_index(std::unordered_map<std::string, std::string> genome,
             const std::string &index_path, int min_support = 3,
             int chain_syncmer_downsample = 1, int n_threads = 0,
             bool reference_bases_needed = true, uint64_t image_offset = 0,
             uint64_t image_bytes = 0, bool packed_reference = false);

  // Builds an engine over an already-loaded index, which may be shared by several engines.
  // The index must embed its reference. Options come only from `config`.
  static std::unique_ptr<LongReadEngine>
  from_shared_index(std::shared_ptr<SeedIndex> index,
                    const ::fa::cpu::options::ResolvedOptions& config);

  // Maps one read with the given Backend (dna::DnaBackend or rna::RnaBackend). This overload
  // uses a fresh scratch; the next one takes a per-worker scratch reused across reads, like
  // minimap2's mm_tbuf_t.
  //
  // `read_name_hash` is tie_name_hash() of the read name and seeds the tie-break between two
  // equally good loci; 0 when there is no name. `self_contig`, where taken, is contig_id()
  // of the read name, read only where the preset leaves a read out of its own vote
  // (DnaLongOptions::skip_self); -1 when it is not known. `name_rank`, where taken, is
  // contig_dual_rank() of the read, read only where the preset prints a pair once
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
  // holds one contig_dual_rank() per read.
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

  std::vector<std::string> chromosome_names() const;
  // The contig named exactly `name`, or -1. Contigs are numbered in name order.
  int contig_id(std::string_view name) const;
  // A read's rank in the --dual=no pair order (bind_dual_order): the number of
  // contigs ordered before a read of this name and length. 0 outside that lane.
  int contig_dual_rank(std::string_view name, std::int64_t length) const;
  // Lengths in chromosome_names() order, available even when the bases were not unpacked.
  std::vector<int64_t> chromosome_lengths() const;
  int get_index_syncmer_s() const;
  // The DNA lane's geometry for a read of `read_length` bp, as map_read
  // derives it (--show-config).
  DnaReadGeometry dna_read_geometry(int read_length) const;
  DnaChainGeometry dna_chain_geometry() const;

private:
  // Sets the vote occurrence cap from the attached index under LongOccPolicy::Platform:
  // max(preset floor, the index's occurrence quantile at the preset's f as minimap2's
  // mm_idx_cal_max_occ), then at most the preset ceiling. Records how it was derived in
  // long_occ_index_resolution. Does nothing for a fixed cap (--max-vote-occ) or without
  // an index. Repeated calls start again from the recorded floor.
  void resolve_index_occ_cap_into(
      ::fa::cpu::options::ResolvedOptions& config);

  static LongOccPolicyConfig
  long_occ_policy_config(const ::fa::cpu::lr::DnaLongOptions& mapping);
  int effective_chain_syncmer_s() const;

  // The seeding parameters the shared seeding code reads.
  LongReadSeedContext
  make_seed_context(
      const ::fa::cpu::options::ResolvedOptions& config) const;
  LongReadSeedContext make_seed_context() const;

  // The --cs / --MD / --eqx request passed to both realizers.
  static ::fa::cpu::output::CigarReplayRequest cigar_replay_request_for(
      const ::fa::cpu::options::CommonOptions& common);

  // Builds the RNA backend's context and seed context from `config` and the reference.
  // Unused in DNA mode.
  void compose_rna_runtime(
                           const ::fa::cpu::options::ResolvedOptions& config,
                           rna::Context &next_rna_ctx,
                           LongReadSeedContext &next_rna_seed_ctx) const;

  void validate_reconfiguration(
      const ::fa::cpu::options::ResolvedOptions& config) const;

  void install_config(
      const ::fa::cpu::options::ResolvedOptions& config);

  void compose_rna_runtime();

private:
  static std::vector<uint8_t>
  reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc);

  void apply_chain_seed_constructor_config(int chain_syncmer_s,
                                           int chain_syncmer_downsample);

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
      if (cfg_.long_read().all_chains && !cfg_.long_read().dual) {
        dctx.dual_rank = name_rank;
        dctx.dual_contig_rank = &dual_contig_rank_;
      }
      LongReadSeedContext seed_ctx = make_seed_context();
      if (cfg_.long_read().skip_self)
        seed_ctx.self_contig = self_contig;
      AlignResult alignment = Backend::map_read(
          dctx, seed_ctx, worker_scratch, read, fwd_enc, rc_enc,
          shared_fwd_syncmer_seeds);
      // Every DNA read passes here, on the thread that mapped it.
      ::fa::cpu::lr::dna::settle_kept_scratch(worker_scratch);
      return alignment;
    }
  }

private:
  void ingest_genome(std::unordered_map<std::string, std::string> genome);

  // Reference names and lengths from the index's name and offset tables.
  void bind_index_reference_metadata();

  // Whether `contig` comes before a read of this name and length in the
  // --dual=no pair order: by length, then name.
  bool dual_before(int contig, std::string_view name,
                   std::int64_t length) const;

  // The all-chains lane's --dual=no pair order: dual_order_ lists the contigs
  // in it, dual_contig_rank_ gives each contig's place. Empty outside that
  // lane.
  void bind_dual_order();

  // The DNA backend's per-read context: the options it reads and views of the reference.
  ::fa::cpu::lr::DnaContext make_dna_context() const;

  ::fa::cpu::options::ResolvedOptions cfg_;
  std::vector<std::string> chr_names_;
  std::vector<std::vector<uint8_t>> chr_encs_;
  // Always filled; chr_encs_ is empty when the run needs no reference bases or reads them
  // packed.
  std::vector<int64_t> chr_lengths_;
  // The constructors' packed_reference: the DNA context reads the bases from the index.
  bool packed_reference_ = false;
  // bind_dual_order(); read-only while mapping.
  std::vector<int> dual_order_;
  std::vector<int> dual_contig_rank_;
  std::shared_ptr<SeedIndex> index_ = std::make_shared<SeedIndex>();
  // Occurrence distribution of the attached index, computed once per index.
  std::shared_ptr<const ::fa::cpu::FaixOccDistribution> occ_dist_;
  const SeedIndex* occ_dist_source_ = nullptr;
  // Built by compose_rna_runtime(); read-only while mapping.
  rna::Context rna_ctx_;
  LongReadSeedContext rna_seed_ctx_;
};

// Maps a stream of read batches on a BatchWindow of threading::kBatchWindow batches, so
// workers move on to the next batch instead of waiting for the slowest read of the current
// one; it takes a read of a batch younger than the oldest only while the bases in flight
// stay within `max_bases`. Results equal align_batch's for the same batch: they are stored
// at each read's original index and batches come back in submission order. The calling
// thread submits and collects; it does not map reads.
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
    // Empty, or one contig_dual_rank() per read.
    std::vector<int> name_ranks;
  };

  // `engine` must outlive the session and must not be reconfigured while it exists.
  // `n_threads` is the number of scratch slots, i.e. threads that may map reads. The window
  // owns tids [0, n_threads - n_helpers); the other `n_helpers` tids belong to threads of
  // another stage that call try_help(). At least one worker is always the window's own.
  WindowedAlignSession(const LongReadEngine &engine, int n_threads, int n_helpers,
                       std::uint64_t max_bases)
      : engine_(&engine), n_threads_(n_threads < 1 ? 1 : n_threads),
        n_helpers_(clamp_helpers(n_threads_, n_helpers)),
        worker_scratch_(static_cast<size_t>(n_threads_)),
        window_(n_threads_ - n_helpers_, ::fa::cpu::threading::kBatchWindow,
                max_bases) {}

  WindowedAlignSession(const WindowedAlignSession &) = delete;
  WindowedAlignSession &operator=(const WindowedAlignSession &) = delete;

  // Queues a batch and returns at once. Precondition: in_flight() < window().
  // `read_name_hashes` is empty or holds one hash per read, `self_contigs` is empty or
  // holds one contig_id() per read, and `name_ranks` is empty or holds one
  // contig_dual_rank() per read.
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
    window_.submit(
        items,
        [this, state](std::int64_t position, int tid) {
          map_one(*state, position, tid);
        },
        [state](std::int64_t position) {
          const int i = state->order[static_cast<size_t>(position)];
          return static_cast<std::uint64_t>(
              state->reads[static_cast<size_t>(i)].size());
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
