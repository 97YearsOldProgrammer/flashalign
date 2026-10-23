// The engine's plain members: loading, reference encoding, the DNA context and the RNA runtime.
#include "engine.h"

#include <limits>

/**** aligner ****/

namespace fa {
namespace cpu {
namespace lr {

void lpt_read_order(const std::vector<std::string> &reads, int n_threads,
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

const ::fa::cpu::options::ResolvedOptions& LongReadEngine::config() const { return cfg_; }

void LongReadEngine::reconfigure(const ::fa::cpu::options::ResolvedOptions& config) {
  install_config(config);
}

LongReadEngine::LongReadEngine(std::unordered_map<std::string, std::string> genome,
                               int k, int min_support, int chain_syncmer_s,
                               int chain_syncmer_downsample, int build_threads) {
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

std::unique_ptr<LongReadEngine>
LongReadEngine::from_index(std::unordered_map<std::string, std::string> genome,
                           const std::string &index_path, int min_support,
                           int chain_syncmer_downsample, int n_threads,
                           bool reference_bases_needed, uint64_t image_offset,
                           uint64_t image_bytes) {
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

std::unique_ptr<LongReadEngine>
LongReadEngine::from_shared_index(std::shared_ptr<SeedIndex> index,
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

std::vector<std::string> LongReadEngine::chromosome_names() const { return chr_names_; }

int LongReadEngine::contig_id(std::string_view name) const {
  const auto it = std::lower_bound(
      chr_names_.begin(), chr_names_.end(), name,
      [](const std::string &contig, std::string_view wanted) {
        return contig < wanted;
      });
  return it != chr_names_.end() && *it == name
             ? static_cast<int>(it - chr_names_.begin())
             : -1;
}

int LongReadEngine::contig_dual_rank(std::string_view name, std::int64_t length) const {
  return static_cast<int>(
      std::partition_point(dual_order_.begin(), dual_order_.end(),
                           [&](int contig) {
                             return dual_before(contig, name, length);
                           }) -
      dual_order_.begin());
}

std::vector<int64_t> LongReadEngine::chromosome_lengths() const { return chr_lengths_; }

int LongReadEngine::get_index_syncmer_s() const { return index_->build_syncmer_s(); }

DnaReadGeometry LongReadEngine::dna_read_geometry(int read_length) const {
  DnaReadGeometry geometry;
  geometry.vote_diag_bin_width =
      effective_vote_diag_bin_width(make_seed_context(), read_length);
  return geometry;
}

DnaChainGeometry LongReadEngine::dna_chain_geometry() const {
  const ::fa::cpu::lr::DnaContext dctx = make_dna_context();
  const auto params = ::fa::cpu::lr::dna_candidate_chain_params(
      dctx, ::fa::cpu::lr::internal::pass_diagonal_band(dctx),
      std::max(1, dctx.opts.k));
  DnaChainGeometry geometry;
  geometry.band = params.bw;
  geometry.reach = params.max_dist_x;
  return geometry;
}

void LongReadEngine::resolve_index_occ_cap_into(
    ::fa::cpu::options::ResolvedOptions& config) {
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
  record.quantile_f = mapping.long_occ_quantile_f;
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

LongOccPolicyConfig
LongReadEngine::long_occ_policy_config(const ::fa::cpu::lr::DnaLongOptions& mapping) {
  LongOccPolicyConfig cfg;
  cfg.policy = mapping.long_occ_policy;
  cfg.primary_occ_cap = mapping.long_primary_occ_cap;
  cfg.platform_occ_cap = mapping.long_occ_cap;
  return cfg;
}

int LongReadEngine::effective_chain_syncmer_s() const {
  return effective_closed_syncmer_s(
      std::max(1, cfg_.index.k), cfg_.index.syncmer_s);
}

LongReadSeedContext
LongReadEngine::make_seed_context(
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

LongReadSeedContext LongReadEngine::make_seed_context() const {
  return make_seed_context(cfg_);
}

::fa::cpu::output::CigarReplayRequest LongReadEngine::cigar_replay_request_for(
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

void LongReadEngine::compose_rna_runtime(
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
  // The fine harvest's cap is the vote's N; none when the cap is off.
  const int occ_cap =
      effective_long_primary_occ_cap(long_occ_policy_config(mapping));
  opts.cigar_local_global_occ =
      occ_cap > 0 ? occ_cap : std::numeric_limits<int>::max();
  opts.fine_chain_query_gap = mapping.cigar_dp_max_gap;
  opts.fine_chain_band = mapping.cigar_dp_bw;
  opts.query_occ_frac = static_cast<float>(mapping.query_occ_frac);
  opts.occ_dist = mapping.occ_dist;
  if (rna_mapping) {
    opts.rival_lifecycle.pri_ratio = rna_mapping->rival_pri_ratio;
    opts.rival_lifecycle.min_diff = rna_mapping->rival_min_diff;
    opts.rival_lifecycle.realize_max = rna_mapping->rival_realize_max;
    opts.mapq_qcov_tau = rna_mapping->rna_mapq_qcov_tau;
    opts.mapq_disjoint_vote_damp = rna_mapping->rna_mapq_disjoint_vote_damp;
    opts.vote_seed_spacing = rna_mapping->vote_seed_spacing;
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

void LongReadEngine::validate_reconfiguration(
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

void LongReadEngine::install_config(
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
  bind_dual_order();
}

void LongReadEngine::compose_rna_runtime() {
  rna::Context next_rna_ctx;
  LongReadSeedContext next_rna_seed_ctx;
  compose_rna_runtime(cfg_, next_rna_ctx, next_rna_seed_ctx);
  rna_ctx_ = std::move(next_rna_ctx);
  rna_seed_ctx_ = std::move(next_rna_seed_ctx);
}

std::vector<uint8_t>
LongReadEngine::reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc) {
  return ::fa::cpu::lr::reverse_complement_encoded_u8(fwd_enc);
}

void LongReadEngine::apply_chain_seed_constructor_config(int chain_syncmer_s,
                                                         int chain_syncmer_downsample) {
  cfg_.index.syncmer_s = std::max(0, chain_syncmer_s);
  cfg_.index.syncmer_downsample =
      std::max(1, chain_syncmer_downsample);
}

void LongReadEngine::ingest_genome(std::unordered_map<std::string, std::string> genome) {
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

void LongReadEngine::bind_index_reference_metadata() {
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

bool LongReadEngine::dual_before(int contig, std::string_view name,
                                 std::int64_t length) const {
  const std::size_t at = static_cast<std::size_t>(contig);
  if (chr_lengths_[at] != length)
    return chr_lengths_[at] < length;
  return chr_names_[at] < name;
}

void LongReadEngine::bind_dual_order() {
  dual_order_.clear();
  dual_contig_rank_.clear();
  if (!cfg_.long_read().all_chains || cfg_.long_read().dual)
    return;
  dual_order_.resize(chr_names_.size());
  std::iota(dual_order_.begin(), dual_order_.end(), 0);
  std::stable_sort(dual_order_.begin(), dual_order_.end(),
                   [&](int a, int b) {
                     return dual_before(
                         a, chr_names_[static_cast<std::size_t>(b)],
                         chr_lengths_[static_cast<std::size_t>(b)]);
                   });
  dual_contig_rank_.resize(dual_order_.size());
  for (std::size_t rank = 0; rank < dual_order_.size(); ++rank)
    dual_contig_rank_[static_cast<std::size_t>(dual_order_[rank])] =
        static_cast<int>(rank);
}

::fa::cpu::lr::DnaContext LongReadEngine::make_dna_context() const {
  // dna/dna.h repeats these defaults rather than include options/options.h.
  static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.dna_dense_diag_min_runs ==
                    ::fa::cpu::lr::kDnaDenseDiagMinRuns,
                "dna/dna.h dna_dense_diag_min_runs mirror drifted");
  static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.cigar_dp_min_ksw_len ==
                    ::fa::cpu::lr::kDnaMinKswLen,
                "dna/dna.h cigar_dp_min_ksw_len mirror drifted");
  static_assert(::fa::cpu::lr::ResolvedDnaOptions{}.dna_pool_gate_occ ==
                    ::fa::cpu::lr::kDnaPoolGateOcc,
                "dna/dna.h dna_pool_gate_occ mirror drifted");
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
  dctx.opts.dna_dense_diag_min_runs = mapping.dna_dense_diag_min_runs;
  dctx.opts.dna_pool_gate_occ = chain_occ.pool_gate_occ;
  dctx.opts.dna_query_occ_frac = static_cast<float>(mapping.query_occ_frac);
  dctx.opts.occ_dist = mapping.occ_dist;
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
  dctx.opts.inversion_gap = mapping.inversion_gap;
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
  dctx.opts.chain_mapq_hifi_margin = mapping.chain_mapq_hifi_margin;
  dctx.opts.inversion_probe_local_gate = mapping.inversion_probe_local_gate;
  dctx.opts.secondary_max = mapping.secondary_max;
  dctx.opts.rival_split_ends = mapping.rival_split_ends;
  dctx.opts.chain_syncmer_s = cfg_.index.syncmer_s;
  dctx.opts.chain_syncmer_downsample = cfg_.index.syncmer_downsample;
  dctx.opts.vote_admission_ratio = mapping.vote_admission_ratio;
  dctx.opts.query_tiles = mapping.query_tiles;
  dctx.opts.all_chains = mapping.all_chains;
  dctx.opts.catalogue_lane_bound =
      ::fa::cpu::lr::dna_chain_max_candidates(mapping);
  return dctx;
}

} // namespace lr
} // namespace cpu
} // namespace fa
