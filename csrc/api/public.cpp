// The public C++ API (include/flashalign/) over the internal api::LongReadAligner.
#include <flashalign/aligner.hpp>

#include "aligner.h"
#include "../core/types.h"
#include "../index/index.h"
#include "../io/fastx.h"
#include "../io/paf.h"
#include "../options/presets.h"
#include "../seeding/tie_hash.h"
#include "../threading/parallel_for.h"

#include <algorithm>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace flashalign {
namespace {

fa::cpu::options::CsMode internal_cs(CsMode mode) {
  switch (mode) {
    case CsMode::Short: return fa::cpu::options::CsMode::Short;
    case CsMode::Long: return fa::cpu::options::CsMode::Long;
    case CsMode::None: break;
  }
  return fa::cpu::options::CsMode::None;
}

CsMode public_cs(fa::cpu::options::CsMode mode) {
  switch (mode) {
    case fa::cpu::options::CsMode::Short: return CsMode::Short;
    case fa::cpu::options::CsMode::Long: return CsMode::Long;
    case fa::cpu::options::CsMode::None: break;
  }
  return CsMode::None;
}

fa::cpu::options::ResolvedOptions resolve_config(
    const Config& source, int index_k = -1, int index_s = -1) {
  if (!fa::cpu::options::preset_is_valid(source.preset)) {
    throw std::invalid_argument(
        "flashalign: unknown preset: " + source.preset);
  }
  fa::cpu::options::ResolveRequest request;
  request.preset = source.preset;
  if (index_k > 0 && index_s > 0) {
    request.index.has_index = true;
    request.index.k = index_k;
    request.index.syncmer_s = index_s;
  } else {
    if (source.k > 0) request.user.k = source.k;
    if (source.syncmer_s > 0)
      request.user.syncmer_s = source.syncmer_s;
  }
  if (source.min_support >= 0)
    request.user.min_support = source.min_support;
  if (source.max_query_seeds >= 0) {
    request.user.max_query_seeds_per_strand = source.max_query_seeds;
  }
  if (source.tile_score_hit)
    request.user.query_tile_supported_reward = *source.tile_score_hit;
  if (source.tile_score_block)
    request.user.query_tile_block_open_cost = *source.tile_score_block;
  if (source.tile_score_null)
    request.user.query_tile_null_cost = *source.tile_score_null;
  if (source.tile_score_miss)
    request.user.query_tile_unsupported_cost = *source.tile_score_miss;
  request.user.enable_full_read_cigar = source.full_read_cigar;
  request.user.syncmer_downsample =
      std::max(1, source.syncmer_downsample);
  if (source.vote_diag_bin_width != -1) {
    if (source.vote_diag_bin_width <= 0) {
      throw std::invalid_argument(
          "flashalign: vote_diag_bin_width must be -1 (inherit preset) or "
          "a positive integer (got " +
          std::to_string(source.vote_diag_bin_width) + ")");
    }
    request.user.vote_diag_bin_width = source.vote_diag_bin_width;
  }
  if (source.long_primary_occ_cap >= 0)
    request.user.primary_occ_cap = source.long_primary_occ_cap;
  if (source.long_occ_cap >= 0)
    request.user.long_occ_cap = source.long_occ_cap;
  request.user.num_threads =
      source.threads > 0
          ? source.threads
          : fa::cpu::threading::default_thread_count();
  if (source.dp_match >= 0)
    request.user.dp_match = source.dp_match;
  if (source.dp_mismatch >= 0)
    request.user.dp_mismatch = source.dp_mismatch;
  if (source.dp_score_n >= 0)
    request.user.dp_ambi = source.dp_score_n;
  if (source.dp_gap_open1 >= 0)
    request.user.dp_gap_open1 = source.dp_gap_open1;
  if (source.dp_gap_extend1 >= 0)
    request.user.dp_gap_extend1 = source.dp_gap_extend1;
  if (source.dp_gap_open2 >= 0)
    request.user.dp_gap_open2 = source.dp_gap_open2;
  if (source.dp_gap_extend2 >= 0)
    request.user.dp_gap_extend2 = source.dp_gap_extend2;
  if (source.dp_zdrop >= 0)
    request.user.dp_tail_zdrop = source.dp_zdrop;
  // dp_zdrop_inv is -z's second value: -1 follows dp_zdrop, as a lone -z does.
  if (source.dp_zdrop_inv >= 0) {
    if (source.dp_zdrop < 0) {
      throw std::invalid_argument(
          "flashalign: dp_zdrop_inv needs dp_zdrop, as -z INT1,INT2 does");
    }
    if (source.dp_zdrop < source.dp_zdrop_inv) {
      throw std::invalid_argument(
          "flashalign: dp_zdrop (" + std::to_string(source.dp_zdrop) +
          ") is less than dp_zdrop_inv (" +
          std::to_string(source.dp_zdrop_inv) +
          "); lower dp_zdrop_inv too, or set it to -1 to follow dp_zdrop");
    }
    request.user.dp_inversion_zdrop = source.dp_zdrop_inv;
  }
  if (source.dp_end_bonus >= -1) {
    request.user.dp_tail_end_bonus = source.dp_end_bonus;
  }
  if (source.dp_min_score >= 0)
    request.user.dp_min_dp_max = source.dp_min_score;
  // Only an explicit bandwidth is passed; the resolver rejects one under a splice preset.
  if (source.dp_bw > 0) request.user.dp_bw = source.dp_bw;
  if (source.dp_bw_long > 0) request.user.dp_bw_long = source.dp_bw_long;
  request.user.cs = internal_cs(source.cs);
  request.user.emit_md = source.emit_md;
  request.user.emit_eqx = source.emit_eqx;
  if (fa::cpu::options::is_rna_preset(source.preset)) {
    request.user.rna_min_intron = std::max(1, source.rna_min_intron);
    request.user.rna_max_intron =
        std::max(*request.user.rna_min_intron, source.rna_max_intron);
    request.user.rna_strand = source.rna_strand_mode;
    if (!source.rna_junction_bed.empty())
      request.user.rna_junction_bed = source.rna_junction_bed;
    if (source.rna_junction_bonus >= 0)
      request.user.rna_junction_bonus = source.rna_junction_bonus;
    if (source.rna_rival_pri_ratio >= 0.0)
      request.user.pri_ratio = source.rna_rival_pri_ratio;
    if (source.rna_max_loci >= 0)
      request.user.rna_max_loci = source.rna_max_loci;
  }
  fa::cpu::options::ResolvedOptions config =
      fa::cpu::options::resolve_options(request).resolved();
  return config;
}

Config to_public_config(const fa::cpu::options::ResolvedOptions& source,
                        std::string preset) {
  const auto& mapping = source.long_read();
  const auto* rna = source.rna();
  Config config;
  config.preset = std::move(preset);
  config.k = source.index.k;
  config.min_support = source.common.min_support;
  config.max_query_seeds = source.common.max_query_seeds_per_strand;
  // A splice preset rejects --tile-score, so report None there, as dp_bw below.
  if (!rna) {
    config.tile_score_hit = mapping.query_partition.supported_tile_reward;
    config.tile_score_block = mapping.query_partition.block_open_cost;
    config.tile_score_null = mapping.query_partition.null_tile_cost;
    config.tile_score_miss =
        mapping.query_partition.unsupported_ownership_cost;
  }
  config.full_read_cigar = source.common.enable_full_read_cigar;
  config.syncmer_s = source.index.syncmer_s;
  config.syncmer_downsample = source.index.syncmer_downsample;
  // An explicit width turns the adaptive width off, so an adaptive one reports -1.
  config.vote_diag_bin_width = mapping.vote_diag_bin_width_adaptive
                                   ? -1
                                   : mapping.vote_diag_bin_width;
  config.long_occ_cap = mapping.long_occ_cap;
  config.long_primary_occ_cap = mapping.long_primary_occ_cap;
  config.threads = source.common.num_threads;
  // DNA: the gap-fill row, which these set; the end row is the preset's, and
  // on an asm preset the two rows are one.
  config.dp_match = rna ? mapping.cigar_dp_match : mapping.fill_dp_match;
  config.dp_mismatch =
      rna ? mapping.cigar_dp_mismatch : mapping.fill_dp_mismatch;
  config.dp_score_n = rna ? mapping.cigar_dp_ambi : mapping.fill_dp_ambi;
  config.dp_gap_open1 =
      rna ? mapping.cigar_dp_gap_open1 : mapping.fill_dp_gap_open1;
  config.dp_gap_extend1 =
      rna ? mapping.cigar_dp_gap_extend1 : mapping.fill_dp_gap_extend1;
  config.dp_gap_open2 =
      rna ? mapping.cigar_dp_gap_open2 : mapping.fill_dp_gap_open2;
  config.dp_gap_extend2 =
      rna ? mapping.cigar_dp_gap_extend2 : mapping.fill_dp_gap_extend2;
  config.dp_zdrop =
      rna ? mapping.cigar_dp_tail_zdrop : mapping.fill_dp_tail_zdrop;
  config.dp_zdrop_inv = rna ? rna->splice_inversion_zdrop
                            : mapping.fill_dp_inversion_zdrop;
  config.dp_end_bonus = mapping.cigar_dp_tail_end_bonus;
  config.dp_min_score = mapping.cigar_dp_min_dp_max;
  // A splice preset fixes the bandwidth, so report -1 and keep config() -> reconfigure() valid.
  config.dp_bw = rna ? -1 : mapping.cigar_dp_bw;
  config.dp_bw_long = rna ? -1 : mapping.cigar_dp_bw_long;
  config.cs = public_cs(source.common.cs);
  config.emit_md = source.common.emit_md;
  config.emit_eqx = source.common.emit_eqx;
  config.rna_min_intron = rna ? rna->min_intron : 20;
  config.rna_max_intron = rna ? rna->max_intron : 200000;
  config.rna_strand_mode = rna ? rna->strand_mode : -1;
  config.rna_junction_bed = rna ? rna->junction_bed : std::string{};
  config.rna_junction_bonus = rna ? rna->splice_junction_bonus : -1;
  config.rna_rival_pri_ratio = rna ? rna->rival_pri_ratio : -1.0;
  config.rna_max_loci = rna ? rna->max_locus_chains : -1;
  return config;
}

std::unordered_map<std::string, std::string> genome_map(
    const std::vector<std::pair<std::string, std::string>>& references) {
  std::unordered_map<std::string, std::string> genome;
  for (const auto& reference : references) {
    if (!reference.second.empty())
      genome.emplace(reference.first, reference.second);
  }
  return genome;
}

}  // namespace

std::pair<int, int> preset_seeding(const std::string& preset) {
  if (!fa::cpu::options::preset_is_valid(preset)) {
    throw std::invalid_argument(
        "flashalign: unknown preset: " + preset);
  }
  const fa::cpu::options::PresetSeeding seeding =
      fa::cpu::options::resolve_preset_seeding(preset);
  return {seeding.k, seeding.syncmer_s};
}

struct Index::Impl {
  explicit Impl(std::shared_ptr<fa::cpu::SeedIndex> value)
      : index(std::move(value)) {}
  std::shared_ptr<fa::cpu::SeedIndex> index;
  // Embedded reference for Index::sequence(), one code per base, unpacked on first use.
  // Shared by every copy of the Index.
  std::mutex reference_mutex;
  bool reference_unpacked = false;
  std::vector<std::vector<std::uint8_t>> reference;
};

Index::Index()
    : impl_(std::make_shared<Impl>(
          std::make_shared<fa::cpu::SeedIndex>())) {}
Index::Index(std::shared_ptr<Impl> impl) : impl_(std::move(impl)) {}
Index::~Index() = default;
Index::Index(const Index&) = default;
Index& Index::operator=(const Index&) = default;
Index::Index(Index&&) noexcept = default;
Index& Index::operator=(Index&&) noexcept = default;

Index Index::load(const std::string& path) {
  fa::cpu::FaixLoadStatus status;
  auto index =
      std::make_shared<fa::cpu::SeedIndex>(
          fa::cpu::load_seed_index(path, &status));
  if (index->empty())
    throw std::runtime_error(
        "flashalign: cannot load index " + path + ": " +
        fa::cpu::faix_load_error_message(status));
  return Index(std::make_shared<Impl>(std::move(index)));
}

Index Index::build(
    const std::vector<std::pair<std::string, std::string>>& references,
    int k, int syncmer_s, int threads) {
  std::vector<std::size_t> order(references.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return references[lhs].first < references[rhs].first;
  });
  std::vector<std::string> names;
  std::vector<std::vector<std::uint8_t>> encoded;
  for (std::size_t i : order) {
    if (references[i].second.empty()) continue;
    names.push_back(references[i].first);
    encoded.push_back(fa::cpu::encode_sequence_u8(references[i].second));
  }
  fa::cpu::FaixBuildConfig config;
  config.k = k;
  config.syncmer_s = syncmer_s;
  config.build_threads =
      threads > 0 ? threads : fa::cpu::threading::default_thread_count();
  fa::cpu::FaixBuildStatus build_status;
  auto index = std::make_shared<fa::cpu::SeedIndex>(
      fa::cpu::build_seed_index(encoded, config, &build_status));
  if (index->empty()) {
    throw std::runtime_error("flashalign: failed to build the index: " +
                             fa::cpu::faix_build_error_message(build_status));
  }
  if (!index->set_reference_payload(std::move(names), encoded)) {
    throw std::runtime_error(
        "flashalign: failed to attach the reference to the index");
  }
  return Index(std::make_shared<Impl>(std::move(index)));
}

Index Index::build_from_fasta(
    const std::string& path, int k, int syncmer_s, int threads) {
  if (fa::cpu::api::reference_kind(path) == fa::cpu::api::ReferenceKind::Index) {
    throw std::runtime_error(
        "flashalign: " + path + " is already a FlashAlign index; load it "
        "rather than build one from it");
  }
  // The CLI's reader: plain or gzip input, or "-", uppercased.
  fa::cpu::io::FastxReader reader(path);
  fa::cpu::io::FastxRecord record;
  std::vector<std::pair<std::string, std::string>> references;
  while (reader.next(record)) {
    if (record.seq.empty()) continue;
    references.emplace_back(std::move(record.name), std::move(record.seq));
  }
  if (references.empty()) {
    throw std::runtime_error(
        "flashalign: reference FASTA has no sequences: " + path);
  }
  // build() sorts contigs by name, as the CLI does.
  return build(references, k, syncmer_s, threads);
}

bool Index::is_index_file(const std::string& path) {
  // Reads only the header; an unopenable path is not an index.
  if (path.empty() || path == "-") return false;
  return fa::cpu::api::LongReadAligner::probe_index_header(path).is_index;
}

void Index::save(const std::string& path) const {
  if (!impl_->index->save(path))
    throw std::runtime_error("flashalign: failed to save seed index: " + path);
}
int Index::k() const { return impl_->index->k(); }
int Index::syncmer_s() const { return impl_->index->build_syncmer_s(); }
bool Index::has_reference() const {
  return impl_->index->has_reference_payload();
}
std::int64_t Index::total_bases() const {
  return static_cast<std::int64_t>(impl_->index->total_bp());
}
double Index::memory_megabytes() const {
  return static_cast<double>(impl_->index->memory_bytes()) / 1e6;
}
std::string Index::preset() const { return impl_->index->build_preset(); }
std::vector<std::string> Index::reference_names() const {
  return impl_->index->chromosome_names();
}

std::vector<std::int64_t> Index::reference_lengths() const {
  std::vector<std::int64_t> lengths;
  const std::uint64_t* offsets = impl_->index->chrom_offsets_data();
  const std::size_t count =
      static_cast<std::size_t>(impl_->index->chrom_count());
  if (offsets == nullptr) return lengths;
  lengths.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    lengths.push_back(
        static_cast<std::int64_t>(offsets[i + 1] - offsets[i]));
  }
  return lengths;
}

std::string Index::sequence(
    const std::string& name, std::int64_t start, std::int64_t end) const {
  if (!impl_->index->has_reference_payload()) {
    throw std::logic_error(
        "flashalign: this index embeds no reference sequence (built with "
        "--idx-no-seq); rebuild it without --idx-no-seq to read sequence from it");
  }
  const std::vector<std::string>& names = impl_->index->chromosome_names();
  const auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end()) {
    throw std::out_of_range(
        "flashalign: unknown reference name: " + name);
  }
  const std::size_t which =
      static_cast<std::size_t>(found - names.begin());
  {
    std::lock_guard<std::mutex> lock(impl_->reference_mutex);
    if (!impl_->reference_unpacked) {
      if (!impl_->index->reference_chromosomes_u8(impl_->reference, 0)) {
        throw std::logic_error(
            "flashalign: failed to unpack the index's embedded reference");
      }
      impl_->reference_unpacked = true;
    }
  }
  if (which >= impl_->reference.size()) return {};
  const std::vector<std::uint8_t>& contig = impl_->reference[which];
  const std::int64_t length = static_cast<std::int64_t>(contig.size());
  // Clamped as mappy's Aligner.seq().
  const std::int64_t lo = start < 0 ? 0 : std::min(start, length);
  const std::int64_t hi = end < 0 ? length : std::min(end, length);
  if (lo >= hi) return {};
  std::string out;
  out.resize(static_cast<std::size_t>(hi - lo));
  for (std::int64_t i = lo; i < hi; ++i) {
    const std::uint8_t code = contig[static_cast<std::size_t>(i)];
    out[static_cast<std::size_t>(i - lo)] =
        code < 4 ? "ACGT"[code] : 'N';
  }
  return out;
}

struct Aligner::Impl {
  // What config() reports: the resolved configuration, unaffected by a MapRequest.
  Config public_config;
  // The Config the caller supplied. A MapRequest re-resolves from it with only cs and md
  // replaced, and never changes it.
  Config requested_config;
  // k and s passed to the resolver at the last install.
  int index_k = -1;
  int index_s = -1;
  fa::cpu::api::LongReadAligner aligner;
  std::unordered_map<std::string, std::int64_t> reference_lengths;
  // Held shared while mapping with the installed cs/MD request, exclusive to change it.
  mutable std::shared_mutex request_mutex;
  CsMode installed_cs = CsMode::None;
  bool installed_md = false;
  // A MapRequest with its unset fields taken from the configured Config.
  struct Demand {
    CsMode cs;
    bool md;
  };

  Impl(Config public_view, Config requested, int k, int s,
       fa::cpu::api::LongReadAligner value)
      : public_config(std::move(public_view)),
        requested_config(std::move(requested)),
        index_k(k),
        index_s(s),
        aligner(std::move(value)) {
    const auto names = aligner.chromosome_names();
    const auto lengths = aligner.chromosome_lengths();
    for (std::size_t i = 0; i < names.size() && i < lengths.size(); ++i)
      reference_lengths.emplace(names[i], lengths[i]);
    installed_cs = public_config.cs;
    installed_md = public_config.emit_md;
  }

  // Sets target_len on every record and is_secondary on secondary hypotheses and their
  // supplementary segments.
  void annotate(Alignment& record, bool secondary) const {
    record.is_secondary = secondary;
    if (!record.chromosome.empty()) {
      const auto found = reference_lengths.find(record.chromosome);
      record.target_len =
          found == reference_lengths.end() ? 0 : found->second;
    }
    for (Alignment& supplementary : record.supplementary)
      annotate(supplementary, secondary);
    for (Alignment& alternative : record.secondary)
      annotate(alternative, true);
  }

  // Caller holds request_mutex, shared or exclusive.
  Demand effective(const MapRequest& request) const {
    return {request.cs.value_or(requested_config.cs),
            request.md.value_or(requested_config.emit_md)};
  }

  bool demand_installed(const Demand& demand) const {
    return installed_cs == demand.cs && installed_md == demand.md;
  }

  // Caller holds request_mutex exclusively. Reconfigures the engine for `demand` without
  // changing requested_config or public_config.
  void install_demand(const Demand& demand) {
    if (demand_installed(demand)) return;
    Config next = requested_config;
    next.cs = demand.cs;
    next.emit_md = demand.md;
    aligner.reconfigure(resolve_config(next, index_k, index_s));
    installed_cs = demand.cs;
    installed_md = demand.md;
  }

  // `read_name_hash` is tie_name_hash() of the read name, or 0.
  Alignment map(const std::string& read, const MapRequest& request,
                std::uint32_t read_name_hash) {
    Alignment result;
    {
      std::shared_lock<std::shared_mutex> lock(request_mutex);
      if (demand_installed(effective(request))) {
        result = aligner.align(read, read_name_hash);
      } else {
        lock.unlock();
        std::unique_lock<std::shared_mutex> write(request_mutex);
        install_demand(effective(request));
        result = aligner.align(read, read_name_hash);
      }
    }
    annotate(result, false);
    return result;
  }
};

Aligner::Aligner(Index index, Config config) {
  if (!index.impl_ || !index.impl_->index || index.impl_->index->empty())
    throw std::invalid_argument("flashalign: Aligner requires a non-empty index");
  const int index_k = index.impl_->index->k();
  const int index_s = index.impl_->index->build_syncmer_s();
  auto resolved = resolve_config(config, index_k, index_s);
  auto aligner = fa::cpu::api::LongReadAligner::from_shared_index(
      index.impl_->index, resolved);
  impl_ = std::make_unique<Impl>(
      to_public_config(resolved, config.preset), config, index_k, index_s,
      std::move(aligner));
}

Aligner::Aligner(
    const std::vector<std::pair<std::string, std::string>>& references,
    Config config) {
  auto resolved = resolve_config(config);
  fa::cpu::api::LongReadAligner aligner(
      genome_map(references), resolved.index.k,
      resolved.common.min_support, resolved.index.syncmer_s,
      resolved.index.syncmer_downsample, resolved.common.num_threads);
  aligner.reconfigure(resolved);
  impl_ = std::make_unique<Impl>(
      to_public_config(resolved, config.preset), config, -1, -1,
      std::move(aligner));
}

Aligner::~Aligner() = default;
Aligner::Aligner(Aligner&&) noexcept = default;
Aligner& Aligner::operator=(Aligner&&) noexcept = default;

Alignment Aligner::map(const std::string& read) const {
  return map(read, MapRequest{});
}
std::vector<Alignment> Aligner::map_batch(
    const std::vector<std::string>& reads) const {
  return map_batch(reads, MapRequest{});
}

Alignment Aligner::map(
    const std::string& read, const MapRequest& request) const {
  // No name, as minimap2 without a qname.
  return impl_->map(read, request, 0);
}

std::vector<Alignment> Aligner::map_batch(
    const std::vector<std::string>& reads,
    const MapRequest& request) const {
  std::vector<Alignment> results;
  {
    std::shared_lock<std::shared_mutex> lock(impl_->request_mutex);
    if (impl_->demand_installed(impl_->effective(request))) {
      results = impl_->aligner.align_batch(reads);
    } else {
      lock.unlock();
      std::unique_lock<std::shared_mutex> write(impl_->request_mutex);
      impl_->install_demand(impl_->effective(request));
      results = impl_->aligner.align_batch(reads);
    }
  }
  for (Alignment& result : results) impl_->annotate(result, false);
  return results;
}

std::string Aligner::paf(
    const std::string& name, const std::string& read,
    bool with_cigar) const {
  // The name seeds the tie-break as in the CLI, so the rows match the CLI's for that read.
  Alignment result =
      impl_->map(read, MapRequest{}, fa::cpu::lr::tie_name_hash(name));
  if (!result.mapped()) return {};
  fa::cpu::io::FastxRecord record;
  record.name = name;
  record.seq = read;
  std::ostringstream output;
  fa::cpu::output::write_paf_record(
      output, record, result, impl_->reference_lengths, with_cigar);
  for (const auto& supplementary : result.supplementary) {
    if (!supplementary.mapped()) continue;
    fa::cpu::output::write_paf_record(
        output, record, supplementary, impl_->reference_lengths, with_cigar,
        fa::cpu::output::EmittedRole::Supplementary);
  }
  // No secondary rows, as in the CLI's default output.
  return output.str();
}
Config Aligner::config() const {
  std::shared_lock<std::shared_mutex> lock(impl_->request_mutex);
  return impl_->public_config;
}
void Aligner::reconfigure(Config config) {
  std::unique_lock<std::shared_mutex> lock(impl_->request_mutex);
  const int index_k = impl_->aligner.config().index.k;
  const int index_s = impl_->aligner.get_index_syncmer_s();
  auto resolved = resolve_config(config, index_k, index_s);
  impl_->aligner.reconfigure(resolved);
  impl_->requested_config = config;
  impl_->index_k = index_k;
  impl_->index_s = index_s;
  impl_->public_config = to_public_config(resolved, config.preset);
  impl_->installed_cs = impl_->public_config.cs;
  impl_->installed_md = impl_->public_config.emit_md;
}
std::vector<std::string> Aligner::reference_names() const {
  return impl_->aligner.chromosome_names();
}
std::vector<std::int64_t> Aligner::reference_lengths() const {
  return impl_->aligner.chromosome_lengths();
}

}  // namespace flashalign
