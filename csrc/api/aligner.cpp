// The only translation unit that includes engine/aligner.h, so the engine templates are
// compiled once, here. Every method forwards to LongReadEngine.
#include "aligner.h"

#include "../core/sequence.h"
#include "../engine/aligner.h"
#include "../index/faix.h"
#include "../options/presets.h"
#include "../rna/backend.h"

#include <cstddef>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace fa {
namespace cpu {
namespace api {

struct LongReadAligner::Impl {
  // Mapping calls hold it shared; reconfiguration holds it exclusively.
  mutable std::shared_mutex state_mutex;

  std::unique_ptr<::fa::cpu::lr::LongReadEngine> a;
  // Stateless backend tag chosen from the config; align calls std::visit it.
  std::variant<::fa::cpu::lr::dna::DnaBackend, ::fa::cpu::lr::rna::RnaBackend>
      backend;

  void select_backend() {
    if (a && a->config().is_rna())
      backend.emplace<::fa::cpu::lr::rna::RnaBackend>();
    else
      backend.emplace<::fa::cpu::lr::dna::DnaBackend>();
  }

  template <class Function>
  decltype(auto) read(Function &&function) const {
    std::shared_lock<std::shared_mutex> lock(state_mutex);
    return std::forward<Function>(function)();
  }

  template <class Function>
  decltype(auto) write(Function &&function) {
    std::unique_lock<std::shared_mutex> lock(state_mutex);
    return std::forward<Function>(function)();
  }
};

PresetSeeding resolve_preset_seeding(std::string_view preset) {
  const ::fa::cpu::options::PresetSeeding s =
      ::fa::cpu::options::resolve_preset_seeding(preset);
  return PresetSeeding{s.k, s.syncmer_s};
}
bool is_hifi_preset(std::string_view preset) {
  return ::fa::cpu::options::is_hifi_preset(preset);
}
bool is_assembly_preset(std::string_view preset) {
  return ::fa::cpu::options::is_assembly_preset(preset);
}
bool is_rna_preset(std::string_view preset) {
  return ::fa::cpu::options::is_rna_preset(preset);
}
bool is_rna_hifi_preset(std::string_view preset) {
  return ::fa::cpu::options::is_rna_hifi_preset(preset);
}
bool preset_is_valid(std::string_view preset) {
  return ::fa::cpu::options::preset_is_valid(preset);
}
std::string accepted_preset_names() {
  return ::fa::cpu::options::accepted_preset_names();
}
LongReadAligner::LongReadAligner() = default;
LongReadAligner::LongReadAligner(LongReadAligner &&) noexcept = default;
LongReadAligner &LongReadAligner::operator=(LongReadAligner &&) noexcept = default;
LongReadAligner::~LongReadAligner() = default;

LongReadAligner::LongReadAligner(
    std::unordered_map<std::string, std::string> genome, int k,
    int min_support, int chain_syncmer_s, int chain_syncmer_downsample,
    int build_threads)
    : impl_(std::make_unique<Impl>()) {
  impl_->a = std::make_unique<::fa::cpu::lr::LongReadEngine>(
      std::move(genome), k, min_support, chain_syncmer_s,
      chain_syncmer_downsample, build_threads);
  impl_->select_backend();
}

LongReadAligner LongReadAligner::from_index(
    std::unordered_map<std::string, std::string> genome,
    const std::string &index_path, int min_support,
    int chain_syncmer_downsample, int n_threads,
    bool reference_bases_needed, IndexImageSpan image) {
  LongReadAligner self;
  self.impl_ = std::make_unique<Impl>();
  self.impl_->a = ::fa::cpu::lr::LongReadEngine::from_index(
      std::move(genome), index_path, min_support,
      chain_syncmer_downsample, n_threads, reference_bases_needed,
      image.offset, image.bytes);
  self.impl_->select_backend();
  return self;
}

LongReadAligner LongReadAligner::from_shared_index(
    std::shared_ptr<SeedIndex> index,
    const ResolvedOptions& config) {
  LongReadAligner self;
  self.impl_ = std::make_unique<Impl>();
  self.impl_->a =
      ::fa::cpu::lr::LongReadEngine::from_shared_index(
          std::move(index), config);
  self.impl_->select_backend();
  return self;
}

IndexProbe LongReadAligner::probe_index_header(const std::string &path) {
  const ::fa::cpu::FaixProbe probe = ::fa::cpu::probe_faix_header(path);
  IndexProbe out;
  out.is_index = probe.is_index;
  if (!probe.is_index)
    return out;
  if (!probe.status) {
    out.error = ::fa::cpu::faix_load_error_message(probe.status);
    return out;
  }
  const ::fa::cpu::FaixHeader &h = probe.header;
  out.k = static_cast<int>(h.k);
  out.syncmer_s = static_cast<int>(h.syncmer_s);
  out.reference_sequences = h.contig_count;
  out.has_reference = (h.flags & ::fa::cpu::kFaixFlagHasReference) != 0;
  out.multipart = (h.flags & ::fa::cpu::kFaixFlagMultipart) != 0;
  out.preset = ::fa::cpu::faix_header_preset(h);
  return out;
}

ResolvedOptions LongReadAligner::config() const {
  return impl_->read([&] { return impl_->a->config(); });
}
void LongReadAligner::reconfigure(const ResolvedOptions& config) {
  impl_->write([&] {
    impl_->a->reconfigure(config);
    impl_->select_backend();
  });
}

int LongReadAligner::get_index_syncmer_s() const {
  return impl_->read([&] { return impl_->a->get_index_syncmer_s(); });
}
Alignment LongReadAligner::align(const std::string &read,
                                 std::uint32_t read_name_hash) const {
  return impl_->read([&] {
    return std::visit(
        [&](auto backend) {
          using Backend = decltype(backend);
          return impl_->a->template align<Backend>(read, read_name_hash);
        },
        impl_->backend);
  });
}

std::vector<Alignment>
LongReadAligner::align_batch(
    const std::vector<std::string> &reads,
    const std::vector<std::uint32_t> *read_name_hashes) const {
  return impl_->read([&] {
    return std::visit(
        [&](auto backend) {
          using Backend = decltype(backend);
          return impl_->a->template align_batch<Backend>(reads,
                                                         read_name_hashes);
        },
        impl_->backend);
  });
}

// Type-erases the backend of engine/aligner.h's WindowedAlignSession: one virtual call per
// batch.
struct WindowedAlignSession::Impl {
  struct Arm {
    virtual ~Arm() = default;
    virtual void submit(std::vector<std::string> reads,
                        std::vector<std::uint32_t> read_name_hashes) = 0;
    virtual AlignedReadBatch collect() = 0;
    virtual std::size_t in_flight() const = 0;
    virtual bool try_help(int tid) = 0;
  };

  template <class Backend> struct BackendArm final : Arm {
    ::fa::cpu::lr::WindowedAlignSession<Backend> session;

    BackendArm(const ::fa::cpu::lr::LongReadEngine &engine, int n_threads,
               int window, int n_helpers)
        : session(engine, n_threads, window, n_helpers) {}

    void submit(std::vector<std::string> reads,
                std::vector<std::uint32_t> read_name_hashes) override {
      session.submit(std::move(reads), std::move(read_name_hashes));
    }
    AlignedReadBatch collect() override {
      auto batch = session.collect();
      AlignedReadBatch out;
      out.reads = std::move(batch->reads);
      out.results = std::move(batch->results);
      return out;
    }
    std::size_t in_flight() const override { return session.in_flight(); }
    bool try_help(int tid) override { return session.try_help(tid); }
  };

  std::unique_ptr<Arm> arm;
};

WindowedAlignSession::WindowedAlignSession(const LongReadAligner &aligner,
                                           int window, int helper_threads)
    : impl_(std::make_unique<Impl>()) {
  const ::fa::cpu::lr::LongReadEngine &engine = *aligner.impl_->a;
  // The thread count align_batch uses; helper threads are part of it.
  const int n_threads = engine.config().common.num_threads;
  std::visit(
      [&](auto backend) {
        using Backend = decltype(backend);
        impl_->arm = std::make_unique<Impl::BackendArm<Backend>>(
            engine, n_threads, window, helper_threads);
      },
      aligner.impl_->backend);
}

// Joins the session's workers.
WindowedAlignSession::~WindowedAlignSession() = default;

void WindowedAlignSession::submit(std::vector<std::string> reads,
                                  std::vector<std::uint32_t> read_name_hashes) {
  impl_->arm->submit(std::move(reads), std::move(read_name_hashes));
}
AlignedReadBatch WindowedAlignSession::collect() {
  return impl_->arm->collect();
}
std::size_t WindowedAlignSession::in_flight() const {
  return impl_->arm->in_flight();
}
bool WindowedAlignSession::try_help(int tid) {
  return impl_->arm->try_help(tid);
}

std::unique_ptr<WindowedAlignSession>
LongReadAligner::make_windowed_session(int window, int helper_threads) const {
  return std::unique_ptr<WindowedAlignSession>(
      new WindowedAlignSession(*this, window, helper_threads));
}

std::vector<std::string> LongReadAligner::chromosome_names() const {
  return impl_->read([&] { return impl_->a->chromosome_names(); });
}
std::vector<int64_t> LongReadAligner::chromosome_lengths() const {
  return impl_->read([&] { return impl_->a->chromosome_lengths(); });
}

} // namespace api
} // namespace cpu
} // namespace fa
