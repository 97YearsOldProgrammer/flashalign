// Instantiates the engine templates for both lanes behind the façade; the reference-kind probe.
#include "api.h"

#include "../engine/engine.h"
#include "../index/format.h"
#include "../index/index.h"
#include "../io/io.h"
#include "../options/options.h"
#include "../rna/rna.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <cerrno>
#include <cstddef>
#include <cstring>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

/**** aligner ****/
// The only translation unit that instantiates the engine templates (engine/engine.h), so
// they are compiled once, here. Every method forwards to LongReadEngine.

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
bool is_assembly_preset(std::string_view preset) {
  return ::fa::cpu::options::is_assembly_preset(preset);
}
bool is_overlap_preset(std::string_view preset) {
  return ::fa::cpu::options::is_overlap_preset(preset);
}
bool lane_runs_stage(LaneStage stage, std::string_view preset) {
  return ::fa::cpu::options::lane_runs_stage(stage, preset);
}
bool is_rna_preset(std::string_view preset) {
  return ::fa::cpu::options::is_rna_preset(preset);
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
    int build_threads, bool packed_reference)
    : impl_(std::make_unique<Impl>()) {
  impl_->a = std::make_unique<::fa::cpu::lr::LongReadEngine>(
      std::move(genome), k, min_support, chain_syncmer_s,
      chain_syncmer_downsample, build_threads, packed_reference);
  impl_->select_backend();
}

LongReadAligner LongReadAligner::from_index(
    std::unordered_map<std::string, std::string> genome,
    const std::string &index_path, int min_support,
    int chain_syncmer_downsample, int n_threads,
    bool reference_bases_needed, IndexImageSpan image, bool packed_reference) {
  LongReadAligner self;
  self.impl_ = std::make_unique<Impl>();
  self.impl_->a = ::fa::cpu::lr::LongReadEngine::from_index(
      std::move(genome), index_path, min_support,
      chain_syncmer_downsample, n_threads, reference_bases_needed,
      image.offset, image.bytes, packed_reference);
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
DnaReadGeometry LongReadAligner::dna_read_geometry(int read_length) const {
  return impl_->read([&] { return impl_->a->dna_read_geometry(read_length); });
}
DnaChainGeometry LongReadAligner::dna_chain_geometry() const {
  return impl_->read([&] { return impl_->a->dna_chain_geometry(); });
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
    const std::vector<std::uint32_t> *read_name_hashes,
    const std::vector<int> *self_contigs,
    const std::vector<int> *name_ranks) const {
  return impl_->read([&] {
    return std::visit(
        [&](auto backend) {
          using Backend = decltype(backend);
          return impl_->a->template align_batch<Backend>(
              reads, read_name_hashes, self_contigs, name_ranks);
        },
        impl_->backend);
  });
}

// Type-erases the backend of engine/engine.h's WindowedAlignSession: one virtual call per
// batch.
struct WindowedAlignSession::Impl {
  struct Arm {
    virtual ~Arm() = default;
    virtual void submit(std::vector<std::string> reads,
                        std::vector<std::uint32_t> read_name_hashes,
                        std::vector<int> self_contigs,
                        std::vector<int> name_ranks) = 0;
    virtual AlignedReadBatch collect() = 0;
    virtual std::size_t in_flight() const = 0;
    virtual std::size_t window() const = 0;
    virtual bool try_help(int tid) = 0;
  };

  template <class Backend> struct BackendArm final : Arm {
    ::fa::cpu::lr::WindowedAlignSession<Backend> session;

    BackendArm(const ::fa::cpu::lr::LongReadEngine &engine, int n_threads,
               int n_helpers, std::uint64_t max_bases)
        : session(engine, n_threads, n_helpers, max_bases) {}

    void submit(std::vector<std::string> reads,
                std::vector<std::uint32_t> read_name_hashes,
                std::vector<int> self_contigs,
                std::vector<int> name_ranks) override {
      session.submit(std::move(reads), std::move(read_name_hashes),
                     std::move(self_contigs), std::move(name_ranks));
    }
    AlignedReadBatch collect() override {
      auto batch = session.collect();
      AlignedReadBatch out;
      out.reads = std::move(batch->reads);
      out.results = std::move(batch->results);
      return out;
    }
    std::size_t in_flight() const override { return session.in_flight(); }
    std::size_t window() const override { return session.window(); }
    bool try_help(int tid) override { return session.try_help(tid); }
  };

  std::unique_ptr<Arm> arm;
};

WindowedAlignSession::WindowedAlignSession(const LongReadAligner &aligner,
                                           int helper_threads,
                                           std::uint64_t max_bases)
    : impl_(std::make_unique<Impl>()) {
  const ::fa::cpu::lr::LongReadEngine &engine = *aligner.impl_->a;
  // The thread count align_batch uses; helper threads are part of it.
  const int n_threads = engine.config().common.num_threads;
  std::visit(
      [&](auto backend) {
        using Backend = decltype(backend);
        impl_->arm = std::make_unique<Impl::BackendArm<Backend>>(
            engine, n_threads, helper_threads, max_bases);
      },
      aligner.impl_->backend);
}

// Joins the session's workers.
WindowedAlignSession::~WindowedAlignSession() = default;

void WindowedAlignSession::submit(std::vector<std::string> reads,
                                  std::vector<std::uint32_t> read_name_hashes,
                                  std::vector<int> self_contigs,
                                  std::vector<int> name_ranks) {
  impl_->arm->submit(std::move(reads), std::move(read_name_hashes),
                     std::move(self_contigs), std::move(name_ranks));
}
AlignedReadBatch WindowedAlignSession::collect() {
  return impl_->arm->collect();
}
std::size_t WindowedAlignSession::in_flight() const {
  return impl_->arm->in_flight();
}
std::size_t WindowedAlignSession::window() const {
  return impl_->arm->window();
}
bool WindowedAlignSession::try_help(int tid) {
  return impl_->arm->try_help(tid);
}

std::unique_ptr<WindowedAlignSession>
LongReadAligner::make_windowed_session(int helper_threads,
                                       std::uint64_t max_bases) const {
  return std::unique_ptr<WindowedAlignSession>(
      new WindowedAlignSession(*this, helper_threads, max_bases));
}

std::vector<std::string> LongReadAligner::chromosome_names() const {
  return impl_->read([&] { return impl_->a->chromosome_names(); });
}
std::vector<int64_t> LongReadAligner::chromosome_lengths() const {
  return impl_->read([&] { return impl_->a->chromosome_lengths(); });
}
int LongReadAligner::contig_id(std::string_view name) const {
  return impl_->read([&] { return impl_->a->contig_id(name); });
}
int LongReadAligner::contig_dual_rank(std::string_view name,
                                      std::int64_t length) const {
  return impl_->read([&] { return impl_->a->contig_dual_rank(name, length); });
}
} // namespace api
} // namespace cpu
} // namespace fa

/**** reference_kind ****/

namespace fa {
namespace cpu {
namespace api {

ReferenceKind reference_kind(const std::string &path) {
  if (path == "-")
    return ReferenceKind::Sequences;
  const auto refuse = [&](const std::string &why) {
    return std::runtime_error(path + ": " + why);
  };
  // stat() first: opening a named pipe would block, and closing it could end its writer.
  struct stat st;
  if (::stat(path.c_str(), &st) != 0)
    throw refuse(std::strerror(errno));
  if (S_ISDIR(st.st_mode))
    throw refuse("is a directory");
  if (!S_ISREG(st.st_mode))
    return ReferenceKind::Sequences;
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw refuse(std::strerror(errno));

  char magic[8] = {};
  if (::pread(fd, magic, sizeof(magic), 0) == static_cast<ssize_t>(sizeof(magic)) &&
      (std::memcmp(magic, kFaixMagic, sizeof(magic)) == 0 ||
       std::memcmp(magic, kFaixDevelopmentMagic, sizeof(magic)) == 0)) {
    ::close(fd);
    return ReferenceKind::Index;
  }

  // zlib reads plain input as is. gzclose() closes the descriptor.
  gzFile in = gzdopen(fd, "r");
  if (in == nullptr) {
    ::close(fd);
    throw refuse("cannot read the file");
  }
  char buffer[4096];
  char first = 0;
  while (first == 0) {
    const int n = gzread(in, buffer, sizeof(buffer));
    if (n < 0) {
      int code = 0;
      const std::string why = gzerror(in, &code);
      gzclose(in);
      throw refuse("cannot decompress: " + why);
    }
    if (n == 0) {
      // A gzip file cut before any content: say so before the refusal below.
      int code = Z_OK;
      (void)gzerror(in, &code);
      if (code == Z_BUF_ERROR)
        ::fa::cpu::io::warn_truncated_gzip(path);
      break;
    }
    for (int i = 0; i < n; ++i) {
      const char c = buffer[i];
      if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
        first = c;
        break;
      }
    }
  }
  gzclose(in);
  if (first == 0)
    throw refuse("the file is empty");
  if (first != '>' && first != '@')
    throw refuse("not a FASTA or FASTQ file, nor a FlashAlign index");
  return ReferenceKind::Sequences;
}

} // namespace api
} // namespace cpu
} // namespace fa
