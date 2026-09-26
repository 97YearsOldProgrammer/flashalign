// Internal façade over the long-read engine for the CLI, the public API and the Python
// binding. It hides the engine and backend headers behind a PIMPL.
#pragma once

#include "../core/types.h"
#include "../options/resolve.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Forward declaration; repeating the alias is valid alongside index/index.h.
namespace fa {
namespace cpu {
class FaixIndex;
using SeedIndex = FaixIndex;
} // namespace cpu
} // namespace fa

namespace fa {
namespace cpu {
namespace api {

using Alignment = ::fa::cpu::AlignResult;

using ::fa::cpu::options::BackendMode;
using ::fa::cpu::options::IndexMetadata;
using ::fa::cpu::options::ResolveRequest;
using ::fa::cpu::options::ResolvedMapOptions;
using ::fa::cpu::options::UserOverrides;
using ::fa::cpu::options::resolve_options;

using ::fa::cpu::options::CommonOptions;
using ::fa::cpu::options::IndexIdentity;
using ::fa::cpu::options::ResolvedOptions;
using ::fa::cpu::lr::DnaLongOptions;
using RnaLongOptions = ::fa::cpu::lr::rna::RnaLongOptions;

// Preset queries.
struct PresetSeeding {
  int k = 0;
  int syncmer_s = 0;
};
PresetSeeding resolve_preset_seeding(std::string_view preset);
bool is_hifi_preset(std::string_view preset);
bool is_rna_preset(std::string_view preset);
bool is_rna_hifi_preset(std::string_view preset);
int native_splice_k();
bool native_splice_k_supported(int k);
bool preset_is_valid(std::string_view preset);
std::string accepted_preset_names();
int preset_cli_chain_max_gap(std::string_view preset);

// What a reference argument holds, judged from its content (reference_kind).
enum class ReferenceKind { Index, Sequences };

// A FlashAlign index by its magic, including one from a development build that the loader
// then refuses by name; otherwise FASTA or FASTQ, plain or gzip, by the first record
// character after decompression. "-" and other inputs that are not regular files are
// Sequences, since they cannot be inspected without consuming them. Throws
// std::runtime_error naming `path` for a file that cannot be read, is empty, or is neither.
ReferenceKind reference_kind(const std::string &path);

// What an index file's header says (LongReadAligner::probe_index_header).
struct IndexProbe {
  // The file starts with an index magic, current or from a development build.
  bool is_index = false;
  // Why this build cannot use the index; empty when it can. The fields below are set only
  // when it can.
  std::string error;
  int k = 0;
  int syncmer_s = 0;
  std::uint64_t reference_sequences = 0;
  bool has_reference = false;
  // A `flashalign index -I` file, mapped one part at a time.
  bool multipart = false;
  // The preset the index was built with, or "" when it records none. A run without -x
  // uses it.
  std::string preset;
};

// Byte span of one part of a multi-part index; both 0 means the whole file.
struct IndexImageSpan {
  std::uint64_t offset = 0;
  std::uint64_t bytes = 0;
};

class WindowedAlignSession;

// A batch returned by WindowedAlignSession::collect(): the submitted reads, given back, and
// one result per read in submission order.
struct AlignedReadBatch {
  std::vector<std::string> reads;
  std::vector<Alignment> results;
};

class LongReadAligner {
public:
  LongReadAligner(LongReadAligner &&) noexcept;
  LongReadAligner &operator=(LongReadAligner &&) noexcept;
  ~LongReadAligner();
  LongReadAligner(const LongReadAligner &) = delete;
  LongReadAligner &operator=(const LongReadAligner &) = delete;

  // Builds the index over `genome` ({name: sequence}).
  explicit LongReadAligner(std::unordered_map<std::string, std::string> genome,
                           int k = 15, int min_support = 3,
                           int chain_syncmer_s = 9,
                           int chain_syncmer_downsample = 1,
                           int build_threads = 0);
  // Loads a prebuilt .faix; k and s are taken from the index (LongReadEngine::from_index).
  // `n_threads` <= 0 uses the default. With `reference_bases_needed` false (DNA map-only
  // PAF) the reference is not unpacked. `image` selects one part of a multi-part index.
  static LongReadAligner
  from_index(std::unordered_map<std::string, std::string> genome,
             const std::string &index_path, int min_support = 3,
             int chain_syncmer_downsample = 1, int n_threads = 0,
             bool reference_bases_needed = true, IndexImageSpan image = {});
  // Aligner over an already-loaded index, which may be shared.
  static LongReadAligner from_shared_index(std::shared_ptr<SeedIndex> index,
                                           const ResolvedOptions& config);
  // Reads only the file header. Never throws; a path that cannot be opened is not an index.
  static IndexProbe probe_index_header(const std::string &path);

  ResolvedOptions config() const;
  void reconfigure(const ResolvedOptions& config);

  int get_index_syncmer_s() const;

  // `read_name_hash` is tie_name_hash() of the read name, used to break exact ties; 0 for no
  // name.
  Alignment align(const std::string &read,
                  std::uint32_t read_name_hash = 0) const;
  // `read_name_hashes` is null or holds one hash per read.
  std::vector<Alignment>
  align_batch(const std::vector<std::string> &reads,
              const std::vector<std::uint32_t> *read_name_hashes = nullptr) const;

  // Session that maps up to `window` batches at once with the configured thread count.
  // `helper_threads` of those threads are not started here; the caller's own threads lend
  // time through WindowedAlignSession::try_help. At least one worker stays the session's.
  // The session borrows this aligner: destroy it first, and do not reconfigure meanwhile.
  // Results equal align_batch for any window size.
  std::unique_ptr<WindowedAlignSession>
  make_windowed_session(int window, int helper_threads = 0) const;

  std::vector<std::string> chromosome_names() const;
  std::vector<int64_t> chromosome_lengths() const;

private:
  friend class WindowedAlignSession;
  LongReadAligner();
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Ordered windowed mapping session, owned by one thread. submit() queues a batch and
// returns; collect() waits for the oldest batch and rethrows its first per-read exception.
// Keep in_flight() below the window before each submit.
class WindowedAlignSession {
public:
  ~WindowedAlignSession();
  WindowedAlignSession(const WindowedAlignSession &) = delete;
  WindowedAlignSession &operator=(const WindowedAlignSession &) = delete;

  // `read_name_hashes` is empty or holds one hash per read.
  void submit(std::vector<std::string> reads,
              std::vector<std::uint32_t> read_name_hashes = {});
  AlignedReadBatch collect();
  std::size_t in_flight() const;

  // For helper threads only: maps at most one read and returns whether one was available.
  // `tid` must lie in [thread_count - helper_threads, thread_count) after clamping. Safe to
  // call concurrently with submit/collect and with other tids.
  bool try_help(int tid);

private:
  friend class LongReadAligner;
  WindowedAlignSession(const LongReadAligner &aligner, int window,
                       int helper_threads);
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace api
} // namespace cpu
} // namespace fa
