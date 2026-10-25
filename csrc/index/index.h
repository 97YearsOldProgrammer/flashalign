// The seed index as the lanes see it: FaixIndex, the reference view, the occurrence distribution, the build and storage entry points.
#pragma once

#include "format.h"
#include "seed.h"                        // htable_valid_k, closed-syncmer reference extraction
#include "../core/core.h"               // detail::splitmix64
#include "../threading/threading.h"     // default_thread_count

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

/**** faix ****/
// FaixIndex, the query-time seed index: load, save and lookup over the format in format.h.
// build_faix_index(), in the build section below, produces one through from_compact_owned().

namespace fa { namespace cpu {

// What the header of a file says: one open and one short read, no mapping. `is_index` is
// true when the file starts with an index magic, current or from a development build;
// `status` then says whether this build can read it, and `header` holds its fields when it
// can.
struct FaixProbe {
    bool is_index = false;
    FaixLoadStatus status;
    FaixHeader header;
};

// Never throws. A path that cannot be opened, or holds fewer than 8 bytes, is not an index.
FaixProbe probe_faix_header(const std::string& path);

// The table of a multi-part container: everything before the part images, enough to name
// every contig and find each part's image. Success is `loaded()`.
struct FaixMultipartTable {
    FaixHeader header;
    std::vector<uint64_t> chr_offsets;      // global prefix sums, contig_count + 1
    std::vector<std::string> chrom_names;   // global, in contig-id order
    std::vector<FaixPartEntry> parts;       // in file and contig order

    bool loaded() const noexcept;
    // Bases in part `i`: the global offsets bracket its contig run.
    std::uint64_t part_bp(std::size_t i) const noexcept;
};

// Reads and checks a container's table, and the header of every part image against it.
FaixMultipartTable load_faix_multipart_table(const std::string& path,
                                             FaixLoadStatus* status = nullptr);

class FaixIndex {
public:
    FaixIndex() = default;
    ~FaixIndex();

    FaixIndex(const FaixIndex&) = delete;
    FaixIndex& operator=(const FaixIndex&) = delete;

    FaixIndex(FaixIndex&& other) noexcept;

    FaixIndex& operator=(FaixIndex&& other) noexcept;

    void swap(FaixIndex& other) noexcept;

    void release();

    bool empty() const;

    int k() const;
    int build_syncmer_s() const;
    // The preset this index was built with, or "" when it records none.
    std::string build_preset() const;
    uint64_t total_bp() const;
    uint64_t distinct_keys() const;
    // Inline: called per seed from dna/seeds.cpp (slice, slice_batch).
    uint64_t chrom_count() const { return header_.contig_count; }
    uint64_t memory_bytes() const;

    // Inline: called per seed from dna/seeds.cpp (slice, slice_batch).
    const uint64_t* chrom_offsets_data() const { return chr_offsets_; }
    // The record table (the shards' slices end to end), for enumerating stored keys.
    const uint32_t* tags_data() const;
    const uint64_t* payloads_data() const;
    const std::vector<std::string>& chromosome_names() const;

    bool has_reference_payload() const;

    // Materializes the embedded reference as one u8-per-base vector per contig, threaded.
    // `n_threads` <= 0 means threading::default_thread_count(). On failure `out` is empty.
    bool reference_chromosomes_u8(
        std::vector<std::vector<uint8_t>>& out,
        int n_threads = 0
    ) const;

    // Bases [lo, hi) of `contig` into `out`, one base code per byte, as
    // reference_chromosomes_u8 unpacks them; minimap2's mm_idx_getseq. Needs
    // has_reference_payload() and lo <= hi <= the contig's length.
    void reference_bases_u8(uint64_t contig, uint64_t lo, uint64_t hi,
                            uint8_t* out) const;

    // Attaches the contig names and the reference, one base code per byte and contig, in
    // contig order.
    bool set_reference_payload(
        std::vector<std::string> chrom_names,
        const std::vector<std::vector<uint8_t>>& contigs
    );

    // Names only (`flashalign index --idx-no-seq`): the index still names and measures every
    // contig, which map-only output needs, but has_reference_payload() is false, so CIGAR
    // output is refused. Any attached reference is dropped.
    bool set_chromosome_names(std::vector<std::string> chrom_names);

    // Attaches the contig names and the reference packed as in the file.
    bool set_packed_reference_payload(
        std::vector<std::string> chrom_names,
        std::vector<uint8_t>&& packed_reference
    );

    // Inline: called per seed from vote/vote.cpp and rna/vote.cpp.
    KmerPostingView lookup(uint64_t key) const {
        if (!tags_ || !shard_dir_ || header_.slot_count == 0) return {};
        return faix_sharded_lookup(
            tags_, payloads_, postings_, chr_offsets_, header_.contig_count,
            occupancy_bits_, rank_index_, shard_dir_, shard_mask_,
            static_cast<int>(header_.shard_bits), 2 * static_cast<int>(header_.k),
            header_.rank_block_bits, key);
    }

    // Batched lookup() with overlapped memory loads. `out` must have room for `n` views; it
    // is fully written, misses as {}.
    void lookup_batch(
        const uint64_t* keys,
        size_t n,
        KmerPostingView* out
    ) const;

    // The narrowing half of lookup_interval, for a caller that already holds a key's full
    // posting view. The span aliases the index's postings and lives as long as the index.
    // Inline: called per seed from rna/harvest.cpp.
    KmerPostingIntervalView narrow_interval(
        const KmerPostingView& full,
        int chr_idx,
        uint32_t chr_lo,
        uint32_t chr_hi
    ) const {
        if (!chr_offsets_ ||
            chr_idx < 0 ||
            chr_idx >= static_cast<int>(header_.contig_count) ||
            chr_hi <= chr_lo) {
            return {};
        }

        if (!full.found()) return {};
        if (full.count == 0) {
            return {full.positions, 0, full.occurrence};
        }

        const uint64_t chr_length =
            chr_offsets_[static_cast<size_t>(chr_idx) + 1] -
            chr_offsets_[static_cast<size_t>(chr_idx)];
        const uint32_t local_lo =
            static_cast<uint32_t>(std::min<uint64_t>(chr_lo, chr_length));
        const uint32_t local_hi =
            static_cast<uint32_t>(std::min<uint64_t>(chr_hi, chr_length));
        if (local_hi <= local_lo) {
          return {full.positions, 0, full.occurrence};
        }
        // Bounds use z = 0, so both orientations at local_lo are in and both at local_hi out.
        const PackedRefPos lo_key =
            pack_ref_pos(static_cast<uint32_t>(chr_idx), local_lo, 0);
        const PackedRefPos hi_key =
            pack_ref_pos(static_cast<uint32_t>(chr_idx), local_hi, 0);
        const uint32_t lo = full.positions.lower_bound_packed(lo_key);
        const uint32_t hi = full.positions.subspan(lo, full.count - lo)
                                .lower_bound_packed(hi_key) +
                            lo;
        return {
            full.positions.subspan(lo, hi - lo),
            hi - lo,
            full.occurrence,
        };
    }

    KmerPostingIntervalView lookup_interval(
        uint64_t key,
        int chr_idx,
        uint32_t chr_lo,
        uint32_t chr_hi
    ) const;

    // Writes this index as one complete image to `out` and reports the byte count the
    // layout implies, computed before writing so a caller can check what it wrote. Refuses
    // an index without its contig names.
    bool write_image(FaixByteSink& out, uint64_t& expected_bytes) const;

    bool save(const std::string& path) const;

    // `n_threads` bounds the loader's parallel reads; <= 0 takes the default.
    static FaixIndex load(
        const std::string& path, FaixLoadStatus* status = nullptr,
        int n_threads = 0);

    // load() over the `image_bytes` starting at byte `image_offset` of `path`, one part of a
    // multi-part container. Both 0 means the whole file.
    static FaixIndex load_image(
        const std::string& path, uint64_t image_offset, uint64_t image_bytes,
        FaixLoadStatus* status = nullptr, int n_threads = 0);

    // `header` carries the build's k, s, shard and rank-block bits, counts and preset; the
    // names and the reference are attached afterwards.
    static FaixIndex
    from_compact_owned(FaixHeader header, std::vector<uint64_t> chr_offsets,
                       std::vector<uint8_t> occupancy_bits,
                       std::vector<uint64_t> rank_index,
                       std::vector<uint32_t> compact_tags,
                       std::vector<uint64_t> compact_payloads,
                       std::vector<PackedRefPos> compact_postings,
                       std::vector<FaixShardEntry> shard_dir);

private:
    void rebind_views_to_owned();

    // Attaches one name per contig and records the name table's size.
    bool set_names(std::vector<std::string> chrom_names);

    FaixHeader header_{};
    std::vector<uint64_t> chr_offsets_owned_;
    std::vector<uint32_t> tags_owned_;
    std::vector<uint64_t> payloads_owned_;
    std::vector<PackedRefPos> postings_owned_;
    std::vector<std::string> chrom_names_owned_;
    std::vector<uint8_t> occupancy_bits_owned_;
    std::vector<uint64_t> rank_index_owned_;
    std::vector<FaixShardEntry> shard_dir_owned_;
    std::vector<uint8_t> packed_reference_owned_;

    const uint64_t* chr_offsets_ = nullptr;
    const uint32_t* tags_ = nullptr;
    const uint64_t* payloads_ = nullptr;
    const PackedRefPos* postings_ = nullptr;
    const uint8_t* occupancy_bits_ = nullptr;
    const uint64_t* rank_index_ = nullptr;
    const FaixShardEntry* shard_dir_ = nullptr;
    uint64_t shard_mask_ = 0;
    const uint8_t* packed_reference_ = nullptr;

    void* mapped_base_ = nullptr;
    size_t mapped_size_ = 0;
};

}}  // namespace fa::cpu

/**** occ_stats ****/
// Occurrence distribution of a built index, and minimap2's mid-occurrence quantile over
// it. Every key's exact occurrence is stored: an inline singleton is 1, otherwise the
// payload's high word is the posting count.

namespace fa {
namespace cpu {

// Occurrences at or below this bound get a dense bucket; the rare tail above it is kept
// exactly in a std::map.
inline constexpr std::uint32_t kFaixOccHistogramBound = 1u << 16;

struct FaixOccDistribution {
  // counts[v - 1] == number of distinct keys whose occurrence is exactly v,
  // for 1 <= v <= kFaixOccHistogramBound.
  std::vector<std::uint64_t> counts;
  // Exact key counts for occurrence values above kFaixOccHistogramBound.
  std::map<std::uint32_t, std::uint64_t> tail;
  std::uint64_t n_distinct = 0;
  // Every key's occurrences summed: all the reference seeds.
  std::uint64_t total_postings = 0;
};

// One streaming pass over the compact record table (the shards laid end to end).
FaixOccDistribution faix_occ_distribution(const FaixIndex& index);

// minimap2's mm_idx_cal_max_occ(mi, f), evaluated from the cumulative counts. With the n
// distinct occurrences sorted ascending, minimap2 returns occ[(size_t)((1 - f) * n)] + 1:
// the smallest v whose cumulative count exceeds that 0-based rank. f <= 0 disables the
// quantile (INT32_MAX), as in minimap2.
std::int64_t faix_cal_max_occ(const FaixOccDistribution& dist,
                              double f);

} // namespace cpu
} // namespace fa

/**** faix_publish ****/

namespace fa {
namespace cpu {

// Sequential transactional file: bytes go to a private partial in the same directory,
// renamed into place only after an exact-size check and fsync. An existing destination
// survives any failure.
class FaixTransactionalFile final : public FaixByteSink {
public:
  explicit FaixTransactionalFile(const std::string& final_path);

  ~FaixTransactionalFile() override;

  FaixTransactionalFile(const FaixTransactionalFile&) = delete;
  FaixTransactionalFile& operator=(const FaixTransactionalFile&) = delete;

  bool opened() const noexcept;

  bool write_bytes(const void* data, uint64_t size) override;

  // Overwrites `size` bytes at `offset` within what has already been written; it never
  // extends the file.
  bool patch_bytes(uint64_t offset, const void* data, uint64_t size);

  uint64_t written() const noexcept;

  bool publish(uint64_t expected_bytes);

private:
  std::string final_path_;
  std::string partial_path_;
  int fd_ = -1;
  uint64_t written_ = 0;
  bool published_ = false;
};

} // namespace cpu
} // namespace fa

/**** faix_multipart_writer ****/
// Writer for the multi-part container of `flashalign index -I` (format.h). One
// transactional file in layout order: the header and global tables, then each part's image
// as the caller builds it (so only one part's index is held at a time), and finally the
// part table patched in at the front. Nothing appears at `path` until finish() publishes; a
// failure leaves the destination as it was.

namespace fa {
namespace cpu {

class FaixMultipartWriter {
public:
  // What the container declares up front. The parts' image fields are filled in as images
  // are appended; everything else is checked against each part as it arrives.
  struct Plan {
    int k = 0;
    int syncmer_s = 0;
    bool has_reference = true; // every part embeds its 4-bit reference
    std::string preset;        // recorded on the container as on the parts
    std::vector<uint64_t> chr_offsets; // global prefix sums, chrom_count + 1
    std::vector<std::string> names;    // global, in contig-id order
    std::vector<FaixPartEntry> parts;  // first_contig / contig_count set
  };

  FaixMultipartWriter(const std::string& path, Plan plan);

  bool opened() const noexcept;

  // The header, the global offsets, a zeroed part table and the global names. finish()
  // patches the part table.
  bool begin();

  // The next part in order: the index of exactly the contigs the plan gave it, with the
  // plan's names, lengths, seeding, reference policy and preset.
  bool append_part(const FaixIndex& part);

  // Every part appended: patch the part table, publish.
  bool finish();

private:
  Plan plan_;
  FaixTransactionalFile out_;
  uint64_t part_table_offset_ = 0;
  size_t next_part_ = 0;
  bool begun_ = false;
};

} // namespace cpu
} // namespace fa

/**** build ****/
// In-memory index builder: emit reference seeds into shards across threads, sort and pack
// each shard, and assemble a FaixIndex.

namespace fa { namespace cpu {

// Why a build returned an empty index; the counterpart of FaixLoadError.
enum class FaixBuildError : std::uint8_t {
    None,
    UnsupportedK,           // k outside [1, kFaixMaxK]
    NoSequences,            // empty reference
    TooManySequences,       // more than kFaixMaxContigCount contigs
    ContigTooLong,          // a contig above kFaixMaxContigBp
    InvalidChrOffsets,      // offsets not from 0 / not monotone / too few
    SeedPositionOverflow,   // seed emitter produced an out-of-range position
    IncompleteStream,       // streamed source did not deliver every contig once
    ReferenceTooLarge,      // a shard or a key needs more than 2^32 - 1 keys or postings
    Internal,               // the packed tables came out inconsistent
    PresetNameTooLong,      // cfg.preset does not fit the header's 16 bytes
};

struct FaixBuildStatus {
    FaixBuildError error = FaixBuildError::None;
    // Zero-based index of the offending reference sequence, when the error
    // names one (ContigTooLong, InvalidChrOffsets); otherwise 0.
    std::uint64_t sequence_index = 0;
    // The offending magnitude, when the error has one: the contig length for
    // ContigTooLong, k for UnsupportedK, the sequence count for TooManySequences.
    // Otherwise 0.
    std::uint64_t value = 0;

    explicit operator bool() const noexcept {
        return error == FaixBuildError::None;
    }
};

std::string faix_build_error_message(const FaixBuildStatus& status);

// Records a refusal in the optional `status` and returns the empty index.
inline FaixIndex faix_build_fail(
    FaixBuildStatus* status,
    FaixBuildError error,
    std::uint64_t sequence_index = 0,
    std::uint64_t value = 0
) {
    if (status) {
        status->error = error;
        status->sequence_index = sequence_index;
        status->value = value;
    }
    return {};
}

// Slots per shard are the next power of two at or above keys / kFaixLoadFactor.
static constexpr double kFaixLoadFactor = 0.50;
// The index has 1 << shard_bits shards, at least this many bits. The per-key tag (2k -
// shard_bits key bits plus the multi bit) must fit a uint32_t, so htable_dispatch_bucket_bits()
// raises it to 2k - 31 when needed, up to 16.
static constexpr int kFaixMinShardBits = 11;

struct FaixBuildConfig {
    int k = 17;
    int syncmer_s = 0;                 // closed-syncmer s-mer length (5 HiFi / 9 ONT)
    int build_threads = 1;              // <= 0: threading::default_thread_count()
    // The canonical preset name to record in the header, or empty for none. It is read back
    // by the aligner and does not affect the table.
    std::string preset;
};

int htable_effective_closed_syncmer_s(const FaixBuildConfig& cfg);

struct KmerPositionRecord {
    uint64_t key = 0;
    PackedRefPos pos = 0;
};
static_assert(sizeof(KmerPositionRecord) == 16,
              "KmerPositionRecord must stay 16 bytes");

bool htable_position_record_less(
    const KmerPositionRecord& a,
    const KmerPositionRecord& b
);

void htable_sort_position_records(
    std::vector<KmerPositionRecord>& records
);

// One key of a build shard while it is packed: a singleton carries its position in
// `payload`, flagged by kInlineKey in `key`; a multi-hit key carries count << 32 | begin
// into the shard's scratch postings.
struct KmerBucket {
    static constexpr uint64_t kInlineKey = 1ULL << 63;
    static constexpr uint64_t kBeginMask = 0xffffffffULL;

    uint64_t key = kEmptyKmerKey;
    uint64_t payload = 0;

    uint64_t key_value() const { return key & ~kInlineKey; }

    uint64_t begin() const {
      return inlined() ? payload : payload & kBeginMask;
    }

    uint32_t count() const {
      return inlined() ? 1u : static_cast<uint32_t>(payload >> 32);
    }

    bool inlined() const { return (key & kInlineKey) != 0; }

    // Records a singleton's position (`inline_position`, count 1) or a multi-hit key's
    // begin and count; false when they do not fit.
    bool set_payload(uint64_t begin, uint64_t count, bool inline_position) {
      if ((key & kInlineKey) != 0)
        return false;
      if (inline_position) {
        if (count != 1)
          return false;
        key |= kInlineKey;
        payload = begin;
        return true;
      }
      if (begin > kBeginMask || count > 0xffffffffULL)
        return false;
      payload = begin | (count << 32);
      return true;
    }
};
static_assert(sizeof(KmerBucket) == 16, "KmerBucket must stay compact");

struct KmerBuildScratchBucket {
    std::vector<KmerPositionRecord> records;
    std::vector<KmerBucket> groups;   // begin is bucket-local posting offset
    std::vector<PackedRefPos> postings;
};

int htable_build_threads(int requested_threads);

int htable_dispatch_bucket_bits(const FaixBuildConfig& cfg);

// Shard id: the low shard_bits of the invertible key mix, as the lookup computes it. The
// mix folds high key bits down, so shards stay balanced despite GC-skewed trailing bases.
uint32_t htable_dispatch_bucket(uint64_t key, uint32_t mask, int twok);

double htable_reference_seed_density_estimate(const FaixBuildConfig& cfg);

// Appends one emit batch to a bucket (under that bucket's lock), growing capacity by 1.25x
// rather than 2x. The per-bucket reservation from the closed-syncmer density 2/(k-s+1)
// underestimates real references, so most buckets outgrow it; a 2x step would leave the
// emit phase holding nearly twice the records' size. The later per-bucket sort makes the
// output independent of append order.
inline void htable_append_bucket_batch(std::vector<KmerPositionRecord>& out,
                                       std::vector<KmerPositionRecord>& batch) {
  const size_t want = out.size() + batch.size();
  if (want > out.capacity()) {
    out.reserve(std::max(want, out.capacity() + out.capacity() / 4));
  }
  out.insert(out.end(), batch.begin(), batch.end());
  batch.clear();
}

// `emit` is called as emit(key, pos, z); see emit_closed_syncmer_reference_seeds for z.
template <class BaseT, class Emit>
inline bool htable_emit_reference_seeds(
    const std::vector<BaseT>& chr,
    uint64_t base_offset,
    const FaixBuildConfig& cfg,
    Emit&& emit
) {
    ClosedSyncmerConfig scfg;
    scfg.k = cfg.k;
    scfg.s = cfg.syncmer_s;
    return emit_closed_syncmer_reference_seeds(
        chr, base_offset, scfg, std::forward<Emit>(emit));
}

template <class BaseT>
inline std::vector<uint64_t> build_chr_offsets_for_encoded_sequences(
    const std::vector<std::vector<BaseT>>& chr_encs
) {
    std::vector<uint64_t> chr_offsets(chr_encs.size() + 1, 0);
    for (size_t i = 0; i < chr_encs.size(); i++) {
        chr_offsets[i + 1] = chr_offsets[i] + static_cast<uint64_t>(chr_encs[i].size());
    }
    return chr_offsets;
}

}}  // namespace fa::cpu

namespace fa { namespace cpu {

FaixIndex htable_finish_faix_index_from_scratch(
    std::vector<KmerBuildScratchBucket>& scratch,
    const FaixBuildConfig& cfg,
    const std::vector<uint64_t>& chr_offsets,
    uint64_t chrom_count,
    int build_threads,
    FaixBuildStatus* status = nullptr
);

template <class BaseT>
inline FaixIndex build_faix_index(
    const std::vector<std::vector<BaseT>>& chr_encs,
    const FaixBuildConfig& cfg = FaixBuildConfig(),
    FaixBuildStatus* status = nullptr
) {
  if (status) *status = FaixBuildStatus();
  if (!htable_valid_k(cfg.k) || cfg.k > kFaixMaxK)
    return faix_build_fail(status, FaixBuildError::UnsupportedK, 0,
                           static_cast<uint64_t>(cfg.k));
  if (chr_encs.empty())
    return faix_build_fail(status, FaixBuildError::NoSequences);
  if (static_cast<uint64_t>(chr_encs.size()) > kFaixMaxContigCount) {
    return faix_build_fail(status, FaixBuildError::TooManySequences, 0,
                           static_cast<uint64_t>(chr_encs.size()));
  }
  for (size_t i = 0; i < chr_encs.size(); ++i) {
    if (static_cast<uint64_t>(chr_encs[i].size()) > kFaixMaxContigBp)
      return faix_build_fail(status, FaixBuildError::ContigTooLong, i,
                             static_cast<uint64_t>(chr_encs[i].size()));
  }

    const std::vector<uint64_t> chr_offsets = build_chr_offsets_for_encoded_sequences(chr_encs);
    const int build_threads = htable_build_threads(cfg.build_threads);
    const int dispatch_bits = htable_dispatch_bucket_bits(cfg);
    const uint32_t dispatch_count = 1u << dispatch_bits;
    const uint32_t dispatch_mask = dispatch_count - 1;
    const int twok = 2 * cfg.k;

    std::vector<KmerBuildScratchBucket> scratch(dispatch_count);

    if (build_threads <= 1 || chr_encs.size() <= 1) {
        for (size_t chr_idx = 0; chr_idx < chr_encs.size(); chr_idx++) {
            const auto& chr = chr_encs[chr_idx];
            const uint64_t base_offset = chr_offsets[chr_idx];
            const bool ok = htable_emit_reference_seeds(
                chr, base_offset, cfg,
                [&](uint64_t key, RefPos pos, uint32_t z) {
                  const uint32_t bucket =
                      htable_dispatch_bucket(key, dispatch_mask, twok);
                  scratch[bucket].records.push_back(
                      {key,
                       pack_ref_pos(static_cast<uint32_t>(chr_idx),
                                    static_cast<uint32_t>(pos - base_offset),
                                    z)});
                });
            if (!ok)
              return faix_build_fail(status,
                                     FaixBuildError::SeedPositionOverflow,
                                     chr_idx);
        }
    } else {
        // Threaded emit: each worker fills small thread-local per-bucket batches and
        // flushes a full batch into the shared bucket under its lock, so records exist only
        // once. Buckets are pre-reserved from the closed-syncmer density. The later
        // per-bucket sort makes the output independent of append order.
        constexpr size_t kEmitBatch = 512;
        const double density = htable_reference_seed_density_estimate(cfg);
        const uint64_t total_bp = chr_offsets.empty() ? 0 : chr_offsets.back();
        const size_t bucket_reserve = static_cast<size_t>(
            static_cast<double>(total_bp) * density /
            static_cast<double>(dispatch_count));
        if (bucket_reserve) {
            for (auto& b : scratch) b.records.reserve(bucket_reserve);
        }

        std::vector<std::mutex> bucket_locks(dispatch_count);
        std::vector<std::vector<std::vector<KmerPositionRecord>>> thread_batch(
            static_cast<size_t>(build_threads),
            std::vector<std::vector<KmerPositionRecord>>(dispatch_count));
        for (auto& tb : thread_batch)
            for (auto& v : tb) v.reserve(kEmitBatch);

        auto flush = [&](uint32_t bucket, std::vector<KmerPositionRecord>& batch) {
            std::lock_guard<std::mutex> lk(bucket_locks[static_cast<size_t>(bucket)]);
            htable_append_bucket_batch(
                scratch[static_cast<size_t>(bucket)].records, batch);
        };

        std::atomic<bool> invalid_pos{false};
        threading::parallel_for(
            build_threads, static_cast<int64_t>(chr_encs.size()),
            [&](int64_t chr_idx, int tid) {
                auto& batches = thread_batch[static_cast<size_t>(tid)];
                const auto& chr = chr_encs[static_cast<size_t>(chr_idx)];
                const uint64_t base_offset = chr_offsets[static_cast<size_t>(chr_idx)];
                const bool ok = htable_emit_reference_seeds(
                    chr, base_offset, cfg,
                    [&](uint64_t key, RefPos pos, uint32_t z) {
                      const uint32_t bucket =
                          htable_dispatch_bucket(key, dispatch_mask, twok);
                      auto& batch = batches[bucket];
                      batch.push_back(
                          {key, pack_ref_pos(
                                    static_cast<uint32_t>(chr_idx),
                                    static_cast<uint32_t>(pos - base_offset),
                                    z)});
                      if (batch.size() >= kEmitBatch)
                        flush(bucket, batch);
                    });
                if (!ok) {
                    invalid_pos.store(true, std::memory_order_relaxed);
                }
            });

        if (invalid_pos.load(std::memory_order_relaxed))
          return faix_build_fail(status, FaixBuildError::SeedPositionOverflow);

        // Drain the leftover (< kEmitBatch) tails. Serial, but cheap.
        for (auto& batches : thread_batch) {
            for (uint32_t bucket = 0; bucket < dispatch_count; ++bucket) {
                if (!batches[bucket].empty()) flush(bucket, batches[bucket]);
            }
        }
    }

    return htable_finish_faix_index_from_scratch(
        scratch, cfg, chr_offsets, chr_encs.size(), build_threads, status);
}

template <class NextEncodedSequence>
inline FaixIndex build_faix_index_streamed(
    const std::vector<uint64_t>& chr_offsets,
    const FaixBuildConfig& cfg,
    NextEncodedSequence&& next_encoded_sequence,
    FaixBuildStatus* status = nullptr
) {
  if (status) *status = FaixBuildStatus();
  if (!htable_valid_k(cfg.k) || cfg.k > kFaixMaxK)
    return faix_build_fail(status, FaixBuildError::UnsupportedK, 0,
                           static_cast<uint64_t>(cfg.k));
  if (chr_offsets.size() < 2 || chr_offsets.front() != 0)
    return faix_build_fail(status, FaixBuildError::InvalidChrOffsets);
  if (static_cast<uint64_t>(chr_offsets.size() - 1) > kFaixMaxContigCount)
    return faix_build_fail(status, FaixBuildError::TooManySequences, 0,
                           static_cast<uint64_t>(chr_offsets.size() - 1));
  for (size_t i = 1; i < chr_offsets.size(); ++i) {
    if (chr_offsets[i] < chr_offsets[i - 1])
      return faix_build_fail(status, FaixBuildError::InvalidChrOffsets, i - 1);
    if (chr_offsets[i] - chr_offsets[i - 1] > kFaixMaxContigBp)
      return faix_build_fail(status, FaixBuildError::ContigTooLong, i - 1,
                             chr_offsets[i] - chr_offsets[i - 1]);
  }

    const uint64_t chrom_count = static_cast<uint64_t>(chr_offsets.size() - 1);
    const int build_threads = htable_build_threads(cfg.build_threads);
    const int dispatch_bits = htable_dispatch_bucket_bits(cfg);
    const uint32_t dispatch_count = 1u << dispatch_bits;
    const uint32_t dispatch_mask = dispatch_count - 1;
    const int twok = 2 * cfg.k;

    std::vector<KmerBuildScratchBucket> scratch(dispatch_count);
    std::vector<uint8_t> seen(static_cast<size_t>(chrom_count), 0);
    std::mutex seen_mutex;
    uint64_t seen_count = 0;
    bool invalid = false;
    // `invalid` covers both a bad stream and a refused seed position; the latter is also
    // tracked on its own so the refusal names the right error.
    std::atomic<bool> bad_position{false};

    auto mark_sequence = [&](size_t chr_idx, const std::vector<uint8_t>& chr,
                             uint64_t& base_offset) -> bool {
        if (chr_idx >= static_cast<size_t>(chrom_count) || seen[chr_idx] != 0) {
            return false;
        }
        const uint64_t expected_len = chr_offsets[chr_idx + 1] - chr_offsets[chr_idx];
        if (expected_len != static_cast<uint64_t>(chr.size())) {
            return false;
        }
        seen[chr_idx] = 1;
        ++seen_count;
        base_offset = chr_offsets[chr_idx];
        return true;
    };

    auto emit_sequence_serial = [&](size_t chr_idx, const std::vector<uint8_t>& chr) -> bool {
        uint64_t base_offset = 0;
        if (!mark_sequence(chr_idx, chr, base_offset)) {
            invalid = true;
            return false;
        }
        const bool ok = htable_emit_reference_seeds(
            chr, base_offset, cfg, [&](uint64_t key, RefPos pos, uint32_t z) {
              const uint32_t bucket =
                  htable_dispatch_bucket(key, dispatch_mask, twok);
              scratch[bucket].records.push_back(
                  {key,
                   pack_ref_pos(static_cast<uint32_t>(chr_idx),
                                static_cast<uint32_t>(pos - base_offset), z)});
            });
        if (!ok) {
            invalid = true;
            bad_position.store(true, std::memory_order_relaxed);
        }
        return ok;
    };

    if (build_threads <= 1) {
        for (;;) {
            size_t chr_idx = 0;
            std::vector<uint8_t> chr;
            if (!next_encoded_sequence(chr_idx, chr)) {
                break;
            }
            if (!emit_sequence_serial(chr_idx, chr)) break;
        }
    } else {
        constexpr size_t kEmitBatch = 512;
        const double density = htable_reference_seed_density_estimate(cfg);
        const uint64_t total_bp = chr_offsets.back();
        const size_t bucket_reserve = static_cast<size_t>(
            static_cast<double>(total_bp) * density /
            static_cast<double>(dispatch_count));
        if (bucket_reserve) {
            for (auto& b : scratch) b.records.reserve(bucket_reserve);
        }

        std::vector<std::mutex> bucket_locks(dispatch_count);
        std::vector<std::vector<std::vector<KmerPositionRecord>>> thread_batch(
            static_cast<size_t>(build_threads),
            std::vector<std::vector<KmerPositionRecord>>(dispatch_count));
        for (auto& tb : thread_batch)
            for (auto& v : tb) v.reserve(kEmitBatch);

        auto flush = [&](uint32_t bucket, std::vector<KmerPositionRecord>& batch) {
            std::lock_guard<std::mutex> lk(bucket_locks[static_cast<size_t>(bucket)]);
            htable_append_bucket_batch(
                scratch[static_cast<size_t>(bucket)].records, batch);
        };

        std::atomic<bool> invalid_pos{false};
        threading::parallel_for(
            build_threads, build_threads,
            [&](int64_t /*worker*/, int tid) {
                auto& batches = thread_batch[static_cast<size_t>(tid)];
                for (;;) {
                    if (invalid_pos.load(std::memory_order_relaxed)) break;

                    size_t chr_idx = 0;
                    std::vector<uint8_t> chr;
                    if (!next_encoded_sequence(chr_idx, chr)) break;

                    uint64_t base_offset = 0;
                    {
                        std::lock_guard<std::mutex> seen_lk(seen_mutex);
                        if (!mark_sequence(chr_idx, chr, base_offset)) {
                            invalid_pos.store(true, std::memory_order_relaxed);
                            break;
                        }
                    }

                    const bool ok = htable_emit_reference_seeds(
                        chr, base_offset, cfg,
                        [&](uint64_t key, RefPos pos, uint32_t z) {
                          const uint32_t bucket =
                              htable_dispatch_bucket(key, dispatch_mask, twok);
                          auto& batch = batches[bucket];
                          batch.push_back(
                              {key, pack_ref_pos(static_cast<uint32_t>(chr_idx),
                                                 static_cast<uint32_t>(
                                                     pos - base_offset),
                                                 z)});
                          if (batch.size() >= kEmitBatch)
                            flush(bucket, batch);
                        });
                    if (!ok) {
                        bad_position.store(true, std::memory_order_relaxed);
                        invalid_pos.store(true, std::memory_order_relaxed);
                        break;
                    }
                }

                for (uint32_t bucket = 0; bucket < dispatch_count; ++bucket) {
                    if (!batches[bucket].empty()) flush(bucket, batches[bucket]);
                }
            });
        if (invalid_pos.load(std::memory_order_relaxed)) invalid = true;
    }

    if (invalid || seen_count != chrom_count) {
      // Either the emitter reported a bad seed position, or the streamed source
      // did not deliver every declared contig exactly once.
      return faix_build_fail(status,
                             bad_position.load(std::memory_order_relaxed)
                                 ? FaixBuildError::SeedPositionOverflow
                                 : FaixBuildError::IncompleteStream);
    }

    return htable_finish_faix_index_from_scratch(
        scratch, cfg, chr_offsets, chrom_count, build_threads, status);
}

}}  // namespace fa::cpu

/**** index ****/
// Seed-index facade: the header DNA, RNA and seeding code include. The runtime (FaixIndex)
// and the builder are above, the file format in format.h.

namespace fa { namespace cpu {

using SeedIndex = FaixIndex;
using IndexView = SeedIndex;

// `image_offset` / `image_bytes` name one part of a multi-part container
// (FaixPartEntry); both 0 is the whole file, i.e. a single-part .faix.
inline SeedIndex load_seed_index(
    const std::string& path, FaixLoadStatus* status = nullptr,
    int n_threads = 0, uint64_t image_offset = 0, uint64_t image_bytes = 0) {
    return SeedIndex::load_image(path, image_offset, image_bytes, status,
                                 n_threads);
}

// Builds an index in memory. `status`, when supplied, names the refusal behind an empty
// result (see FaixBuildError).
inline SeedIndex build_seed_index(
    const std::vector<std::vector<uint8_t>>& chr_encs,
    const FaixBuildConfig& cfg,
    FaixBuildStatus* status = nullptr
) {
    return build_faix_index(chr_encs, cfg, status);
}

template <class EmitEncodedSequences>
inline SeedIndex build_seed_index_streamed(
    const std::vector<uint64_t>& chr_offsets,
    const FaixBuildConfig& cfg,
    EmitEncodedSequences&& emit_encoded_sequences,
    FaixBuildStatus* status = nullptr
) {
    return build_faix_index_streamed(
        chr_offsets,
        cfg,
        std::forward<EmitEncodedSequences>(emit_encoded_sequences), status);
}

}}  // namespace fa::cpu

/**** reference_context ****/
// Non-owning view of a loaded reference. The aligner owns every pointed-to object; a
// ReferenceContext is valid only while the aligner lives and must not escape a mapping call.

namespace fa {
namespace cpu {
namespace engine {

// A view of a contig's bases [lo, hi), addressed in contig coordinates.
struct ReferenceWindow {
    const uint8_t* origin = nullptr; // the base at lo
    int64_t lo = 0;
    int64_t hi = 0;

    const uint8_t* ptr(int64_t pos) const noexcept {
        assert(pos >= lo && pos <= hi);
        return origin + (pos - lo);
    }
    uint8_t at(int64_t pos) const noexcept {
        assert(pos >= lo && pos < hi);
        return origin[pos - lo];
    }
};

struct ReferenceContext {
    const SeedIndex* index = nullptr;
    const std::vector<std::string>* names = nullptr;
    // One length per contig, in contig-id order, aligned with `names`. Read from the
    // index's offsets table, so it is present whether or not the bases are.
    const std::vector<int64_t>* lengths = nullptr;
    // The unpacked reference, one byte per base. Empty when the mode reads no bases (DNA
    // map-only) or reads them packed; check has_bases() first, and take contig counts and
    // lengths from `lengths`.
    const std::vector<std::vector<uint8_t>>* encoded = nullptr;
    // The DNA presets: `encoded` is empty and the bases stay packed in `index`, read
    // through window().
    bool packed = false;

    bool valid() const noexcept {
        return index != nullptr && names != nullptr && encoded != nullptr;
    }
    // The reference bases can be read: packed (DNA with a CIGAR or a FASTA target) or
    // materialized (RNA).
    bool has_bases() const noexcept {
        return packed || (encoded != nullptr && !encoded->empty());
    }
    // The bases of `contig` over [lo, hi), clamped to the contig. Unpacked, the view is the
    // whole contig (lo 0, hi its length) and nothing is copied; packed, [lo, hi) is unpacked
    // into `buffer`, as minimap2's mm_idx_getseq. Needs has_bases().
    ReferenceWindow window(int contig, int64_t lo, int64_t hi,
                           std::vector<uint8_t>& buffer) const {
        if (!packed) {
            const std::vector<uint8_t>& bases =
                (*encoded)[static_cast<size_t>(contig)];
            return {bases.data(), 0, static_cast<int64_t>(bases.size())};
        }
        const int64_t length = contig_length(contig);
        lo = std::min(std::max<int64_t>(lo, 0), length);
        hi = std::min(std::max(hi, lo), length);
        // Never empty, so the origin is never null.
        buffer.resize(static_cast<size_t>(std::max<int64_t>(hi - lo, 1)));
        index->reference_bases_u8(static_cast<uint64_t>(contig),
                                  static_cast<uint64_t>(lo),
                                  static_cast<uint64_t>(hi), buffer.data());
        return {buffer.data(), lo, hi};
    }
    int contig_count() const noexcept {
        return lengths == nullptr ? 0 : static_cast<int>(lengths->size());
    }
    // 0 for an id outside the table (a zero-length contig admits no seed span).
    int64_t contig_length(int contig) const noexcept {
        if (lengths == nullptr || contig < 0 || contig >= contig_count())
            return 0;
        return (*lengths)[static_cast<size_t>(contig)];
    }
};

} // namespace engine
} // namespace cpu
} // namespace fa
