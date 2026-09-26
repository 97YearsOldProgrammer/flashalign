// In-memory index builder: emit reference seeds into shards across threads, sort and pack
// each shard, and assemble a FaixIndex.
#pragma once

#include "format.h"
#include "seed.h"                        // htable_valid_k, closed-syncmer reference extraction
#include "faix.h"                        // FaixIndex runtime
#include "../core/hash.h"               // detail::splitmix64
#include "../threading/parallel_for.h"  // default_thread_count

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

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

// Turns the packed scratch into the final compact index, consuming `scratch`.
// `chrom_count` is the number of reference sequences.
FaixIndex assemble_faix_index(
    std::vector<KmerBuildScratchBucket>& scratch,
    const FaixBuildConfig& cfg,
    const std::vector<uint64_t>& chr_offsets,
    uint64_t chrom_count,
    int shard_bits,
    int build_threads,
    FaixBuildStatus* status = nullptr
);

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
