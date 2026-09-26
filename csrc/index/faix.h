// FaixIndex, the query-time seed index: load, save and lookup over the format in format.h.
// build_faix_index() in build.h produces one through from_compact_owned().
#pragma once

#include "format.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>


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

    bool loaded() const noexcept { return !parts.empty(); }
    // Bases in part `i`: the global offsets bracket its contig run.
    std::uint64_t part_bp(std::size_t i) const noexcept {
        const FaixPartEntry& e = parts[i];
        return chr_offsets[static_cast<std::size_t>(e.first_contig +
                                                    e.contig_count)] -
               chr_offsets[static_cast<std::size_t>(e.first_contig)];
    }
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

    bool empty() const { return header_.slot_count == 0; }

    int k() const { return static_cast<int>(header_.k); }
    int build_syncmer_s() const { return static_cast<int>(header_.syncmer_s); }
    // The preset this index was built with, or "" when it records none.
    std::string build_preset() const { return faix_header_preset(header_); }
    uint64_t total_bp() const { return header_.total_bp; }
    uint64_t distinct_keys() const { return header_.key_count; }
    uint64_t chrom_count() const { return header_.contig_count; }
    uint64_t memory_bytes() const {
        if (empty() || !shard_dir_) return 0;
        const uint64_t n_shards = 1ULL << header_.shard_bits;
        uint64_t occ_bytes = 0, rank_entries = 0;
        for (uint64_t s = 0; s < n_shards; ++s) {
            occ_bytes += occupancy_storage_bytes(shard_dir_[s].bucket_count);
            rank_entries += rank_index_count(shard_dir_[s].bucket_count,
                                             header_.rank_block_bits);
        }
        const uint64_t chrom_bytes =
            (header_.contig_count + 1) * sizeof(uint64_t);
        const uint64_t shard_bytes = n_shards * sizeof(FaixShardEntry);
        const uint64_t rank_bytes = rank_entries * sizeof(uint64_t);
        const uint64_t bucket_bytes =
            header_.key_count * (sizeof(uint32_t) + sizeof(uint64_t));
        const uint64_t posting_bytes =
            header_.multi_posting_count * sizeof(PackedRefPos);
        const uint64_t ref_bytes =
            has_reference_payload() ? packed_reference_bytes(header_.total_bp) : 0;
        return chrom_bytes + shard_bytes + occ_bytes + rank_bytes + bucket_bytes
            + posting_bytes + header_.name_bytes + ref_bytes;
    }

    const uint64_t* chrom_offsets_data() const { return chr_offsets_; }
    // The record table (the shards' slices end to end), for enumerating stored keys.
    const uint32_t* tags_data() const { return tags_; }
    const uint64_t* payloads_data() const { return payloads_; }
    const std::vector<std::string>& chromosome_names() const { return chrom_names_owned_; }

    bool has_reference_payload() const {
        return (header_.flags & kFaixFlagHasReference) != 0 &&
               packed_reference_ != nullptr;
    }

    // Materializes the embedded reference as one u8-per-base vector per contig, threaded.
    // `n_threads` <= 0 means threading::default_thread_count(). On failure `out` is empty.
    bool reference_chromosomes_u8(
        std::vector<std::vector<uint8_t>>& out,
        int n_threads = 0
    ) const;

    // Attaches the contig names and the reference, one base code per byte and contig, in
    // contig order.
    bool set_reference_payload(
        std::vector<std::string> chrom_names,
        const std::vector<std::vector<uint8_t>>& contigs
    ) {
        if (contigs.size() != static_cast<size_t>(header_.contig_count) ||
            !chr_offsets_) {
            return false;
        }
        std::vector<uint8_t> packed(
            static_cast<size_t>(packed_reference_bytes(header_.total_bp)), 0);
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (chr_offsets_[i + 1] - chr_offsets_[i] !=
                static_cast<uint64_t>(contigs[i].size())) {
                return false;
            }
            pack_reference_4bit_into(contigs[i].data(),
                                     static_cast<uint64_t>(contigs[i].size()),
                                     chr_offsets_[i], packed);
        }
        return set_packed_reference_payload(std::move(chrom_names), std::move(packed));
    }

    // Names only (`flashalign index --idx-no-seq`): the index still names and measures every
    // contig, which map-only output needs, but has_reference_payload() is false, so CIGAR
    // output is refused. Any attached reference is dropped.
    bool set_chromosome_names(std::vector<std::string> chrom_names) {
        if (!set_names(std::move(chrom_names)))
            return false;
        packed_reference_owned_.clear();
        packed_reference_owned_.shrink_to_fit();
        packed_reference_ = nullptr;
        header_.flags &= ~kFaixFlagHasReference;
        return true;
    }

    // Attaches the contig names and the reference packed as in the file.
    bool set_packed_reference_payload(
        std::vector<std::string> chrom_names,
        std::vector<uint8_t>&& packed_reference
    ) {
        if (!chr_offsets_ ||
            chr_offsets_[static_cast<size_t>(header_.contig_count)] != header_.total_bp ||
            packed_reference.size() !=
                static_cast<size_t>(packed_reference_bytes(header_.total_bp))) {
            return false;
        }
        if (!set_names(std::move(chrom_names)))
            return false;
        packed_reference_owned_ = std::move(packed_reference);
        packed_reference_ = packed_reference_owned_.empty()
            ? nullptr
            : packed_reference_owned_.data();
        header_.flags |= kFaixFlagHasReference;
        return true;
    }

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
    ) const {
        if (n == 0 || !keys || !out) return;
        if (!tags_ || !shard_dir_ || header_.slot_count == 0) {
            for (size_t i = 0; i < n; ++i) out[i] = KmerPostingView{};
            return;
        }
        faix_sharded_lookup_batch(
            tags_, payloads_, postings_, chr_offsets_, header_.contig_count,
            occupancy_bits_, rank_index_, shard_dir_, shard_mask_,
            static_cast<int>(header_.shard_bits), 2 * static_cast<int>(header_.k),
            header_.rank_block_bits, keys, n, out);
    }

    // The narrowing half of lookup_interval, for a caller that already holds a key's full
    // posting view. The span aliases the index's postings and lives as long as the index.
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
    ) const {
        // Checked before the probe, so a degenerate interval costs no lookup.
        if (!chr_offsets_ ||
            chr_idx < 0 ||
            chr_idx >= static_cast<int>(header_.contig_count) ||
            chr_hi <= chr_lo) {
            return {};
        }
        return narrow_interval(lookup(key), chr_idx, chr_lo, chr_hi);
    }

    // Writes this index as one complete image to `out` and reports the byte count the
    // layout implies, computed before writing so a caller can check what it wrote. Refuses
    // an index without its contig names.
    bool write_image(FaixByteSink& out, uint64_t& expected_bytes) const;

    bool save(const std::string& path) const;

    // `n_threads` bounds the loader's parallel reads; <= 0 takes the default.
    static FaixIndex load(
        const std::string& path, FaixLoadStatus* status = nullptr,
        int n_threads = 0) {
        return load_image(path, 0, 0, status, n_threads);
    }

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
                       std::vector<FaixShardEntry> shard_dir) {
      FaixIndex idx;
      idx.header_ = header;
      idx.chr_offsets_owned_ = std::move(chr_offsets);
      idx.occupancy_bits_owned_ = std::move(occupancy_bits);
      idx.rank_index_owned_ = std::move(rank_index);
      idx.tags_owned_ = std::move(compact_tags);
      idx.payloads_owned_ = std::move(compact_payloads);
      idx.postings_owned_ = std::move(compact_postings);
      idx.shard_dir_owned_ = std::move(shard_dir);
      idx.shard_mask_ = (1ULL << idx.header_.shard_bits) - 1ULL;
      idx.rebind_views_to_owned();
      return idx;
    }

private:
    void rebind_views_to_owned() {
        chr_offsets_ = chr_offsets_owned_.empty() ? nullptr : chr_offsets_owned_.data();
        tags_ = tags_owned_.empty() ? nullptr : tags_owned_.data();
        payloads_ = payloads_owned_.empty() ? nullptr : payloads_owned_.data();
        postings_ = postings_owned_.empty() ? nullptr : postings_owned_.data();
        occupancy_bits_ = occupancy_bits_owned_.empty()
            ? nullptr
            : occupancy_bits_owned_.data();
        rank_index_ = rank_index_owned_.empty()
            ? nullptr
            : rank_index_owned_.data();
        shard_dir_ = shard_dir_owned_.empty() ? nullptr : shard_dir_owned_.data();
        packed_reference_ = packed_reference_owned_.empty()
            ? nullptr
            : packed_reference_owned_.data();
    }

    // Attaches one name per contig and records the name table's size.
    bool set_names(std::vector<std::string> chrom_names) {
        if (chrom_names.size() != static_cast<size_t>(header_.contig_count))
            return false;
        uint64_t name_bytes = 0;
        for (const auto& name : chrom_names) {
            if (name.size() >
                static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
                return false;
            }
            name_bytes += sizeof(uint32_t) + static_cast<uint64_t>(name.size());
        }
        chrom_names_owned_ = std::move(chrom_names);
        header_.name_bytes = name_bytes;
        return true;
    }

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
