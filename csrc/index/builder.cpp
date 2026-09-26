#include "build.h"

namespace fa { namespace cpu {

namespace {

uint64_t next_power_of_two_u64(uint64_t x) {
    if (x <= 1) return 1;
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    return x + 1;
}

}  // namespace

std::string faix_build_error_message(const FaixBuildStatus& status) {
    const std::string value = std::to_string(status.value);
    switch (status.error) {
      case FaixBuildError::None:
        return "the reference has no seed k-mers";
      case FaixBuildError::UnsupportedK:
        return "k=" + value + " is not supported; k must be between 1 and " +
               std::to_string(kFaixMaxK);
      case FaixBuildError::NoSequences:
        return "the reference has no sequences";
      case FaixBuildError::TooManySequences:
        return "the reference has " + value + " sequences; the most an index holds is " +
               std::to_string(kFaixMaxContigCount);
      case FaixBuildError::ContigTooLong:
        return "reference sequence #" + std::to_string(status.sequence_index + 1) +
               " is " + value + " bp, longer than the " +
               std::to_string(kFaixMaxContigBp) + " bp SAM and BAM can address";
      case FaixBuildError::InvalidChrOffsets:
        return "the reference offsets do not start at 0 and increase";
      case FaixBuildError::SeedPositionOverflow:
        return "a reference seed position is outside the index range";
      case FaixBuildError::IncompleteStream:
        return "the reference file changed while it was being indexed";
      case FaixBuildError::ReferenceTooLarge:
        return "the reference is too large for one index; split it with "
               "`flashalign index -I`";
      case FaixBuildError::Internal:
        return "internal error while packing the index";
      case FaixBuildError::PresetNameTooLong:
        return "the preset name is longer than " +
               std::to_string(kFaixPresetBytes - 1) + " characters";
    }
    return "unknown build error";
}

int htable_effective_closed_syncmer_s(const FaixBuildConfig& cfg) {
    return effective_closed_syncmer_s(cfg.k, cfg.syncmer_s);
}

bool htable_position_record_less(
    const KmerPositionRecord& a,
    const KmerPositionRecord& b
) {
    if (a.key != b.key) return a.key < b.key;
    return a.pos < b.pos;
}

void htable_sort_position_records(
    std::vector<KmerPositionRecord>& records
) {
    if (records.size() < 4096) {
        std::sort(records.begin(), records.end(), htable_position_record_less);
        return;
    }

    constexpr size_t kRadixBits = 16;
    constexpr size_t kRadix = 1ULL << kRadixBits;
    constexpr uint64_t kRadixMask = kRadix - 1ULL;

    auto free_records = [](KmerPositionRecord* p) { std::free(p); };
    std::unique_ptr<KmerPositionRecord, decltype(free_records)> tmp_owner(
        static_cast<KmerPositionRecord*>(
            std::malloc(records.size() * sizeof(KmerPositionRecord))),
        free_records);
    if (!tmp_owner) {
        std::sort(records.begin(), records.end(), htable_position_record_less);
        return;
    }

    std::array<size_t, kRadix> count{};
    KmerPositionRecord* in = records.data();
    KmerPositionRecord* out = tmp_owner.get();

    auto radix_pass = [&](bool by_pos, unsigned shift) {
        count.fill(0);
        if (by_pos) {
            for (const KmerPositionRecord* p = in; p != in + records.size(); ++p) {
                ++count[static_cast<size_t>((p->pos >> shift) & kRadixMask)];
            }
            size_t running = 0;
            for (size_t i = 0; i < kRadix; ++i) {
                const size_t n = count[i];
                count[i] = running;
                running += n;
            }
            for (const KmerPositionRecord* p = in; p != in + records.size(); ++p) {
                const size_t digit =
                    static_cast<size_t>((p->pos >> shift) & kRadixMask);
                out[count[digit]++] = *p;
            }
        } else {
            for (const KmerPositionRecord* p = in; p != in + records.size(); ++p) {
                ++count[static_cast<size_t>((p->key >> shift) & kRadixMask)];
            }
            size_t running = 0;
            for (size_t i = 0; i < kRadix; ++i) {
                const size_t n = count[i];
                count[i] = running;
                running += n;
            }
            for (const KmerPositionRecord* p = in; p != in + records.size(); ++p) {
                const size_t digit =
                    static_cast<size_t>((p->key >> shift) & kRadixMask);
                out[count[digit]++] = *p;
            }
        }
        std::swap(in, out);
    };

    // Stable LSD radix sort by (key, pos), the same order as the comparator: interval
    // lookups need every key's postings sorted by packed (contig, local, z).
    radix_pass(/*by_pos=*/true, 0);
    radix_pass(/*by_pos=*/true, 16);
    radix_pass(/*by_pos=*/true, 32);
    radix_pass(/*by_pos=*/true, 48);
    radix_pass(/*by_pos=*/false, 0);
    radix_pass(/*by_pos=*/false, 16);
    radix_pass(/*by_pos=*/false, 32);
    radix_pass(/*by_pos=*/false, 48);

    if (in != records.data()) {
        std::copy(in, in + records.size(), records.data());
    }
}

int htable_build_threads(int requested_threads) {
    if (requested_threads > 0) {
        return requested_threads;
    }
    return threading::default_thread_count();
}

int htable_dispatch_bucket_bits(const FaixBuildConfig& cfg) {
    // The tag word packs 2k - shard_bits key bits and the posting bit into a uint32_t, so
    // shard_bits must be at least 2k - 31; capped at 16, which covers k <= 23.
    const int need = 2 * cfg.k - 31;
    return std::min(16, std::max(kFaixMinShardBits, need));
}

uint32_t htable_dispatch_bucket(uint64_t key, uint32_t mask, int twok) {
    return static_cast<uint32_t>(detail::shard_mix_2k(key, twok) & mask);
}

double htable_reference_seed_density_estimate(const FaixBuildConfig& cfg) {
    // Expected closed-syncmer density: about 2 / (k - s + 1).
    const int s = htable_effective_closed_syncmer_s(cfg);
    if (s <= 0 || s > cfg.k) return 0.0;
    return std::min(1.0, 2.0 / static_cast<double>(cfg.k - s + 1));
}

FaixIndex assemble_faix_index(
    std::vector<KmerBuildScratchBucket>& scratch,
    const FaixBuildConfig& cfg,
    const std::vector<uint64_t>& chr_offsets,
    uint64_t chrom_count,
    int shard_bits,
    int build_threads,
    FaixBuildStatus* status
) {
    const size_t n_shards = scratch.size();  // == 1u << shard_bits
    const uint32_t rank_block_bits = kFaixRankBlockBits;
    const int twok = 2 * cfg.k;

    // Pass A (parallel): per-shard sizes from the grouped scratch. Each shard is its own
    // compact directory.
    std::vector<FaixShardEntry> shard_dir(n_shards);
    std::vector<uint64_t> occ_bytes(n_shards, 0);
    std::vector<uint64_t> rank_entries(n_shards, 0);
    threading::parallel_for(
        build_threads, static_cast<int64_t>(n_shards),
        [&](int64_t s, int /*tid*/) {
            const auto& groups = scratch[static_cast<size_t>(s)].groups;
            const uint64_t n_keys = static_cast<uint64_t>(groups.size());
            uint64_t bc = 0;
            if (n_keys > 0) {
                const double need_d =
                    static_cast<double>(n_keys) / kFaixLoadFactor;
                const uint64_t need = static_cast<uint64_t>(need_d + 0.999999);
                bc = next_power_of_two_u64(std::max<uint64_t>(2, need));
            }
            uint64_t pcount = 0;
            for (const auto& g : groups) {
                if (g.count() > 1)
                  pcount += g.count();
            }
            FaixShardEntry e;
            e.bucket_count = bc;
            e.unique_keys = n_keys;
            e.posting_count = pcount;
            shard_dir[static_cast<size_t>(s)] = e;
            occ_bytes[static_cast<size_t>(s)] = occupancy_storage_bytes(bc);
            rank_entries[static_cast<size_t>(s)] =
                rank_index_count(bc, rank_block_bits);
        });

    // Cumulative offsets and totals.
    uint64_t occ_total = 0, rank_total = 0, bucket_total = 0, post_total = 0;
    uint64_t slot_total = 0;
    bool shard_payload_overflow = false;
    for (size_t s = 0; s < n_shards; ++s) {
        FaixShardEntry& e = shard_dir[s];
        if (e.unique_keys > std::numeric_limits<uint32_t>::max() ||
            e.posting_count > std::numeric_limits<uint32_t>::max()) {
          shard_payload_overflow = true;
        }
        e.occ_off = occ_total;
        e.rank_off = rank_total;
        e.bucket_off = bucket_total;
        e.post_off = post_total;
        occ_total += occ_bytes[s];
        rank_total += rank_entries[s];
        bucket_total += e.unique_keys;
        post_total += e.posting_count;
        slot_total += e.bucket_count;
    }
    if (shard_payload_overflow)
      return faix_build_fail(status, FaixBuildError::ReferenceTooLarge);

    FaixHeader header;
    header.k = static_cast<uint32_t>(cfg.k);
    header.syncmer_s = static_cast<uint32_t>(htable_effective_closed_syncmer_s(cfg));
    header.shard_bits = static_cast<uint32_t>(shard_bits);
    header.rank_block_bits = rank_block_bits;
    header.contig_count = chrom_count;
    header.total_bp = chr_offsets.empty() ? 0 : chr_offsets.back();
    header.key_count = bucket_total;
    header.slot_count = slot_total;
    header.multi_posting_count = post_total;
    // Recorded so `flashalign align` can map under the preset the index was built for.
    if (!faix_header_set_preset(header, cfg.preset))
      return faix_build_fail(status, FaixBuildError::PresetNameTooLong);

    // The four concatenated blobs, allocated once and owned by the index.
    std::vector<uint8_t> occupancy(static_cast<size_t>(occ_total), 0);
    std::vector<uint64_t> rank(static_cast<size_t>(rank_total), 0);
    std::vector<uint32_t> tags(static_cast<size_t>(bucket_total));      // (tag<<1)|posting_bit
    std::vector<uint64_t> payloads(static_cast<size_t>(bucket_total));
    std::vector<PackedRefPos> postings(static_cast<size_t>(post_total));

    // Pass B (parallel): place, rank and scatter each shard into its slice, releasing the
    // shard's scratch as soon as it is read.
    std::atomic<bool> failed{false};
    threading::parallel_for(
        build_threads, static_cast<int64_t>(n_shards),
        [&](int64_t s, int /*tid*/) {
            auto& src = scratch[static_cast<size_t>(s)];
            const FaixShardEntry& e = shard_dir[static_cast<size_t>(s)];
            const uint64_t n_keys = e.unique_keys;
            if (n_keys == 0) {
                std::vector<KmerBucket>().swap(src.groups);
                std::vector<PackedRefPos>().swap(src.postings);
                return;
            }
            const uint64_t bc = e.bucket_count;
            const uint64_t mask = bc - 1;
            uint8_t* occ = occupancy.data() + e.occ_off;
            uint64_t* rk = rank.data() + e.rank_off;
            uint32_t* tg = tags.data() + e.bucket_off;
            uint64_t* pl = payloads.data() + e.bucket_off;
            PackedRefPos* pk =
                post_total ? postings.data() + e.post_off : nullptr;

            // place: probe by the key tag (= mixed_key >> shard_bits) so the slot
            // order matches the read path; record per-group slot + tag.
            std::vector<uint64_t> group_slot(static_cast<size_t>(n_keys));
            std::vector<uint32_t> group_tag(static_cast<size_t>(n_keys));
            for (size_t gi = 0; gi < n_keys; ++gi) {
              const uint64_t mixed =
                  detail::shard_mix_2k(src.groups[gi].key_value(), twok);
              const uint32_t tag = static_cast<uint32_t>(mixed >> shard_bits);
              group_tag[gi] = tag;
              uint64_t slot = detail::splitmix64(tag) & mask;
              while (occupancy_bit(occ, slot))
                slot = (slot + 1) & mask;
              occ[static_cast<size_t>(slot >> 3)] |=
                  static_cast<uint8_t>(1u << (slot & 7ULL));
              group_slot[gi] = slot;
            }

            // rank index straight from the shard's occupancy bitmap.
            const uint64_t rank_n = rank_index_count(bc, rank_block_bits);
            {
                const uint64_t block_size = 1ULL << rank_block_bits;
                uint64_t running = 0;
                for (uint64_t b = 0; b + 1 < rank_n; ++b) {
                    rk[b] = running;
                    const uint64_t lo = b * block_size;
                    const uint64_t hi = std::min(bc, (b + 1) * block_size);
                    running += popcount_bit_range(occ, lo, hi);
                }
                if (rank_n) rk[rank_n - 1] = running;
                if (running != n_keys) {
                    failed.store(true, std::memory_order_relaxed);
                    return;
                }
            }

            FaixDirectory dir;
            dir.occupancy_bits = occ;
            dir.rank_index = rk;
            dir.rank_entries = rank_n;
            dir.rank_block_bits = rank_block_bits;
            dir.unique_key_count = n_keys;

            // compact position per group + prefix-summed shard-LOCAL post offsets.
            std::vector<uint32_t> group_pos(static_cast<size_t>(n_keys));
            std::vector<uint32_t> posting_begin(static_cast<size_t>(n_keys) + 1, 0);
            for (size_t gi = 0; gi < n_keys; ++gi) {
                const uint64_t pos = directory_rank_before(dir, group_slot[gi]);
                group_pos[gi] = static_cast<uint32_t>(pos);
                const uint32_t count = src.groups[gi].count();
                posting_begin[static_cast<size_t>(pos) + 1] =
                    (count > 1) ? count : 0u;
            }
            std::vector<uint64_t>().swap(group_slot);
            for (size_t i = 1; i < posting_begin.size(); ++i) {
                posting_begin[i] += posting_begin[i - 1];
            }

            // scatter into the shard's tag/payload/posting slice. The tag word
            // is (tag << 1) | posting_bit: bit 0 = 1 for multi, 0 for singleton.
            for (size_t gi = 0; gi < n_keys; ++gi) {
                const KmerBucket& group = src.groups[gi];
                const size_t pos = static_cast<size_t>(group_pos[gi]);
                const uint32_t count = group.count();
                uint64_t payload = 0;
                uint32_t posting_bit = 0;
                if (count <= 1) {
                    // singleton: the payload is the packed (contig, local, z) word.
                    const uint64_t begin = group.begin();
                    payload =
                        group.inlined()
                            ? begin
                            : (begin < src.postings.size()
                                   ? src.postings[static_cast<size_t>(begin)]
                                   : 0u);
                    if (!group.inlined() && begin >= src.postings.size()) {
                        failed.store(true, std::memory_order_relaxed);
                    }
                } else {
                  // Multi: payload owns count and shard-local begin; the
                  // posting array contains positions only.
                  const uint64_t begin = group.begin();
                  const uint64_t out_begin = posting_begin[pos];
                  if (begin > src.postings.size() ||
                      static_cast<uint64_t>(count) >
                          src.postings.size() - begin) {
                    failed.store(true, std::memory_order_relaxed);
                  } else {
                    posting_bit = kFaixTagMultiBit;
                    if (!faix_pack_multi_payload(out_begin, count, payload)) {
                      failed.store(true, std::memory_order_relaxed);
                      continue;
                    }
                    std::copy(src.postings.begin() +
                                  static_cast<std::ptrdiff_t>(begin),
                              src.postings.begin() +
                                  static_cast<std::ptrdiff_t>(begin + count),
                              pk + out_begin);
                  }
                }
                tg[pos] = (group_tag[gi] << 1) | posting_bit;
                pl[pos] = payload;
            }

            std::vector<KmerBucket>().swap(src.groups);
            std::vector<PackedRefPos>().swap(src.postings);
        });

    if (failed.load(std::memory_order_relaxed))
      return faix_build_fail(status, FaixBuildError::Internal);

    return FaixIndex::from_compact_owned(
        header,
        chr_offsets,
        std::move(occupancy),
        std::move(rank),
        std::move(tags),
        std::move(payloads),
        std::move(postings),
        std::move(shard_dir));
}

FaixIndex htable_finish_faix_index_from_scratch(
    std::vector<KmerBuildScratchBucket>& scratch,
    const FaixBuildConfig& cfg,
    const std::vector<uint64_t>& chr_offsets,
    uint64_t chrom_count,
    int build_threads,
    FaixBuildStatus* status
) {
    std::atomic<bool> packed_bucket_overflow{false};

    // Sort + pack each bucket independently; one worker per bucket.
    threading::parallel_for(
        build_threads, static_cast<int64_t>(scratch.size()),
        [&](int64_t bucket, int /*tid*/) {
        auto& dst = scratch[static_cast<size_t>(bucket)];
        if (dst.records.empty()) return;

        htable_sort_position_records(dst.records);

        size_t group_count = 0;
        uint64_t multi_posting_count = 0;
        for (size_t i = 0; i < dst.records.size();) {
            size_t j = i + 1;
            while (j < dst.records.size() && dst.records[j].key == dst.records[i].key) {
                ++j;
            }
            const uint64_t total = static_cast<uint64_t>(j - i);
            ++group_count;
            if (total > 1) multi_posting_count += total;
            i = j;
        }
        dst.groups.reserve(group_count);
        dst.postings.reserve(static_cast<size_t>(multi_posting_count));

        for (size_t i = 0; i < dst.records.size();) {
            size_t j = i + 1;
            while (j < dst.records.size() && dst.records[j].key == dst.records[i].key) {
                ++j;
            }

            const uint64_t total = static_cast<uint64_t>(j - i);

            // Every key is kept in full. A singleton is inlined with its position in begin; a
            // multi-hit key's begin indexes the bucket postings.
            KmerBucket group;
            group.key = dst.records[i].key;
            const uint64_t begin = (total == 1)
                ? static_cast<uint64_t>(dst.records[i].pos)
                : static_cast<uint64_t>(dst.postings.size());
            if (!group.set_payload(begin, total, total == 1)) {
                packed_bucket_overflow.store(true, std::memory_order_relaxed);
                i = j;
                continue;
            }
            dst.groups.push_back(group);

            if (total > 1) {
                for (uint64_t t = 0; t < total; ++t) {
                    dst.postings.push_back(dst.records[i + static_cast<size_t>(t)].pos);
                }
            }
            i = j;
        }

        std::vector<KmerPositionRecord>().swap(dst.records);
        });

    if (packed_bucket_overflow.load(std::memory_order_relaxed))
      return faix_build_fail(status, FaixBuildError::ReferenceTooLarge);

    // scratch.size() == 1u << shard_bits: the build buckets are the index shards.
    const int shard_bits = htable_dispatch_bucket_bits(cfg);
    return assemble_faix_index(scratch, cfg, chr_offsets, chrom_count,
                               shard_bits, build_threads, status);
}

}}  // namespace fa::cpu
