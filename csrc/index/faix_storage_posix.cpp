// POSIX storage for .faix files: save, load, probe and the container table.

#include "faix.h"
#include "faix_publish.h"
#include "threading/parallel_for.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

namespace fa { namespace cpu {

namespace {

bool checked_add(uint64_t& total, uint64_t value) {
  if (value > std::numeric_limits<uint64_t>::max() - total)
    return false;
  total += value;
  return true;
}

bool checked_product(uint64_t left, uint64_t right, uint64_t& out) {
  if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left)
    return false;
  out = left * right;
  return true;
}

// Reads the index image into the anonymous mapping load() hands out, with large parallel
// preads. A file mapping would page in the table one fault at a time on first use, and its
// clean pages can be reclaimed under memory pressure and faulted back mid-run; anonymous
// memory is neither (minimap2 likewise reads its .mmi into malloc'd memory). The caller
// then makes the image read-only and drops the page-cache copy. The reader count is bounded
// by `n_threads` when > 0, else the default thread count. The image is the `image_size`
// bytes at `file_offset`: the whole file, or one part of a multi-part container.
bool read_file_image(int fd, void* base, uint64_t file_offset,
                     uint64_t image_size, int n_threads) {
  constexpr uint64_t kSlice = 64ull << 20;
  constexpr size_t kRead = 8u << 20;
  if (image_size == 0) return true;
  const int64_t slices =
      static_cast<int64_t>((image_size + kSlice - 1) / kSlice);
  const int threads = std::max(
      1, std::min(8, n_threads > 0
                         ? n_threads
                         : ::fa::cpu::threading::default_thread_count()));
  std::atomic<bool> ok{true};
  char* const out = static_cast<char*>(base);
  ::fa::cpu::threading::parallel_for(
      threads, slices, [&](int64_t slice, int) {
        uint64_t off = static_cast<uint64_t>(slice) * kSlice;
        const uint64_t end = std::min(off + kSlice, image_size);
        while (off < end) {
          const size_t want =
              static_cast<size_t>(std::min<uint64_t>(kRead, end - off));
          const ssize_t n = ::pread(fd, out + off, want,
                                    static_cast<off_t>(file_offset + off));
          if (n <= 0) {
            ok.store(false, std::memory_order_relaxed);
            return;
          }
          off += static_cast<uint64_t>(n);
        }
      });
  return ok.load(std::memory_order_relaxed);
}

// Exactly `size` bytes at `offset`, or false.
bool pread_exact(int fd, void* out, uint64_t size, uint64_t offset) {
  char* cursor = static_cast<char*>(out);
  while (size > 0) {
    const ssize_t n = ::pread(fd, cursor, static_cast<size_t>(size),
                              static_cast<off_t>(offset));
    if (n < 0 && errno == EINTR)
      continue;
    if (n <= 0)
      return false;
    cursor += n;
    offset += static_cast<uint64_t>(n);
    size -= static_cast<uint64_t>(n);
  }
  return true;
}

struct FdCloser {
  int fd;
  ~FdCloser() {
    if (fd >= 0)
      close(fd);
  }
};

// Reads and decodes the header of the image of `bytes` bytes at `offset`.
FaixLoadError read_header(int fd, uint64_t offset, uint64_t bytes, FaixHeader& h) {
  char head[kFaixHeaderBytes] = {};
  const uint64_t got = std::min<uint64_t>(bytes, kFaixHeaderBytes);
  if (!pread_exact(fd, head, got, offset))
    return FaixLoadError::ReadFailed;
  if (got < 8)
    return FaixLoadError::NotIndex;
  if (got < kFaixHeaderBytes) {
    if (std::memcmp(head, kFaixDevelopmentMagic, 8) == 0)
      return FaixLoadError::DevelopmentBuild;
    return std::memcmp(head, kFaixMagic, 8) == 0 ? FaixLoadError::Truncated
                                                 : FaixLoadError::NotIndex;
  }
  return faix_header_decode(head, h);
}

} // namespace

FaixIndex::~FaixIndex() { release(); }

FaixIndex::FaixIndex(FaixIndex&& other) noexcept { swap(other); }

FaixIndex& FaixIndex::operator=(FaixIndex&& other) noexcept {
        if (this != &other) {
            release();
            swap(other);
        }
        return *this;
    }

void FaixIndex::swap(FaixIndex& other) noexcept {
        using std::swap;
        swap(header_, other.header_);
        swap(chr_offsets_owned_, other.chr_offsets_owned_);
        swap(tags_owned_, other.tags_owned_);
        swap(payloads_owned_, other.payloads_owned_);
        swap(postings_owned_, other.postings_owned_);
        swap(chrom_names_owned_, other.chrom_names_owned_);
        swap(occupancy_bits_owned_, other.occupancy_bits_owned_);
        swap(rank_index_owned_, other.rank_index_owned_);
        swap(shard_dir_owned_, other.shard_dir_owned_);
        swap(packed_reference_owned_, other.packed_reference_owned_);
        swap(chr_offsets_, other.chr_offsets_);
        swap(tags_, other.tags_);
        swap(payloads_, other.payloads_);
        swap(postings_, other.postings_);
        swap(occupancy_bits_, other.occupancy_bits_);
        swap(rank_index_, other.rank_index_);
        swap(shard_dir_, other.shard_dir_);
        swap(shard_mask_, other.shard_mask_);
        swap(packed_reference_, other.packed_reference_);
        swap(mapped_base_, other.mapped_base_);
        swap(mapped_size_, other.mapped_size_);
    }

void FaixIndex::release() {
        if (mapped_base_) {
            munmap(mapped_base_, mapped_size_);
            mapped_base_ = nullptr;
            mapped_size_ = 0;
        }
        header_ = {};
        chr_offsets_owned_.clear();
        tags_owned_.clear();
        payloads_owned_.clear();
        postings_owned_.clear();
        chrom_names_owned_.clear();
        occupancy_bits_owned_.clear();
        rank_index_owned_.clear();
        shard_dir_owned_.clear();
        packed_reference_owned_.clear();
        chr_offsets_ = nullptr;
        tags_ = nullptr;
        payloads_ = nullptr;
        postings_ = nullptr;
        occupancy_bits_ = nullptr;
        rank_index_ = nullptr;
        shard_dir_ = nullptr;
        shard_mask_ = 0;
        packed_reference_ = nullptr;
    }

bool FaixIndex::save(const std::string& path) const {
        FaixTransactionalFile out(path);
        if (!out.opened()) return false;
        uint64_t expected_bytes = 0;
        if (!write_image(out, expected_bytes)) return false;
        return out.publish(expected_bytes);
    }

bool FaixIndex::write_image(FaixByteSink& out,
                            uint64_t& expected_bytes) const {
        expected_bytes = 0;
        if (!chr_offsets_ || !tags_ || !payloads_ || !occupancy_bits_ ||
            !rank_index_ || !shard_dir_ ||
            chrom_names_owned_.size() != static_cast<size_t>(header_.contig_count)) {
            return false;
        }

        FaixHeader out_header = header_;
        out_header.version = kFaixFormatVersion;
        const bool write_reference = has_reference_payload();
        if (write_reference)
            out_header.flags |= kFaixFlagHasReference;
        else
            out_header.flags &= ~kFaixFlagHasReference;
        const uint64_t reference_bytes =
            write_reference ? packed_reference_bytes(out_header.total_bp) : 0;

        char header[kFaixHeaderBytes] = {};
        faix_header_to_bytes(out_header, header);
        const uint64_t chrom_entries = out_header.contig_count + 1;
        const uint64_t n_shards = 1ULL << out_header.shard_bits;
        uint64_t occ_total = 0, rank_total = 0;
        for (uint64_t s = 0; s < n_shards; ++s) {
            occ_total += occupancy_storage_bytes(shard_dir_[s].bucket_count);
            rank_total += rank_index_count(shard_dir_[s].bucket_count,
                                           out_header.rank_block_bits);
        }
        uint64_t chrom_bytes = 0, shard_bytes = 0, rank_bytes = 0;
        uint64_t payload_bytes = 0, posting_bytes = 0, tag_bytes = 0;
        if (!checked_product(chrom_entries, sizeof(uint64_t), chrom_bytes) ||
            !checked_product(n_shards, sizeof(FaixShardEntry), shard_bytes) ||
            !checked_product(rank_total, sizeof(uint64_t), rank_bytes) ||
            !checked_product(out_header.key_count, sizeof(uint64_t),
                             payload_bytes) ||
            !checked_product(out_header.multi_posting_count, sizeof(PackedRefPos),
                             posting_bytes) ||
            !checked_product(out_header.key_count, sizeof(uint32_t), tag_bytes)) {
          return false;
        }
        expected_bytes = kFaixHeaderBytes;
        if (!checked_add(expected_bytes, chrom_bytes) ||
            !checked_add(expected_bytes, shard_bytes) ||
            !checked_add(expected_bytes, occ_total) ||
            !checked_add(expected_bytes, rank_bytes) ||
            !checked_add(expected_bytes, payload_bytes) ||
            !checked_add(expected_bytes, posting_bytes) ||
            !checked_add(expected_bytes, tag_bytes) ||
            !checked_add(expected_bytes, out_header.name_bytes) ||
            !checked_add(expected_bytes, reference_bytes)) {
          expected_bytes = 0;
          return false;
        }

        if (!out.write_bytes(header, kFaixHeaderBytes) ||
            !out.write_bytes(chr_offsets_, chrom_bytes) ||
            !out.write_bytes(shard_dir_, shard_bytes) ||
            !out.write_bytes(occupancy_bits_, occ_total) ||
            !out.write_bytes(rank_index_, rank_bytes) ||
            !out.write_bytes(payloads_, payload_bytes)) {
          return false;
        }
        if (out_header.multi_posting_count > 0) {
            if (!postings_ || !out.write_bytes(postings_, posting_bytes))
              return false;
        }
        if (!out.write_bytes(tags_, tag_bytes))
          return false;
        for (const auto& name : chrom_names_owned_) {
            const uint32_t len = static_cast<uint32_t>(name.size());
            if (!out.write_bytes(&len, sizeof(len)) ||
                !out.write_bytes(name.data(), name.size()))
              return false;
        }
        if (write_reference &&
            !out.write_bytes(packed_reference_, reference_bytes)) {
          return false;
        }
        return true;
    }

FaixProbe probe_faix_header(const std::string& path) {
    FaixProbe probe;
    // Only a regular file can be an index; a named pipe is not even opened.
    struct stat st;
    if (path.empty() || stat(path.c_str(), &st) != 0) {
        probe.status.error = FaixLoadError::OpenFailed;
        return probe;
    }
    if (!S_ISREG(st.st_mode)) {
        probe.status.error = FaixLoadError::NotIndex;
        return probe;
    }
    const FdCloser file{open(path.c_str(), O_RDONLY)};
    if (file.fd < 0) {
        probe.status.error = FaixLoadError::OpenFailed;
        return probe;
    }
    FaixHeader h;
    h.version = 0;
    probe.status.error =
        read_header(file.fd, 0, static_cast<uint64_t>(st.st_size), h);
    probe.status.observed_version = h.version;
    probe.is_index = probe.status.error != FaixLoadError::OpenFailed &&
                     probe.status.error != FaixLoadError::ReadFailed &&
                     probe.status.error != FaixLoadError::NotIndex;
    if (probe.status)
        probe.header = h;
    return probe;
}

FaixMultipartTable load_faix_multipart_table(const std::string& path,
                                             FaixLoadStatus* status) {
    FaixMultipartTable table;
    if (status) *status = {};
    const auto fail = [&](FaixLoadError error) {
      if (status) status->error = error;
      return FaixMultipartTable{};
    };

    const FdCloser file{open(path.c_str(), O_RDONLY)};
    if (file.fd < 0) return fail(FaixLoadError::OpenFailed);
    struct stat st;
    if (fstat(file.fd, &st) != 0) return fail(FaixLoadError::ReadFailed);
    const uint64_t file_size = static_cast<uint64_t>(st.st_size);
    FaixHeader& h = table.header;
    h.version = 0;
    const FaixLoadError header_error = read_header(file.fd, 0, file_size, h);
    if (status) status->observed_version = h.version;
    if (header_error != FaixLoadError::None) return fail(header_error);
    if ((h.flags & kFaixFlagMultipart) == 0) return fail(FaixLoadError::InvalidHeader);

    const uint64_t chrom_entries = h.contig_count + 1;
    const uint64_t chrom_bytes = chrom_entries * sizeof(uint64_t);
    const uint64_t part_bytes = h.part_count * sizeof(FaixPartEntry);
    uint64_t table_bytes = kFaixHeaderBytes + chrom_bytes + part_bytes;
    if (!checked_add(table_bytes, h.name_bytes)) return fail(FaixLoadError::InvalidHeader);
    if (table_bytes > file_size) return fail(FaixLoadError::Truncated);
    // One read for everything before the images.
    std::vector<char> body(static_cast<size_t>(table_bytes - kFaixHeaderBytes));
    if (!pread_exact(file.fd, body.data(), body.size(), kFaixHeaderBytes))
        return fail(FaixLoadError::ReadFailed);

    const char* cursor = body.data();
    table.chr_offsets.resize(static_cast<size_t>(chrom_entries));
    std::memcpy(table.chr_offsets.data(), cursor, static_cast<size_t>(chrom_bytes));
    cursor += chrom_bytes;
    if (!validate_chrom_offsets(table.chr_offsets, h.contig_count, h.total_bp))
        return fail(FaixLoadError::InvalidLayout);
    for (uint64_t c = 0; c < h.contig_count; ++c) {
        if (table.chr_offsets[static_cast<size_t>(c + 1)] -
                table.chr_offsets[static_cast<size_t>(c)] >
            kFaixMaxContigBp) {
            return fail(FaixLoadError::ContigTooLong);
        }
    }
    std::vector<FaixPartEntry> parts(static_cast<size_t>(h.part_count));
    std::memcpy(parts.data(), cursor, static_cast<size_t>(part_bytes));
    cursor += part_bytes;
    // The parts tile [0, contig_count) in order and their images tile
    // [table_bytes, file_size) in the same order, with no gap or overlap.
    uint64_t next_contig = 0;
    uint64_t next_image = table_bytes;
    for (const FaixPartEntry& e : parts) {
        if (e.contig_count == 0 || e.first_contig != next_contig ||
            e.contig_count > h.contig_count - next_contig ||
            e.image_offset != next_image ||
            e.image_bytes < kFaixHeaderBytes ||
            e.image_bytes > file_size - next_image) {
            return fail(FaixLoadError::InvalidLayout);
        }
        next_contig += e.contig_count;
        next_image += e.image_bytes;
    }
    if (next_contig != h.contig_count || next_image != file_size)
        return fail(FaixLoadError::InvalidLayout);
    if (!parse_chrom_names(cursor, h.name_bytes, h.contig_count, table.chrom_names))
        return fail(FaixLoadError::InvalidLayout);
    table.parts = std::move(parts);

    // Every part is an image of the same format and seeding, for exactly its contigs.
    for (std::size_t i = 0; i < table.parts.size(); ++i) {
        const FaixPartEntry& e = table.parts[i];
        FaixHeader part;
        part.version = 0;
        const FaixLoadError part_error =
            read_header(file.fd, e.image_offset, e.image_bytes, part);
        if (part_error != FaixLoadError::None) {
            if (status) status->observed_version = part.version;
            return fail(part_error == FaixLoadError::NotIndex
                            ? FaixLoadError::InvalidLayout
                            : part_error);
        }
        if ((part.flags & kFaixFlagMultipart) != 0 ||
            (part.flags & kFaixFlagHasReference) !=
                (h.flags & kFaixFlagHasReference) ||
            part.k != h.k || part.syncmer_s != h.syncmer_s ||
            part.contig_count != e.contig_count ||
            part.total_bp != table.part_bp(i) ||
            std::memcmp(part.preset, h.preset, kFaixPresetBytes) != 0) {
            return fail(FaixLoadError::InvalidLayout);
        }
    }
    return table;
}

FaixIndex FaixIndex::load_image(
    const std::string& path, uint64_t image_offset, uint64_t image_bytes,
    FaixLoadStatus* status, int n_threads) {
        FaixIndex idx;
        if (status) *status = {};
        const auto fail = [&](FaixLoadError error) {
          if (status) status->error = error;
          return FaixIndex{};
        };

        const FdCloser file{open(path.c_str(), O_RDONLY)};
        if (file.fd < 0) return fail(FaixLoadError::OpenFailed);
        struct stat st;
        if (fstat(file.fd, &st) != 0) return fail(FaixLoadError::ReadFailed);
        // The image is a span of the file: the whole of it unless a container
        // part was named. A span that runs past the file is a truncated image.
        const uint64_t whole_file = static_cast<uint64_t>(st.st_size);
        if (image_offset == 0 && image_bytes == 0) image_bytes = whole_file;
        if (image_offset > whole_file || image_bytes > whole_file - image_offset)
            return fail(FaixLoadError::Truncated);

        FaixHeader& h = idx.header_;
        h.version = 0;
        const FaixLoadError header_error =
            read_header(file.fd, image_offset, image_bytes, h);
        if (status) status->observed_version = h.version;
        if (header_error != FaixLoadError::None) return fail(header_error);
        if ((h.flags & kFaixFlagMultipart) != 0) return fail(FaixLoadError::Multipart);
        const bool has_reference = (h.flags & kFaixFlagHasReference) != 0;

        // The layout's size follows from the header and the shard directory, so a truncated
        // or damaged file is refused before the image is read. The header bounds keep these
        // products far from overflow.
        const uint64_t chrom_entries = h.contig_count + 1;
        const uint64_t chrom_bytes = chrom_entries * sizeof(uint64_t);
        const uint64_t n_shards = 1ULL << h.shard_bits;
        const uint64_t shard_bytes = n_shards * sizeof(FaixShardEntry);
        if (kFaixHeaderBytes + chrom_bytes + shard_bytes > image_bytes)
            return fail(FaixLoadError::Truncated);
        std::vector<FaixShardEntry> dir(static_cast<size_t>(n_shards));
        if (!pread_exact(file.fd, dir.data(), shard_bytes,
                         image_offset + kFaixHeaderBytes + chrom_bytes)) {
            return fail(FaixLoadError::ReadFailed);
        }
        uint64_t occ_total = 0, rank_total = 0, bucket_total = 0;
        uint64_t post_total = 0, slot_total = 0;
        for (const FaixShardEntry& e : dir) {
            if ((e.bucket_count & (e.bucket_count - 1)) != 0 ||
                e.unique_keys > e.bucket_count ||
                e.unique_keys > std::numeric_limits<uint32_t>::max() ||
                e.posting_count > std::numeric_limits<uint32_t>::max() ||
                (e.bucket_count + 7ULL) / 8ULL > image_bytes ||
                e.occ_off != occ_total || e.rank_off != rank_total ||
                e.bucket_off != bucket_total || e.post_off != post_total) {
                return fail(FaixLoadError::InvalidLayout);
            }
            if (!checked_add(occ_total, occupancy_storage_bytes(e.bucket_count)) ||
                !checked_add(rank_total,
                             rank_index_count(e.bucket_count, h.rank_block_bits)) ||
                !checked_add(bucket_total, e.unique_keys) ||
                !checked_add(post_total, e.posting_count) ||
                !checked_add(slot_total, e.bucket_count)) {
                return fail(FaixLoadError::InvalidLayout);
            }
        }
        if (bucket_total != h.key_count || post_total != h.multi_posting_count ||
            slot_total != h.slot_count) {
            return fail(FaixLoadError::InvalidLayout);
        }
        uint64_t rank_bytes = 0, tag_bytes = 0, payload_bytes = 0;
        uint64_t posting_bytes = 0;
        if (!checked_product(rank_total, sizeof(uint64_t), rank_bytes) ||
            !checked_product(bucket_total, sizeof(uint32_t), tag_bytes) ||
            !checked_product(bucket_total, sizeof(uint64_t), payload_bytes) ||
            !checked_product(post_total, sizeof(PackedRefPos), posting_bytes)) {
            return fail(FaixLoadError::InvalidLayout);
        }
        const uint64_t reference_bytes =
            has_reference ? packed_reference_bytes(h.total_bp) : 0;
        uint64_t expected_bytes = kFaixHeaderBytes;
        if (!checked_add(expected_bytes, chrom_bytes) ||
            !checked_add(expected_bytes, shard_bytes) ||
            !checked_add(expected_bytes, occ_total) ||
            !checked_add(expected_bytes, rank_bytes) ||
            !checked_add(expected_bytes, payload_bytes) ||
            !checked_add(expected_bytes, posting_bytes) ||
            !checked_add(expected_bytes, tag_bytes) ||
            !checked_add(expected_bytes, h.name_bytes) ||
            !checked_add(expected_bytes, reference_bytes)) {
            return fail(FaixLoadError::InvalidLayout);
        }
        if (expected_bytes > image_bytes) return fail(FaixLoadError::Truncated);
        if (expected_bytes != image_bytes) return fail(FaixLoadError::InvalidLayout);

        // An anonymous image of the file, not a file mapping (see read_file_image()).
        const size_t file_size = static_cast<size_t>(image_bytes);
        void* base = mmap(nullptr, file_size, PROT_READ | PROT_WRITE,
                          MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (base == MAP_FAILED) return fail(FaixLoadError::OutOfMemory);
        if (!read_file_image(file.fd, base, image_offset, image_bytes, n_threads)) {
            munmap(base, file_size);
            return fail(FaixLoadError::ReadFailed);
        }
        // Read-only from here on, so a stray write faults.
        (void)mprotect(base, file_size, PROT_READ);
#ifdef POSIX_FADV_DONTNEED
        // The page cache holds a second copy of what was just read; drop it.
        (void)posix_fadvise(file.fd, static_cast<off_t>(image_offset),
                            static_cast<off_t>(image_bytes),
                            POSIX_FADV_DONTNEED);
#endif
        const auto unmap_fail = [&](FaixLoadError error) {
          munmap(base, file_size);
          return fail(error);
        };

        const char* body = static_cast<const char*>(base) + kFaixHeaderBytes;
        const uint64_t* disk_chr_offsets = reinterpret_cast<const uint64_t*>(body);
        idx.chr_offsets_owned_.assign(disk_chr_offsets,
                                      disk_chr_offsets + static_cast<size_t>(chrom_entries));
        if (!validate_chrom_offsets(idx.chr_offsets_owned_, h.contig_count, h.total_bp))
            return unmap_fail(FaixLoadError::InvalidLayout);
        for (uint64_t c = 0; c < h.contig_count; ++c) {
            if (idx.chr_offsets_owned_[c + 1] - idx.chr_offsets_owned_[c] >
                kFaixMaxContigBp) {
                return unmap_fail(FaixLoadError::ContigTooLong);
            }
        }
        idx.chr_offsets_ = idx.chr_offsets_owned_.data();

        const char* cursor = body + chrom_bytes;
        idx.shard_dir_ = reinterpret_cast<const FaixShardEntry*>(cursor);
        cursor += shard_bytes;
        idx.occupancy_bits_ = reinterpret_cast<const uint8_t*>(cursor);
        cursor += occ_total;
        idx.rank_index_ = reinterpret_cast<const uint64_t*>(cursor);
        cursor += rank_bytes;
        idx.payloads_ = reinterpret_cast<const uint64_t*>(cursor);
        cursor += payload_bytes;
        idx.postings_ =
            post_total > 0 ? reinterpret_cast<const PackedRefPos*>(cursor) : nullptr;
        cursor += posting_bytes;
        idx.tags_ = reinterpret_cast<const uint32_t*>(cursor);
        cursor += tag_bytes;
        if (!parse_chrom_names(cursor, h.name_bytes, h.contig_count,
                               idx.chrom_names_owned_)) {
            return unmap_fail(FaixLoadError::InvalidLayout);
        }
        cursor += h.name_bytes;
        if (has_reference)
            idx.packed_reference_ = reinterpret_cast<const uint8_t*>(cursor);
        idx.shard_mask_ = n_shards - 1;
        idx.mapped_base_ = base;
        idx.mapped_size_ = file_size;
        return idx;
    }

}}  // namespace fa::cpu
