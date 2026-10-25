// Using a .faix file: load and probe, header decoding, the 4-bit reference codec, the name table, the occurrence distribution, the posting lookup.
#include "index.h"
#include "../threading/threading.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

/**** faix_storage_posix ****/
// POSIX storage for .faix files: save, load, probe and the container table.

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

/**** faix_format ****/
// The header decoder, the load errors, the 4-bit reference codec and the name table.

namespace fa { namespace cpu {

namespace {

// One unpack task: [lo, hi) in global base coordinates within contig `chrom`. Large enough
// to amortize task overhead, small enough that one long contig spreads over all threads.
constexpr uint64_t kUnpackSliceBases = 32ULL << 20;

struct UnpackSlice {
    uint64_t chrom;
    uint64_t lo;
    uint64_t hi;
};

}  // namespace

FaixLoadError faix_header_decode(const char* p, FaixHeader& h) {
    if (std::memcmp(p, kFaixDevelopmentMagic, 8) == 0)
        return FaixLoadError::DevelopmentBuild;
    if (std::memcmp(p, kFaixMagic, 8) != 0)
        return FaixLoadError::NotIndex;
    std::memcpy(&h.version, p + 8, 4);
    if (h.version == 0)
        return FaixLoadError::InvalidHeader;
    if (h.version > kFaixFormatVersion)
        return FaixLoadError::NewerFormat;
    std::memcpy(&h.flags, p + 12, 4);
    std::memcpy(&h.k, p + 16, 4);
    std::memcpy(&h.syncmer_s, p + 20, 4);
    std::memcpy(&h.shard_bits, p + 24, 4);
    std::memcpy(&h.rank_block_bits, p + 28, 4);
    std::memcpy(&h.contig_count, p + 32, 8);
    std::memcpy(&h.total_bp, p + 40, 8);
    std::memcpy(&h.name_bytes, p + 48, 8);
    std::memcpy(&h.key_count, p + 56, 8);
    std::memcpy(&h.slot_count, p + 64, 8);
    std::memcpy(&h.multi_posting_count, p + 72, 8);
    std::memcpy(&h.part_count, p + 80, 8);
    std::memcpy(h.preset, p + kFaixPresetOffset, kFaixPresetBytes);
    // A later writer may use a flag or the zero bytes without a new version.
    if ((h.flags & ~kFaixKnownFlags) != 0)
        return FaixLoadError::NewerFormat;
    for (size_t i = 88; i < kFaixPresetOffset; ++i) {
        if (p[i] != 0)
            return FaixLoadError::NewerFormat;
    }

    const int k = static_cast<int>(h.k);
    if (h.k < 1 || h.k > static_cast<uint32_t>(kFaixMaxK) || h.syncmer_s < 1 ||
        h.syncmer_s > h.k || h.contig_count == 0 ||
        h.contig_count > kFaixMaxContigCount) {
        return FaixLoadError::InvalidHeader;
    }
    if ((h.flags & kFaixFlagMultipart) != 0) {
        if (h.part_count == 0 || h.part_count > h.contig_count ||
            h.shard_bits != 0 || h.rank_block_bits != 0 || h.key_count != 0 ||
            h.slot_count != 0 || h.multi_posting_count != 0) {
            return FaixLoadError::InvalidHeader;
        }
    } else if (h.part_count != 0 || h.shard_bits < 1 || h.shard_bits > 24 ||
               2 * k - static_cast<int>(h.shard_bits) > 31 ||
               h.rank_block_bits < 1 || h.rank_block_bits > 24 ||
               h.key_count > h.slot_count) {
        return FaixLoadError::InvalidHeader;
    }
    return FaixLoadError::None;
}

std::string faix_load_error_message(const FaixLoadStatus& status) {
    switch (status.error) {
      case FaixLoadError::None: return "no error";
      case FaixLoadError::OpenFailed: return "cannot open the file";
      case FaixLoadError::ReadFailed: return "cannot read the file";
      case FaixLoadError::OutOfMemory: return "not enough memory";
      case FaixLoadError::NotIndex: return "not a FlashAlign index";
      case FaixLoadError::DevelopmentBuild:
        return "an index from a FlashAlign development build; rebuild it with "
               "`flashalign index`";
      case FaixLoadError::NewerFormat:
        if (status.observed_version > kFaixFormatVersion) {
          return "written by a newer FlashAlign (index format " +
                 std::to_string(status.observed_version) +
                 "; this build reads format " +
                 std::to_string(kFaixFormatVersion) + ")";
        }
        return "written by a newer FlashAlign (index format " +
               std::to_string(kFaixFormatVersion) +
               " with features this build does not know)";
      case FaixLoadError::Multipart:
        return "a multi-part index (flashalign index -I), which only "
               "`flashalign align` reads";
      case FaixLoadError::Truncated: return "the file is truncated";
      case FaixLoadError::InvalidHeader: return "the header is damaged";
      case FaixLoadError::InvalidLayout: return "the file is damaged";
      case FaixLoadError::ContigTooLong:
        return "a reference sequence is longer than 2147483647 bp, the most SAM "
               "and BAM can address";
    }
    return "unknown error";
}

void pack_reference_4bit_into(
    const uint8_t* reference,
    uint64_t decoded_bytes,
    uint64_t decoded_offset,
    std::vector<uint8_t>& out
) {
    if (decoded_bytes == 0) return;
    if (!reference) return;
    const uint64_t needed = packed_reference_bytes(decoded_offset + decoded_bytes);
    if (out.size() < needed) out.resize(static_cast<size_t>(needed), 0);
    for (uint64_t i = 0; i < decoded_bytes; ++i) {
        const uint64_t pos = decoded_offset + i;
        const uint8_t base = reference[i] & 0x0fU;
        uint8_t& byte = out[static_cast<size_t>(pos >> 1)];
        if ((pos & 1ULL) == 0) {
            byte = static_cast<uint8_t>((byte & 0xf0U) | base);
        } else {
            byte = static_cast<uint8_t>((byte & 0x0fU) | static_cast<uint8_t>(base << 4));
        }
    }
}

void unpack_reference_4bit_to(
    const uint8_t* packed,
    uint64_t lo,
    uint64_t hi,
    uint8_t* out
) {
    if (hi <= lo) return;
    uint64_t i = lo;
    // An odd first base is the high nibble of its byte; peel it so the loop starts even.
    if ((i & 1ULL) != 0) {
        out[0] = static_cast<uint8_t>(
            (packed[static_cast<size_t>(i >> 1)] >> 4) & 0x0fU);
        ++i;
    }
    // Each source byte is two output bases, low nibble first; this loop vectorizes.
    const uint64_t pairs = (hi - i) >> 1;
    const uint8_t* src = packed + static_cast<size_t>(i >> 1);
    uint8_t* dst = out + static_cast<size_t>(i - lo);
    for (uint64_t p = 0; p < pairs; ++p) {
        const uint8_t byte = src[static_cast<size_t>(p)];
        dst[static_cast<size_t>(2 * p)] = static_cast<uint8_t>(byte & 0x0fU);
        dst[static_cast<size_t>(2 * p + 1)] = static_cast<uint8_t>(byte >> 4);
    }
    i += pairs * 2;
    // A leftover last base, at an even index, is the low nibble of its byte.
    if (i < hi) {
        out[static_cast<size_t>(i - lo)] = static_cast<uint8_t>(
            packed[static_cast<size_t>(i >> 1)] & 0x0fU);
    }
}

// Contigs are independent and a contig's output bases split into disjoint ranges, so the
// unpack runs in parallel: one pass sizes the contigs, one unpacks the slices.
bool FaixIndex::reference_chromosomes_u8(
    std::vector<std::vector<uint8_t>>& out,
    int n_threads
) const {
    out.clear();
    if (!has_reference_payload() || !chr_offsets_) return false;
    const size_t n_chrom = static_cast<size_t>(header_.contig_count);

    // Validate the whole offset table first, so a refusal leaves `out` empty.
    for (size_t i = 0; i < n_chrom; ++i) {
        if (chr_offsets_[i + 1] < chr_offsets_[i] ||
            chr_offsets_[i + 1] > header_.total_bp) {
            return false;
        }
    }

    std::vector<UnpackSlice> slices;
    for (size_t i = 0; i < n_chrom; ++i) {
        const uint64_t end = chr_offsets_[i + 1];
        for (uint64_t lo = chr_offsets_[i]; lo < end;
             lo += kUnpackSliceBases) {
            slices.push_back({static_cast<uint64_t>(i), lo,
                              std::min(lo + kUnpackSliceBases, end)});
        }
    }

    int threads = n_threads > 0 ? n_threads
                                : threading::default_thread_count();
    if (threads < 1) threads = 1;

    out.resize(n_chrom);
    // vector::resize zero-fills every base; spread that over the workers, one per contig.
    threading::parallel_for(
        threads, static_cast<int64_t>(n_chrom), [&](int64_t i, int) {
            const size_t c = static_cast<size_t>(i);
            out[c].resize(
                static_cast<size_t>(chr_offsets_[c + 1] - chr_offsets_[c]));
        });

    // Slices partition the output, so workers never write the same byte; neighbours may only
    // read the same packed byte.
    threading::parallel_for(
        threads, static_cast<int64_t>(slices.size()), [&](int64_t s, int) {
            const UnpackSlice& slice = slices[static_cast<size_t>(s)];
            const size_t c = static_cast<size_t>(slice.chrom);
            uint8_t* dst = out[c].data()
                + static_cast<size_t>(slice.lo - chr_offsets_[c]);
            unpack_reference_4bit_to(packed_reference_, slice.lo, slice.hi, dst);
        });
    return true;
}

void FaixIndex::reference_bases_u8(
    uint64_t contig,
    uint64_t lo,
    uint64_t hi,
    uint8_t* out
) const {
    const uint64_t base = chr_offsets_[contig];
    unpack_reference_4bit_to(packed_reference_, base + lo, base + hi, out);
}

bool parse_chrom_names(
    const char* data,
    uint64_t bytes,
    uint64_t count,
    std::vector<std::string>& out
) {
    if (!data && bytes > 0) return false;
    out.clear();
    out.reserve(static_cast<size_t>(count));
    uint64_t offset = 0;
    for (uint64_t i = 0; i < count; ++i) {
        if (offset + sizeof(uint32_t) > bytes) return false;
        uint32_t len = 0;
        std::memcpy(&len, data + offset, sizeof(len));
        offset += sizeof(uint32_t);
        if (offset + static_cast<uint64_t>(len) > bytes) return false;
        out.emplace_back(data + offset, data + offset + len);
        offset += static_cast<uint64_t>(len);
    }
    return offset == bytes;
}

bool validate_chrom_offsets(
    const std::vector<uint64_t>& offsets,
    uint64_t chrom_count,
    uint64_t total_bp
) {
    if (offsets.size() != static_cast<size_t>(chrom_count + 1)) return false;
    if (offsets.empty() || offsets.front() != 0) return false;
    for (size_t i = 1; i < offsets.size(); ++i) {
        if (offsets[i] < offsets[i - 1]) return false;
    }
    return offsets.back() == total_bp;
}

}}  // namespace fa::cpu

/**** occ_stats ****/

namespace fa {
namespace cpu {

FaixOccDistribution faix_occ_distribution(const FaixIndex& index) {
  FaixOccDistribution dist;
  dist.counts.assign(kFaixOccHistogramBound, 0);
  const std::uint32_t* tags = index.tags_data();
  const std::uint64_t* payloads = index.payloads_data();
  if (index.empty() || tags == nullptr || payloads == nullptr)
    return dist;
  const std::uint64_t records = index.distinct_keys();
  for (std::uint64_t r = 0; r < records; ++r) {
    const std::uint32_t occ =
        (tags[r] & kFaixTagMultiBit) == 0
            ? 1u
            : static_cast<std::uint32_t>(payloads[r] >> 32);
    if (occ == 0)
      continue; // structurally impossible; never counted
    ++dist.n_distinct;
    dist.total_postings += occ;
    if (occ <= kFaixOccHistogramBound) {
      ++dist.counts[occ - 1];
    } else {
      ++dist.tail[occ];
    }
  }
  return dist;
}

std::int64_t faix_cal_max_occ(const FaixOccDistribution& dist,
                              double f) {
  if (!(f > 0.0))
    return static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max());
  if (dist.n_distinct == 0)
    return 0;
  const double position = (1.0 - f) * static_cast<double>(dist.n_distinct);
  std::uint64_t rank =
      position <= 0.0 ? 0 : static_cast<std::uint64_t>(position);
  if (rank >= dist.n_distinct)
    rank = dist.n_distinct - 1;
  std::uint64_t cumulative = 0;
  for (std::uint32_t v = 1; v <= kFaixOccHistogramBound; ++v) {
    cumulative += dist.counts[v - 1];
    if (cumulative > rank)
      return static_cast<std::int64_t>(v) + 1;
  }
  for (const auto& entry : dist.tail) {
    cumulative += entry.second;
    if (cumulative > rank)
      return static_cast<std::int64_t>(entry.first) + 1;
  }
  return 0; // unreachable: the ranks are bounded by n_distinct
}

} // namespace cpu
} // namespace fa

/**** faix ****/

namespace fa { namespace cpu {

bool FaixMultipartTable::loaded() const noexcept { return !parts.empty(); }

std::uint64_t FaixMultipartTable::part_bp(std::size_t i) const noexcept {
    const FaixPartEntry& e = parts[i];
    return chr_offsets[static_cast<std::size_t>(e.first_contig +
                                                e.contig_count)] -
           chr_offsets[static_cast<std::size_t>(e.first_contig)];
}

bool FaixIndex::empty() const { return header_.slot_count == 0; }

int FaixIndex::k() const { return static_cast<int>(header_.k); }

int FaixIndex::build_syncmer_s() const { return static_cast<int>(header_.syncmer_s); }

std::string FaixIndex::build_preset() const { return faix_header_preset(header_); }

uint64_t FaixIndex::total_bp() const { return header_.total_bp; }

uint64_t FaixIndex::distinct_keys() const { return header_.key_count; }

uint64_t FaixIndex::memory_bytes() const {
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

const uint32_t* FaixIndex::tags_data() const { return tags_; }

const uint64_t* FaixIndex::payloads_data() const { return payloads_; }

const std::vector<std::string>&
FaixIndex::chromosome_names() const { return chrom_names_owned_; }

bool FaixIndex::has_reference_payload() const {
    return (header_.flags & kFaixFlagHasReference) != 0 &&
           packed_reference_ != nullptr;
}

bool FaixIndex::set_reference_payload(
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

bool FaixIndex::set_chromosome_names(std::vector<std::string> chrom_names) {
    if (!set_names(std::move(chrom_names)))
        return false;
    packed_reference_owned_.clear();
    packed_reference_owned_.shrink_to_fit();
    packed_reference_ = nullptr;
    header_.flags &= ~kFaixFlagHasReference;
    return true;
}

bool FaixIndex::set_packed_reference_payload(
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

void FaixIndex::lookup_batch(
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

KmerPostingIntervalView FaixIndex::lookup_interval(
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

FaixIndex FaixIndex::load(
    const std::string& path, FaixLoadStatus* status,
    int n_threads) {
    return load_image(path, 0, 0, status, n_threads);
}

FaixIndex
FaixIndex::from_compact_owned(FaixHeader header, std::vector<uint64_t> chr_offsets,
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

void FaixIndex::rebind_views_to_owned() {
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

bool FaixIndex::set_names(std::vector<std::string> chrom_names) {
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

}}  // namespace fa::cpu
