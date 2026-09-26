// The seed-index file (.faix), format version 1.
//
// Every integer is little-endian. A file is one index image, or a multi-part container
// (flashalign index -I) that holds one image per part. Both start with a 128-byte header:
//
//   offset  size  field
//        0     8  magic "FAIX\r\n\x1a\n"
//        8     4  version: 1
//       12     4  flags: bit 0 HAS_REFERENCE, bit 1 MULTIPART; every other bit is 0
//       16     4  k, 1 to 23
//       20     4  s, the closed-syncmer s-mer length, 1 to k
//       24     4  shard_bits, 1 to 24 with 2k - shard_bits <= 31 (0 in a container)
//       28     4  rank_block_bits, 1 to 24 (0 in a container)
//       32     8  contig_count, 1 to 2^31 - 1
//       40     8  total_bp, the sum of the contig lengths
//       48     8  name_bytes
//       56     8  key_count (0 in a container)
//       64     8  slot_count (0 in a container)
//       72     8  multi_posting_count (0 in a container)
//       80     8  part_count (0 in an image)
//       88    24  zero
//      112    16  preset name, NUL-padded; all zero when none is recorded
//
// An image continues with:
//   offsets    contig_count + 1 u64: 0, then the running sum of the contig lengths, each
//              at most 2^31 - 1
//   shards     2^shard_bits FaixShardEntry
//   occupancy  per shard, one bit per slot (bit i & 7 of byte i >> 3), padded to a multiple
//              of 8 bytes
//   rank       per shard of n > 0 slots, ceil(n / 2^rank_block_bits) + 1 u64: entry j counts
//              the occupied slots below j * 2^rank_block_bits, the last one all of them
//   payloads   key_count u64
//   postings   multi_posting_count u64
//   tags       key_count u32
//   names      contig_count names, each a u32 length and its bytes, name_bytes in all
//   reference  with HAS_REFERENCE, (total_bp + 1) / 2 bytes: base i is the low nibble of
//              byte i / 2 when i is even, else the high one; A C G T are 0 1 2 3, any other
//              base is 4, and a last unused nibble is 0
// Each shard's slices follow the previous shard's, at the offsets its FaixShardEntry holds.
//
// A container continues with:
//   offsets    as in an image, for the whole reference
//   parts      part_count FaixPartEntry, in contig order
//   names      as in an image, for the whole reference
//   images     one image per part, back to back; contig i of a part is contig
//              first_contig + i of the reference
//
// A k-mer is coded 2 bits a base, A 0, C 1, G 2, T 3, first base highest; a k-mer with any
// other base is never a seed. Its key is min(fwd, rc) of its code and its reverse
// complement's, and its orientation bit z is 1 when rc < fwd. A palindrome (fwd == rc) is
// never a seed. Reference seeds are closed syncmers: of the k - s + 1 s-mers, one with the
// least sketch_order_hash(min(fwd, rc)) is the first or the last.
//
// For x = shard_mix_2k(key, 2k), a key lives in shard x & (2^shard_bits - 1) under tag
// x >> shard_bits. A shard of n slots, a power of two, is an open-addressed table probed
// linearly from slot splitmix64(tag) & (n - 1); a lookup stops at an unoccupied slot. The
// record of an occupied slot is the number of occupied slots before it in its shard, and
// record r of a shard has tag word tag << 1 | multi at tags[r] of the shard's slice and its
// payload at payloads[r]. With multi 0 the payload is the key's one position; with multi 1
// it is count << 32 | begin, and the key's count >= 2 positions are
// postings[begin, begin + count) of the shard's slice, ascending. A position is
// contig << 33 | offset << 1 | z, z being the orientation bit of the reference k-mer there.
//
// The hash functions are in core/hash.h. A reader refuses another version, a flag it does
// not know and non-zero zero bytes. Changing any rule above, or the layout, needs format
// version 2.
#pragma once

#include "../core/hash.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace fa { namespace cpu {

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
              "the index file is little-endian and is read and written in place");
#endif

static constexpr char kFaixMagic[8] = {'F', 'A', 'I', 'X', '\r', '\n', '\x1a', '\n'};
// The magic of indexes written by FlashAlign development builds, refused by name.
static constexpr char kFaixDevelopmentMagic[8] = {'F', 'A', 'H', 'T', 'B', 'L', '0', '2'};
static constexpr uint32_t kFaixFormatVersion = 1;
static constexpr size_t kFaixHeaderBytes = 128;
static constexpr size_t kFaixPresetOffset = 112;
static constexpr size_t kFaixPresetBytes = 16;

static constexpr uint32_t kFaixFlagHasReference = 1u << 0;
static constexpr uint32_t kFaixFlagMultipart = 1u << 1;
static constexpr uint32_t kFaixKnownFlags = kFaixFlagHasReference | kFaixFlagMultipart;

static constexpr int kFaixMaxK = 23;
// SAM and BAM positions are signed 32-bit.
static constexpr uint64_t kFaixMaxContigBp =
    static_cast<uint64_t>(std::numeric_limits<int32_t>::max());
// A position word has 31 bits for the contig, which must also hold the one-past-the-end
// sentinel (contig_count, 0).
static constexpr uint64_t kFaixMaxContigCount = (static_cast<uint64_t>(1) << 31) - 1;
// The rank block this build writes; readers take it from the header.
static constexpr uint32_t kFaixRankBlockBits = 9;
static constexpr uint32_t kFaixTagMultiBit = 1u;

// Why an index file cannot be read.
enum class FaixLoadError : std::uint8_t {
    None,
    OpenFailed,
    ReadFailed,
    OutOfMemory,
    NotIndex,
    DevelopmentBuild,
    NewerFormat,
    Multipart,  // a container given to the single-image loader
    Truncated,
    InvalidHeader,
    InvalidLayout,
    ContigTooLong,
};

struct FaixLoadStatus {
    FaixLoadError error = FaixLoadError::None;
    std::uint32_t observed_version = 0;

    explicit operator bool() const noexcept {
        return error == FaixLoadError::None;
    }
};

std::string faix_load_error_message(const FaixLoadStatus& status);

// The header fields; the magic and the zero bytes are implied.
struct FaixHeader {
    uint32_t version = kFaixFormatVersion;
    uint32_t flags = 0;
    uint32_t k = 0;
    uint32_t syncmer_s = 0;
    uint32_t shard_bits = 0;
    uint32_t rank_block_bits = 0;
    uint64_t contig_count = 0;
    uint64_t total_bp = 0;
    uint64_t name_bytes = 0;
    uint64_t key_count = 0;
    uint64_t slot_count = 0;
    uint64_t multi_posting_count = 0;
    uint64_t part_count = 0;
    char preset[kFaixPresetBytes] = {};
};

// `out` has room for kFaixHeaderBytes.
inline void faix_header_to_bytes(const FaixHeader& h, char* out) {
    std::memset(out, 0, kFaixHeaderBytes);
    std::memcpy(out + 0, kFaixMagic, 8);
    std::memcpy(out + 8, &h.version, 4);
    std::memcpy(out + 12, &h.flags, 4);
    std::memcpy(out + 16, &h.k, 4);
    std::memcpy(out + 20, &h.syncmer_s, 4);
    std::memcpy(out + 24, &h.shard_bits, 4);
    std::memcpy(out + 28, &h.rank_block_bits, 4);
    std::memcpy(out + 32, &h.contig_count, 8);
    std::memcpy(out + 40, &h.total_bp, 8);
    std::memcpy(out + 48, &h.name_bytes, 8);
    std::memcpy(out + 56, &h.key_count, 8);
    std::memcpy(out + 64, &h.slot_count, 8);
    std::memcpy(out + 72, &h.multi_posting_count, 8);
    std::memcpy(out + 80, &h.part_count, 8);
    std::memcpy(out + kFaixPresetOffset, h.preset, kFaixPresetBytes);
}

// Decodes kFaixHeaderBytes at `p` and checks everything the header alone can show: the
// magic, the version, the flags, the zero bytes and the ranges of the fields. Anything but
// None leaves `h` unspecified, except that `h.version` is set once the magic matches.
FaixLoadError faix_header_decode(const char* p, FaixHeader& h);

// The recorded preset name, or "" when none is recorded. A field with no NUL reads as none.
inline std::string faix_header_preset(const FaixHeader& header) {
    for (size_t i = 0; i < kFaixPresetBytes; ++i) {
        if (header.preset[i] == '\0')
            return std::string(header.preset, i);
    }
    return std::string();
}

// False when `name` does not fit with its terminator.
inline bool faix_header_set_preset(FaixHeader& header, std::string_view name) {
    if (name.size() >= kFaixPresetBytes)
        return false;
    std::memset(header.preset, 0, kFaixPresetBytes);
    std::memcpy(header.preset, name.data(), name.size());
    return true;
}

// One part of a container.
struct FaixPartEntry {
    uint64_t image_offset = 0;  // byte offset of the part's image in the file
    uint64_t image_bytes = 0;
    uint64_t first_contig = 0;  // the reference's id of the part's first contig
    uint64_t contig_count = 0;  // at least 1
};
static_assert(sizeof(FaixPartEntry) == 32, "FaixPartEntry must stay packed");

// Reference positions are flattened-global uint64 at run time and packed position words
// (see above) in the index. With z in bit 0, unsigned order is (contig, offset, z); bounds
// use z = 0, so a lower bound at (contig, p) admits both orientations at p and an exclusive
// upper bound excludes both.
using RefPos = uint64_t;
using PackedRefPos = uint64_t;

inline PackedRefPos pack_ref_pos(uint32_t contig_id, uint32_t local_pos,
                                 uint32_t z = 0) {
  return (static_cast<uint64_t>(contig_id) << 33) |
         (static_cast<uint64_t>(local_pos) << 1) |
         static_cast<uint64_t>(z & 1u);
}

inline uint32_t packed_ref_contig(PackedRefPos pos) {
  return static_cast<uint32_t>(pos >> 33);
}

inline uint32_t packed_ref_local(PackedRefPos pos) {
  return static_cast<uint32_t>(pos >> 1);
}

// The only accessor that yields the orientation bit; every positional accessor
// above masks it out.
inline uint32_t packed_ref_strand(PackedRefPos pos) {
  return static_cast<uint32_t>(pos & 1ULL);
}

// Strand compatibility. A posting is reached through the query seed's canonical key, so
// the query and reference k-mers share it, and each one's forward bases are the key
// (z = 0) or its reverse complement (z = 1). Their forward bases match on the forward lane
// exactly when the bits agree, and on the reverse lane when they differ:
// (z_query ^ z_ref) == lane. Palindromic keys, which z cannot orient, are never seeds.
// `query_z` is the seed's forward-occurrence bit; it is not flipped for the reverse lane.
inline bool packed_ref_orientation_compatible(uint32_t query_z,
                                              PackedRefPos packed,
                                              bool reverse_lane) {
  return (query_z ^ packed_ref_strand(packed)) == (reverse_lane ? 1u : 0u);
}

inline bool faix_pack_multi_payload(uint64_t posting_begin,
                                    uint64_t posting_count, uint64_t& out) {
  if (posting_begin > std::numeric_limits<uint32_t>::max() ||
      posting_count < 2 ||
      posting_count > std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  out = (posting_count << 32) | posting_begin;
  return true;
}

// Converts a flattened half-open bound to its packed ordering key. total_bp maps to the
// synthetic (chromosome_count, 0) sentinel, so chromosome_count must fit the 31-bit contig
// field. Bounds use z = 0 and are orientation-blind.
inline bool packed_ref_bound_from_global(const uint64_t* chromosome_offsets,
                                         uint64_t chromosome_count,
                                         RefPos global, PackedRefPos& out) {
  if (!chromosome_offsets || chromosome_count == 0 ||
      chromosome_count > kFaixMaxContigCount ||
      global > chromosome_offsets[chromosome_count]) {
    return false;
  }
  if (global == chromosome_offsets[chromosome_count]) {
    out = pack_ref_pos(static_cast<uint32_t>(chromosome_count), 0, 0);
    return true;
  }
  const uint64_t* first_gt =
      std::upper_bound(chromosome_offsets + 1,
                       chromosome_offsets + chromosome_count + 1, global);
  const uint64_t contig =
      static_cast<uint64_t>(first_gt - chromosome_offsets - 1);
  const uint64_t local = global - chromosome_offsets[contig];
  if (contig > kFaixMaxContigCount ||
      local > std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  out = pack_ref_pos(static_cast<uint32_t>(contig),
                     static_cast<uint32_t>(local), 0);
  return true;
}

class RefPosSpan {
public:
  RefPosSpan() = default;
  RefPosSpan(const PackedRefPos* data, uint32_t count,
             const uint64_t* chromosome_offsets, uint64_t chromosome_count)
      : data_(data), count_(count), chromosome_offsets_(chromosome_offsets),
        chromosome_count_(chromosome_count) {}

  bool valid() const noexcept {
    return data_ != nullptr && chromosome_offsets_ != nullptr;
  }
  uint32_t size() const noexcept { return count_; }
  const PackedRefPos* data() const noexcept { return data_; }
  const uint64_t* chromosome_offsets() const noexcept {
    return chromosome_offsets_;
  }
  uint64_t chromosome_count() const noexcept { return chromosome_count_; }
  bool operator==(const RefPosSpan& other) const noexcept {
    return data_ == other.data_ && count_ == other.count_ &&
           chromosome_offsets_ == other.chromosome_offsets_ &&
           chromosome_count_ == other.chromosome_count_;
  }
  bool operator!=(const RefPosSpan& other) const noexcept {
    return !(*this == other);
  }

  PackedRefPos packed_at(uint32_t i) const noexcept { return data_[i]; }
  RefPos operator[](uint32_t i) const noexcept {
    const PackedRefPos packed = data_[i];
    const uint32_t contig = packed_ref_contig(packed);
    return chromosome_offsets_[contig] + packed_ref_local(packed);
  }

  RefPosSpan subspan(uint32_t offset, uint32_t count) const noexcept {
    return data_ ? RefPosSpan(data_ + offset, count, chromosome_offsets_,
                              chromosome_count_)
                 : RefPosSpan{};
  }

  uint32_t lower_bound_packed(PackedRefPos value) const noexcept {
    const PackedRefPos* found = std::lower_bound(data_, data_ + count_, value);
    return static_cast<uint32_t>(found - data_);
  }

  uint32_t lower_bound(RefPos value) const noexcept {
    PackedRefPos packed = 0;
    if (!packed_ref_bound_from_global(chromosome_offsets_, chromosome_count_,
                                      value, packed)) {
      return count_;
    }
    return lower_bound_packed(packed);
  }

private:
  const PackedRefPos* data_ = nullptr;
  uint32_t count_ = 0;
  const uint64_t* chromosome_offsets_ = nullptr;
  uint64_t chromosome_count_ = 0;
};

struct KmerPostingView {
  RefPosSpan positions;
  uint32_t count = 0;
  uint32_t occurrence = 0; // Exact reference occurrence count (never capped).

  // Success always has a non-null pointer; failure returns {}.
  bool found() const { return positions.valid(); }
};

struct KmerPostingIntervalView {
  RefPosSpan positions;
  uint32_t count = 0;
  uint32_t global_count = 0;

  bool found() const { return positions.valid(); }
};

}}  // namespace fa::cpu

namespace fa { namespace cpu {

// chr_bounds: cumulative-bp prefix sums of length n_chr+1, with chr_bounds[0]==0
// and chr_bounds[n_chr]==total_bp. Returns the index of the chromosome whose
// half-open span [chr_bounds[c], chr_bounds[c+1]) contains global_pos, clamped
// to [0, n_chr-1]. Chromosome-local offset is then global_pos - chr_bounds[c].
inline int chromosome_index_for_global_pos(const uint64_t* chr_bounds,
                                           int n_chr, RefPos global_pos) {
  if (!chr_bounds || n_chr <= 1)
    return 0;
  const uint64_t pos = static_cast<uint64_t>(global_pos);
  const uint64_t* const begin = chr_bounds + 1;
  const uint64_t* const end = chr_bounds + n_chr + 1;
  const uint64_t* const first_gt = std::upper_bound(begin, end, pos);
  int chr_idx = static_cast<int>(first_gt - chr_bounds) - 1;
  if (chr_idx < 0)
    chr_idx = 0;
  if (chr_idx >= n_chr)
    chr_idx = n_chr - 1;
  return chr_idx;
}

// Advances a (chr_idx, chr_lo, chr_hi) cursor to the contig holding global position `g`,
// for ascending positions. One step covers the common case and a longer jump binary-
// searches the rest. The cursor never moves backwards (a g below chr_lo is left for the
// caller to reject) and never passes contig n_chr - 1.
inline void advance_contig_cursor(const uint64_t* chr_bounds, int n_chr,
                                  uint64_t g, int& chr_idx, uint64_t& chr_lo,
                                  uint64_t& chr_hi) {
  if (chr_idx + 1 >= n_chr || g < chr_hi)
    return;
  // The common case: the next contig.
  ++chr_idx;
  chr_lo = chr_hi;
  chr_hi = chr_bounds[static_cast<size_t>(chr_idx + 1)];
  if (chr_idx + 1 >= n_chr || g < chr_hi)
    return;
  // Still short: find the first contig c > chr_idx with g < chr_bounds[c + 1], or stop at
  // the last contig.
  const uint64_t* const first = chr_bounds + chr_idx + 2;
  const uint64_t* const last = chr_bounds + n_chr + 1;
  const uint64_t* const it = std::upper_bound(first, last, g);
  chr_idx = (it == last) ? n_chr - 1 : static_cast<int>(it - chr_bounds) - 1;
  chr_lo = chr_bounds[static_cast<size_t>(chr_idx)];
  chr_hi = chr_bounds[static_cast<size_t>(chr_idx + 1)];
}
}}  // namespace fa::cpu

namespace fa { namespace cpu {

inline uint64_t packed_reference_bytes(uint64_t decoded_bytes) {
    return (decoded_bytes + 1ULL) / 2ULL;
}

// Packs `decoded_bytes` bases, one code per byte, into `out` from base `decoded_offset` on,
// growing `out` as needed.
void pack_reference_4bit_into(
    const uint8_t* reference,
    uint64_t decoded_bytes,
    uint64_t decoded_offset,
    std::vector<uint8_t>& out
);

struct FaixByteSink {
  virtual ~FaixByteSink() = default;
  virtual bool write_bytes(const void* data, uint64_t size) = 0;
};

// Unpack the 4-bit reference range [lo, hi) into `out`, which must have room for hi - lo
// bytes. Neither allocates nor bounds-checks.
void unpack_reference_4bit_to(
    const uint8_t* packed,
    uint64_t lo,
    uint64_t hi,
    uint8_t* out
);

// Parses `count` names from the `bytes` bytes at `data`; false unless the names fill them
// exactly.
bool parse_chrom_names(
    const char* data,
    uint64_t bytes,
    uint64_t count,
    std::vector<std::string>& out
);

// Whether `offsets` holds chrom_count + 1 non-decreasing values from 0 to total_bp.
bool validate_chrom_offsets(
    const std::vector<uint64_t>& offsets,
    uint64_t chrom_count,
    uint64_t total_bp
);

}}  // namespace fa::cpu

namespace fa { namespace cpu {

// Compact-directory sizing: occupancy bitmap and rank index.

inline uint64_t directory_align_up_u64(uint64_t value, uint64_t alignment) {
    if (alignment <= 1) return value;
    const uint64_t rem = value % alignment;
    return rem == 0 ? value : value + (alignment - rem);
}

inline uint64_t occupancy_storage_bytes(uint64_t bucket_count) {
    return directory_align_up_u64((bucket_count + 7ULL) / 8ULL, 8);
}

inline uint64_t rank_index_count(uint64_t bucket_count, uint32_t block_bits) {
    if (bucket_count == 0 || block_bits >= 63) return 0;
    const uint64_t block_size = 1ULL << block_bits;
    return ((bucket_count + block_size - 1ULL) >> block_bits) + 1ULL;
}

inline bool occupancy_bit(const uint8_t* bits, uint64_t slot) {
    return (bits[static_cast<size_t>(slot >> 3)] &
            static_cast<uint8_t>(1u << (slot & 7ULL))) != 0;
}

inline uint64_t popcount_bit_range(
    const uint8_t* bits,
    uint64_t begin_bit,
    uint64_t end_bit
) {
    if (!bits || end_bit <= begin_bit) return 0;
    uint64_t count = 0;
    uint64_t byte_idx = begin_bit >> 3;
    uint32_t bit_idx = static_cast<uint32_t>(begin_bit & 7ULL);
    uint64_t remaining = end_bit - begin_bit;

    if (bit_idx != 0) {
        const uint32_t take =
            static_cast<uint32_t>(std::min<uint64_t>(8ULL - bit_idx, remaining));
        const uint8_t mask = static_cast<uint8_t>(
            ((1u << take) - 1u) << bit_idx);
        count += static_cast<uint64_t>(
            __builtin_popcount(static_cast<unsigned>(bits[byte_idx] & mask)));
        ++byte_idx;
        remaining -= take;
    }

    // Eight bytes per popcount; the memcpy is an unaligned-safe load.
    while (remaining >= 64) {
        uint64_t word;
        std::memcpy(&word, bits + byte_idx, sizeof(word));
        count += static_cast<uint64_t>(__builtin_popcountll(word));
        byte_idx += 8;
        remaining -= 64;
    }

    while (remaining >= 8) {
        count += static_cast<uint64_t>(
            __builtin_popcount(static_cast<unsigned>(bits[byte_idx])));
        ++byte_idx;
        remaining -= 8;
    }

    if (remaining > 0) {
        const uint8_t mask =
            static_cast<uint8_t>((1u << static_cast<uint32_t>(remaining)) - 1u);
        count += static_cast<uint64_t>(
            __builtin_popcount(static_cast<unsigned>(bits[byte_idx] & mask)));
    }
    return count;
}

// Non-owning view of one directory: the pointers alias owned vectors (built index) or
// mapped file bytes (loaded index). tags[] and payloads[] are parallel, one entry per
// stored key in slot order. A tag is mixed_key >> shard_bits; the low bits are the shard
// id, so (shard, tag) identifies the key exactly.
struct FaixDirectory {
    const uint32_t* tags = nullptr;       // (key_tag << 1) | multi bit
    const uint64_t* payloads =
        nullptr; // single: packed (contig, local, z) | multi: count/begin
    const PackedRefPos* postings = nullptr; // multi positions only (no prefix)
    const uint64_t* chromosome_offsets = nullptr;
    uint64_t chromosome_count = 0;
    const uint8_t* occupancy_bits = nullptr;
    const uint64_t* rank_index = nullptr;
    uint64_t rank_entries = 0;
    uint32_t rank_block_bits = kFaixRankBlockBits;
    uint64_t bucket_count = 0;
    uint64_t unique_key_count = 0;
    uint64_t compact_posting_count = 0;
};

// Dense record index of the first occupied slot < `slot`: rank-sample at the
// slot's block + a popcount of the occupancy bitmap from the block start.
inline uint64_t directory_rank_before(const FaixDirectory& d, uint64_t slot) {
    if (!d.rank_index || d.rank_entries == 0) return 0;
    const uint64_t block = slot >> d.rank_block_bits;
    if (block + 1 >= d.rank_entries) return d.unique_key_count;
    const uint64_t block_begin = block << d.rank_block_bits;
    return d.rank_index[block] +
        popcount_bit_range(d.occupancy_bits, block_begin, slot);
}

// Decode one dense record into a posting view via the tag word's multi bit:
// inline singleton (packed position in payloads[record]), or a direct slice of
// the shard posting array (payload low/high uint32 words are begin/count).
inline KmerPostingView directory_bucket_view(
    const FaixDirectory& d,
    uint64_t record
) {
    if ((d.tags[static_cast<size_t>(record)] & kFaixTagMultiBit) == 0) {
        // singleton: the payload is the packed (contig, local, z) coordinate
        if (packed_ref_contig(d.payloads[static_cast<size_t>(record)]) >=
            d.chromosome_count) {
          return {};
        }
        return {RefPosSpan(&d.payloads[static_cast<size_t>(record)], 1,
                           d.chromosome_offsets, d.chromosome_count),
                1, 1};
    }
    const uint64_t payload = d.payloads[static_cast<size_t>(record)];
    const uint64_t begin = static_cast<uint32_t>(payload);
    const uint32_t count = static_cast<uint32_t>(payload >> 32);
    if (!d.postings || count < 2 || begin > d.compact_posting_count ||
        static_cast<uint64_t>(count) > d.compact_posting_count - begin) {
      return {};
    }
    // Postings are sorted, so the last carries the largest contig id: a constant-time
    // corruption guard.
    if (packed_ref_contig(d.postings[begin + count - 1]) >=
        d.chromosome_count) {
      return {};
    }
    return {RefPosSpan(d.postings + begin, count, d.chromosome_offsets,
                       d.chromosome_count),
            count, count};
}

// Probes one shard's open-addressed table for `tag` through the occupancy bitmap and rank
// directory; the probe hash is splitmix64(tag), as at build time. The dense record is
// carried across collisions: stepping past an occupied slot adds one to the rank, and only
// a wrap to slot 0 needs a fresh rank.
inline KmerPostingView directory_lookup(
    const FaixDirectory& d,
    uint32_t tag
) {
    if (!d.occupancy_bits || !d.rank_index || !d.tags) return {};
    const uint64_t mask = d.bucket_count - 1;
    uint64_t slot = detail::splitmix64(tag) & mask;
    uint64_t record = 0;
    bool have_record = false;
    for (uint64_t probes = 0; probes < d.bucket_count; ++probes) {
        if (!occupancy_bit(d.occupancy_bits, slot)) {
            return {};
        }
        if (!have_record) {
            record = directory_rank_before(d, slot);
            have_record = true;
        }
        if (record >= d.unique_key_count) {
            return {};
        }
        const uint32_t tagword = d.tags[static_cast<size_t>(record)];
        if ((tagword >> 1) == tag) {  // high 31 bits = key tag
            return directory_bucket_view(d, record);
        }
        const uint64_t next = (slot + 1) & mask;
        if (next < slot) {  // wrapped to slot 0: rank discontinuity -> recompute
            have_record = false;
        } else {  // advanced past an occupied slot: next rank = record + 1
            ++record;
        }
        slot = next;
    }
    return {};
}

// One entry per shard. occ_off is in bytes, the other offsets in elements. Payload begin
// fields are local to the shard's posting slice. bucket_count == 0 is an empty shard.
struct FaixShardEntry {
    uint64_t bucket_count = 0;   // open-addressed slots in this shard (power of two)
    uint64_t unique_keys = 0;    // dense compact buckets in this shard
    uint64_t posting_count = 0;  // compact (multi-hit) postings in this shard
    uint64_t occ_off = 0;        // byte offset into the occupancy blob
    uint64_t rank_off = 0;       // element offset into the rank blob
    uint64_t bucket_off = 0;     // element offset into the buckets blob
    uint64_t post_off = 0;       // element offset into the postings blob
};
static_assert(sizeof(FaixShardEntry) == 56, "FaixShardEntry must stay packed");

// Build the per-shard FaixDirectory view from the concatenated blob bases.
// tags[] and payloads[] share the same per-shard offset (e.bucket_off).
inline FaixDirectory
faix_shard_directory(const uint32_t* tags, const uint64_t* payloads,
                     const PackedRefPos* postings,
                     const uint64_t* chromosome_offsets,
                     uint64_t chromosome_count, const uint8_t* occupancy_bits,
                     const uint64_t* rank_index, const FaixShardEntry& e,
                     uint32_t rank_block_bits) {
  FaixDirectory d;
  d.tags = tags + e.bucket_off;
  d.payloads = payloads + e.bucket_off;
  d.postings = postings ? postings + e.post_off : nullptr;
  d.chromosome_offsets = chromosome_offsets;
  d.chromosome_count = chromosome_count;
  d.occupancy_bits = occupancy_bits + e.occ_off;
  d.rank_index = rank_index + e.rank_off;
  d.rank_entries = rank_index_count(e.bucket_count, rank_block_bits);
  d.rank_block_bits = rank_block_bits;
  d.bucket_count = e.bucket_count;
  d.unique_key_count = e.unique_keys;
  d.compact_posting_count = e.posting_count;
  return d;
}

// Two-level lookup: shard from the low mixed bits, then probe that shard's
// compact directory by the high-bit tag. Returns the same KmerPostingView
// contract as the global directory_lookup.
inline KmerPostingView faix_sharded_lookup(
    const uint32_t* tags, const uint64_t* payloads,
    const PackedRefPos* postings, const uint64_t* chromosome_offsets,
    uint64_t chromosome_count, const uint8_t* occupancy_bits,
    const uint64_t* rank_index, const FaixShardEntry* shard_dir,
    uint64_t shard_mask, int shard_bits, int twok, uint32_t rank_block_bits,
    uint64_t key) {
  if (!shard_dir || !tags)
    return {};
  const uint64_t mixed = detail::shard_mix_2k(key, twok);
  const FaixShardEntry& e = shard_dir[static_cast<size_t>(mixed & shard_mask)];
  if (e.bucket_count == 0)
    return {};
  const FaixDirectory d = faix_shard_directory(
      tags, payloads, postings, chromosome_offsets, chromosome_count,
      occupancy_bits, rank_index, e, rank_block_bits);
  return directory_lookup(d, static_cast<uint32_t>(mixed >> shard_bits));
}

// Batched faix_sharded_lookup. One lookup is a chain of dependent loads (occupancy, rank,
// tag, payload, posting tail), so the batch runs the same probe in stages and prefetches
// each stage's loads for several keys before any is consumed. Every key gets the same view
// as the scalar path.

// Keys per pipeline stage step (see faix_sharded_lookup_batch): a lane's
// prefetches are consumed three steps later, and four lanes are in flight.
static constexpr size_t kFaixPipelineLane = 8;

namespace detail {

// One in-flight scalar probe, frozen between stages. `probe` is the index of
// the scalar loop iteration currently being executed, so the loop bound
// `probes < d.bucket_count` can be reproduced exactly across the stage split.
struct FaixLookupProbeState {
    FaixDirectory d;
    uint64_t mask = 0;
    uint64_t slot = 0;
    uint64_t record = 0;
    uint64_t probe = 0;
    uint32_t tag = 0;
    bool have_record = false;
    bool live = false;     // still probing: stages 1/2 must run it
    bool matched = false;  // stage 2 found the tag; stage 3 decodes the record
};

}  // namespace detail

inline void faix_sharded_lookup_batch(
    const uint32_t* tags, const uint64_t* payloads,
    const PackedRefPos* postings, const uint64_t* chromosome_offsets,
    uint64_t chromosome_count, const uint8_t* occupancy_bits,
    const uint64_t* rank_index, const FaixShardEntry* shard_dir,
    uint64_t shard_mask, int shard_bits, int twok, uint32_t rank_block_bits,
    const uint64_t* keys, size_t key_count, KmerPostingView* out) {
  if (!keys || !out || key_count == 0)
    return;
  for (size_t i = 0; i < key_count; ++i)
    out[i] = KmerPostingView{};
  // Same refusal as the scalar's `if (!shard_dir || !tags)`.
  if (!shard_dir || !tags)
    return;

  // The scalar probe split at its three dependent loads into four stages. Keys move through
  // them in lanes of kFaixPipelineLane: each step runs stage 0 on lane t, stage 1 on t-1,
  // stage 2 on t-2 and stage 3 on t-3, so other lanes' work overlaps each prefetch.
  auto stage0 = [&](detail::FaixLookupProbeState* st, size_t base, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      detail::FaixLookupProbeState& s = st[i];
      s.live = false;
      s.matched = false;
      const uint64_t mixed = detail::shard_mix_2k(keys[base + i], twok);
      const FaixShardEntry& e =
          shard_dir[static_cast<size_t>(mixed & shard_mask)];
      if (e.bucket_count == 0)
        continue;
      s.d = faix_shard_directory(tags, payloads, postings, chromosome_offsets,
                                 chromosome_count, occupancy_bits, rank_index,
                                 e, rank_block_bits);
      // Mirrors directory_lookup()'s first line.
      if (!s.d.occupancy_bits || !s.d.rank_index || !s.d.tags)
        continue;
      s.tag = static_cast<uint32_t>(mixed >> shard_bits);
      s.mask = s.d.bucket_count - 1;
      s.slot = detail::splitmix64(s.tag) & s.mask;
      s.record = 0;
      s.probe = 0;
      s.have_record = false;
      s.live = true;
      const uint64_t block = s.slot >> s.d.rank_block_bits;
      __builtin_prefetch(&s.d.occupancy_bits[s.slot >> 3], 0, 3);
      __builtin_prefetch(
          &s.d.occupancy_bits[(block << s.d.rank_block_bits) >> 3], 0, 3);
      if (s.d.rank_entries != 0 && block + 1 < s.d.rank_entries)
        __builtin_prefetch(&s.d.rank_index[block], 0, 3);
    }
  };
  auto stage1 = [&](detail::FaixLookupProbeState* st, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      detail::FaixLookupProbeState& s = st[i];
      if (!s.live)
        continue;
      if (!occupancy_bit(s.d.occupancy_bits, s.slot)) {
        s.live = false;
        continue;
      }
      s.record = directory_rank_before(s.d, s.slot);
      s.have_record = true;
      if (s.record >= s.d.unique_key_count) {
        s.live = false;
        continue;
      }
      __builtin_prefetch(&s.d.tags[static_cast<size_t>(s.record)], 0, 3);
      __builtin_prefetch(&s.d.payloads[static_cast<size_t>(s.record)], 0, 3);
    }
  };
  auto stage2 = [&](detail::FaixLookupProbeState* st, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      detail::FaixLookupProbeState& s = st[i];
      if (!s.live)
        continue;
      s.live = false;
      for (;;) {
        const uint32_t tagword = s.d.tags[static_cast<size_t>(s.record)];
        if ((tagword >> 1) == s.tag) {  // high 31 bits = key tag
          s.matched = true;
          if ((tagword & kFaixTagMultiBit) != 0) {
            // Multi-hit: prefetch the posting tail directory_bucket_view reads in
            // stage 3, within the bounds it checks.
            const uint64_t payload =
                s.d.payloads[static_cast<size_t>(s.record)];
            const uint64_t begin = static_cast<uint32_t>(payload);
            const uint64_t count = static_cast<uint32_t>(payload >> 32);
            if (s.d.postings && count >= 2 &&
                begin <= s.d.compact_posting_count &&
                count <= s.d.compact_posting_count - begin) {
              __builtin_prefetch(&s.d.postings[begin + count - 1], 0, 3);
            }
          }
          break;
        }
        const uint64_t next = (s.slot + 1) & s.mask;
        if (next < s.slot) {  // wrapped to slot 0: rank discontinuity
          s.have_record = false;
        } else {
          ++s.record;
        }
        s.slot = next;
        ++s.probe;
        if (s.probe >= s.d.bucket_count) {  // scalar loop bound exhausted
          break;
        }
        if (!occupancy_bit(s.d.occupancy_bits, s.slot)) {
          break;
        }
        if (!s.have_record) {
          s.record = directory_rank_before(s.d, s.slot);
          s.have_record = true;
        }
        if (s.record >= s.d.unique_key_count) {
          break;
        }
      }
    }
  };
  auto stage3 = [&](detail::FaixLookupProbeState* st, size_t base, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      detail::FaixLookupProbeState& s = st[i];
      if (!s.matched)
        continue;
      out[base + i] = directory_bucket_view(s.d, s.record);
    }
  };

  detail::FaixLookupProbeState st[4][kFaixPipelineLane];
  const size_t lanes = (key_count + kFaixPipelineLane - 1) / kFaixPipelineLane;
  const auto lane_base = [](size_t lane) { return lane * kFaixPipelineLane; };
  const auto lane_size = [&](size_t lane) {
    return std::min<size_t>(kFaixPipelineLane, key_count - lane_base(lane));
  };
  for (size_t t = 0; t < lanes + 3; ++t) {
    if (t < lanes)
      stage0(st[t & 3], lane_base(t), lane_size(t));
    if (t >= 1 && t - 1 < lanes)
      stage1(st[(t - 1) & 3], lane_size(t - 1));
    if (t >= 2 && t - 2 < lanes)
      stage2(st[(t - 2) & 3], lane_size(t - 2));
    if (t >= 3 && t - 3 < lanes)
      stage3(st[(t - 3) & 3], lane_base(t - 3), lane_size(t - 3));
  }
}
}}  // namespace fa::cpu
