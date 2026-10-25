// Nucleotide encoding, AlignResult, the CIGAR, sequence, hash and checked-range helpers every directory uses, the flat int64 map and the radix sorts.
#pragma once

#include <flashalign/alignment.hpp>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

/**** types ****/
// Nucleotide encoding and shared alignment types.

namespace fa { namespace cpu {

// A C G T/U -> 0 1 2 3, anything else -> 4.
inline int32_t nuc_encode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return 4;
    }
}

inline std::vector<uint8_t> encode_sequence_u8(const std::string& seq) {
    std::vector<uint8_t> enc(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        enc[i] = static_cast<uint8_t>(nuc_encode(seq[i]));
    }
    return enc;
}

using AlignmentOrigin = ::flashalign::AlignmentOrigin;
using AlignResult = ::flashalign::Alignment;

// The heap bytes of vectors at capacity.
template <class... Vectors>
std::size_t capacity_bytes(const Vectors&... vectors) {
    return (std::size_t{0} + ... +
            (vectors.capacity() * sizeof(typename Vectors::value_type)));
}

}}  // namespace fa::cpu

/**** cigar ****/
// CIGAR building, parsing, clipping and replay, shared by the DP and output code.

namespace fa { namespace cpu { namespace output {

// Appends `count` of `op`, merged into the last run when the operator matches; a count <= 0
// writes nothing.
void append_cigar_run(std::string& cigar, int count, char op);

// Appends every run of `suffix` through append_cigar_run.
void append_cigar_string(std::string& cigar, const std::string& suffix);

// Packed form: one uint32 per run, `count << 4 | op` with op indexing "MIDNSHP=XB", as in BAM
// and the ksw2 traceback. The packed helpers below produce the same runs as their text
// counterparts, so the DNA realizer can build CIGARs as runs and render text once.
using PackedCigar = std::vector<std::uint32_t>;

constexpr char kPackedCigarOps[] = "MIDNSHP=XB";

constexpr int packed_run_count(std::uint32_t run) {
    return static_cast<int>(run >> 4);
}
constexpr int packed_run_op(std::uint32_t run) {
    return static_cast<int>(run & 0xfu);
}

// Packed append_cigar_string.
void append_packed_cigar(PackedCigar& cigar, const PackedCigar& suffix);

// Puts `prefix` in front of `cigar`, as `append_packed_cigar(prefix, cigar)` would, but in
// time linear in `prefix`. Assumes `cigar` has no adjacent same-op runs, so only its leading
// run can merge; the rest is copied as is.
void prepend_packed_cigar(PackedCigar& cigar, PackedCigar prefix);

// Query bases consumed (M/I/S/=/X), or -1 when a run has a zero count.
int packed_cigar_query_consumed(const PackedCigar& cigar);

// Text rendering of a packed CIGAR.
std::string packed_cigar_text(const PackedCigar& cigar);

// (count, op) pairs; empty when the string is malformed or has a zero count.
std::vector<std::pair<int, char>> parse_cigar_ops(const std::string& cigar);

struct TerminalHardClip {
    std::string cigar;  // a copy with leading/trailing terminal 'S' rewritten to 'H'
    int lead = 0;       // length of the leading terminal clip op (0 if none)
    int trail = 0;      // length of the trailing terminal clip op (0 if none)
};

// Rewrites the terminal soft clips of a supplementary record's CIGAR to hard clips, as
// minimap2 does, and reports their lengths so the caller can trim SEQ/QUAL to match.
TerminalHardClip hard_clip_terminal(const std::string& cigar);

// q_lo/q_hi are on the forward read, while a CIGAR is in reference orientation: on the reverse
// strand the leading clip covers the read's 3' end. The helpers below do this conversion.
struct ForwardQueryInterval {
    int q_lo = 0;
    int q_hi = 0;
    bool valid = false;
};

struct ReferenceOrientedTerminalClips {
    int lead = 0;
    int trail = 0;
};

constexpr int clamp_query_coord(int x, int read_len) {
    return x < 0 ? 0 : (x > read_len ? read_len : x);
}

constexpr ReferenceOrientedTerminalClips
reference_oriented_terminal_clips(int read_len, int forward_q_lo,
                                  int forward_q_hi, bool is_reverse) {
    const int rl = read_len < 0 ? 0 : read_len;
    const int q0 = clamp_query_coord(forward_q_lo, rl);
    const int q1raw = clamp_query_coord(forward_q_hi, rl);
    const int q1 = q1raw < q0 ? q0 : q1raw;
    return is_reverse
        ? ReferenceOrientedTerminalClips{rl - q1, q0}
        : ReferenceOrientedTerminalClips{q0, rl - q1};
}

// Forward-read interval covered by a reference-oriented CIGAR; S and H clips both count.
ForwardQueryInterval forward_query_interval_from_cigar(
    const std::string& cigar, int read_len, bool is_reverse);

// Replaces the terminal clips of a reference-oriented CIGAR with the soft clips for the
// forward interval [forward_q_lo, forward_q_hi). Returns an empty string when the body does
// not consume exactly that interval.
std::string clipped_reference_oriented_cigar(
    const std::string& reference_cigar, int read_len, int forward_q_lo,
    int forward_q_hi, bool is_reverse);

static_assert(reference_oriented_terminal_clips(1000, 100, 600, false).lead == 100,
              "forward segment leading clip");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, false).trail == 400,
              "forward segment trailing clip");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, true).lead == 400,
              "reverse segment leading clip is mirrored");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, true).trail == 100,
              "reverse segment trailing clip is mirrored");

// Result of replaying a CIGAR against the query (reference orientation) and one reference
// sequence. `valid` is set only when every operation stays in bounds, the clips and body
// cover the whole query, and there is no interior clip.
//
// The counts follow minimap2's mm_update_extra: ambiguous (N) M/I/D residues are left out of
// block_len and matches and counted once in ambiguities, and
// edit_distance = block_len - matches + ambiguities. Consumed lengths include them.
struct CigarReplay {
    bool valid = false;
    int query_consumed = 0;         // M/I/S/H/=/X in full-query coordinates
    int aligned_query_consumed = 0; // M/I/=/X only
    int target_consumed = 0;        // M/D/N/=/X
    int target_start = -1;
    int target_end = -1;
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    int ambiguities = 0;
    int edit_distance = 0;
    int block_len = 0;
    int lead_clip = 0;
    int trail_clip = 0;
    // cs and MD strings without their tag prefixes; empty unless requested.
    std::string cs;
    std::string md;
    // The CIGAR with M/=/X split into =/X runs (--eqx); empty unless requested.
    std::string eqx_cigar;
};

// Optional outputs of replay_cigar; the default asks for none. They follow minimap2's
// write_cs_ds_core, write_MD_core and mm_update_cigar_eqx and never change the counts.
//
//   cs, over the aligned core (clips excluded): ":<len>" (Short) or "=<BASES>" (Long) for an
//   identical run, "*<ref><qry>" for a mismatch, "+<bases>" for I, "-<bases>" for D and
//   "~<t0><t1><len><t[len-2]><t[len-1]>" for N, on the forward reference. An identical run
//   ends at every M/=/X op boundary, so "5M5M" gives ":5:5".
//   MD: running match count; a mismatch writes "<count><REF>", a D writes "<count>^<REFS>",
//   and I and N neither write nor reset the count. Unlike minimap2, the final count is always
//   written, even "0", as the SAM MD grammar requires.
//   =/X: each M/=/X op becomes alternating = and X runs; other ops are copied.
//
// cs, MD and =/X compare base codes, so N against N is a match there although the counts
// treat it as an ambiguity.
struct CigarReplayRequest {
  enum class Cs { None, Short, Long };
  Cs cs = Cs::None;
  bool md = false;
  bool eqx = false;

  bool any() const { return cs != Cs::None || md || eqx; }
};

CigarReplay replay_cigar(const std::string& cigar, const uint8_t* qseq,
                         int qlen, const uint8_t* tseq, int tlen, int ref_start,
                         const CigarReplayRequest& request);

// Counts only.
CigarReplay replay_cigar(const std::string& cigar,
                                const uint8_t* qseq, int qlen,
                                const uint8_t* tseq, int tlen,
                                int ref_start);

}}} // namespace fa::cpu::output

/**** sequence ****/

namespace fa {
namespace cpu {
namespace lr {

// Reverse complement of a 0-4 encoded sequence (A C G T N); N stays N.
std::vector<uint8_t>
reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc);

} // namespace lr
} // namespace cpu
} // namespace fa

/**** hash ****/
// Integer hashes shared by the index and seeding.

namespace fa { namespace cpu {

// Seed-ordering hash for syncmer sketches.
inline uint64_t sketch_order_hash(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

namespace detail {

// Open-addressing probe hash for the k-mer table and the seeding lookup caches.
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Bijection on the low `twok` bits of a 2-bit-packed k-mer (twok = 2k): xorshifts and odd
// multipliers mod 2^twok. The sharded index takes the shard from the low bits of the mix and
// stores the rest as the tag, so (shard, tag) encodes the key losslessly. Mixing folds high
// bits into the shard bits, which keeps shards balanced despite GC-skewed trailing bases.
inline uint64_t shard_mix_2k(uint64_t key, int twok) {
    const uint64_t mask = (twok >= 64) ? ~0ULL : ((1ULL << twok) - 1ULL);
    const int s1 = (twok >= 3) ? twok / 3 : 1;
    const int s2 = (twok >= 2) ? twok / 2 : 1;
    uint64_t x = key & mask;
    x = (x ^ (x >> s1)) & mask;
    x = (x * 0x9E3779B97F4A7C15ULL) & mask;  // odd: invertible mod 2^twok
    x = (x ^ (x >> s2)) & mask;
    x = (x * 0xD6E8FEB86659FD93ULL) & mask;  // odd
    x = (x ^ (x >> s1)) & mask;
    return x;
}

}  // namespace detail

}}  // namespace fa::cpu

/**** checked_range ****/
// Checked conversions for `int` sequence coordinates and counts. The supported domain is
// [0, INT_MAX], and a half-open end may equal INT_MAX. Check values while still 64-bit;
// narrowing first and validating afterwards is wrong.

namespace fa {
namespace cpu {
namespace mapping {

constexpr std::int64_t kSignedCoordinateMinimum = 0;
constexpr std::int64_t kSignedCoordinateMaximum =
    static_cast<std::int64_t>(std::numeric_limits<int>::max());

enum class CoordinateDomainRefusal : std::uint8_t {
  None = 0,
  UnrepresentableCount,
  NegativeCoordinate,
  InvalidRange,
  EndpointOverflow,
  ProductOverflow,
  ReferenceShapeMismatch,
};

inline bool checked_size_to_int(std::size_t value, int& converted) noexcept {
  if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return false;
  converted = static_cast<int>(value);
  return true;
}

inline bool checked_nonnegative_to_int(std::int64_t value,
                                       int& converted) noexcept {
  if (value < kSignedCoordinateMinimum ||
      value > kSignedCoordinateMaximum)
    return false;
  converted = static_cast<int>(value);
  return true;
}

inline bool checked_half_open_range(std::int64_t begin, std::int64_t end,
                                    std::uint64_t length) noexcept {
  return begin >= 0 && end >= begin &&
      static_cast<std::uint64_t>(end) <= length;
}

inline bool checked_count_range(std::int64_t begin, std::uint64_t count,
                                std::uint64_t size) noexcept {
  if (begin < 0) return false;
  const std::uint64_t unsigned_begin = static_cast<std::uint64_t>(begin);
  return unsigned_begin <= size && count <= size - unsigned_begin;
}

inline bool checked_add(std::int64_t left, std::int64_t right,
                        std::int64_t& sum) noexcept {
  const std::int64_t maximum = std::numeric_limits<std::int64_t>::max();
  const std::int64_t minimum = std::numeric_limits<std::int64_t>::min();
  if ((right > 0 && left > maximum - right) ||
      (right < 0 && left < minimum - right))
    return false;
  sum = left + right;
  return true;
}

inline bool checked_signed_endpoint(std::int64_t begin, std::int64_t span,
                                    std::int64_t length,
                                    int& endpoint) noexcept {
  if (begin < kSignedCoordinateMinimum || span < 0 ||
      length < kSignedCoordinateMinimum ||
      begin > kSignedCoordinateMaximum ||
      length > kSignedCoordinateMaximum)
    return false;
  if (span > kSignedCoordinateMaximum - begin) return false;
  const std::int64_t widened_endpoint = begin + span;
  if (widened_endpoint > length) return false;
  endpoint = static_cast<int>(widened_endpoint);
  return true;
}

inline bool checked_u64_product(std::uint64_t left, std::uint64_t right,
                                std::uint64_t& product) noexcept {
  const std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
  if (left != 0 && right > maximum / left) return false;
  product = left * right;
  return true;
}

}  // namespace mapping
}  // namespace cpu
}  // namespace fa

namespace fa::cpu::lr {

// C's (int) conversion, truncating toward zero, saturated so an out-of-range
// value is not undefined behaviour.
inline int truncate_to_int(double value) noexcept {
  if (!(value == value))
    return 0; // NaN
  if (value >= 2147483647.0)
    return 2147483647;
  if (value <= -2147483648.0)
    return -2147483648;
  return static_cast<int>(value);
}

}  // namespace fa::cpu::lr

/**** flat_int64_map ****/
// Open-addressed int64 -> Value map for the diagonal-bin vote, covering the subset of
// std::unordered_map it uses. Entries sit in a dense array in insertion order; the index
// holds generation-tagged entry numbers, so clear() is O(1). Iteration is in insertion order.

namespace fa {
namespace cpu {
namespace lr {

template <class Value> class FlatInt64Map {
public:
  struct value_type {
    int64_t first;
    Value second;
  };
  using iterator = value_type *;
  using const_iterator = const value_type *;

  FlatInt64Map() = default;

  iterator begin() { return entries_.data(); }
  iterator end() { return entries_.data() + entries_.size(); }
  const_iterator begin() const { return entries_.data(); }
  const_iterator end() const { return entries_.data() + entries_.size(); }

  size_t size() const { return entries_.size(); }
  bool empty() const { return entries_.empty(); }
  // Index capacity, as unordered_map::bucket_count.
  size_t bucket_count() const { return slots_.size(); }

  // Bumps the generation instead of clearing the index, except on wrap-around.
  void clear() {
    entries_.clear();
    if (++generation_ == 0) {
      std::fill(slots_.begin(), slots_.end(), uint64_t{0});
      generation_ = 1;
    }
  }

  // Room for n entries at load <= 0.5 without a rehash.
  void reserve(size_t n) {
    const size_t want = pow2_at_least(n == 0 ? 1 : n * 2);
    if (want > slots_.size())
      rebuild_index(want);
    entries_.reserve(n);
  }

  // Get-or-insert a value-initialized Value. The reference is valid until the next insert.
  Value &operator[](int64_t key) {
    ensure_index();
    size_t mask = slots_.size() - 1;
    size_t slot = hash(key) & mask;
    for (;;) {
      const uint64_t s = slots_[slot];
      if (!slot_live(s))
        break; // empty slot -> insert below
      if (entries_[slot_index(s)].first == key)
        return entries_[slot_index(s)].second;
      slot = (slot + 1) & mask;
    }
    if ((entries_.size() + 1) * 2 > slots_.size()) {
      rebuild_index(slots_.size() * 2);
      mask = slots_.size() - 1;
      slot = hash(key) & mask;
      while (slot_live(slots_[slot]))
        slot = (slot + 1) & mask;
    }
    const uint32_t new_idx = static_cast<uint32_t>(entries_.size());
    entries_.push_back(value_type{key, Value{}});
    slots_[slot] = pack(generation_, new_idx);
    return entries_.back().second;
  }

  iterator find(int64_t key) {
    const int64_t i = find_index(key);
    return i < 0 ? end() : entries_.data() + i;
  }
  const_iterator find(int64_t key) const {
    const int64_t i = find_index(key);
    return i < 0 ? end() : entries_.data() + i;
  }

  void swap(FlatInt64Map &o) {
    entries_.swap(o.entries_);
    slots_.swap(o.slots_);
    std::swap(generation_, o.generation_);
  }

  // The heap bytes of the entries and the index.
  size_t heap_bytes() const { return capacity_bytes(entries_, slots_); }

private:
  static size_t pow2_at_least(size_t x) {
    size_t c = 1;
    while (c < x)
      c <<= 1;
    return c;
  }
  // Keys are (chr << 32) | bin. The default hash folds chr into the low bits, so contiguous
  // bins land in consecutive slots while equal bins on different chromosomes do not all
  // collide. FA_DNA_LONG_VOTE_HASH selects splitmix64 (0) or the identity (2) instead.
#ifndef FA_DNA_LONG_VOTE_HASH
#define FA_DNA_LONG_VOTE_HASH 1
#endif
  static uint64_t hash(int64_t key) {
#if FA_DNA_LONG_VOTE_HASH == 0
    return fa::cpu::detail::splitmix64(static_cast<uint64_t>(key));
#elif FA_DNA_LONG_VOTE_HASH == 2
    return static_cast<uint64_t>(key);
#else
    const uint64_t k = static_cast<uint64_t>(key);
    return k ^ (k >> 32);
#endif
  }
  static uint64_t pack(uint32_t gen, uint32_t idx) {
    return (static_cast<uint64_t>(gen) << 32) | idx;
  }
  bool slot_live(uint64_t s) const {
    return static_cast<uint32_t>(s >> 32) == generation_;
  }
  static uint32_t slot_index(uint64_t s) {
    return static_cast<uint32_t>(s & 0xFFFFFFFFu);
  }
  void ensure_index() {
    if (slots_.empty())
      rebuild_index(16);
  }

  int64_t find_index(int64_t key) const {
    if (slots_.empty())
      return -1;
    const size_t mask = slots_.size() - 1;
    size_t slot = hash(key) & mask;
    for (;;) {
      const uint64_t s = slots_[slot];
      if (!slot_live(s))
        return -1;
      if (entries_[slot_index(s)].first == key)
        return static_cast<int64_t>(slot_index(s));
      slot = (slot + 1) & mask;
    }
  }

  // Rebuilds the index at new_cap (a power of two, at least 16) from the live entries.
  void rebuild_index(size_t new_cap) {
    new_cap = pow2_at_least(new_cap < 16 ? 16 : new_cap);
    slots_.assign(new_cap, uint64_t{0});
    if (generation_ == 0)
      generation_ = 1; // 0 marks an empty slot
    const size_t mask = new_cap - 1;
    for (uint32_t i = 0; i < entries_.size(); ++i) {
      size_t slot = hash(entries_[i].first) & mask;
      while (slot_live(slots_[slot]))
        slot = (slot + 1) & mask;
      slots_[slot] = pack(generation_, i);
    }
  }

  std::vector<value_type> entries_; // dense, insertion order
  std::vector<uint64_t> slots_;     // gen<<32 | entry index; live iff gen match
  uint32_t generation_ = 1;
};

template <class Value>
inline void swap(FlatInt64Map<Value> &a, FlatInt64Map<Value> &b) {
  a.swap(b);
}

} // namespace lr
} // namespace cpu
} // namespace fa

/**** radix_sort ****/
// LSD radix sorts with the same output as std::stable_sort on the same key, without its
// temporary buffer or comparator calls. Sort composite keys least significant key first,
// one call per key.

namespace fa {
namespace cpu {
namespace radix {

// The byte order that is the signed order.
inline std::uint32_t key_u32(std::int32_t value) {
  return static_cast<std::uint32_t>(value) ^ 0x80000000u;
}
inline std::uint32_t key_u32(std::uint32_t value) { return value; }

// Caller-owned scratch (ping-pong buffer and the four 256-bucket histograms), reused
// across calls.
template <class T>
struct Scratch {
  std::vector<T> tmp;
  std::vector<std::uint32_t> count;
};

// Stable ascending sort of `items` by the 32-bit key `key(item)` (int32 or uint32). A byte
// position shared by every key is skipped.
template <class T, class Key>
void stable_sort_u32(std::vector<T>& items, Scratch<T>& scratch, Key key) {
  const std::size_t n = items.size();
  if (n < 2)
    return;
  std::vector<T>& tmp = scratch.tmp;
  std::vector<std::uint32_t>& count = scratch.count;
  tmp.resize(n);
  count.assign(4 * 256, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::uint32_t k = key_u32(key(items[i]));
    ++count[k & 255u];
    ++count[256 + ((k >> 8) & 255u)];
    ++count[512 + ((k >> 16) & 255u)];
    ++count[768 + (k >> 24)];
  }
  std::vector<T>* src = &items;
  std::vector<T>* dst = &tmp;
  for (int pass = 0; pass < 4; ++pass) {
    std::uint32_t* bucket = count.data() + 256 * pass;
    bool trivial = false;
    for (int b = 0; b < 256; ++b) {
      if (bucket[b] == n) {
        trivial = true;
        break;
      }
    }
    if (trivial)
      continue;
    std::uint32_t sum = 0;
    for (int b = 0; b < 256; ++b) {
      const std::uint32_t c = bucket[b];
      bucket[b] = sum;
      sum += c;
    }
    const int shift = 8 * pass;
    const std::vector<T>& from = *src;
    std::vector<T>& to = *dst;
    for (std::size_t i = 0; i < n; ++i) {
      const std::uint32_t k = key_u32(key(from[i]));
      to[bucket[(k >> shift) & 255u]++] = from[i];
    }
    std::swap(src, dst);
  }
  if (src != &items)
    items.swap(tmp);
}

} // namespace radix
} // namespace cpu
} // namespace fa
