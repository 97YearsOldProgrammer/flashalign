// k-mer encoding and closed-syncmer seed extraction, shared by the index builder
// (reference seeds) and the per-read query path (query seeds).
#pragma once

#include "../core/core.h"   // sketch_order_hash, detail::splitmix64, nuc_encode
#include "format.h"         // RefPos

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

namespace fa { namespace cpu {

inline bool htable_valid_k(int k) {
    return k > 0 && k <= 31;
}

// No k-mer key takes this value (keys have at most 62 bits).
static constexpr uint64_t kEmptyKmerKey = std::numeric_limits<uint64_t>::max();

inline uint64_t htable_kmer_mask(int k) {
    return (k <= 0) ? 0ULL : ((1ULL << (2 * k)) - 1ULL);
}

// The key is min(forward, reverse complement), as the index stores it. `out_is_rc` reports
// the orientation bit z: true iff the key is the reverse complement of the forward k-mer.
// `out_is_pal` reports a palindrome (key == its reverse complement, even k only), which is
// never a seed.
template <class BaseT>
inline bool encode_kmer_key(
    const BaseT* enc,
    int k,
    uint64_t* out_key,
    bool* out_is_rc = nullptr,
    bool* out_is_pal = nullptr
) {
    if (out_is_rc) *out_is_rc = false;
    if (out_is_pal) *out_is_pal = false;
    if (!out_key || !enc || !htable_valid_k(k)) return false;
    uint64_t fwd = 0;
    uint64_t rev = 0;
    for (int i = 0; i < k; i++) {
        const int base = static_cast<int>(enc[i]);
        if (base < 0 || base > 3) return false;
        fwd = (fwd << 2) | static_cast<uint64_t>(base);
        rev |= static_cast<uint64_t>(3 - base) << (2 * i);
    }
    *out_key = std::min(fwd, rev);
    if (out_is_rc) *out_is_rc = rev < fwd;
    if (out_is_pal) *out_is_pal = fwd == rev;
    return true;
}

class RollingKmerEncoder {
public:
    explicit RollingKmerEncoder(int k = 0)
        : k_(k)
        , mask_(htable_valid_k(k) ? htable_kmer_mask(k) : 0ULL)
    {}

    void reset() {
        forward_ = 0;
        reverse_ = 0;
        valid_bases_ = 0;
    }

    bool push_base(int32_t base) {
        if (!htable_valid_k(k_)) return false;
        if (base < 0 || base > 3) {
            reset();
            return false;
        }
        return push_base_unchecked(base);
    }

    bool push_base_unchecked(int32_t base) {
        forward_ = ((forward_ << 2) | static_cast<uint64_t>(base)) & mask_;
        reverse_ = (reverse_ >> 2)
            | (static_cast<uint64_t>(3 - base) << (2 * (k_ - 1)));
        if (valid_bases_ < k_) valid_bases_++;
        return valid_bases_ >= k_;
    }

    uint64_t key() const {
        return std::min(forward_, reverse_);
    }

    // The orientation bit z: true iff key() is the reverse complement of the forward k-mer.
    // False for a palindrome, matching min()'s choice.
    bool key_is_rc() const { return reverse_ < forward_; }

    // True iff the key is its own reverse complement (even k only); never a seed.
    bool key_is_palindrome() const {
        return forward_ == reverse_;
    }

    int valid_bases() const { return valid_bases_; }

    // Canonical rightmost s-mer of the window, from the k-state: the low 2s bits of forward_
    // are the last s bases and the high 2s bits of reverse_ their reverse complement. Valid
    // once valid_bases_ >= s; equals RollingKmerEncoder(s).key().
    uint64_t suffix_smer_key(int s) const {
        const uint64_t s_mask = (1ULL << (2 * s)) - 1ULL;
        const uint64_t f = forward_ & s_mask;
        const uint64_t r = reverse_ >> (2 * (k_ - s));
        return std::min(f, r);
    }

private:
    int k_ = 0;
    uint64_t mask_ = 0;
    uint64_t forward_ = 0;
    uint64_t reverse_ = 0;
    int valid_bases_ = 0;
};

}}  // namespace fa::cpu

namespace fa { namespace cpu {

// z is the forward occurrence's orientation bit: 1 iff the forward k-mer at read_pos is the
// reverse complement of `key`. Reverse-complement projection copies it unchanged; the
// lane belongs to the consumer. z fits in the struct's padding.
struct QuerySeed {
    uint64_t key = 0;
    int read_pos = 0;
    uint8_t z = 0;
};
static_assert(sizeof(QuerySeed) == 16,
              "QuerySeed must stay 16 bytes: z lives in the existing padding");

// As minimap2's mm_seed_collect_all: a query seed is tandem when an adjacent seed in query
// order has the same key (MM_SEED_TANDEM).
inline bool query_seed_has_tandem_neighbor(
    const std::vector<QuerySeed>& seeds,
    size_t index
) {
    if (index >= seeds.size()) return false;
    return (index > 0 && seeds[index - 1].key == seeds[index].key)
        || (index + 1 < seeds.size()
            && seeds[index + 1].key == seeds[index].key);
}

}}  // namespace fa::cpu

namespace fa { namespace cpu {

struct ClosedSyncmerConfig {
    int k = 17;
    int s = 0;              // closed-syncmer s-mer length; set explicitly (s>0)
    int downsample = 1;
};

// Clamps s to [1, k]; the CLI already refuses an out-of-range s.
inline int effective_closed_syncmer_s(int k, int explicit_s) {
    if (!htable_valid_k(k)) return 0;
    return std::max(1, std::min(k, explicit_s));
}

// A k-mer is a closed syncmer iff its minimum-order s-mer sits at its left or right end.
// Used only for windows k - s + 1 > 64, which k <= 31 never reaches.
template <class BaseT>
inline bool closed_syncmer_selected_generic(const BaseT* kmer, int k, int s) {
    if (!kmer || !htable_valid_k(k) || !htable_valid_k(s) || s > k) {
        return false;
    }
    uint64_t best_order = std::numeric_limits<uint64_t>::max();
    bool best_at_left = false;
    bool best_at_right = false;
    const int last = k - s;
    for (int offset = 0; offset <= last; ++offset) {
        uint64_t subkey = 0;
        if (!encode_kmer_key(kmer + offset, s, &subkey)) {
            return false;
        }
        const uint64_t order = sketch_order_hash(subkey);
        if (order < best_order) {
            best_order = order;
            best_at_left = offset == 0;
            best_at_right = offset == last;
        } else if (order == best_order) {
            best_at_left = best_at_left || offset == 0;
            best_at_right = best_at_right || offset == last;
        }
    }
    return best_at_left || best_at_right;
}

// Order-rank table for the s-mer ordering. sketch_order_hash (Murmur3 fmix64) is a
// bijection on uint64, so order_rank[x] = the rank of sketch_order_hash(x) among all 4^s
// s-mers preserves every comparison, and selection is unchanged while the per-base hash
// becomes one table load.

// Built once per S, thread-safely; S <= 16 keeps every rank in uint32_t.
template <int S>
inline const uint32_t* closed_syncmer_order_rank_table_impl() {
    static_assert(S >= 1 && S <= 16, "order-rank table requires 1 <= S <= 16");
    static const std::vector<uint32_t> table = [] {
        const size_t n = static_cast<size_t>(1) << (2 * S);
        std::vector<uint32_t> idx(n);
        for (size_t i = 0; i < n; ++i) idx[i] = static_cast<uint32_t>(i);
        std::sort(idx.begin(), idx.end(), [](uint32_t a, uint32_t b) {
            return sketch_order_hash(a) < sketch_order_hash(b);
        });
        std::vector<uint32_t> rank(n);
        for (size_t r = 0; r < n; ++r) rank[idx[r]] = static_cast<uint32_t>(r);
        return rank;
    }();
    return table.data();
}

// The order-rank table to use for `s`, or nullptr to hash directly with sketch_order_hash.
// Only the s = 5 preset has one; larger tables do not stay cached, so s = 9 hashes.
inline const uint32_t* closed_syncmer_order_rank_table_for(int s) {
    return s == 5 ? closed_syncmer_order_rank_table_impl<5>() : nullptr;
}

// Sliding window of the last W = k - s + 1 s-mer order values, tracking the window minimum
// and its multiplicity: closed-syncmer selection only asks whether the oldest or newest
// value equals the minimum.
struct SyncmerMinCountWindow {
    static constexpr int CAP = 64;
    static constexpr int MASK = CAP - 1;
    static_assert((CAP & MASK) == 0, "CAP must be a power of two");
    // k <= 31 gives W <= 31 < CAP; callers guard W <= CAP.

    uint64_t ord[CAP];
    int w = 0;        // logical window size W
    int head = 0;     // index of the oldest (left endpoint)
    int count = 0;
    uint64_t min_ord = std::numeric_limits<uint64_t>::max();
    int min_cnt = 0;

    void reset(int window) {
        w = window;
        head = 0;
        count = 0;
        min_ord = std::numeric_limits<uint64_t>::max();
        min_cnt = 0;
    }

    uint64_t oldest() const { return ord[head]; }

    void recompute_min() {
        uint64_t m = std::numeric_limits<uint64_t>::max();
        int c = 0;
        for (int i = 0; i < count; ++i) {
            const uint64_t o = ord[(head + i) & MASK];
            if (o < m) {
                m = o;
                c = 1;
            } else if (o == m) {
                ++c;
            }
        }
        min_ord = m;
        min_cnt = c;
    }

    void push(uint64_t o) {
        if (count == w) {  // window full: evict the oldest
            const uint64_t old = ord[head];
            head = (head + 1) & MASK;
            --count;
            if (old == min_ord && --min_cnt == 0) recompute_min();
        }
        ord[(head + count) & MASK] = o;
        ++count;
        if (o < min_ord) {
            min_ord = o;
            min_cnt = 1;
        } else if (o == min_ord) {
            ++min_cnt;
        }
    }
};

}}  // namespace fa::cpu

namespace fa { namespace cpu {

// Reference closed-syncmer seeds. A palindromic k-mer (even k only) is never a seed, on the
// reference or the query: its strand is unknowable, so the orientation bit could not
// resolve a hit. minimap2's sketching skips the same k-mers.
//
// `emit` is called as emit(key, pos, z): z is the orientation bit of the reference k-mer at
// `pos` (0 when its forward bases are the key itself).
template <class BaseT, class Emit>
inline bool emit_closed_syncmer_reference_seeds(
    const std::vector<BaseT>& chr,
    uint64_t base_offset,
    const ClosedSyncmerConfig& cfg,
    Emit&& emit
) {
    const int s = effective_closed_syncmer_s(cfg.k, cfg.s);
    if (!htable_valid_k(cfg.k) || !htable_valid_k(s) || s > cfg.k) return false;

    RollingKmerEncoder k_encoder(cfg.k);
    const int syncmer_window = cfg.k - s + 1;
    if (syncmer_window > 64) {
        for (size_t i = 0; i + static_cast<size_t>(cfg.k) <= chr.size(); ++i) {
            const BaseT* kmer = chr.data() + i;
            if (!closed_syncmer_selected_generic(kmer, cfg.k, s)) continue;
            uint64_t key = 0;
            bool is_rc = false;
            bool is_pal = false;
            if (!encode_kmer_key(kmer, cfg.k, &key, &is_rc, &is_pal) || is_pal)
                continue;
            emit(key,
                 static_cast<RefPos>(base_offset + static_cast<uint64_t>(i)),
                 static_cast<uint32_t>(is_rc ? 1u : 0u));
        }
        return true;
    }

    // One rolling k-state yields each rightmost s-mer; the min/count window tracks the
    // window minimum, and the order-rank table (when active) replaces the hash.
    const int W = syncmer_window;
    const uint32_t* otab = closed_syncmer_order_rank_table_for(s);
    SyncmerMinCountWindow window;
    window.reset(W);
    uint64_t current_order = std::numeric_limits<uint64_t>::max();

    for (size_t i = 0; i < chr.size(); ++i) {
        const int32_t base = static_cast<int32_t>(chr[i]);
        if (base < 0 || base > 3) {
            k_encoder.reset();
            window.reset(W);
            current_order = std::numeric_limits<uint64_t>::max();
            continue;
        }
        const bool full = k_encoder.push_base_unchecked(base);
        if (k_encoder.valid_bases() >= s) {
            const uint64_t sk = k_encoder.suffix_smer_key(s);
            current_order = otab
                ? static_cast<uint64_t>(otab[static_cast<size_t>(sk)])
                : sketch_order_hash(sk);
            window.push(current_order);
        } else {
            current_order = std::numeric_limits<uint64_t>::max();
        }

        if (!full) continue;
        const uint64_t best_order = window.min_ord;
        if (window.oldest() != best_order && current_order != best_order) {
            continue;
        }

        if (k_encoder.key_is_palindrome()) continue;
        const int pos = static_cast<int>(i + 1 - cfg.k);
        emit(k_encoder.key(),
             static_cast<RefPos>(base_offset + static_cast<uint64_t>(pos)),
             static_cast<uint32_t>(k_encoder.key_is_rc() ? 1u : 0u));
    }
    return true;
}

}}  // namespace fa::cpu

namespace fa { namespace cpu {

// What a read's closed-syncmer stream is reserved for: its expected count, 2 in
// every k - s + 1 k-mers and one in cfg.downsample of those, and an eighth
// more, so a read up to that much denser fills it without a reallocation.
inline size_t closed_syncmer_stream_reserve(int read_len,
                                            const ClosedSyncmerConfig& cfg,
                                            int s) {
    const size_t window = static_cast<size_t>(std::max(1, cfg.k - s + 1));
    const size_t thinning = static_cast<size_t>(std::max(1, cfg.downsample));
    const size_t expected =
        2 * static_cast<size_t>(std::max(0, read_len)) / (window * thinning);
    return std::max<size_t>(1, expected + expected / 8);
}

// The query stream for a window W = k - s + 1 <= 64, a chunk of positions at a time,
// with no branch on the bases: (1) each position's k-mer state, whether it has k valid
// bases, and the order of its rightmost s-mer; (2) the minimum of the last W orders, from
// blocks of W positions (the suffix minimum of the block holding the window's left end and
// the prefix minimum of the block holding its right end); (3) the positions whose leftmost
// or rightmost s-mer holds that minimum; (4) their k-mers, palindromes excepted. The state
// rolls through a non-ACGT base unreset: k valid bases in a row rewrite all of it, and the
// window of a k-mer with k valid bases holds no older order. With kEncode, enc[i] =
// nuc_encode(seq[i]) is written first. Seeds go to sink(seeds, n) in read order.
template <bool kEncode, class BaseT, class Sink>
inline void closed_syncmer_query_seed_chunks(
    const char* seq,
    BaseT* enc,
    int read_len,
    const ClosedSyncmerConfig& cfg,
    int s,
    Sink&& sink
) {
    constexpr int kChunk = 512;
    constexpr int kHist = 64;  // orders kept from the previous chunk, >= W
    const int k = cfg.k;
    const int W = k - s + 1;
    const uint64_t k_mask = htable_kmer_mask(k);
    const int rc_shift = 2 * (k - 1);
    const uint64_t s_mask = (1ULL << (2 * s)) - 1ULL;
    const int s_shift = 2 * (k - s);
    const uint32_t* otab = closed_syncmer_order_rank_table_for(s);
    const uint64_t thinning = static_cast<uint64_t>(std::max(1, cfg.downsample));

    uint64_t ord[kHist + kChunk];  // ord[kHist + i]: the order at chunk position i
    uint64_t suf[kHist + kChunk];
    uint64_t fwd_at[kChunk];
    uint64_t rev_at[kChunk];
    uint8_t full_at[kChunk];
    int sel_at[kChunk];
    QuerySeed picked[kChunk];
    std::fill(ord, ord + kHist, 0);
    uint64_t fwd = 0;
    uint64_t rev = 0;
    int run = 0;  // valid bases ending here

    for (int c0 = 0; c0 < read_len; c0 += kChunk) {
        const int n = std::min(kChunk, read_len - c0);
        uint64_t* o = ord + kHist;
        auto roll = [&](auto order_of) {
            for (int i = 0; i < n; ++i) {
                uint32_t code;
                if constexpr (kEncode) {
                    code = static_cast<uint32_t>(nuc_encode(seq[c0 + i]));
                    enc[c0 + i] = static_cast<uint8_t>(code);
                } else {
                    code = static_cast<uint32_t>(static_cast<int32_t>(enc[c0 + i]));
                }
                const uint64_t b = code & 3u;
                run = code <= 3u ? run + 1 : 0;
                fwd = ((fwd << 2) | b) & k_mask;
                rev = (rev >> 2) | ((3u - b) << rc_shift);
                fwd_at[i] = fwd;
                rev_at[i] = rev;
                full_at[i] = run >= k;
                o[i] = order_of(std::min(fwd & s_mask, rev >> s_shift));
            }
        };
        if (otab) {
            roll([otab](uint64_t sk) { return static_cast<uint64_t>(otab[sk]); });
        } else {
            roll([](uint64_t sk) { return sketch_order_hash(sk); });
        }

        // Blocks of W positions laid from the chunk's start, one more before it in the
        // history: suf[j] is the minimum from j to its block's end, `lo` from its block's
        // start to j. The window ending at j is one block, or the tail of one and the head
        // of the next, so its minimum is min(suf[j - W + 1], lo).
        const int end = kHist + n;
        uint64_t lo = ord[kHist - 1];
        for (int j = kHist - 1; j >= kHist - W; --j) suf[j] = lo = std::min(lo, ord[j]);
        int sel_n = 0;
        for (int b = kHist; b < end; b += W) {
            const int e = std::min(b + W, end);
            lo = ord[e - 1];
            for (int j = e - 1; j >= b; --j) suf[j] = lo = std::min(lo, ord[j]);
            lo = ord[b];
            for (int j = b; j < e; ++j) {
                lo = std::min(lo, ord[j]);
                const uint64_t m = std::min(suf[j - W + 1], lo);
                sel_at[sel_n] = j - kHist;
                sel_n += full_at[j - kHist] & ((m == ord[j - W + 1]) | (m == ord[j]));
            }
        }

        int picked_n = 0;
        for (int q = 0; q < sel_n; ++q) {
            const int i = sel_at[q];
            const uint64_t f = fwd_at[i];
            const uint64_t r = rev_at[i];
            picked[picked_n] =
                QuerySeed{std::min(f, r), c0 + i + 1 - k, static_cast<uint8_t>(r < f)};
            picked_n += f != r;
        }
        if (thinning > 1) {
            int thinned = 0;
            for (int q = 0; q < picked_n; ++q) {
                picked[thinned] = picked[q];
                thinned += sketch_order_hash(picked[q].key) % thinning == 0;
            }
            picked_n = thinned;
        }
        sink(static_cast<const QuerySeed*>(picked), picked_n);
        if (c0 + kChunk < read_len) std::copy(ord + n, ord + n + kHist, ord);
    }
}

// Query closed-syncmer seeds in read order, thinned by key hash when cfg.downsample > 1,
// passed to emit(seeds, n) in batches.
template <class BaseT, class Emit>
inline void for_each_closed_syncmer_query_seed(
    const BaseT* enc,
    int read_len,
    const ClosedSyncmerConfig& cfg,
    Emit&& emit
) {
    const int s = effective_closed_syncmer_s(cfg.k, cfg.s);
    if (!enc || !htable_valid_k(cfg.k) || !htable_valid_k(s) || s > cfg.k ||
        read_len < cfg.k) {
        return;
    }

    const int syncmer_window = cfg.k - s + 1;
    const int thinning = std::max(1, cfg.downsample);
    if (syncmer_window > 64) {
        for (int pos = 0; pos <= read_len - cfg.k; ++pos) {
            if (!closed_syncmer_selected_generic(enc + pos, cfg.k, s)) continue;
            uint64_t key = 0;
            bool is_rc = false;
            bool is_pal = false;
            if (!encode_kmer_key(enc + pos, cfg.k, &key, &is_rc, &is_pal) ||
                is_pal)
                continue;
            if (thinning > 1 &&
                sketch_order_hash(key) % static_cast<uint64_t>(thinning) != 0) {
                continue;
            }
            const QuerySeed seed{key, pos, static_cast<uint8_t>(is_rc ? 1 : 0)};
            emit(&seed, 1);
        }
        return;
    }

    closed_syncmer_query_seed_chunks<false>(nullptr, enc, read_len, cfg, s, emit);
}

template <class BaseT>
inline void extract_closed_syncmer_query_seeds_into(
    const BaseT* enc,
    int read_len,
    const ClosedSyncmerConfig& cfg,
    std::vector<QuerySeed>& out
) {
    out.clear();
    const int s = effective_closed_syncmer_s(cfg.k, cfg.s);
    if (!enc || !htable_valid_k(cfg.k) || !htable_valid_k(s) || s > cfg.k ||
        read_len < cfg.k) {
        return;
    }
    out.reserve(closed_syncmer_stream_reserve(read_len, cfg, s));
    for_each_closed_syncmer_query_seed(
        enc,
        read_len,
        cfg,
        [&](const QuerySeed* seeds, int n) { out.insert(out.end(), seeds, seeds + n); });
}

template <class Emit>
inline void encode_and_for_each_closed_syncmer_query_seed(
    const char* seq,
    int read_len,
    const ClosedSyncmerConfig& cfg,
    std::vector<uint8_t>& enc,
    Emit&& emit
) {
    enc.clear();
    if (!seq || read_len <= 0) return;
    enc.resize(static_cast<size_t>(read_len));

    const int s = effective_closed_syncmer_s(cfg.k, cfg.s);
    const bool can_extract =
        htable_valid_k(cfg.k) && htable_valid_k(s) && s <= cfg.k &&
        read_len >= cfg.k;
    const int syncmer_window = can_extract ? (cfg.k - s + 1) : 0;
    if (!can_extract || syncmer_window > 64) {
        for (int i = 0; i < read_len; ++i) {
            enc[static_cast<size_t>(i)] =
                static_cast<uint8_t>(nuc_encode(seq[static_cast<size_t>(i)]));
        }
        if (can_extract) {
            for_each_closed_syncmer_query_seed(
                enc.data(), read_len, cfg, std::forward<Emit>(emit));
        }
        return;
    }

    // As for_each_closed_syncmer_query_seed, with the 2-bit encoding folded into the pass.
    closed_syncmer_query_seed_chunks<true>(seq, enc.data(), read_len, cfg, s, emit);
}

inline void encode_and_extract_closed_syncmer_query_seeds_into(
    const char* seq,
    int read_len,
    const ClosedSyncmerConfig& cfg,
    std::vector<uint8_t>& enc,
    std::vector<QuerySeed>& out
) {
    out.clear();
    const int s = effective_closed_syncmer_s(cfg.k, cfg.s);
    if (seq && read_len >= cfg.k && htable_valid_k(cfg.k) &&
        htable_valid_k(s) && s <= cfg.k) {
        out.reserve(closed_syncmer_stream_reserve(read_len, cfg, s));
    }
    encode_and_for_each_closed_syncmer_query_seed(
        seq,
        read_len,
        cfg,
        enc,
        [&](const QuerySeed* seeds, int n) { out.insert(out.end(), seeds, seeds + n); });
}

// The reverse-complement lane is the forward stream with read_pos mirrored: same key, same
// postings. z is copied unchanged (see QuerySeed).
inline void project_query_seeds_to_rc_coordinates(
    const std::vector<QuerySeed>& fwd_seeds,
    int read_len,
    int k,
    std::vector<QuerySeed>& rc_out
) {
    rc_out.clear();
    if (read_len < k) return;
    rc_out.reserve(fwd_seeds.size());
    for (auto it = fwd_seeds.rbegin(); it != fwd_seeds.rend(); ++it) {
        const int rc_pos = read_len - k - it->read_pos;
        if (rc_pos >= 0) {
            rc_out.push_back(QuerySeed{it->key, rc_pos, it->z});
        }
    }
}

}}  // namespace fa::cpu
