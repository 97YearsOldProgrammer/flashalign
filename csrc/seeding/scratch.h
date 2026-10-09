// Per-thread scratch storage for long-read seeding and voting.
#pragma once

#include "../index/format.h"      // KmerPostingView
#include "../core/hash.h"         // detail::splitmix64
#include "../index/index.h"       // IndexView
#include "../index/seed.h"        // QuerySeed, kEmptyKmerKey
#include "types.h"                // ChainWindow* records, VotePeak, DnaLongSeed*
#include "../core/flat_int64_map.h"      // FlatInt64Map

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace lr {

// Diagonal-bin vote table. The vote comparators are total orders over the (chr, bin)
// buckets, so output does not depend on the table's iteration order.
#if FA_DNA_LONG_VOTE_FLAT
using ChainWindowVoteTable = FlatInt64Map<ChainWindowBucketVote>;
#else
using ChainWindowVoteTable = std::unordered_map<int64_t, ChainWindowBucketVote>;
#endif

struct ChainWindowPeakScratch {
    std::vector<QuerySeed> seeds;
    std::vector<QuerySeed> extracted_seeds;
    DnaLongSeedBundle seed_bundle;
    std::vector<KmerPostingView> seed_views;
    std::vector<uint8_t> seed_used;
    std::vector<ChainWindowRetainedSeed> retained_seeds;
    // Caller-prepared DNA views still owed the vote's exact refine pass. clear_for_window
    // keeps them; window_anchor_peaks consumes and clears them.
    std::vector<DnaLongSeedView> pending_exact_refine_views;
    ChainWindowVoteTable buckets;
    std::vector<ChainWindowRankedBucket> ranked;
    std::vector<int> ref_starts;
    std::vector<uint32_t> occurrence_values;

    void clear_for_window() {
        seeds.clear();
        extracted_seeds.clear();
        seed_bundle.clear();
        seed_views.clear();
        seed_used.clear();
        retained_seeds.clear();
        buckets.clear();
        ranked.clear();
        ref_starts.clear();
        occurrence_values.clear();
    }

    template <class T>
    static void release_vec_if_excessive(
        std::vector<T>& v,
        size_t target_retain
    ) {
        target_retain = std::max<size_t>(target_retain, 1024);
        if (v.capacity() > target_retain * 4) {
            std::vector<T>().swap(v);
        }
    }

    void release_excess_for_read(size_t seed_retain) {
        release_vec_if_excessive(seeds, seed_retain);
        release_vec_if_excessive(extracted_seeds, seed_retain);
        release_vec_if_excessive(seed_views, seed_retain);
        release_vec_if_excessive(seed_used, seed_retain);
        release_vec_if_excessive(retained_seeds, seed_retain);
        release_vec_if_excessive(pending_exact_refine_views, seed_retain * 2);
        release_vec_if_excessive(ranked, seed_retain * 2);
        release_vec_if_excessive(ref_starts, seed_retain * 4);
        release_vec_if_excessive(occurrence_values, seed_retain);
        const size_t bucket_retain = std::max<size_t>(2048, seed_retain * 2);
        if (buckets.bucket_count() > bucket_retain * 4) {
            ChainWindowVoteTable().swap(buckets);
        }
    }
};

struct ChainSeedLookupCache {
    // Hash buckets hold dense entry IDs; captured IDs survive a rehash.
    std::vector<uint32_t> buckets;
    std::vector<uint32_t> stamps;
    std::vector<uint64_t> keys;
    std::vector<KmerPostingView> views;
    size_t mask = 0;
    size_t used = 0;
    uint32_t epoch = 1;
    std::vector<uint32_t> warm_slots_scratch;

    static size_t next_power_of_two(size_t x) {
        size_t cap = 1;
        while (cap < x) cap <<= 1;
        return cap;
    }

    void allocate(size_t cap) {
        cap = std::max<size_t>(16, next_power_of_two(cap));
        buckets.resize(cap);
        stamps.assign(cap, 0);
        keys.clear();
        views.clear();
        if (keys.capacity() > cap / 2) {
            std::vector<uint64_t>().swap(keys);
            std::vector<KmerPostingView>().swap(views);
        }
        keys.reserve(cap / 2);
        views.reserve(cap / 2);
        mask = cap - 1;
        used = 0;
        epoch = 1;
    }

    void reset(size_t expected_entries) {
        const size_t target = std::max<size_t>(16, expected_entries * 2 + 1);
        const size_t target_cap = next_power_of_two(target);
        const size_t max_retain = std::max<size_t>(1024, target_cap * 8);
        if (buckets.size() > max_retain || buckets.size() < target) {
            allocate(target);
        } else {
            ++epoch;
            if (epoch == 0) {
                std::fill(stamps.begin(), stamps.end(), 0);
                epoch = 1;
            }
            keys.clear();
            views.clear();
            used = 0;
        }
    }

    void rehash(size_t new_cap) {
        const size_t cap = std::max<size_t>(16, next_power_of_two(new_cap));
        buckets.resize(cap);
        stamps.assign(cap, 0);
        mask = cap - 1;
        epoch = 1;
        for (size_t entry = 0; entry < used; ++entry) {
            size_t slot = static_cast<size_t>(
                fa::cpu::detail::splitmix64(keys[entry])) & mask;
            while (stamps[slot] == epoch) slot = (slot + 1) & mask;
            stamps[slot] = epoch;
            buckets[slot] = static_cast<uint32_t>(entry);
        }
    }

    struct ProbeResult {
        size_t slot;
        bool hit;
    };

    ProbeResult probe(uint64_t key) {
        if (buckets.empty()) allocate(16);
        size_t slot = static_cast<size_t>(fa::cpu::detail::splitmix64(key)) & mask;
        for (;;) {
            if (stamps[slot] != epoch) return {slot, false};
            if (keys[buckets[slot]] == key) return {slot, true};
            slot = (slot + 1) & mask;
        }
    }

    KmerPostingView lookup(const IndexView& index, uint64_t key) {
        ProbeResult p = probe(key);
        if (p.hit) return views[buckets[p.slot]];
        KmerPostingView view = index.lookup(key);
        size_t slot = p.slot;
        if ((used + 1) * 2 >= buckets.size()) {
            rehash(buckets.size() * 2);
            slot = probe(key).slot;
        }
        stamps[slot] = epoch;
        buckets[slot] = static_cast<uint32_t>(used);
        keys.push_back(key);
        views.push_back(view);
        ++used;
        return view;
    }

    // Resolve misses in batches, then return each position's view and dense ID.
    void warm(const IndexView& index, const QuerySeed* seeds, size_t n,
              KmerPostingView* out_views = nullptr,
              uint32_t* out_slots = nullptr) {
        if (seeds == nullptr || n == 0) return;
        constexpr size_t kWarmPending = 256;
        uint64_t pending_keys[kWarmPending];
        size_t pending_slots[kWarmPending];
        KmerPostingView pending_views[kWarmPending];
        size_t pending = 0;
        auto flush = [&]() {
            if (pending == 0) return;
            index.lookup_batch(pending_keys, pending, pending_views);
            for (size_t i = 0; i < pending; ++i)
                views[pending_slots[i]] = pending_views[i];
            pending = 0;
        };
        uint32_t* slots = out_slots;
        if (slots == nullptr && out_views != nullptr) {
            warm_slots_scratch.resize(n);
            slots = warm_slots_scratch.data();
        }
        for (size_t i = 0; i < n; ++i) {
            const uint64_t key = seeds[i].key;
            ProbeResult p = probe(key);
            if (p.hit) {
                if (slots != nullptr) slots[i] = buckets[p.slot];
                continue;
            }
            size_t slot = p.slot;
            if ((used + 1) * 2 >= buckets.size()) {
                flush();
                rehash(buckets.size() * 2);
                slot = probe(key).slot;
            }
            const auto entry = static_cast<uint32_t>(used);
            stamps[slot] = epoch;
            buckets[slot] = entry;
            keys.push_back(key);
            views.push_back(KmerPostingView{});
            ++used;
            if (slots != nullptr) slots[i] = entry;
            pending_keys[pending] = key;
            pending_slots[pending] = entry;
            ++pending;
            if (pending == kWarmPending) flush();
        }
        flush();
        if (out_views != nullptr)
            for (size_t i = 0; i < n; ++i) out_views[i] = views[slots[i]];
    }

    static constexpr size_t npos = static_cast<size_t>(-1);
    size_t slot_of(uint64_t key) const {
        if (buckets.empty()) return npos;
        size_t slot = static_cast<size_t>(fa::cpu::detail::splitmix64(key)) & mask;
        for (;;) {
            if (stamps[slot] != epoch) return npos;
            const uint32_t entry = buckets[slot];
            if (keys[entry] == key) return entry;
            slot = (slot + 1) & mask;
        }
    }
    bool slot_holds(size_t entry, uint64_t key) const {
        return entry < used && keys[entry] == key;
    }
    const KmerPostingView& view_at(size_t entry) const { return views[entry]; }
    uint64_t key_at(size_t entry) const { return keys[entry]; }
    size_t entry_count() const { return used; }
    size_t capacity() const { return buckets.size(); }
    void prefetch_slot(size_t entry) const {
        if (entry >= used) return;
        __builtin_prefetch(&keys[entry]);
        const char* view = reinterpret_cast<const char*>(&views[entry]);
        __builtin_prefetch(view);
        __builtin_prefetch(view + sizeof(KmerPostingView) - 1);
    }

    bool find_cached_view(uint64_t key, KmerPostingView& view) const {
        const size_t entry = slot_of(key);
        if (entry == npos) return false;
        view = views[entry];
        return true;
    }
};

struct ChainAnchorScratch {
    std::vector<VotePeak> raw;
    ChainWindowPeakScratch fwd;
    ChainWindowPeakScratch rc;
    ChainWindowPeakScratch shadow_rc;
    ChainSeedLookupCache lookup_cache;
    std::vector<QuerySeed> full_syncmer_rc_seeds;
    std::vector<QuerySeed> shared_window_fwd_seeds;
    std::vector<QuerySeed> shared_window_rc_seeds;
    // The read's capture, by seed position (dna/query_seed_pool.h).
    std::vector<KmerPostingView> captured_fwd_views;
    std::vector<uint32_t> captured_fwd_slots;
    std::vector<KmerPostingView> captured_rc_views;
    std::vector<uint32_t> captured_rc_slots;
    // Each lane's exact-refine views ([0] forward, [1] reverse), kept for the winner's
    // line fit after the vote and the reverse prepare have cleared their own copies.
    std::vector<DnaLongSeedView> slope_refine_views[2];

    void clear_for_read() {
        raw.clear();
        fwd.clear_for_window();
        rc.clear_for_window();
        shadow_rc.clear_for_window();
        full_syncmer_rc_seeds.clear();
        shared_window_fwd_seeds.clear();
        shared_window_rc_seeds.clear();
        captured_fwd_views.clear();
        captured_fwd_slots.clear();
        captured_rc_views.clear();
        captured_rc_slots.clear();
        for (int lane = 0; lane < 2; ++lane)
            slope_refine_views[lane].clear();
    }

    template <class T>
    static void release_vec_if_excessive(
        std::vector<T>& v,
        size_t target_retain
    ) {
        target_retain = std::max<size_t>(target_retain, 1024);
        if (v.capacity() > target_retain * 4) {
            std::vector<T>().swap(v);
        }
    }

    void release_excess_for_read(size_t read_len) {
        const size_t seed_retain = std::max<size_t>(4096, read_len / 2);
        const size_t window_retain = std::max<size_t>(64, read_len / 256);
        release_vec_if_excessive(raw, window_retain * 4);
        release_vec_if_excessive(full_syncmer_rc_seeds, seed_retain);
        release_vec_if_excessive(shared_window_fwd_seeds, seed_retain);
        release_vec_if_excessive(shared_window_rc_seeds, seed_retain);
        release_vec_if_excessive(captured_fwd_views, seed_retain);
        release_vec_if_excessive(captured_fwd_slots, seed_retain);
        release_vec_if_excessive(captured_rc_views, seed_retain);
        release_vec_if_excessive(captured_rc_slots, seed_retain);
        for (int lane = 0; lane < 2; ++lane)
            release_vec_if_excessive(slope_refine_views[lane], seed_retain * 2);
        fwd.release_excess_for_read(seed_retain);
        rc.release_excess_for_read(seed_retain);
        shadow_rc.release_excess_for_read(seed_retain);
    }
};

}}}  // namespace fa::cpu::lr
