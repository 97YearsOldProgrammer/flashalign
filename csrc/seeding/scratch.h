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
    std::vector<uint64_t> keys;
    std::vector<KmerPostingView> views;
    std::vector<uint32_t> stamps;
    size_t mask = 0;
    size_t used = 0;
    uint32_t epoch = 1;
    // warm()'s per-position slots when the caller wants views but not slots,
    // and the old-slot -> new-slot map of a rehash inside warm().
    std::vector<uint32_t> warm_slots_scratch;
    std::vector<uint32_t> rehash_map_scratch;

    static size_t next_power_of_two(size_t x) {
        size_t cap = 1;
        while (cap < x) cap <<= 1;
        return cap;
    }

    void allocate(size_t cap) {
        cap = std::max<size_t>(16, next_power_of_two(cap));
        keys.assign(cap, kEmptyKmerKey);
        views.assign(cap, KmerPostingView{});
        stamps.assign(cap, 0);
        mask = cap - 1;
        used = 0;
        epoch = 1;
    }

    void reset(size_t expected_entries) {
        const size_t target =
            std::max<size_t>(16, expected_entries * 2 + 1);
        const size_t target_cap = next_power_of_two(target);
        const size_t max_retain =
            std::max<size_t>(1024, target_cap * 8);
        if (!keys.empty() && keys.size() > max_retain) {
            allocate(target);
        } else if (keys.size() < target || keys.empty()) {
            allocate(target);
        } else {
            ++epoch;
            if (epoch == 0) {
                std::fill(stamps.begin(), stamps.end(), 0);
                epoch = 1;
            }
            used = 0;
        }
    }

    // `old_to_new`, when given, receives the new slot of every old slot that
    // held an entry (UINT32_MAX elsewhere), so a caller holding slot numbers
    // can translate them without probing the table again.
    void rehash(size_t new_cap, std::vector<uint32_t>* old_to_new = nullptr) {
        struct Moved {
            uint64_t key;
            KmerPostingView view;
            uint32_t old_slot;
        };
        std::vector<Moved> entries;
        entries.reserve(used);
        for (size_t i = 0; i < keys.size(); ++i) {
            if (stamps[i] == epoch) {
                entries.push_back({keys[i], views[i], static_cast<uint32_t>(i)});
            }
        }
        if (old_to_new != nullptr)
            old_to_new->assign(keys.size(), UINT32_MAX);
        allocate(new_cap);
        for (const auto& entry : entries) {
            size_t slot =
                static_cast<size_t>(fa::cpu::detail::splitmix64(entry.key)) &
                mask;
            while (stamps[slot] == epoch) slot = (slot + 1) & mask;
            stamps[slot] = epoch;
            keys[slot] = entry.key;
            views[slot] = entry.view;
            ++used;
            if (old_to_new != nullptr)
                (*old_to_new)[entry.old_slot] = static_cast<uint32_t>(slot);
        }
    }

    struct ProbeResult {
        size_t slot;  // hit: matching slot; miss: first empty slot (insert point)
        bool hit;
    };

    // One open-addressed walk: the matching slot on a hit, or the first empty slot (the
    // insertion point) on a miss.
    ProbeResult probe(uint64_t key) {
        if (keys.empty()) allocate(16);
        size_t slot =
            static_cast<size_t>(fa::cpu::detail::splitmix64(key)) & mask;
        for (;;) {
            if (stamps[slot] != epoch) return {slot, false};
            if (keys[slot] == key) return {slot, true};
            slot = (slot + 1) & mask;
        }
    }

    KmerPostingView lookup(const IndexView& index, uint64_t key) {
        ProbeResult p = probe(key);
        if (p.hit) {
            return views[p.slot];
        }
        KmerPostingView view = index.lookup(key);

        // Reuse the empty slot probe() found; only a table grow forces a re-probe.
        size_t slot = p.slot;
        if ((used + 1) * 2 >= keys.size()) {
            rehash(keys.size() * 2);
            slot = static_cast<size_t>(
                       fa::cpu::detail::splitmix64(key)) & mask;
            while (stamps[slot] == epoch) slot = (slot + 1) & mask;
        }
        stamps[slot] = epoch;
        keys[slot] = key;
        views[slot] = view;
        ++used;
        return view;
    }

    // Batched cache warm-up: the same end state as calling lookup() for every seed key in
    // order, but index lookups are deferred into batches whose memory loads overlap. A
    // deferred key is inserted into its slot at once, so repeats hit and the grow trigger
    // fires after the same keys; pending lookups are resolved before any rehash and at the
    // end, so no placeholder view survives. `out_views` / `out_slots`, when given, receive
    // per seed position the cached view and its slot, valid when warm() returns.
    void warm(const IndexView& index, const QuerySeed* seeds, size_t n,
              KmerPostingView* out_views = nullptr,
              uint32_t* out_slots = nullptr) {
        if (seeds == nullptr || n == 0) return;
        // Misses go to the batched lookup in runs of kWarmPending, long enough for its
        // prefetch pipeline to fill.
        constexpr size_t kWarmPending = 256;
        uint64_t pending_keys[kWarmPending];
        size_t pending_slots[kWarmPending];
        KmerPostingView pending_views[kWarmPending];
        size_t pending = 0;
        auto flush = [&]() {
            if (pending == 0) return;
            index.lookup_batch(pending_keys, pending, pending_views);
            for (size_t i = 0; i < pending; ++i) {
                views[pending_slots[i]] = pending_views[i];
            }
            pending = 0;
        };
        // Positional views are copied from the slots after the final flush, so no
        // unresolved placeholder is copied.
        uint32_t* slots = out_slots;
        if (slots == nullptr && out_views != nullptr) {
            warm_slots_scratch.resize(n);
            slots = warm_slots_scratch.data();
        }
        for (size_t i = 0; i < n; ++i) {
            const uint64_t key = seeds[i].key;
            ProbeResult p = probe(key);
            if (p.hit) {
                if (slots != nullptr) slots[i] = static_cast<uint32_t>(p.slot);
                continue;
            }
            size_t slot = p.slot;
            if ((used + 1) * 2 >= keys.size()) {
                flush();  // resolve placeholders before rehash() moves slots
                rehash(keys.size() * 2,
                       slots != nullptr ? &rehash_map_scratch : nullptr);
                slot = static_cast<size_t>(
                           fa::cpu::detail::splitmix64(key)) & mask;
                while (stamps[slot] == epoch) slot = (slot + 1) & mask;
                // Translate the slots recorded so far through the rehash's map.
                if (slots != nullptr)
                    for (size_t j = 0; j < i; ++j)
                        slots[j] = rehash_map_scratch[slots[j]];
            }
            stamps[slot] = epoch;
            keys[slot] = key;
            views[slot] = KmerPostingView{};
            ++used;
            if (slots != nullptr) slots[i] = static_cast<uint32_t>(slot);
            pending_keys[pending] = key;
            pending_slots[pending] = slot;
            ++pending;
            if (pending == kWarmPending) flush();
        }
        flush();
        if (out_views != nullptr)
            for (size_t i = 0; i < n; ++i) out_views[i] = views[slots[i]];
    }

    // The slot holding `key`, or npos. Const: nothing is resolved.
    static constexpr size_t npos = static_cast<size_t>(-1);
    size_t slot_of(uint64_t key) const {
        if (keys.empty()) return npos;
        size_t slot =
            static_cast<size_t>(fa::cpu::detail::splitmix64(key)) & mask;
        for (;;) {
            if (stamps[slot] != epoch) return npos;
            if (keys[slot] == key) return slot;
            slot = (slot + 1) & mask;
        }
    }
    // True when `slot` currently holds `key`; validates a slot recorded earlier.
    bool slot_holds(size_t slot, uint64_t key) const {
        return slot < keys.size() && stamps[slot] == epoch && keys[slot] == key;
    }
    const KmerPostingView& view_at(size_t slot) const { return views[slot]; }
    size_t capacity() const { return keys.size(); }

    // Read-only access: unlike lookup(), never resolves a missing key or mutates the cache.
    bool find_cached_view(uint64_t key, KmerPostingView& view) const {
        if (keys.empty()) return false;
        size_t slot =
            static_cast<size_t>(fa::cpu::detail::splitmix64(key)) & mask;
        for (;;) {
            if (stamps[slot] != epoch) return false;
            if (keys[slot] == key) {
                view = views[slot];
                return true;
            }
            slot = (slot + 1) & mask;
        }
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
        fwd.release_excess_for_read(seed_retain);
        rc.release_excess_for_read(seed_retain);
        shadow_rc.release_excess_for_read(seed_retain);
    }
};

}}}  // namespace fa::cpu::lr
