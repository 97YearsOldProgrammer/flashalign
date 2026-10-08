#include "retained_seed_density.h"

#include "../core/flat_int64_map.h"
#include "../seeding/posting_density.h"

#include <algorithm>
#include <limits>
#include <utility>

namespace fa::cpu::lr {
namespace {

bool compatible_views(const KmerPostingView& left,
                      const KmerPostingView& right) noexcept {
  if (!left.found() || !right.found()) return left.found() == right.found();
  return left.positions.data() == right.positions.data() &&
         left.positions.chromosome_offsets() ==
             right.positions.chromosome_offsets() &&
         left.count == right.count && left.occurrence == right.occurrence;
}

// lower_bound of one value over up to kSliceBatch sorted ranges in lock step:
// each round issues one independent load per live range, so the misses
// overlap. A null range yields 0.
void batched_lower_bound_packed(const PackedRefPos* const* data,
                                const std::uint32_t* count,
                                std::size_t lanes,
                                PackedRefPos value,
                                std::uint32_t* out) {
  const PackedRefPos* base[RetainedSeedDensity::kSliceBatch];
  std::uint32_t length[RetainedSeedDensity::kSliceBatch];
  PackedRefPos probe[RetainedSeedDensity::kSliceBatch] = {};
  for (std::size_t lane = 0; lane < lanes; ++lane) {
    base[lane] = data[lane];
    length[lane] = data[lane] != nullptr ? count[lane] : 0;
  }
  for (;;) {
    bool live = false;
    for (std::size_t lane = 0; lane < lanes; ++lane) {
      if (length[lane] == 0) continue;
      probe[lane] = base[lane][length[lane] >> 1];
      live = true;
    }
    if (!live) break;
    for (std::size_t lane = 0; lane < lanes; ++lane) {
      if (length[lane] == 0) continue;
      const std::uint32_t half = length[lane] >> 1;
      if (probe[lane] < value) {
        base[lane] += half + 1;
        length[lane] -= half + 1;
      } else {
        length[lane] = half;
      }
    }
  }
  for (std::size_t lane = 0; lane < lanes; ++lane) {
    out[lane] = data[lane] != nullptr
                    ? static_cast<std::uint32_t>(base[lane] - data[lane])
                    : 0;
  }
}

}  // namespace

bool RetainedSeedDensity::build(
    const SeedIndex& index,
    const std::vector<ChainWindowRetainedSeed>* forward,
    const std::vector<ChainWindowRetainedSeed>* reverse,
    const std::vector<QuerySeed>* fine_forward,
    const std::vector<QuerySeed>* fine_reverse,
    const ChainSeedLookupCache* lookup_cache,
    const std::vector<std::uint32_t>* fine_forward_slots,
    const std::vector<std::uint32_t>* fine_reverse_slots) {
  entries_.clear();
  forward_.clear();
  reverse_.clear();
  fine_forward_.clear();
  fine_reverse_.clear();
  validated_offsets_ = nullptr;
  validated_chrom_count_ = 0;
  if (index.empty() ||
      !posting_density_boundaries_valid(
          index.chrom_offsets_data(), index.chrom_count())) {
    return false;
  }
  validated_offsets_ = index.chrom_offsets_data();
  validated_chrom_count_ = index.chrom_count();

  const std::size_t forward_count = forward ? forward->size() : 0;
  const std::size_t reverse_count = reverse ? reverse->size() : 0;
  const std::size_t input_seed_records = forward_count + reverse_count;
  entries_.reserve(input_seed_records);
  forward_.reserve(forward_count);
  reverse_.reserve(reverse_count);

  // The key -> entry map. When the capture handed over each fine seed's cache
  // slot, the cache's slot space is the map: one int32 per slot and no second
  // hash table. A slot is trusted only after checking it still holds its key.
  // Otherwise a flat map is used; both paths produce the same entries.
  const bool slot_keyed =
      lookup_cache != nullptr && fine_forward_slots != nullptr &&
      fine_reverse_slots != nullptr && fine_forward != nullptr &&
      fine_reverse != nullptr &&
      fine_forward_slots->size() == fine_forward->size() &&
      fine_reverse_slots->size() == fine_reverse->size() &&
      lookup_cache->capacity() != 0;
  // On the slot-keyed path `by_key` holds only keys the cache does not.
  FlatInt64Map<std::uint32_t> by_key;
  if (slot_keyed) {
    entry_of_slot_.assign(lookup_cache->capacity(), -1);
    by_key.reserve(16);
  } else {
    by_key.reserve(input_seed_records + 1);
  }
  // The entry holding `key`, or -1. `slot_hint` is the seed's recorded slot
  // (npos for none); `slot_out` receives the key's slot, npos if uncached.
  const auto entry_of = [&](std::uint64_t key, std::size_t slot_hint,
                            std::size_t& slot_out) -> std::int32_t {
    slot_out = ChainSeedLookupCache::npos;
    if (slot_keyed) {
      std::size_t slot = slot_hint;
      if (!lookup_cache->slot_holds(slot, key)) slot = lookup_cache->slot_of(key);
      slot_out = slot;
      if (slot != ChainSeedLookupCache::npos) return entry_of_slot_[slot];
    }
    auto found = by_key.find(static_cast<std::int64_t>(key));
    return found == by_key.end() ? -1 : static_cast<std::int32_t>(found->second);
  };
  const auto remember = [&](std::uint64_t key, std::size_t slot,
                            std::uint32_t entry_index) {
    if (slot_keyed && slot != ChainSeedLookupCache::npos)
      entry_of_slot_[slot] = static_cast<std::int32_t>(entry_index);
    else
      by_key[static_cast<std::int64_t>(key)] = entry_index;
  };
  bool valid = true;
  auto ingest = [&](const std::vector<ChainWindowRetainedSeed>* input,
                    std::vector<RetainedSeedRef>& output) {
    if (input == nullptr) return;
    for (const ChainWindowRetainedSeed& retained : *input) {
      if (!retained.view.found()) {
        valid = false;
        continue;
      }
      std::size_t slot = ChainSeedLookupCache::npos;
      const std::int32_t held =
          entry_of(retained.seed.key, ChainSeedLookupCache::npos, slot);
      std::uint32_t entry_index = 0;
      if (held < 0) {
        if (entries_.size() >= static_cast<std::size_t>(
                                   std::numeric_limits<std::uint32_t>::max())) {
          valid = false;
          continue;
        }
        entry_index = static_cast<std::uint32_t>(entries_.size());
        remember(retained.seed.key, slot, entry_index);
        entries_.push_back(
            RetainedSeedDensityEntry{retained.seed.key, retained.view});
      } else {
        entry_index = static_cast<std::uint32_t>(held);
        if (!compatible_views(entries_[entry_index].view, retained.view))
          valid = false;
      }
      output.push_back(
          RetainedSeedRef{retained.seed, entry_index, retained.rescued});
      if (retained.rescued) entries_[entry_index].rescued = true;
    }
  };
  ingest(forward, forward_);
  ingest(reverse, reverse_);
  if (!valid) return false;

  // Per strand, the rescued vote seeds by (read_pos, key).
  std::vector<std::pair<int, std::uint64_t>> rescued[2];
  const std::vector<ChainWindowRetainedSeed>* votes[2] = {forward, reverse};
  for (int strand = 0; strand < 2; ++strand) {
    if (votes[strand] == nullptr) continue;
    for (const ChainWindowRetainedSeed& retained : *votes[strand])
      if (retained.rescued)
        rescued[strand].emplace_back(retained.seed.read_pos,
                                     retained.seed.key);
    std::sort(rescued[strand].begin(), rescued[strand].end());
  }

  // Per strand, the tandem fine seeds by (read_pos, key).
  std::vector<std::pair<int, std::uint64_t>> tandem_seeds[2];
  auto ingest_fine = [&](const std::vector<QuerySeed>* input,
                         const std::vector<std::uint32_t>* slots, int strand,
                         std::vector<RetainedSeedRef>& output) {
    if (input == nullptr) return;
    if (lookup_cache == nullptr) {
      valid = false;
      return;
    }
    output.reserve(input->size());
    for (std::size_t position = 0; position < input->size(); ++position) {
      const QuerySeed& seed = (*input)[position];
      std::size_t slot = ChainSeedLookupCache::npos;
      const std::int32_t held = entry_of(
          seed.key,
          slot_keyed ? static_cast<std::size_t>((*slots)[position])
                     : ChainSeedLookupCache::npos,
          slot);
      std::uint32_t entry_index = 0;
      if (held < 0) {
        KmerPostingView view;
        if (slot != ChainSeedLookupCache::npos) {
          view = lookup_cache->view_at(slot);
        } else if (!lookup_cache->find_cached_view(seed.key, view)) {
          valid = false;
          continue;
        }
        if (!view.found()) continue;
        if (entries_.size() >= static_cast<std::size_t>(
                                   std::numeric_limits<std::uint32_t>::max())) {
          valid = false;
          continue;
        }
        entry_index = static_cast<std::uint32_t>(entries_.size());
        remember(seed.key, slot, entry_index);
        entries_.push_back(RetainedSeedDensityEntry{seed.key, view});
      } else {
        entry_index = static_cast<std::uint32_t>(held);
      }
      const std::vector<std::pair<int, std::uint64_t>>& twins = rescued[strand];
      const bool twin =
          !twins.empty() &&
          std::binary_search(twins.begin(), twins.end(),
                             std::make_pair(seed.read_pos, seed.key));
      const bool tandem = query_seed_has_tandem_neighbor(*input, position);
      output.push_back(RetainedSeedRef{seed, entry_index, twin, tandem});
      if (twin) entries_[entry_index].rescued = true;
      if (tandem) tandem_seeds[strand].emplace_back(seed.read_pos, seed.key);
    }
  };
  ingest_fine(fine_forward, fine_forward_slots, 0, fine_forward_);
  ingest_fine(fine_reverse, fine_reverse_slots, 1, fine_reverse_);
  if (!valid) return false;

  std::vector<RetainedSeedRef>* vote_refs[2] = {&forward_, &reverse_};
  for (int strand = 0; strand < 2; ++strand) {
    std::vector<std::pair<int, std::uint64_t>>& tandem = tandem_seeds[strand];
    if (tandem.empty()) continue;
    std::sort(tandem.begin(), tandem.end());
    for (RetainedSeedRef& ref : *vote_refs[strand])
      ref.tandem = std::binary_search(
          tandem.begin(), tandem.end(),
          std::make_pair(ref.seed.read_pos, ref.seed.key));
  }

  return true;
}

KmerPostingIntervalView RetainedSeedDensity::slice(
    const SeedIndex& index,
    std::uint32_t entry_index,
    int chromosome,
    std::uint32_t low,
    std::uint32_t high) const {
  if (entry_index >= entries_.size() || chromosome < 0 ||
      chromosome >= static_cast<int>(index.chrom_count()) || high <= low) {
    return {};
  }
  const std::uint64_t* offsets = index.chrom_offsets_data();
  if ((offsets != validated_offsets_ ||
       index.chrom_count() != validated_chrom_count_) &&
      !posting_density_boundaries_valid(offsets, index.chrom_count()))
    return {};
  const KmerPostingView full = entries_[entry_index].view;
  if (!full.found()) return {};

  const std::uint64_t chromosome_base =
      offsets[static_cast<std::size_t>(chromosome)];
  const std::uint64_t chromosome_end =
      offsets[static_cast<std::size_t>(chromosome) + 1];
  std::uint64_t global_low =
      chromosome_base + static_cast<std::uint64_t>(low);
  std::uint64_t global_high =
      chromosome_base + static_cast<std::uint64_t>(high);
  global_low = std::min(global_low, chromosome_end);
  global_high = std::min(global_high, chromosome_end);
  if (global_high <= global_low) {
    return {full.positions, 0, full.occurrence};
  }

  const std::uint32_t interval_begin = full.positions.lower_bound(global_low);
  const std::uint32_t interval_end =
      interval_begin +
      full.positions.subspan(interval_begin, full.count - interval_begin)
          .lower_bound(global_high);
  return {full.positions.subspan(interval_begin, interval_end - interval_begin),
          interval_end - interval_begin, full.occurrence};
}

void RetainedSeedDensity::slice_batch(const SeedIndex& index,
                                      const std::uint32_t* entry_indices,
                                      std::size_t n,
                                      int chromosome,
                                      std::uint32_t low,
                                      std::uint32_t high,
                                      KmerPostingIntervalView* out) {
  if (n == 0 || entry_indices == nullptr || out == nullptr) return;

  // The per-window part of slice(), computed once.
  const std::uint64_t* offsets = index.chrom_offsets_data();
  const bool window_valid =
      chromosome >= 0 &&
      chromosome < static_cast<int>(index.chrom_count()) && high > low &&
      ((offsets == validated_offsets_ &&
        index.chrom_count() == validated_chrom_count_) ||
       posting_density_boundaries_valid(offsets, index.chrom_count()));
  bool empty_window = true;
  PackedRefPos packed_low = 0;
  PackedRefPos packed_high = 0;
  bool packed_low_valid = false;
  bool packed_high_valid = false;
  if (window_valid) {
    const std::uint64_t chromosome_base =
        offsets[static_cast<std::size_t>(chromosome)];
    const std::uint64_t chromosome_end =
        offsets[static_cast<std::size_t>(chromosome) + 1];
    const std::uint64_t global_low = std::min(
        chromosome_base + static_cast<std::uint64_t>(low), chromosome_end);
    const std::uint64_t global_high = std::min(
        chromosome_base + static_cast<std::uint64_t>(high), chromosome_end);
    empty_window = global_high <= global_low;
    if (!empty_window) {
      packed_low_valid = packed_ref_bound_from_global(
          offsets, index.chrom_count(), global_low, packed_low);
      packed_high_valid = packed_ref_bound_from_global(
          offsets, index.chrom_count(), global_high, packed_high);
    }
  }

  const KmerPostingView* views[kSliceBatch];
  const PackedRefPos* range_data[kSliceBatch];
  std::uint32_t range_count[kSliceBatch];
  std::uint32_t bound[kSliceBatch];
  std::uint32_t interval_begin[kSliceBatch];
  std::uint32_t interval_length[kSliceBatch];

  for (std::size_t block = 0; block < n; block += kSliceBatch) {
    const std::size_t lanes = std::min(kSliceBatch, n - block);
    for (std::size_t lane = 0; lane < lanes; ++lane) {
      out[block + lane] = {};
      views[lane] = nullptr;
      range_data[lane] = nullptr;
      range_count[lane] = 0;
      interval_begin[lane] = 0;
      interval_length[lane] = 0;
      if (!window_valid) continue;
      const std::uint32_t entry_index = entry_indices[block + lane];
      if (entry_index >= entries_.size()) continue;
      const KmerPostingView& full = entries_[entry_index].view;
      if (!full.found()) continue;
      if (empty_window) {
        out[block + lane] = {full.positions, 0, full.occurrence};
        continue;
      }
      if (full.positions.chromosome_offsets() != offsets ||
          full.positions.chromosome_count() != index.chrom_count()) {
        // A view with another chromosome directory packs the bounds
        // differently; take the scalar path.
        out[block + lane] = slice(index, entry_index, chromosome, low, high);
        continue;
      }
      views[lane] = &full;
      if (packed_low_valid) {
        range_data[lane] = full.positions.data();
        range_count[lane] = full.positions.size();
      } else {
        // An unpackable bound makes lower_bound return the span size.
        interval_begin[lane] = full.positions.size();
      }
    }

    if (packed_low_valid) {
      batched_lower_bound_packed(range_data, range_count, lanes, packed_low,
                                 bound);
      for (std::size_t lane = 0; lane < lanes; ++lane)
        if (views[lane] != nullptr) interval_begin[lane] = bound[lane];
    }

    for (std::size_t lane = 0; lane < lanes; ++lane) {
      range_data[lane] = nullptr;
      range_count[lane] = 0;
      if (views[lane] == nullptr) continue;
      const RefPosSpan tail = views[lane]->positions.subspan(
          interval_begin[lane], views[lane]->count - interval_begin[lane]);
      if (packed_high_valid) {
        range_data[lane] = tail.data();
        range_count[lane] = tail.size();
      } else {
        // Same fallback over the trailing subspan.
        interval_length[lane] = tail.size();
      }
    }

    if (packed_high_valid) {
      batched_lower_bound_packed(range_data, range_count, lanes, packed_high,
                                 bound);
      for (std::size_t lane = 0; lane < lanes; ++lane)
        if (views[lane] != nullptr) interval_length[lane] = bound[lane];
    }

    for (std::size_t lane = 0; lane < lanes; ++lane) {
      if (views[lane] == nullptr) continue;
      out[block + lane] = {
          views[lane]->positions.subspan(interval_begin[lane],
                                         interval_length[lane]),
          interval_length[lane], views[lane]->occurrence};
    }
  }
}

}  // namespace fa::cpu::lr
