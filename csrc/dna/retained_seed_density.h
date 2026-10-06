// Read-local posting views of the retained and fine seeds, deduplicated by
// canonical key while each strand keeps its own query positions. A candidate
// slices a view to its reference window with two lower bounds and no new key
// lookup.
#pragma once

#include "../index/index.h"
#include "../seeding/scratch.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

struct RetainedSeedDensityEntry {
  std::uint64_t key = 0;
  KmerPostingView view;
  // Some RetainedSeedRef of this key is rescued.
  bool rescued = false;
};

struct RetainedSeedRef {
  QuerySeed seed;
  std::uint32_t entry = 0;
  // A vote seed the empty-tile rescue admitted past the cap
  // (options/dna_profile.h kDnaTileRescueOcc), or its fine twin, which
  // placement's harvest gate admits too.
  bool rescued = false;
};

class RetainedSeedDensity {
 public:
  // A rescued vote seed's fine twin is the fine seed of its strand with its
  // read_pos and key: both come from the strand's one seed stream.
  bool build(const SeedIndex& index,
             const std::vector<ChainWindowRetainedSeed>* forward,
             const std::vector<ChainWindowRetainedSeed>* reverse,
             const std::vector<QuerySeed>* fine_forward = nullptr,
             const std::vector<QuerySeed>* fine_reverse = nullptr,
             const ChainSeedLookupCache* lookup_cache = nullptr,
             const std::vector<std::uint32_t>* fine_forward_slots = nullptr,
             const std::vector<std::uint32_t>* fine_reverse_slots = nullptr);

  const std::vector<RetainedSeedRef>& forward() const { return forward_; }
  const std::vector<RetainedSeedRef>& reverse() const { return reverse_; }
  const std::vector<RetainedSeedRef>& fine_forward() const {
    return fine_forward_;
  }
  const std::vector<RetainedSeedRef>& fine_reverse() const {
    return fine_reverse_;
  }

  // Entries, one per distinct key; every RetainedSeedRef::entry is below it.
  std::size_t entry_count() const { return entries_.size(); }

  // The genome-wide occurrence of an entry's key: the global_count of every
  // interval slice() finds for it.
  std::uint32_t occurrence(std::uint32_t entry_index) const {
    return entries_[entry_index].view.occurrence;
  }

  // Whether a rescued seed refers to the entry.
  bool holds_rescued(std::uint32_t entry_index) const {
    return entries_[entry_index].rescued;
  }

  // The entry's whole posting list.
  const KmerPostingView& view(std::uint32_t entry_index) const {
    return entries_[entry_index].view;
  }

  // Intervals resolved per lock-step round of slice_batch, which accepts any
  // n and chunks internally.
  static constexpr std::size_t kSliceBatch = 16;

  KmerPostingIntervalView slice(const SeedIndex& index,
                                std::uint32_t entry_index,
                                int chromosome,
                                std::uint32_t low,
                                std::uint32_t high) const;

  // Resolves n intervals for the same (chromosome, low, high) window; result i
  // equals slice(index, entry_indices[i], chromosome, low, high). The binary
  // searches of up to kSliceBatch entries are interleaved so their cache
  // misses overlap.
  void slice_batch(const SeedIndex& index,
                   const std::uint32_t* entry_indices,
                   std::size_t n,
                   int chromosome,
                   std::uint32_t low,
                   std::uint32_t high,
                   KmerPostingIntervalView* out);

 private:
  std::vector<RetainedSeedDensityEntry> entries_;
  // Entry index per lookup-cache slot for the slot-keyed build (-1: none).
  std::vector<std::int32_t> entry_of_slot_;
  std::vector<RetainedSeedRef> forward_;
  std::vector<RetainedSeedRef> reverse_;
  std::vector<RetainedSeedRef> fine_forward_;
  std::vector<RetainedSeedRef> fine_reverse_;
  // The chromosome directory build() validated; slicing with the same
  // directory skips re-validation.
  const std::uint64_t* validated_offsets_ = nullptr;
  std::uint64_t validated_chrom_count_ = 0;
};

}  // namespace fa::cpu::lr
