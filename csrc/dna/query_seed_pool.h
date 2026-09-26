// Read-local access to the query seeds of both strands. Non-owning: the
// syncmer streams stay with the caller and the posting views in the
// ChainSeedLookupCache.
#pragma once

#include "../seeding/context.h"
#include "../seeding/scratch.h"
#include "../seeding/syncmer.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace dna {

enum class DnaQuerySeedStrand : unsigned char { Forward, Reverse };

class DnaQuerySeedPool {
public:
  DnaQuerySeedPool(const LongReadSeedContext& context, int read_length,
                   ChainSeedLookupCache& lookup_cache)
      : context_(context), read_length_(read_length),
        lookup_cache_(lookup_cache) {}

  void bind(DnaQuerySeedStrand strand,
            const std::vector<QuerySeed>& full_stream) {
    stream_slot(strand) = &full_stream;
    prepared_slot(strand) = false;
    captured_slot(strand) = false;
  }

  // Optional caller-owned stores for the capture's views and cache slots, one
  // entry per seed of the bound stream. Without them the capture only warms
  // the cache.
  void bind_capture_stores(std::vector<KmerPostingView>* forward_views,
                           std::vector<std::uint32_t>* forward_slots,
                           std::vector<KmerPostingView>* reverse_views,
                           std::vector<std::uint32_t>* reverse_slots) {
    forward_views_ = forward_views;
    forward_slots_ = forward_slots;
    reverse_views_ = reverse_views;
    reverse_slots_ = reverse_slots;
  }

  // The captured views / slots of a strand, or null when the strand was not
  // captured into a store.
  const std::vector<KmerPostingView>* captured_views(
      DnaQuerySeedStrand strand) const {
    const std::vector<KmerPostingView>* store =
        strand == DnaQuerySeedStrand::Forward ? forward_views_ : reverse_views_;
    const std::vector<QuerySeed>* seeds = stream(strand);
    return store != nullptr && seeds != nullptr && store->size() == seeds->size()
               ? store
               : nullptr;
  }
  const std::vector<std::uint32_t>* captured_slots(
      DnaQuerySeedStrand strand) const {
    const std::vector<std::uint32_t>* store =
        strand == DnaQuerySeedStrand::Forward ? forward_slots_ : reverse_slots_;
    const std::vector<QuerySeed>* seeds = stream(strand);
    return store != nullptr && seeds != nullptr && store->size() == seeds->size()
               ? store
               : nullptr;
  }

  bool bound(DnaQuerySeedStrand strand) const {
    const auto* seeds = stream(strand);
    return seeds != nullptr && !seeds->empty();
  }

  // Resolves every query key into the read-local cache. Later exact chains
  // read only these cached views.
  bool capture_posting_views(DnaQuerySeedStrand strand) {
    const std::vector<QuerySeed>* seeds = stream(strand);
    bool& captured = captured_slot(strand);
    if (captured) return true;
    if (seeds == nullptr || context_.index == nullptr) return false;
    if (strand == DnaQuerySeedStrand::Reverse && forward_captured_) {
      // The rc stream holds the forward stream's keys in reverse order, so
      // its views and slots are the forward stores reversed.
      if (reverse_views_ != nullptr && forward_views_ != nullptr &&
          forward_ != nullptr && forward_views_->size() == forward_->size() &&
          seeds->size() == forward_->size()) {
        mirror(*forward_views_, *reverse_views_);
      } else if (reverse_views_ != nullptr) {
        reverse_views_->clear();
      }
      if (reverse_slots_ != nullptr && forward_slots_ != nullptr &&
          forward_ != nullptr && forward_slots_->size() == forward_->size() &&
          seeds->size() == forward_->size()) {
        mirror(*forward_slots_, *reverse_slots_);
      } else if (reverse_slots_ != nullptr) {
        reverse_slots_->clear();
      }
      captured = true;
      return true;
    }
    // warm() batches the lookups so their cache misses overlap.
    std::vector<KmerPostingView>* views =
        strand == DnaQuerySeedStrand::Forward ? forward_views_ : reverse_views_;
    std::vector<std::uint32_t>* slots =
        strand == DnaQuerySeedStrand::Forward ? forward_slots_ : reverse_slots_;
    if (views != nullptr) views->resize(seeds->size());
    if (slots != nullptr) slots->resize(seeds->size());
    lookup_cache_.warm(*context_.index, seeds->data(), seeds->size(),
                       views != nullptr ? views->data() : nullptr,
                       slots != nullptr ? slots->data() : nullptr);
    captured = true;
    return true;
  }

  // Selects the representative vote seeds of the whole read, once per strand;
  // a second call returns false.
  bool prepare_representatives(
      DnaQuerySeedStrand strand, int start, int end, int downsample,
      std::vector<QuerySeed>& rebased_seeds, DnaLongSeedBundle& bundle,
      ChainWindowPeakScratch& selection_scratch) {
    const std::vector<QuerySeed>* seeds = stream(strand);
    bool& prepared = prepared_slot(strand);
    if (prepared || seeds == nullptr || start != 0 || end != read_length_)
      return false;
    prepared = true;
    const int window_downsample = std::max(1, downsample);
    // At downsample 1 the whole-read window is the bound stream itself, so it
    // is used in place instead of copied.
    const bool identity = syncmer_stream_window_is_identity(
        *seeds, start, end, context_.k, window_downsample);
    if (identity) {
      rebased_seeds.clear();
    } else {
      rebase_syncmer_stream_window(*seeds, start, end, context_.k,
                                   window_downsample, rebased_seeds);
    }
    // Captured views are positional, so they apply only to the stream in
    // place.
    const std::vector<KmerPostingView>* views =
        identity ? captured_views(strand) : nullptr;
    return select_chain_closed_syncmer_occ_aware_seed_bundle_from_seeds_into(
        context_, identity ? *seeds : rebased_seeds, std::max(0, end - start),
        &lookup_cache_, bundle, selection_scratch,
        views != nullptr ? views->data() : nullptr);
  }

private:
  const std::vector<QuerySeed>* stream(DnaQuerySeedStrand strand) const {
    return strand == DnaQuerySeedStrand::Forward ? forward_ : reverse_;
  }

  const std::vector<QuerySeed>*& stream_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_ : reverse_;
  }

  bool& prepared_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_prepared_
                                                 : reverse_prepared_;
  }

  bool& captured_slot(DnaQuerySeedStrand strand) {
    return strand == DnaQuerySeedStrand::Forward ? forward_captured_
                                                 : reverse_captured_;
  }

  template <class T>
  static void mirror(const std::vector<T>& from, std::vector<T>& to) {
    to.resize(from.size());
    for (std::size_t i = 0; i < from.size(); ++i)
      to[i] = from[from.size() - 1 - i];
  }

  const LongReadSeedContext& context_;
  int read_length_ = 0;
  ChainSeedLookupCache& lookup_cache_;
  const std::vector<QuerySeed>* forward_ = nullptr;
  const std::vector<QuerySeed>* reverse_ = nullptr;
  bool forward_prepared_ = false;
  bool reverse_prepared_ = false;
  bool forward_captured_ = false;
  bool reverse_captured_ = false;
  std::vector<KmerPostingView>* forward_views_ = nullptr;
  std::vector<std::uint32_t>* forward_slots_ = nullptr;
  std::vector<KmerPostingView>* reverse_views_ = nullptr;
  std::vector<std::uint32_t>* reverse_slots_ = nullptr;
};

} // namespace dna
} // namespace lr
} // namespace cpu
} // namespace fa
