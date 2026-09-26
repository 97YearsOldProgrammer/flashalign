// One-pass RNA fine-anchor harvest: each seed is visited once, its regions' reference
// intervals are merged before posting lookup, postings are routed into the two strand
// lanes, and the deterministic caps are applied.
#include "skeleton_harvest.h"

#include "../../core/checked_range.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <map>
#include <tuple>
#include <utility>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {
namespace {

void sort_unique(std::vector<uint32_t>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

uint32_t evidence_weight(uint32_t occurrence) {
  uint32_t log = 0;
  uint32_t value = std::max<uint32_t>(1, occurrence);
  while (value > 1) {
    value = (value + 1) >> 1;
    ++log;
  }
  return std::max<uint32_t>(1, 16u > log ? 16u - log : 1u);
}

bool anchor_overlaps_node(const FineAnchor& anchor, const SkeletonNode& node,
                          int64_t diagonal_tolerance,
                          uint32_t minimum_evidence) {
  const int64_t diagonal =
      static_cast<int64_t>(anchor.reference_begin) - anchor.query_begin;
  return anchor.query_end > node.query_begin &&
         anchor.query_begin < node.query_end &&
         anchor.reference_end > node.reference_begin &&
         anchor.reference_begin < node.reference_end &&
         std::llabs(diagonal - node.median_diagonal) <= diagonal_tolerance &&
         anchor.evidence_weight >= minimum_evidence;
}

// Cap priority: anchors on a skeleton node first, then the rarest (genome-wide
// occurrence).
struct PriorityLess {
  bool operator()(const FineAnchor* left, const FineAnchor* right) const {
    const bool left_core = !left->compatible_core_nodes.empty();
    const bool right_core = !right->compatible_core_nodes.empty();
    if (left_core != right_core)
      return left_core > right_core;
    if (left->source_occurrence != right->source_occurrence)
      return left->source_occurrence < right->source_occurrence;
    if (left->evidence_weight != right->evidence_weight)
      return left->evidence_weight > right->evidence_weight;
    if (left->region_ids.size() != right->region_ids.size())
      return left->region_ids.size() > right->region_ids.size();
    if (left->query_begin != right->query_begin)
      return left->query_begin < right->query_begin;
    if (left->reference_begin != right->reference_begin)
      return left->reference_begin < right->reference_begin;
    return left->anchor_id < right->anchor_id;
  }
};

uint64_t fnv_mix_u64(uint64_t hash, uint64_t value) {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (byte * 8)) & 0xffu;
    hash *= 1099511628211ull;
  }
  return hash;
}

uint64_t pool_hash(const CandidatePool& pool) {
  uint64_t hash = 1469598103934665603ull;
  hash = fnv_mix_u64(hash, static_cast<uint64_t>(pool.kind));
  for (const FineAnchor& anchor : pool.anchors) {
    if (!anchor.final_pool)
      continue;
    hash = fnv_mix_u64(hash, anchor.anchor_id);
    hash = fnv_mix_u64(hash, static_cast<uint32_t>(anchor.query_begin));
    hash = fnv_mix_u64(hash, static_cast<uint32_t>(anchor.reference_begin));
    hash = fnv_mix_u64(hash, static_cast<uint32_t>(anchor.span));
    hash = fnv_mix_u64(hash, anchor.source_seed_id);
    hash = fnv_mix_u64(hash, anchor.source_occurrence);
    hash = fnv_mix_u64(hash, anchor.ignore_flag ? 1u : 0u);
    hash = fnv_mix_u64(hash, anchor.tandem_flag ? 1u : 0u);
    hash = fnv_mix_u64(hash, anchor.long_join_flag ? 1u : 0u);
    hash = fnv_mix_u64(hash, anchor.region_kind_mask);
    for (uint32_t value : anchor.region_ids)
      hash = fnv_mix_u64(hash, value);
    hash = fnv_mix_u64(hash, std::numeric_limits<uint64_t>::max());
    for (uint32_t value : anchor.compatible_core_nodes)
      hash = fnv_mix_u64(hash, value);
  }
  return hash;
}

void rebuild_ledger(CandidatePool& pool, uint64_t skipped_occurrence) {
  // `skipped_occurrence` postings were not enumerated because their key is over the
  // occurrence cap; full enumeration would have made that many occurrence-dropped
  // anchors, so they are folded in here.
  CandidatePoolLedger ledger;
  ledger.raw_seen = skipped_occurrence;
  ledger.occurrence_filtered = skipped_occurrence;
  for (const FineAnchor& anchor : pool.anchors) {
    ledger.raw_seen += anchor.raw_multiplicity;
    if (anchor.drop_reason == FineAnchorDropReason::kOccurrence) {
      ledger.occurrence_filtered += anchor.raw_multiplicity;
      continue;
    }
    if (anchor.drop_reason == FineAnchorDropReason::kGeometry) {
      ledger.geometry_filtered += anchor.raw_multiplicity;
      continue;
    }
    switch (anchor.drop_reason) {
    case FineAnchorDropReason::kFeatureCap:
      ++ledger.feature_capped;
      break;
    case FineAnchorDropReason::kRegionCap:
      ++ledger.region_capped;
      break;
    case FineAnchorDropReason::kGlobalCap:
      ++ledger.global_capped;
      break;
    default:
      break;
    }
    if (anchor.final_pool)
      ++ledger.final_pool;
  }
  pool.ledger = ledger;
  pool.content_hash = pool_hash(pool);
}

void apply_budgets(CandidatePool& pool,
                   const std::vector<SearchRegion>& regions,
                   const SkeletonHarvestParams& params,
                   uint64_t skipped_occurrence) {
  const PriorityLess priority_less;
  std::map<uint32_t, std::vector<FineAnchor*>> by_feature;
  for (FineAnchor& anchor : pool.anchors) {
    anchor.survived_feature_cap = false;
    anchor.survived_region_cap = false;
    anchor.survived_global_cap = false;
    anchor.final_pool = false;
    if (!anchor.survived_occurrence || !anchor.survived_geometry)
      continue;
    anchor.drop_reason = FineAnchorDropReason::kNone;
    by_feature[anchor.source_seed_id].push_back(&anchor);
  }
  size_t feature_survivors = 0;
  for (auto& item : by_feature) {
    std::vector<FineAnchor*>& values = item.second;
    std::sort(values.begin(), values.end(), priority_less);
    for (size_t i = 0; i < values.size(); ++i) {
      values[i]->survived_feature_cap = i < params.per_feature_budget;
      if (!values[i]->survived_feature_cap)
        values[i]->drop_reason = FineAnchorDropReason::kFeatureCap;
      else
        ++feature_survivors;
    }
  }

  // A single region whose budget holds every feature-cap survivor keeps all of them, so
  // the region pass is skipped. This is the usual case: one whole-query locus region.
  if (regions.size() == 1 &&
      static_cast<uint64_t>(regions[0].per_region_budget) >=
          static_cast<uint64_t>(feature_survivors)) {
    for (FineAnchor& anchor : pool.anchors)
      anchor.survived_region_cap = anchor.survived_feature_cap;
  } else {
    // Indexed by position in pool.anchors: identifiers may leave gaps for skipped
    // evidence, so they are not dense.
    std::vector<char> kept_by_region(pool.anchors.size(), 0);
    const FineAnchor* anchor_base = pool.anchors.data();
    for (const SearchRegion& region : regions) {
      std::vector<FineAnchor*> values;
      for (FineAnchor& anchor : pool.anchors) {
        if (!anchor.survived_feature_cap ||
            !std::binary_search(anchor.region_ids.begin(),
                                anchor.region_ids.end(), region.region_id)) {
          continue;
        }
        values.push_back(&anchor);
      }
      std::sort(values.begin(), values.end(), priority_less);
      const size_t keep =
          std::min<size_t>(region.per_region_budget, values.size());
      for (size_t i = 0; i < keep; ++i)
        kept_by_region[static_cast<size_t>(values[i] - anchor_base)] = 1;
    }
    for (FineAnchor& anchor : pool.anchors) {
      if (!anchor.survived_feature_cap)
        continue;
      anchor.survived_region_cap =
          kept_by_region[static_cast<size_t>(&anchor - anchor_base)] != 0;
      if (!anchor.survived_region_cap)
        anchor.drop_reason = FineAnchorDropReason::kRegionCap;
    }
  }

  std::vector<FineAnchor*> survivors;
  for (FineAnchor& anchor : pool.anchors)
    if (anchor.survived_region_cap)
      survivors.push_back(&anchor);
  // Under the budget every survivor is final, so the sort is skipped.
  if (survivors.size() <= static_cast<size_t>(params.global_budget)) {
    for (FineAnchor* anchor : survivors) {
      anchor->survived_global_cap = true;
      anchor->final_pool = true;
    }
  } else {
    std::sort(survivors.begin(), survivors.end(), priority_less);
    for (size_t i = 0; i < survivors.size(); ++i) {
      FineAnchor& anchor = *survivors[i];
      anchor.survived_global_cap = i < params.global_budget;
      anchor.final_pool = anchor.survived_global_cap;
      if (!anchor.final_pool)
        anchor.drop_reason = FineAnchorDropReason::kGlobalCap;
    }
  }
  rebuild_ledger(pool, skipped_occurrence);
}

// Assigns identifiers over records in ascending query order, advancing past each skipped
// block by the anchors it would have contributed, so anchor_id and the content hash
// match full enumeration. Shared with the frame projection.
void assign_exact_control_ids(std::vector<FineAnchor>& anchors,
                              const std::vector<SkippedBlock>& skipped) {
  size_t block = 0;
  uint32_t next_anchor_id = 0;
  for (FineAnchor& anchor : anchors) {
    while (block < skipped.size() &&
           skipped[block].query_begin < anchor.query_begin) {
      next_anchor_id += static_cast<uint32_t>(skipped[block].count);
      ++block;
    }
    anchor.anchor_id = next_anchor_id++;
  }
}

// Canonical records get their region tags deduplicated, their core-node tags
// from the skeleton, and their exact-control identifiers.
void canonicalize_raw(std::vector<FineAnchor>&& canonical,
                      const std::vector<SkeletonNode>& skeleton,
                      const SkeletonHarvestParams& params,
                      const std::vector<SkippedBlock>& skipped,
                      CandidatePool& pool) {
  pool.anchors = std::move(canonical);
  for (FineAnchor& anchor : pool.anchors) {
    sort_unique(anchor.region_ids);
    for (const SkeletonNode& node : skeleton) {
      if (anchor_overlaps_node(anchor, node, params.core_diagonal_tolerance,
                               params.minimum_core_evidence_weight)) {
        anchor.compatible_core_nodes.push_back(node.node_id);
      }
    }
    sort_unique(anchor.compatible_core_nodes);
  }
  assign_exact_control_ids(pool.anchors, skipped);
}

// Merges into `merged`, reusing the caller's buffers (called once per seed).
void merge_reference_intervals(
    const std::vector<const SearchRegion*>& regions,
    std::vector<std::pair<int64_t, int64_t>>& scratch,
    std::vector<std::pair<int64_t, int64_t>>& merged) {
  merged.clear();
  if (regions.empty())
    return;
  if (regions.size() == 1) {
    merged.push_back({regions[0]->reference_begin, regions[0]->reference_end});
    return;
  }
  scratch.clear();
  scratch.reserve(regions.size());
  for (const SearchRegion* region : regions)
    scratch.push_back({region->reference_begin, region->reference_end});
  std::sort(scratch.begin(), scratch.end());
  for (const auto& value : scratch) {
    if (merged.empty() || value.first > merged.back().second) {
      merged.push_back(value);
    } else {
      merged.back().second = std::max(merged.back().second, value.second);
    }
  }
}

// Low-complexity seed mask. QuerySeed::key is the 2-bit packed canonical k-mer, so the
// bases are read off the key with no sequence access. Both tests are invariant under
// reverse complement, so canonicality does not matter:
//   * the most frequent base covers >= 11 of 15 positions (a homopolymer or poly-A);
//   * <= 4 distinct 3-mers among the 13 overlapping ones (a short-period repeat).
// The thresholds are k = 15 ratios, applied proportionally for other k.
constexpr int kMaskHomopolymerNumerator = 11;
constexpr int kMaskHomopolymerDenominator = 15;
constexpr int kMaskDistinctTrimerNumerator = 4;
constexpr int kMaskDistinctTrimerDenominator = 13;

bool kmer_key_is_low_complexity(uint64_t key, int k) {
  if (k < 3 || k > 32)
    return false;
  int counts[4] = {0, 0, 0, 0};
  uint32_t bases[32];
  for (int i = 0; i < k; ++i) {
    const uint32_t base =
        static_cast<uint32_t>((key >> (2 * (k - 1 - i))) & 0x3u);
    bases[i] = base;
    ++counts[base];
  }
  int most = counts[0];
  for (int i = 1; i < 4; ++i)
    if (counts[i] > most)
      most = counts[i];
  if (most * kMaskHomopolymerDenominator >= kMaskHomopolymerNumerator * k)
    return true;
  // 64 possible 3-mers, so a bitset is the whole distinct-count.
  uint64_t seen = 0;
  for (int i = 0; i + 2 < k; ++i) {
    const uint32_t code = (bases[i] << 4) | (bases[i + 1] << 2) | bases[i + 2];
    seen |= (1ull << code);
  }
  int distinct = 0;
  for (uint64_t bits = seen; bits; bits &= bits - 1)
    ++distinct;
  return distinct * kMaskDistinctTrimerDenominator <=
         kMaskDistinctTrimerNumerator * (k - 2);
}

} // namespace

SkeletonCandidatePools harvest_skeleton_candidate_pools(
    const SeedIndex& index, const std::vector<QuerySeed>& strand_seeds, int chr,
    int32_t query_length, const std::vector<SkeletonNode>& skeleton,
    const std::vector<SearchRegion>& regions,
    const SkeletonHarvestParams& params, SeedPostingMemo* posting_memo) {
  SkeletonCandidatePools result;
  result.uncapped.kind = CandidatePoolKind::kSkeletonUnionUncapped;
  result.bounded.kind = CandidatePoolKind::kSkeletonUnionBounded;
  result.nominated.kind = CandidatePoolKind::kSkeletonUnionBounded;
  result.opposite.kind = CandidatePoolKind::kSkeletonUnionBounded;
  int seed_count = 0;
  int skeleton_count = 0;
  int region_count = 0;
  const bool counts_representable =
      mapping::checked_size_to_int(strand_seeds.size(), seed_count) &&
      mapping::checked_size_to_int(skeleton.size(), skeleton_count) &&
      mapping::checked_size_to_int(regions.size(), region_count) &&
      params.raw_candidate_safety_limit <=
          static_cast<uint64_t>(mapping::kSignedCoordinateMaximum);
  if (!counts_representable) {
    result.uncapped.refused = result.bounded.refused = true;
    result.uncapped.refusal = result.bounded.refusal =
        AnchoringRefusal::HarvestCountDomain;
    return result;
  }
  if (index.empty() || chr < 0 ||
      static_cast<uint64_t>(chr) >=
          static_cast<uint64_t>(index.chrom_count()) ||
      strand_seeds.empty() || query_length <= 0 || regions.empty()) {
    result.uncapped.refused = result.bounded.refused = true;
    result.uncapped.refusal = result.bounded.refusal =
        AnchoringRefusal::InvalidHarvestInput;
    return result;
  }
  const int span = index.k();
  const uint64_t* offsets = index.chrom_offsets_data();
  if (span <= 0 || !offsets) {
    result.uncapped.refused = result.bounded.refused = true;
    result.uncapped.refusal = result.bounded.refusal =
        AnchoringRefusal::InvalidIndexGeometry;
    return result;
  }
  const uint64_t chromosome_base = offsets[static_cast<size_t>(chr)];
  const uint32_t occurrence_cap = params.occurrence_cap;

  // Canonical (query_begin, reference_begin) records, appended already sorted and
  // unique: seeds ascend in query with distinct starts (checked below) and a seed's
  // in-window postings ascend in reference.
  std::vector<FineAnchor> canonical;
  // Routing flags per canonical record: bit 0 compatible with this frame's lane, bit 1
  // with the other lane.
  std::vector<uint8_t> route_flags;
  uint64_t raw_seen = 0;
  HarvestWorkCounters& work = result.work;
  RoutingCounters& routing = result.routing;
  // Per-seed scratch, hoisted out of the loop.
  std::vector<const SearchRegion*> applicable;
  std::vector<std::pair<int64_t, int64_t>> interval_scratch;
  std::vector<std::pair<int64_t, int64_t>> intervals;
  applicable.reserve(regions.size());

  // One seed per read position: a syncmer stream holds at most one seed per position and
  // the reverse projection is injective, so a duplicate means a stream this harvest is
  // not defined over. The identifier-preserving skip, the sorted append and the frame
  // projection all rely on it, so it is refused.
  for (size_t i = 1; i < strand_seeds.size(); ++i) {
    if (strand_seeds[i].read_pos == strand_seeds[i - 1].read_pos) {
      result.uncapped.refused = result.bounded.refused = true;
      result.uncapped.refusal = result.bounded.refusal =
          AnchoringRefusal::DuplicateQueryStart;
      return result;
    }
  }
  const bool skip_enabled = params.pre_enumeration_occurrence_skip;
  const bool mask_enabled = params.mask_low_complexity;
  std::vector<SkippedBlock> skipped_blocks;

  for (size_t seed_id = 0; seed_id < strand_seeds.size(); ++seed_id) {
    const QuerySeed& seed = strand_seeds[seed_id];
    const int32_t query_begin = seed.read_pos;
    ++work.seeds_visited;
    // The mask is read first, so a masked seed costs one key inspection.
    if (mask_enabled && kmer_key_is_low_complexity(seed.key, span)) {
      ++work.seeds_masked_low_complexity;
      continue;
    }
    int query_end = 0;
    if (!mapping::checked_signed_endpoint(query_begin, span, query_length,
                                          query_end))
      continue;
    applicable.clear();
    for (const SearchRegion& region : regions) {
      if (query_end > region.query_begin && query_begin < region.query_end) {
        applicable.push_back(&region);
      }
    }
    if (applicable.empty()) {
      ++work.seeds_without_region;
      continue;
    }
    work.regions_queried += applicable.size();
    merge_reference_intervals(applicable, interval_scratch, intervals);
    work.merged_intervals += intervals.size();
    bool key_present = false;
    bool key_over_occurrence = false;
    uint64_t seed_skipped = 0;
    // The seed's full posting view, resolved lazily at the first interval, at most once.
    KmerPostingView full_view;
    bool full_resolved = false;
    for (const auto& interval : intervals) {
      ++work.interval_queries;
      if (!full_resolved) {
        full_view = posting_memo
                        ? posting_memo->resolve(index, strand_seeds, seed_id)
                        : index.lookup(seed.key);
        full_resolved = true;
      }
      const KmerPostingIntervalView postings = index.narrow_interval(
          full_view, chr, static_cast<uint32_t>(interval.first),
          static_cast<uint32_t>(interval.second));
      const uint32_t occurrence = postings.global_count;
      key_present = key_present || postings.found();
      key_over_occurrence = key_over_occurrence ||
                            (postings.found() && occurrence > occurrence_cap);
      if (postings.count == 0 || !postings.positions.valid()) {
        ++work.interval_queries_empty;
        continue;
      }
      work.postings_available += postings.count;
      const bool over_occurrence = occurrence > occurrence_cap;
      if (over_occurrence)
        work.postings_over_occurrence += postings.count;
      // Every posting of an over-cap interval would fail survived_occurrence (the
      // verdict belongs to the key), so skipping is exact as long as each posting would
      // have made one record; keeping the interval's far end inside the checked
      // coordinate domain guarantees that.
      const bool interval_domain_safe =
          interval.second + span <= mapping::kSignedCoordinateMaximum;
      if (over_occurrence && skip_enabled && interval_domain_safe) {
        if (raw_seen + postings.count > params.raw_candidate_safety_limit) {
          // Full enumeration trips the limit on the posting that first exceeds it and
          // reports limit + 1; reproduce that value.
          result.uncapped.refused = result.bounded.refused = true;
          result.uncapped.refusal = result.bounded.refusal =
              AnchoringRefusal::RawCandidateSafetyLimit;
          result.uncapped.ledger.raw_seen = result.bounded.ledger.raw_seen =
              params.raw_candidate_safety_limit + 1;
          return result;
        }
        raw_seen += postings.count;
        seed_skipped += postings.count;
        work.postings_skipped_occurrence += postings.count;
        routing.postings_bulk_skipped += postings.count;
        continue;
      }
      for (uint32_t posting = 0; posting < postings.count; ++posting) {
        ++work.postings_enumerated;
        ++raw_seen;
        if (raw_seen > params.raw_candidate_safety_limit) {
          result.uncapped.refused = result.bounded.refused = true;
          result.uncapped.refusal = result.bounded.refusal =
              AnchoringRefusal::RawCandidateSafetyLimit;
          result.uncapped.ledger.raw_seen = result.bounded.ledger.raw_seen =
              raw_seen;
          return result;
        }
        // The posting word, read once: the position accessor masks the orientation bit
        // out.
        const PackedRefPos word = postings.positions.packed_at(posting);
        const bool compatible_nominated = packed_ref_orientation_compatible(
            seed.z, word, params.frame_reverse_lane);
        const bool compatible_opposite = packed_ref_orientation_compatible(
            seed.z, word, !params.frame_reverse_lane);
        ++routing.postings_inspected;
        if (compatible_nominated)
          ++routing.routed_nominated;
        if (compatible_opposite)
          ++routing.routed_opposite;
        const uint64_t global_position = postings.positions[posting];
        if (global_position < chromosome_base)
          continue;
        const uint64_t reference_offset = global_position - chromosome_base;
        if (reference_offset >
            static_cast<uint64_t>(mapping::kSignedCoordinateMaximum))
          continue;
        const int64_t reference64 = static_cast<int64_t>(reference_offset);
        int reference_end = 0;
        if (!mapping::checked_signed_endpoint(reference64, span,
                                              mapping::kSignedCoordinateMaximum,
                                              reference_end))
          continue;
        const int32_t reference_begin = static_cast<int32_t>(reference64);
        canonical.emplace_back();
        FineAnchor& anchor = canonical.back();
        anchor.query_begin = query_begin;
        anchor.query_end = query_end;
        anchor.reference_begin = reference_begin;
        anchor.reference_end = reference_end;
        anchor.span = span;
        anchor.source_seed_id = static_cast<uint32_t>(seed_id);
        anchor.source_occurrence = postings.global_count;
        anchor.evidence_weight = evidence_weight(occurrence);
        anchor.raw_seen = true;
        anchor.tandem_flag =
            query_seed_has_tandem_neighbor(strand_seeds, seed_id);
        anchor.raw_multiplicity = 1;
        anchor.survived_occurrence = occurrence <= occurrence_cap;
        for (const SearchRegion* region : applicable) {
          if (!search_region_overlaps_anchor(*region, query_begin,
                                             reference_begin, span)) {
            continue;
          }
          anchor.region_ids.push_back(region->region_id);
          anchor.region_kind_mask |= 1u << static_cast<uint32_t>(region->kind);
        }
        anchor.survived_geometry = !anchor.region_ids.empty();
        route_flags.push_back(
            static_cast<uint8_t>((compatible_nominated ? 1u : 0u) |
                                 (compatible_opposite ? 2u : 0u)));
      }
    }
    if (!key_present)
      ++work.seeds_key_absent;
    else if (key_over_occurrence)
      ++work.seeds_over_occurrence;
    if (seed_skipped > 0) {
      skipped_blocks.push_back(
          SkippedBlock{query_begin, query_end, seed_skipped});
    }
  }
  work.canonical_anchors = canonical.size();

  for (FineAnchor& anchor : canonical) {
    if (!anchor.survived_occurrence)
      anchor.drop_reason = FineAnchorDropReason::kOccurrence;
    else if (!anchor.survived_geometry)
      anchor.drop_reason = FineAnchorDropReason::kGeometry;
  }
  // The two routed pools, from the same canonical records after the drop-reason pass,
  // so each record keeps its occurrence and geometry verdict; every posting lands in
  // exactly one. Nothing reads an unrouted pool, so `result.bounded` stays empty and
  // unrefused.
  std::vector<FineAnchor> nominated_raw;
  std::vector<FineAnchor> opposite_raw;
  nominated_raw.reserve(routing.routed_nominated);
  opposite_raw.reserve(routing.routed_opposite);
  for (size_t i = 0; i < canonical.size() && i < route_flags.size(); ++i) {
    const bool to_nominated = (route_flags[i] & 1u) != 0u;
    const bool to_opposite = (route_flags[i] & 2u) != 0u;
    // `canonical` is dead after this loop, so records are moved.
    if (to_nominated)
      nominated_raw.push_back(std::move(canonical[i]));
    else if (to_opposite)
      opposite_raw.push_back(std::move(canonical[i]));
  }
  // The skipped blocks are already sorted by query. The nominated pool is this frame's:
  // canonicalized against the skeleton, with the same skipped blocks and caps.
  canonicalize_raw(std::move(nominated_raw), skeleton, params, skipped_blocks,
                   result.nominated);
  apply_budgets(result.nominated, regions, params,
                work.postings_skipped_occurrence);
  // The opposite pool is collected in this frame's coordinates and then mirrored, which
  // applies the coordinate and seed-id maps, the order repair, the identifier walk and
  // the caps. No skeleton: the voted peaks' query coordinates do not apply there.
  const std::vector<SkeletonNode> no_skeleton;
  canonicalize_raw(std::move(opposite_raw), no_skeleton, params,
                   skipped_blocks, result.opposite);
  if (regions.size() != 1) {
    // The mirror is defined only for a single whole-query region.
    result.opposite.refused = true;
    result.opposite.refusal = AnchoringRefusal::MirrorProjectionDomain;
  } else {
    result.opposite_skipped_blocks = skipped_blocks;
    mirror_candidate_pool_query_frame(
        result.opposite, result.opposite_skipped_blocks, query_length,
        seed_count, regions, params, work.postings_skipped_occurrence);
  }
  routing.admitted_nominated = result.nominated.ledger.final_pool;
  routing.admitted_opposite =
      result.opposite.refused ? 0 : result.opposite.ledger.final_pool;
  result.routed = true;
  result.skipped_blocks = std::move(skipped_blocks);
  return result;
}

void mirror_candidate_pool_query_frame(CandidatePool& pool,
                                       std::vector<SkippedBlock>& skipped,
                                       int32_t query_length, int32_t seed_count,
                                       const std::vector<SearchRegion>& regions,
                                       const SkeletonHarvestParams& params,
                                       uint64_t skipped_occurrence) {
  if (pool.refused)
    return;
  for (FineAnchor& anchor : pool.anchors) {
    // A seed at q covering [q, q + span) is the mirrored stream's seed at
    // length - span - q.
    anchor.query_begin = query_length - anchor.span - anchor.query_begin;
    anchor.query_end = anchor.query_begin + anchor.span;
    anchor.source_seed_id =
        static_cast<uint32_t>(seed_count - 1) - anchor.source_seed_id;
    // The mirrored frame has no skeleton, so no core tags.
    anchor.compatible_core_nodes.clear();
  }
  // Back into mirror-query-ascending order: reversing the vector fixes the seed order
  // and breaks the within-seed reference order, so each seed's run is reversed back.
  std::reverse(pool.anchors.begin(), pool.anchors.end());
  for (size_t begin = 0; begin < pool.anchors.size();) {
    size_t end = begin + 1;
    while (end < pool.anchors.size() &&
           pool.anchors[end].source_seed_id ==
               pool.anchors[begin].source_seed_id) {
      ++end;
    }
    std::reverse(pool.anchors.begin() + static_cast<std::ptrdiff_t>(begin),
                 pool.anchors.begin() + static_cast<std::ptrdiff_t>(end));
    begin = end;
  }
  // The un-enumerated blocks move with the seeds they belong to, so the
  // identifier walk below still finds them between the same neighbours.
  for (SkippedBlock& block : skipped) {
    const int32_t mirrored_begin = query_length - block.query_end;
    block.query_end = query_length - block.query_begin;
    block.query_begin = mirrored_begin;
  }
  std::reverse(skipped.begin(), skipped.end());

  assign_exact_control_ids(pool.anchors, skipped);
  // Re-cap from scratch: apply_budgets resets every survival flag first, so this is the
  // verdict the mirrored frame would reach. It also rebuilds the ledger and content
  // hash.
  apply_budgets(pool, regions, params, skipped_occurrence);
}

void refinalize_candidate_pool(CandidatePool& pool,
                               const std::vector<SearchRegion>& regions,
                               const SkeletonHarvestParams& params,
                               uint64_t skipped_occurrence) {
  if (pool.refused)
    return;
  apply_budgets(pool, regions, params, skipped_occurrence);
}

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
