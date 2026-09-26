// Fine-anchor candidate pools for the RNA locus harvest: value types, the posting memo
// and the harvest entry points.
#pragma once

#include "../../index/index.h"
#include "../../index/seed.h"
#include "skeleton_regions.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

enum class CandidatePoolKind : uint8_t {
  kSkeletonUnionUncapped = 1,
  kSkeletonUnionBounded = 2,
};

// 3 is unused; the values are part of the pool's provenance record, so they keep their
// numbers.
enum class FineAnchorDropReason : uint8_t {
  kNone = 0,
  kOccurrence = 1,
  kGeometry = 2,
  kFeatureCap = 4,
  kRegionCap = 5,
  kGlobalCap = 6,
};

struct FineAnchor {
  uint32_t anchor_id = 0;
  int32_t query_begin = 0;
  int32_t query_end = 0;
  int32_t reference_begin = 0;
  int32_t reference_end = 0;
  int32_t span = 0;
  uint32_t source_seed_id = 0;
  uint32_t source_occurrence = 0;
  uint32_t evidence_weight = 0;
  std::vector<uint32_t> region_ids;
  std::vector<uint32_t> compatible_core_nodes;
  uint32_t region_kind_mask = 0;
  uint32_t raw_multiplicity = 0;
  // Realization flags carried by the pool; region tags never set them.
  bool ignore_flag = false;
  bool tandem_flag = false;
  bool long_join_flag = false;
  bool raw_seen = false;
  bool survived_occurrence = false;
  bool survived_geometry = false;
  bool survived_feature_cap = false;
  bool survived_region_cap = false;
  bool survived_global_cap = false;
  bool final_pool = false;
  FineAnchorDropReason drop_reason = FineAnchorDropReason::kNone;
};

struct CandidatePoolLedger {
  uint64_t raw_seen = 0;
  uint64_t occurrence_filtered = 0;
  uint64_t geometry_filtered = 0;
  uint64_t feature_capped = 0;
  uint64_t region_capped = 0;
  uint64_t global_capped = 0;
  uint64_t final_pool = 0;
};

struct CandidatePool {
  CandidatePoolKind kind = CandidatePoolKind::kSkeletonUnionUncapped;
  std::vector<FineAnchor> anchors; // includes filtered/capped records
  CandidatePoolLedger ledger;
  uint64_t content_hash = 0;
  bool refused = false;
  AnchoringRefusal refusal = AnchoringRefusal::None;
};

struct SkeletonHarvestParams {
  // Skip enumerating the postings of a key the occurrence cap will reject anyway. The cap
  // reads the key's genome-wide count, which lookup_interval already reports, so those
  // postings would all fail survived_occurrence. They stay counted in the ledger, the
  // safety limit and the anchor identifiers, so the pool, ledger, hash and refusals match
  // full enumeration exactly; false enumerates every posting.
  bool pre_enumeration_occurrence_skip = true;
  // Read against the key's genome-wide count, so a key common anywhere is rejected even
  // where it is unique.
  uint32_t occurrence_cap = 200;
  uint32_t per_feature_budget = 16;
  uint32_t global_budget = 32768;
  // Refuse a seed whose canonical k-mer is low complexity before its postings are read:
  // >= 11/15 of one base, or <= 4 distinct 3-mers among the 13 overlapping ones (scaled
  // from k = 15). Under the log-priced intron transition a lone homopolymer anchor
  // anywhere within max_intron pays for itself, and such k-mers occur everywhere, so the
  // occurrence cap alone does not stop them. A masked seed is bookkept like an absent
  // key; the mirrored frame carries the same keys and masks the same seeds.
  bool mask_low_complexity = true;
  uint64_t raw_candidate_safety_limit = 250000;
  int64_t core_diagonal_tolerance = 64;
  uint32_t minimum_core_evidence_weight = 1;

  // Strand routing: every inspected posting goes to the lane its stored orientation bit
  // is compatible with (packed_ref_orientation_compatible against frame_reverse_lane or
  // its negation); exactly one holds, since no palindromic key is a seed. One raw
  // enumeration yields the nominated pool, in this frame, and the opposite pool,
  // mirrored into the other lane's coordinates. `frame_reverse_lane` is the lane
  // `strand_seeds` belongs to: false forward, true reverse complement.
  bool frame_reverse_lane = false;
};

// Strand-routing work counters. `postings_bulk_skipped` restates
// HarvestWorkCounters::postings_skipped_occurrence: those postings are never read, so
// their orientation is unknown and they enter neither pool.
struct RoutingCounters {
  uint64_t postings_inspected = 0;
  uint64_t postings_bulk_skipped = 0;
  uint64_t routed_nominated = 0;
  uint64_t routed_opposite = 0;
  uint64_t admitted_nominated = 0;
  uint64_t admitted_opposite = 0;
};

// Harvest work counters, kept out of CandidatePool so they never affect the pool's
// ledger or content hash. Postings the harvest declines to enumerate still
// count as available:
// postings_available == postings_enumerated + postings_skipped_occurrence.
struct HarvestWorkCounters {
  uint64_t seeds_visited = 0;
  uint64_t seeds_without_region = 0;  // no applicable search region
  uint64_t seeds_key_absent = 0;      // key not in the index at all
  uint64_t seeds_over_occurrence = 0; // key present, over the occurrence cap
  // Seeds refused by the low-complexity mask; disjoint from the counters below.
  uint64_t seeds_masked_low_complexity = 0;
  uint64_t interval_queries = 0;       // lookup_interval() calls
  uint64_t interval_queries_empty = 0; // calls returning no posting in-window
  uint64_t postings_available = 0;     // postings inside the queried intervals
  uint64_t postings_enumerated = 0;    // postings whose position was read
  // `over_occurrence`: available postings in intervals the occurrence cap rejects;
  // `skipped_occurrence`: those the harvest declined to enumerate (equal when the
  // pre-enumeration skip is on).
  uint64_t postings_over_occurrence = 0;
  uint64_t postings_skipped_occurrence = 0;
  uint64_t canonical_anchors = 0; // distinct (query, reference) records
  uint64_t regions_queried = 0;   // sum over seeds of applicable regions
  uint64_t merged_intervals = 0;  // sum over seeds of merged intervals
};

// Per-read memo of full posting views, shared by both strand streams. One read runs
// several harvests over its two seed vectors, each probing every key; lookup() is pure
// over an immutable index, so caching it changes only CPU. The reverse stream is the
// forward one walked backwards with the same canonical keys, so reverse seed j is
// forward seed n-1-j and storage is in the forward frame. `views` borrow the index's
// posting array; a memo is valid for one index.
//
// Two guards, because a wrong view is silent corruption: bind_reverse refuses a reverse
// stream whose size differs from the forward one's (a dropped seed breaks n-1-j), and
// every hit re-checks the stored key, falling back to a real probe and counting
// `mirror_key_mismatches`, which must stay zero. Freshness is a generation stamp that
// new_read() bumps in O(1); the buffers are never re-zeroed per read.
struct SeedPostingMemo {
  static constexpr std::size_t kNoSlot = static_cast<std::size_t>(-1);

  // Storage frame: the forward stream.
  const void* forward_identity = nullptr;
  std::size_t forward_size = 0;
  // Mirrored frame, bound only when it is provably a mirror of the forward one.
  const void* reverse_identity = nullptr;
  std::size_t reverse_size = 0;

  // 0 means "no read in progress": nothing hits and nothing is stored, so a
  // memo that never saw new_read() degrades to plain lookups.
  std::uint32_t generation = 0;
  std::vector<KmerPostingView> views;
  std::vector<std::uint64_t> keys;
  std::vector<std::uint32_t> stamps;

  // Guard 2's tripwire. Always zero unless the mirroring premise is false.
  std::uint64_t mirror_key_mismatches = 0;

  // Invalidate every slot in O(1). Called once per read by the caller.
  void new_read() {
    if (generation == (std::numeric_limits<std::uint32_t>::max)()) {
      // Wraparound: a stale stamp could alias the new generation, so this one
      // read pays for a real clear.
      std::fill(stamps.begin(), stamps.end(), std::uint32_t{0});
      generation = 1;
      return;
    }
    ++generation;
  }

  // Install the forward stream and size the buffers. Idempotent. Buffers only
  // ever grow: shrinking would cost the O(n) write the stamps exist to avoid.
  void bind_forward(const std::vector<QuerySeed>& seeds) {
    forward_identity = seeds.data();
    forward_size = seeds.size();
    if (views.size() < seeds.size()) {
      views.resize(seeds.size());
      keys.resize(seeds.size());
      stamps.resize(seeds.size(), std::uint32_t{0});
    }
  }

  // Install the reverse stream, but only if it is provably forward's mirror.
  // Idempotent and optional; an unbound reverse stream is simply not memoized.
  void bind_reverse(const std::vector<QuerySeed>& seeds) {
    if (forward_identity == nullptr || seeds.data() == forward_identity ||
        seeds.size() != forward_size) {
      reverse_identity = nullptr;
      reverse_size = 0;
      return;
    }
    reverse_identity = seeds.data();
    reverse_size = seeds.size();
  }

  // Which forward-frame slot backs `stream`'s seed `i`, or kNoSlot when this
  // stream is not a bound frame.
  std::size_t slot_for(const std::vector<QuerySeed>& stream,
                       std::size_t i) const {
    if (generation == 0)
      return kNoSlot;
    const void* identity = stream.data();
    if (identity == forward_identity && stream.size() == forward_size)
      return i < forward_size ? i : kNoSlot;
    if (reverse_identity != nullptr && identity == reverse_identity &&
        stream.size() == reverse_size)
      return i < reverse_size ? reverse_size - 1 - i : kNoSlot;
    return kNoSlot;
  }

  // The seed's full posting view, probing `index` only on a miss.
  KmerPostingView resolve(const SeedIndex& index,
                          const std::vector<QuerySeed>& stream, std::size_t i) {
    const QuerySeed& seed = stream[i];
    const std::size_t slot = slot_for(stream, i);
    if (slot == kNoSlot)
      return index.lookup(seed.key);
    if (stamps[slot] == generation) {
      if (keys[slot] == seed.key)
        return views[slot];
      // Guard 2. Loud where asserts live, safe everywhere.
      assert(false && "SeedPostingMemo mirrored slot key mismatch");
      ++mirror_key_mismatches;
    }
    const KmerPostingView view = index.lookup(seed.key);
    views[slot] = view;
    keys[slot] = seed.key;
    stamps[slot] = generation;
    return view;
  }
};

// A block of anchors full enumeration would have created and dropped on occurrence. The
// harvest visits seeds in ascending query order with distinct query starts, so a skipped
// seed's anchors occupy one contiguous identifier run. The block keeps the seed's query
// interval so the frame projection can mirror it (q -> length - q_end).
struct SkippedBlock {
  int32_t query_begin = 0;
  int32_t query_end = 0;
  uint64_t count = 0;
};

struct SkeletonCandidatePools {
  CandidatePool uncapped;
  CandidatePool bounded;
  // `nominated` is in the harvested frame's coordinates and carries the skeleton's core
  // tags. `opposite` is in the other lane's coordinates (collected from the raw postings,
  // then mirrored) and carries none. Each is capped on its own. `opposite` is refused
  // with MirrorProjectionDomain when the window has more than one search region.
  CandidatePool nominated;
  CandidatePool opposite;
  bool routed = false;
  RoutingCounters routing;
  HarvestWorkCounters work;
  // Ascending in query; kept so the frame projection can replay the identifier walk.
  std::vector<SkippedBlock> skipped_blocks;
  // The same blocks in the opposite lane's coordinates, produced by the same
  // mirror that produced `opposite`.
  std::vector<SkippedBlock> opposite_skipped_blocks;
};

// `posting_memo`, when non-null, serves resolved posting views for a bound stream (see
// SeedPostingMemo); the result is the same without it.
SkeletonCandidatePools harvest_skeleton_candidate_pools(
    const SeedIndex& index, const std::vector<QuerySeed>& strand_seeds, int chr,
    int32_t query_length, const std::vector<SkeletonNode>& skeleton,
    const std::vector<SearchRegion>& regions,
    const SkeletonHarvestParams& params,
    SeedPostingMemo* posting_memo = nullptr);

// Rewrites a bounded pool from one query frame into the other, in place, instead of
// harvesting the mirrored frame again. The reverse seed stream is the forward one
// mirrored (reverse seed j is forward seed n-1-j at length - span - q, with the same
// canonical key), so with the same reference window, one whole-query search region and
// the same caps the mirrored harvest reads exactly the same postings, and every field is
// frame-invariant or a known function of the frame. Preconditions: one search region
// covering the whole query, and `seed_count` equal to both streams' size; the caller
// refuses a window that violates one (MirrorProjectionDomain). Applying it twice restores
// everything but `compatible_core_nodes`, which it clears. A refused pool is returned
// untouched.
void mirror_candidate_pool_query_frame(CandidatePool& pool,
                                       std::vector<SkippedBlock>& skipped,
                                       int32_t query_length, int32_t seed_count,
                                       const std::vector<SearchRegion>& regions,
                                       const SkeletonHarvestParams& params,
                                       uint64_t skipped_occurrence);

// Re-decides the caps on a pool whose admissible evidence the caller narrowed in place
// (by clearing survived_geometry on anchors it will not chain) and rebuilds the ledger
// and content hash, so they describe the pool the chain is taken on. A refused pool is
// returned untouched.
void refinalize_candidate_pool(CandidatePool& pool,
                               const std::vector<SearchRegion>& regions,
                               const SkeletonHarvestParams& params,
                               uint64_t skipped_occurrence);

} // namespace rna
} // namespace lr
} // namespace cpu
} // namespace fa
