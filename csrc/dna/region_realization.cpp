#include "region_realization.h"

#include "family_realization_internal.h"
#include "ordered_anchor_path.h"
#include "placement_chaining.h"
#include "postdp_scoring.h"

#include "../seeding/tie_hash.h"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace fa::cpu::lr {
namespace {

namespace ordered = ::fa::cpu::lr::ordered_anchor;
using family_internal::BlockPlan;
using family_internal::BlockSplit;
using family_internal::Unit;

// minimap2's best_n (-N), mask_level, mask_len, rank_min_len and rank_frac
// (options.c).
constexpr int kBestN = 5;
constexpr float kMaskLevel = 0.5f;
constexpr int kMaskLen = INT_MAX;
constexpr int kRankMinLen = 500;
constexpr float kRankFrac = 0.9f;

// One chain of the selection: its anchors in chain order, deduplicated, and
// its forward-query and reference spans (mm_reg_set_coor).
struct Chain {
  const DnaSelectionItem* item = nullptr;
  std::vector<chaining::Anchor> anchors;
  int q_begin = 0;
  int q_end = 0;
  int ref_begin = 0;
  int ref_end = 0;
};

// A region mm_select_sub kept: an owner with its chain-level subsc and n_sub,
// or an attached chain.
struct Region {
  int item = 0;
  bool owner = false;
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
};

struct Record {
  AlignResult alignment;
  int region = 0;
  int contig = -1;
  bool inversion = false;
  int score = 0;
  int cnt = 0;
  int dp_max = 0;
  int dp_max2 = 0;
  // The record dp_max2 came from, -1 for none.
  int dp_max2_from = -1;
  double subsc = 0.0;
  int n_sub = 0;
  double pool_subsc = 0.0;
  int pool_n_sub = 0;
  int parent = -1;
  // Its rank in an exact dp_max tie (step 4).
  std::uint64_t tie = 0;
};

bool identical_chains(const Chain& a, const Chain& b) {
  return a.q_begin == b.q_begin && a.q_end == b.q_end &&
         a.item->contig == b.item->contig && a.ref_begin == b.ref_begin &&
         a.ref_end == b.ref_end;
}

// hit.c 255-281 with check_strand: every owner, and an attached chain within
// pri_ratio or 2k of its owner, not identical to it, or on the other strand
// above 0.8 * max_gap, at most best_n per read. The roles are in the parent
// walk's order, minimap2's chain order.
std::vector<Region> retain_regions(const DnaContext& context,
                                   const DnaKeptSelection& selection,
                                   const std::vector<Chain>& chains,
                                   int seed_length) {
  std::vector<Region> regions;
  const float pri_ratio = static_cast<float>(context.opts.pri_ratio);
  const int min_diff = 2 * seed_length;
  const int min_strand_sc =
      static_cast<int>(context.opts.cigar_dp_max_gap * 0.8);
  int n_2nd = 0;
  for (const DnaOwnershipRole& role : selection.roles) {
    if (role.item < 0 || static_cast<std::size_t>(role.item) >= chains.size() ||
        role.owner_item < 0 ||
        static_cast<std::size_t>(role.owner_item) >= chains.size())
      continue;
    const Chain& chain = chains[static_cast<std::size_t>(role.item)];
    if (chain.anchors.empty()) continue;
    if (role.owner_item == role.item) {
      regions.push_back({role.item, true, role.subsc, role.n_sub,
                         role.pool_subsc, role.pool_n_sub});
      continue;
    }
    const Chain& parent = chains[static_cast<std::size_t>(role.owner_item)];
    const int score = chain.item->score;
    const int parent_score = parent.item->score;
    if ((static_cast<float>(score) >=
             static_cast<float>(parent_score) * pri_ratio ||
         score + min_diff >= parent_score) &&
        n_2nd < kBestN) {
      if (!identical_chains(chain, parent)) {
        regions.push_back({role.item, false});
        ++n_2nd;
      }
    } else if (n_2nd < kBestN && score > min_strand_sc &&
               chain.item->reverse != parent.item->reverse) {
      regions.push_back({role.item, false});
      ++n_2nd;
    }
  }
  return regions;
}

// minimap2's cap at nearby seeds (align.c 706-768). minimap2 lays the
// regions' anchors out per contig and strand in order of their first anchor's
// reference position (lchain.c compact_a); a segment holding [first, end) of
// its region's anchors looks back from its first anchor and forward from its
// last, and the (min_cnt + 1)-th anchor lying beyond it on both axes caps the
// window at its distance.
class NearbySeeds {
 public:
  NearbySeeds(const std::vector<Region>& regions,
              const std::vector<Chain>& chains)
      : regions_(regions), chains_(chains), order_(regions.size()),
        position_(regions.size()) {
    for (std::size_t r = 0; r < regions.size(); ++r) order_[r] = r;
    std::sort(order_.begin(), order_.end(),
              [this](std::size_t left, std::size_t right) {
                return key(left) < key(right);
              });
    for (std::size_t p = 0; p < order_.size(); ++p) position_[order_[p]] = p;
  }

  ordered::TerminalWindowControl caps(std::size_t region, std::size_t first,
                                      int min_cnt) const {
    ordered::TerminalWindowControl out;
    const std::vector<chaining::Anchor>& own = anchors(region);
    const chaining::Anchor& head = own[first];
    const chaining::Anchor& tail = own.back();
    const int rs0 = head.r;
    const int qs0 = head.q;
    const int re0 = tail.r + tail.span;
    const int qe0 = tail.q + tail.span;
    int count = 0;
    const auto before = [&](const chaining::Anchor& anchor) {
      if (anchor.r >= rs0 || anchor.q >= qs0 || ++count <= min_cnt)
        return false;
      const int l = std::max(rs0 - anchor.r, qs0 - anchor.q);
      out.left_cap = true;
      out.left_cap_query = qs0 - l;
      out.left_cap_target = std::max(0, rs0 - l);
      return true;
    };
    bool found = false;
    for (std::size_t i = first; !found && i-- > 0;) found = before(own[i]);
    for (std::size_t p = position_[region]; !found && p-- > 0;) {
      if (!same_group(order_[p], region)) break;
      const std::vector<chaining::Anchor>& other = anchors(order_[p]);
      for (std::size_t i = other.size(); !found && i-- > 0;)
        found = before(other[i]);
    }
    count = 0;
    const auto after = [&](const chaining::Anchor& anchor) {
      const int x = anchor.r + anchor.span;
      const int y = anchor.q + anchor.span;
      if (x <= re0 || y <= qe0 || ++count <= min_cnt) return false;
      const int l = std::max(x - re0, y - qe0);
      out.right_cap = true;
      out.right_cap_query = qe0 + l;
      out.right_cap_target = re0 + l;
      return true;
    };
    found = false;
    for (std::size_t p = position_[region] + 1;
         !found && p < order_.size(); ++p) {
      if (!same_group(order_[p], region)) break;
      for (const chaining::Anchor& anchor : anchors(order_[p]))
        if ((found = after(anchor))) break;
    }
    return out;
  }

 private:
  const std::vector<chaining::Anchor>& anchors(std::size_t region) const {
    return chains_[static_cast<std::size_t>(regions_[region].item)].anchors;
  }
  std::tuple<int, bool, int, std::size_t> key(std::size_t region) const {
    const DnaSelectionItem& item =
        *chains_[static_cast<std::size_t>(regions_[region].item)].item;
    return {item.contig, item.reverse, anchors(region).front().r, region};
  }
  bool same_group(std::size_t a, std::size_t b) const {
    const DnaSelectionItem& x =
        *chains_[static_cast<std::size_t>(regions_[a].item)].item;
    const DnaSelectionItem& y =
        *chains_[static_cast<std::size_t>(regions_[b].item)].item;
    return x.contig == y.contig && x.reverse == y.reverse;
  }

  const std::vector<Region>& regions_;
  const std::vector<Chain>& chains_;
  std::vector<std::size_t> order_;
  std::vector<std::size_t> position_;
};

// mm_filter_regs (hit.c 301-320); max_clip_ratio is 1, which never fires.
// Returns the records dropped.
int filter_records(const DnaContext& context, std::vector<Record>& records) {
  const int min_cnt = context.opts.cigar_dp_split_min_anchors;
  const std::size_t before = records.size();
  records.erase(
      std::remove_if(records.begin(), records.end(),
                     [&](const Record& record) {
                       return (!record.inversion && record.cnt < min_cnt) ||
                              record.alignment.matches <
                                  context.opts.min_chain_score ||
                              record.dp_max < context.opts.cigar_dp_min_dp_max;
                     }),
      records.end());
  return static_cast<int>(before - records.size());
}

// mm_update_dp_max (align.c 1022-1046).
void update_dp_max(const DnaContext& context,
                   const DnaFamilyRealizationRequest& request, int read_length,
                   std::vector<Record>& records) {
  if (records.size() < 2) return;
  int max = -1;
  int max2 = -1;
  int max_i = -1;
  for (std::size_t i = 0; i < records.size(); ++i) {
    const int dp = records[i].dp_max;
    if (dp > max) {
      max2 = max;
      max = dp;
      max_i = static_cast<int>(i);
    } else if (dp > max2) {
      max2 = dp;
    }
  }
  if (max_i < 0 || max < 0 || max2 < 0) return;
  const AlignResult& best = records[static_cast<std::size_t>(max_i)].alignment;
  if (best.query_end - best.query_start <
      static_cast<double>(read_length) * static_cast<double>(kRankFrac))
    return;
  if (max2 < static_cast<double>(max) * static_cast<double>(kRankFrac)) return;
  double identity = -1.0;
  dna_record_event_identity(context, best, *request.forward_query,
                            *request.reverse_query, identity);
  const double b2 = dna_rank_b2(identity, context.opts.cigar_dp_match,
                                context.opts.cigar_dp_mismatch);
  for (Record& record : records)
    record.dp_max = std::max(
        0, dna_record_recal_dp_max(context, record.alignment,
                                   *request.forward_query,
                                   *request.reverse_query, b2));
}

// mm_set_parent (hit.c 125-186) on the aligned forward-read spans: parent,
// subsc, n_sub and dp_max2. A record's subsc carries over from its chain, and
// the first record keeps its chain's n_sub, as in minimap2.
void set_parents(const DnaContext& context,
                 const std::vector<Chain>& chains,
                 const std::vector<Region>& regions,
                 std::vector<Record>& records) {
  const int sub_diff =
      2 * context.opts.cigar_dp_match + context.opts.cigar_dp_mismatch;
  const auto pool_of = [&](const Record& record) {
    const DnaSelectionItem& item =
        *chains[static_cast<std::size_t>(
                    regions[static_cast<std::size_t>(record.region)].item)]
             .item;
    return !record.inversion ? item.candidate : -1;
  };
  if (records.empty()) return;
  std::vector<std::size_t> w{0};
  std::vector<std::pair<int, int>> cover;
  records[0].parent = 0;
  for (std::size_t i = 1; i < records.size(); ++i) {
    Record& ri = records[i];
    const int si = ri.alignment.query_start;
    const int ei = ri.alignment.query_end;
    int uncov_len = 0;
    cover.clear();
    for (const std::size_t j : w) {
      int sj = records[j].alignment.query_start;
      int ej = records[j].alignment.query_end;
      if (ej <= si || sj >= ei) continue;
      if (sj < si) sj = si;
      if (ej > ei) ej = ei;
      cover.emplace_back(sj, ej);
    }
    bool secondary = false;
    if (!cover.empty()) {
      std::sort(cover.begin(), cover.end());
      int x = si;
      for (const auto& [begin, end] : cover) {
        if (begin > x) uncov_len += begin - x;
        x = end > x ? end : x;
      }
      if (ei > x) uncov_len += ei - x;
      for (const std::size_t j : w) {
        Record& rp = records[j];
        const int sj = rp.alignment.query_start;
        const int ej = rp.alignment.query_end;
        if (ej <= si || sj >= ei) continue;
        const int min = std::min(ej - sj, ei - si);
        const int max = std::max(ej - sj, ei - si);
        const int ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj)
                               : (ej < si ? 0 : ej < ei ? ej - si : ei - si);
        if (static_cast<float>(ol) / static_cast<float>(min) -
                    static_cast<float>(uncov_len) / static_cast<float>(max) <=
                kMaskLevel ||
            uncov_len > kMaskLen)
          continue;
        ri.parent = rp.parent;
        const double sci = ri.score;
        rp.subsc = std::max(rp.subsc, sci);
        bool cnt_sub = ri.cnt >= rp.cnt;
        if (rp.contig != ri.contig || rp.alignment.pos != ri.alignment.pos ||
            rp.alignment.target_end != ri.alignment.target_end || ol != min) {
          if (ri.dp_max > rp.dp_max2) {
            rp.dp_max2 = ri.dp_max;
            rp.dp_max2_from = static_cast<int>(i);
          }
          if (rp.dp_max - ri.dp_max <= sub_diff) cnt_sub = true;
        }
        if (cnt_sub) ++rp.n_sub;
        const int pool = pool_of(rp);
        if (pool >= 0 && pool == pool_of(ri)) {
          rp.pool_subsc = std::max(rp.pool_subsc, sci);
          if (cnt_sub) ++rp.pool_n_sub;
        }
        secondary = true;
        break;
      }
    }
    if (!secondary) {
      w.push_back(i);
      ri.parent = static_cast<int>(i);
      ri.n_sub = 0;
      ri.pool_n_sub = 0;
    }
  }
}

// mm_select_sub (hit.c 255-281) without check_strand, on chain scores, then
// mm_sync_regs.
void select_records(const DnaContext& context, int seed_length,
                    std::vector<Record>& records) {
  const float pri_ratio = static_cast<float>(context.opts.pri_ratio);
  const int min_diff = 2 * seed_length;
  std::vector<int> index(records.size(), -1);
  std::vector<bool> keep(records.size(), false);
  int n_2nd = 0;
  for (std::size_t i = 0; i < records.size(); ++i) {
    const Record& ri = records[i];
    const Record& rp = records[static_cast<std::size_t>(ri.parent)];
    if (ri.parent == static_cast<int>(i) || ri.inversion) {
      keep[i] = true;
    } else if ((static_cast<float>(ri.score) >=
                    static_cast<float>(rp.score) * pri_ratio ||
                ri.score + min_diff >= rp.score) &&
               n_2nd < kBestN) {
      const AlignResult& a = ri.alignment;
      const AlignResult& b = rp.alignment;
      if (!(a.query_start == b.query_start && a.query_end == b.query_end &&
            ri.contig == rp.contig && a.pos == b.pos &&
            a.target_end == b.target_end)) {
        keep[i] = true;
        ++n_2nd;
      }
    }
  }
  std::vector<Record> kept;
  for (std::size_t i = 0; i < records.size(); ++i) {
    if (!keep[i]) continue;
    index[i] = static_cast<int>(kept.size());
    kept.push_back(std::move(records[i]));
  }
  for (Record& record : kept) {
    record.parent = index[static_cast<std::size_t>(record.parent)];
    if (record.dp_max2_from >= 0)
      record.dp_max2_from = index[static_cast<std::size_t>(record.dp_max2_from)];
  }
  records = std::move(kept);
}

}  // namespace

DnaRegionOutcome realize_dna_regions(
    const DnaContext& context, const DnaFamilyRealizationRequest& request) {
  DnaRegionOutcome out;
  DnaFamilyRealizationOutcome& result = out.family;
  if (request.placement == nullptr || !request.placement->accepted) {
    result.failure = request.placement != nullptr &&
                             request.placement->no_owner_chain
                         ? DnaFamilyFailure::NoOwnerChain
                         : DnaFamilyFailure::PlacementChainRefused;
    return out;
  }
  if (request.family == nullptr || request.forward_query == nullptr ||
      request.reverse_query == nullptr || context.ref.index == nullptr ||
      context.ref.names == nullptr || context.ref.encoded == nullptr ||
      !request.family->valid || request.family->read_length <= 0 ||
      static_cast<int>(request.forward_query->size()) !=
          request.family->read_length ||
      static_cast<int>(request.reverse_query->size()) !=
          request.family->read_length)
    return out;
  const DnaPlacementChainingResult& placement = *request.placement;
  const DnaKeptSelection& selection = placement.kept_selection;
  const int read_length = request.family->read_length;
  const int seed_length = request.family->seed_length;

  std::vector<Chain> chains(selection.items.size());
  for (std::size_t i = 0; i < selection.items.size(); ++i) {
    Chain& chain = chains[i];
    chain.item = &selection.items[i];
    const std::vector<chaining::Anchor>* anchors =
        placement.selection_anchors(selection.items[i]);
    if (anchors == nullptr || anchors->empty() || chain.item->contig < 0 ||
        static_cast<std::size_t>(chain.item->contig) >=
            context.ref.encoded->size())
      continue;
    chain.anchors = *anchors;
    family_internal::deduplicate_exact_anchors(chain.anchors);
    const chaining::Anchor& first = chain.anchors.front();
    const chaining::Anchor& last = chain.anchors.back();
    chain.q_begin = chain.item->reverse ? read_length - last.q - last.span
                                        : first.q;
    chain.q_end = chain.item->reverse ? read_length - first.q
                                      : last.q + last.span;
    chain.ref_begin = first.r;
    chain.ref_end = last.r + last.span;
  }

  // 1. Retain.
  const std::vector<Region> regions =
      retain_regions(context, selection, chains, seed_length);
  result.block_count = static_cast<int>(regions.size());
  if (regions.empty()) {
    result.failure = DnaFamilyFailure::NoSelectedFamily;
    return out;
  }

  // 2. Realize each region from its whole anchor vector. A refused segment
  // is dropped; the segments before it stand.
  const NearbySeeds nearby(regions, chains);
  const int min_cnt = context.opts.cigar_dp_split_min_anchors;
  const ordered::NormalizationControl normalization{
      std::max(1, context.opts.cigar_dp_bw), context.opts.min_chain_score,
      std::max(1, context.opts.cigar_dp_max_gap)};
  DnaFamilyFailure first_failure = DnaFamilyFailure::None;
  std::vector<Record> records;
  for (std::size_t r = 0; r < regions.size(); ++r) {
    const Region& region = regions[r];
    const Chain& chain = chains[static_cast<std::size_t>(region.item)];
    const DnaSelectionItem& item = *chain.item;
    const std::vector<std::uint8_t>& reference =
        (*context.ref.encoded)[static_cast<std::size_t>(item.contig)];
    const auto refuse = [&](DnaFamilyFailure failure) {
      if (first_failure == DnaFamilyFailure::None) {
        first_failure = failure;
        result.failed_block = static_cast<int>(r);
      }
    };
    BlockPlan current;
    current.chromosome = item.contig;
    current.reverse = item.reverse;
    current.split_group = static_cast<int>(r);
    std::vector<chaining::Anchor> anchors = chain.anchors;
    int score = item.score;
    int cnt = static_cast<int>(anchors.size());
    bool previous_committed = false;
    Unit previous_unit;
    for (;;) {
      ordered::OrderedAnchorPath raw;
      raw.reference_id = item.contig;
      raw.reverse_complemented = item.reverse;
      raw.query_length = read_length;
      raw.target_length = static_cast<int>(reference.size());
      raw.query_bound_begin = 0;
      raw.query_bound_end = read_length;
      raw.target_bound_begin = 0;
      raw.target_bound_end = static_cast<int>(reference.size());
      raw.minimizer_k = context.ref.index->k();
      raw.selected = std::move(anchors);
      raw.realization_end = raw.selected.size();
      const ordered::NormalizedPath normalized =
          ordered::normalize_for_realization(raw, normalization);
      if (!normalized) {
        refuse(DnaFamilyFailure::InvalidOrderedPath);
        break;
      }
      current.path = normalized.path;
      BlockSplit split;
      DnaFamilyFailure failure = DnaFamilyFailure::None;
      if (!family_internal::plan_and_realize_packets(context, request, current,
                                                     failure, result, split)) {
        refuse(failure);
        break;
      }
      ++result.segment_count;
      // mm_split_reg (hit.c 105-123): the rest takes its share of the score.
      int segment_score = score;
      int segment_cnt = cnt;
      int next_score = 0;
      int next_cnt = 0;
      if (split.has_continuation) {
        next_cnt =
            static_cast<int>(split.continuation.continuation_anchors.size());
        next_score = static_cast<int>(
            static_cast<float>(score) *
                (static_cast<float>(next_cnt) / static_cast<float>(cnt)) +
            .499);
        segment_cnt = cnt - next_cnt;
        segment_score = score - next_score;
      }
      if (split.truncated) ++result.zdrop_truncations;
      bool committed = false;
      Unit unit;
      // A segment cut at its start realized nothing (mm_filter_regs).
      if (!(split.truncated &&
            current.query_end_cut == current.query_begin_cut)) {
        unit = family_internal::unit_from_block(current);
        const ordered::TerminalWindowControl caps = nearby.caps(
            r, chain.anchors.size() - current.path.selected.size(), min_cnt);
        AlignResult alignment;
        if (!family_internal::extend_unit_terminals(context, request, current,
                                                    unit, result, caps)) {
          refuse(DnaFamilyFailure::TerminalRefused);
        } else if (!family_internal::validate_and_commit_unit(
                       context, request, unit, alignment)) {
          refuse(DnaFamilyFailure::FamilyValidation);
        } else {
          // mm_align1_inv between two pieces of a chain-level primary
          // (align.c 1102-1108, 925-926).
          if (region.owner && current.split_continuation &&
              current.split_inv && previous_committed) {
            Record middle;
            if (family_internal::realize_inversion_middle(
                    context, request, previous_unit, unit, result,
                    middle.alignment)) {
              middle.region = static_cast<int>(r);
              middle.contig = item.contig;
              middle.inversion = true;
              records.push_back(std::move(middle));
            }
          }
          Record record;
          record.alignment = std::move(alignment);
          record.alignment.chain_anchors = segment_cnt;
          record.alignment.chain_score = segment_score;
          record.region = static_cast<int>(r);
          record.contig = item.contig;
          record.score = segment_score;
          record.cnt = segment_cnt;
          record.subsc = region.subsc;
          record.n_sub = region.n_sub;
          record.pool_subsc = region.pool_subsc;
          record.pool_n_sub = region.pool_n_sub;
          records.push_back(std::move(record));
          committed = true;
        }
      }
      previous_committed = committed;
      previous_unit = std::move(unit);
      if (!split.has_continuation) break;
      ++result.zdrop_splits;
      current = std::move(split.continuation);
      anchors = std::move(current.continuation_anchors);
      current.continuation_anchors.clear();
      score = next_score;
      cnt = next_cnt;
    }
  }
  result.unit_count = static_cast<int>(records.size());

  // 3. Filter, and for a read of rank_min_len or more rescale and filter
  // again (align.c 1113-1117).
  for (Record& record : records)
    record.dp_max = dna_record_dp_max_segment(
        context, record.alignment, *request.forward_query,
        *request.reverse_query);
  filter_records(context, records);
  if (read_length >= kRankMinLen) {
    update_dp_max(context, request, read_length, records);
    filter_records(context, records);
  }
  if (records.empty()) {
    out.floored = result.unit_count > 0;
    result.failure = out.floored ? DnaFamilyFailure::None
                     : first_failure != DnaFamilyFailure::None
                         ? first_failure
                         : DnaFamilyFailure::NoAnchors;
    return out;
  }

  // 4. Re-select on the aligned spans, in dp_max order (mm_hit_sort). An
  // exact tie goes to the read-seeded hash of the record's contig, strand
  // and reference start (seeding/tie_hash.h), then to contig and start.
  for (Record& record : records)
    record.tie = tie_locus_hash(context.vote_tie_seed, record.contig,
                                record.alignment.is_reverse,
                                record.alignment.pos);
  std::stable_sort(records.begin(), records.end(),
                   [](const Record& left, const Record& right) {
                     if (left.dp_max != right.dp_max)
                       return left.dp_max > right.dp_max;
                     return std::tie(left.tie, left.contig,
                                     left.alignment.pos) <
                            std::tie(right.tie, right.contig,
                                     right.alignment.pos);
                   });
  set_parents(context, chains, regions, records);
  select_records(context, seed_length, records);

  // The first primary is the SAM primary (mm_set_sam_pri), the other
  // primaries its supplementaries, each secondary a hypothesis of its own.
  dna::Result& output = result.output;
  bool head = true;
  for (std::size_t i = 0; i < records.size(); ++i) {
    Record& record = records[i];
    const Chain& chain = chains[static_cast<std::size_t>(
        regions[static_cast<std::size_t>(record.region)].item)];
    const DnaSelectionItem& item = *chain.item;
    if (record.parent != static_cast<int>(i)) {
      record.alignment.mapq = 0;
      output.secondary.push_back(std::move(record.alignment));
      continue;
    }
    DnaRegionPrice price;
    price.candidate = record.inversion ? -1 : item.candidate;
    price.inversion = record.inversion;
    price.score = record.score;
    price.cnt = record.cnt;
    price.score0 = item.score;
    price.anchors0 = static_cast<int>(chain.anchors.size());
    price.chain_q_begin = chain.q_begin;
    price.chain_q_end = chain.q_end;
    price.chain_ref_begin = chain.ref_begin;
    price.chain_ref_end = chain.ref_end;
    price.dp_max = record.dp_max;
    price.dp_max2 = record.dp_max2;
    if (record.dp_max2_from >= 0) {
      const Record& from = records[static_cast<std::size_t>(record.dp_max2_from)];
      const DnaSelectionItem& source =
          *chains[static_cast<std::size_t>(
                      regions[static_cast<std::size_t>(from.region)].item)]
               .item;
      if (!from.inversion)
        price.dp_max2_candidate = source.candidate;
    }
    price.subsc = record.subsc;
    price.n_sub = record.n_sub;
    price.pool_subsc = record.pool_subsc;
    price.pool_n_sub = record.pool_n_sub;
    const ::fa::cpu::voting::CandidateId candidate =
        price.candidate >= 0
            ? static_cast<::fa::cpu::voting::CandidateId>(price.candidate)
            : ::fa::cpu::voting::kNullCandidate;
    if (head) {
      static_cast<AlignResult&>(output) = std::move(record.alignment);
      result.primary_candidate = candidate;
      head = false;
    } else {
      output.supplementary.push_back(std::move(record.alignment));
      output.supplementary_candidates.push_back(candidate);
    }
    out.prices.push_back(price);
  }
  output.primary_part = -1;
  result.failure = DnaFamilyFailure::None;
  return out;
}

}  // namespace fa::cpu::lr
