// Ownership of the read on the dense chains' anchors. The chains are taken
// best score first. A chain whose query extent overlaps an owner's by more
// than half of the shorter, less the share of its own extent no owner covers,
// attaches to the first such owner (minimap2's mm_set_parent at its default
// mask level); otherwise it is a new owner. A chain sharing half the anchors of
// the smaller with an owner, or with a chain already attached to the owner it
// would attach to, is that alignment again and is dropped. The owners are then
// made disjoint: one containing an earlier owner gets a hole there, and two
// that overlap split at the boundary keeping the most anchors on each side.
// Each owner's stretches keep the anchors that start inside them, and the
// bound between neighbouring pieces moves so that every anchor stays whole.
// The blocks tile the read.
#include "chain_ownership.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace fa::cpu::lr {
namespace {

struct AnchorKey {
  int contig;
  int strand;
  std::int32_t q;
  std::int32_t r;
  bool operator==(const AnchorKey& other) const noexcept {
    return contig == other.contig && strand == other.strand && q == other.q &&
           r == other.r;
  }
};

struct AnchorKeyHash {
  std::size_t operator()(const AnchorKey& key) const noexcept {
    std::uint64_t h =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.q)) << 32) |
        static_cast<std::uint32_t>(key.r);
    h ^= ((static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.contig))
           << 1) |
          static_cast<std::uint64_t>(key.strand)) *
         0x9E3779B97F4A7C15ull;
    h ^= h >> 29;
    h *= 0xBF58476D1CE4E5B9ull;
    h ^= h >> 32;
    return static_cast<std::size_t>(h);
  }
};

// One item in the parent walk.
struct Unit {
  int item = 0;
  int n = 0;
  int score = 0;
  int fq0 = 0;
  int fq1 = 0;
  std::vector<int> starts;  // forward starts, ascending
};

// Forward k-mer start of an anchor.
int forward_start(const chaining::Anchor& anchor, bool reverse, int read_length,
                  int seed_length) {
  return reverse ? read_length - seed_length - anchor.q : anchor.q;
}

// p in [ob, oe] maximizing #(fa < p) + #(fb >= p); ties to the p nearest
// (ob + oe) / 2, then the lower p. fa and fb ascending.
int boundary(const std::vector<int>& fa, const std::vector<int>& fb, int ob,
             int oe) {
  std::vector<int> events{ob};
  for (const std::vector<int>* starts : {&fa, &fb}) {
    for (auto it = std::lower_bound(starts->begin(), starts->end(), ob);
         it != starts->end() && *it < oe; ++it)
      events.push_back(*it + 1);
  }
  std::sort(events.begin(), events.end());
  events.erase(std::unique(events.begin(), events.end()), events.end());
  const std::int64_t nb = static_cast<std::int64_t>(fb.size());
  const std::int64_t mid2 = static_cast<std::int64_t>(ob) + oe;
  bool have = false;
  std::tuple<std::int64_t, std::int64_t, std::int64_t> best;
  for (std::size_t index = 0; index < events.size(); ++index) {
    const std::int64_t x = events[index];
    const std::int64_t y = std::min<std::int64_t>(
        index + 1 < events.size() ? events[index + 1] - 1 : oe, oe);
    if (y < x) continue;
    const std::int64_t f =
        (std::lower_bound(fa.begin(), fa.end(), static_cast<int>(x)) -
         fa.begin()) +
        nb -
        (std::lower_bound(fb.begin(), fb.end(), static_cast<int>(x)) -
         fb.begin());
    std::int64_t p;
    if (2 * x >= mid2)
      p = x;
    else if (2 * y <= mid2)
      p = y;
    else
      p = mid2 / 2;
    const std::int64_t off = 2 * p - mid2;
    const auto key = std::make_tuple(-f, off < 0 ? -off : off, p);
    if (!have || key < best) {
      best = key;
      have = true;
    }
  }
  return static_cast<int>(std::get<2>(best));
}

}  // namespace

DnaOwnershipSelection select_chain_owners(
    int read_length, int seed_length, int min_chain_score,
    const std::vector<DnaOwnershipPath>& paths,
    std::vector<DnaOwnershipRole>& roles) {
  DnaOwnershipSelection out;
  const int L = read_length;
  const int k = seed_length;
  for (const DnaOwnershipPath& path : paths) {
    if (path.anchors == nullptr) continue;
    if (static_cast<int>(path.anchors->size()) >= kDnaOwnerMinAnchors &&
        path.score >= min_chain_score)
      out.items.push_back(path);
  }
  const int count = static_cast<int>(out.items.size());
  if (count == 0) return out;

  // Each item's forward query extent, anchors whole.
  std::vector<Unit> units(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const DnaOwnershipPath& item = out.items[static_cast<std::size_t>(i)];
    Unit& unit = units[static_cast<std::size_t>(i)];
    unit.item = i;
    unit.n = static_cast<int>(item.anchors->size());
    unit.score = item.score;
    bool first = true;
    for (const chaining::Anchor& anchor : *item.anchors) {
      const int begin = item.reverse ? L - anchor.q - anchor.span : anchor.q;
      const int end = item.reverse ? L - anchor.q : anchor.q + anchor.span;
      unit.fq0 = first ? begin : std::min(unit.fq0, begin);
      unit.fq1 = first ? end : std::max(unit.fq1, end);
      first = false;
    }
    unit.starts.reserve(static_cast<std::size_t>(unit.n));
    for (const chaining::Anchor& anchor : *item.anchors)
      unit.starts.push_back(forward_start(anchor, item.reverse, L, k));
    std::sort(unit.starts.begin(), unit.starts.end());
  }
  std::vector<int> pool(static_cast<std::size_t>(count));
  std::iota(pool.begin(), pool.end(), 0);
  std::sort(pool.begin(), pool.end(), [&](int left, int right) {
    const Unit& a = units[static_cast<std::size_t>(left)];
    const Unit& b = units[static_cast<std::size_t>(right)];
    const DnaOwnershipPath& pa = out.items[static_cast<std::size_t>(left)];
    const DnaOwnershipPath& pb = out.items[static_cast<std::size_t>(right)];
    return std::make_tuple(-a.score, -a.n, pa.candidate, pa.path) <
           std::make_tuple(-b.score, -b.n, pb.candidate, pb.path);
  });

  // The parent walk. pool_subsc and pool_n_sub count only the attached
  // chains of the owner's own candidate.
  struct Owner {
    int unit = 0;
    double pool_subsc = 0.0;
    int pool_n_sub = 0;
    // For `roles` only.
    double subsc = 0.0;
    int n_sub = 0;
    std::size_t role = 0;
  };
  std::vector<Owner> owners;
  std::vector<std::pair<int, int>> cover;
  // A holder is an owner or an attached chain; slots follow creation order.
  struct Holder {
    int unit = 0;
    int owner = 0;  // its own index, or its parent's
    bool is_owner = false;
  };
  std::vector<Holder> holder_of;
  std::unordered_map<AnchorKey, std::vector<int>, AnchorKeyHash> holders;
  std::vector<int> shared;
  std::vector<int> touched;
  const auto anchor_key = [&](const Unit& unit, int c) {
    const DnaOwnershipPath& item =
        out.items[static_cast<std::size_t>(unit.item)];
    const chaining::Anchor& anchor =
        (*item.anchors)[static_cast<std::size_t>(c)];
    return AnchorKey{item.contig, item.reverse ? 1 : 0, anchor.q, anchor.r};
  };
  const auto hold = [&](int u, int owner, bool is_owner) {
    const int slot = static_cast<int>(holder_of.size());
    holder_of.push_back({u, owner, is_owner});
    shared.push_back(0);
    for (int c = 0; c < units[static_cast<std::size_t>(u)].n; ++c)
      holders[anchor_key(units[static_cast<std::size_t>(u)], c)].push_back(
          slot);
  };
  // The first holder (lowest slot) passing `accept` with half the anchors of
  // the smaller shared, or -1.
  const auto first_holder = [&](const Unit& unit, auto&& accept) {
    int found = -1;
    for (const int slot : touched) {
      const Holder& holder = holder_of[static_cast<std::size_t>(slot)];
      const std::int64_t common = shared[static_cast<std::size_t>(slot)];
      if ((found < 0 || slot < found) && accept(holder) &&
          2 * common >=
              std::min(unit.n, units[static_cast<std::size_t>(holder.unit)].n))
        found = slot;
    }
    return found;
  };
  for (const int u : pool) {
    const Unit& unit = units[static_cast<std::size_t>(u)];
    touched.clear();
    for (int c = 0; c < unit.n; ++c) {
      const auto it = holders.find(anchor_key(unit, c));
      if (it == holders.end()) continue;
      for (const int slot : it->second)
        if (shared[static_cast<std::size_t>(slot)]++ == 0)
          touched.push_back(slot);
    }
    const auto settle = [&] {
      for (const int slot : touched) shared[static_cast<std::size_t>(slot)] = 0;
    };
    // The same alignment as an owner.
    if (first_holder(unit, [](const Holder& holder) {
          return holder.is_owner;
        }) >= 0) {
      settle();
      continue;
    }
    const int si = unit.fq0;
    const int ei = unit.fq1;
    const int len_i = ei - si;
    int parent = -1;
    cover.clear();
    for (const Owner& owner : owners) {
      int sj = units[static_cast<std::size_t>(owner.unit)].fq0;
      int ej = units[static_cast<std::size_t>(owner.unit)].fq1;
      if (ej <= si || sj >= ei) continue;
      if (sj < si) sj = si;
      if (ej > ei) ej = ei;
      cover.emplace_back(sj, ej);
    }
    if (!cover.empty()) {
      std::sort(cover.begin(), cover.end());
      int uncovered = 0;
      int x = si;
      for (const auto& [begin, end] : cover) {
        if (begin > x) uncovered += begin - x;
        x = end > x ? end : x;
      }
      if (ei > x) uncovered += ei - x;
      for (std::size_t j = 0; j < owners.size(); ++j) {
        const int sj = units[static_cast<std::size_t>(owners[j].unit)].fq0;
        const int ej = units[static_cast<std::size_t>(owners[j].unit)].fq1;
        if (ej <= si || sj >= ei) continue;
        const int len_p = ej - sj;
        const int mn = std::min(len_p, len_i);
        const int mx = std::max(len_p, len_i);
        const int ol = std::min(ei, ej) - std::max(si, sj);
        const float covered = static_cast<float>(ol) / static_cast<float>(mn);
        const float spill =
            static_cast<float>(uncovered) / static_cast<float>(mx);
        if (covered - spill > 0.5f) {
          parent = static_cast<int>(j);
          break;
        }
      }
    }
    if (parent < 0) {
      settle();
      owners.push_back({u});
      owners.back().role = roles.size();
      roles.push_back({unit.item, unit.item});
      hold(u, static_cast<int>(owners.size()) - 1, true);
      continue;
    }
    // The same alignment as a chain already attached to that owner.
    const int same_rival =
        first_holder(unit, [parent](const Holder& holder) {
          return !holder.is_owner && holder.owner == parent;
        });
    settle();
    if (same_rival >= 0) continue;
    hold(u, parent, false);
    Owner& owner = owners[static_cast<std::size_t>(parent)];
    const Unit& head = units[static_cast<std::size_t>(owner.unit)];
    if (out.items[static_cast<std::size_t>(unit.item)].candidate ==
        out.items[static_cast<std::size_t>(head.item)].candidate) {
      owner.pool_subsc =
          std::max(owner.pool_subsc, static_cast<double>(unit.score));
      if (unit.n >= head.n) ++owner.pool_n_sub;
    }
    owner.subsc = std::max(owner.subsc, static_cast<double>(unit.score));
    if (unit.n >= head.n) ++owner.n_sub;
    roles.push_back({unit.item, head.item});
  }
  for (const Owner& owner : owners) {
    DnaOwnershipRole& role = roles[owner.role];
    role.subsc = owner.subsc;
    role.n_sub = owner.n_sub;
    role.pool_subsc = owner.pool_subsc;
    role.pool_n_sub = owner.pool_n_sub;
  }

  // Disjoint owners: a later owner containing an earlier one gets a hole
  // there; a partial overlap is split at the boundary that keeps the most
  // anchors.
  const int owner_count = static_cast<int>(owners.size());
  std::vector<int> lo(owner_count), hi(owner_count);
  std::vector<std::vector<std::pair<int, int>>> holes(owner_count);
  // (earlier, later) -> the bound between them: the boundary, or -1 when the
  // later contains the earlier.
  std::map<std::pair<int, int>, int> relation;
  for (int o = 0; o < owner_count; ++o) {
    lo[o] = units[static_cast<std::size_t>(owners[o].unit)].fq0;
    hi[o] = units[static_cast<std::size_t>(owners[o].unit)].fq1;
  }
  for (int i = 0; i < owner_count; ++i) {
    const Unit& a = units[static_cast<std::size_t>(owners[i].unit)];
    for (int j = i + 1; j < owner_count; ++j) {
      const Unit& b = units[static_cast<std::size_t>(owners[j].unit)];
      if (std::min(a.fq1, b.fq1) - std::max(a.fq0, b.fq0) <= 0) continue;
      if (b.fq0 <= a.fq0 && a.fq1 <= b.fq1) {
        holes[j].emplace_back(a.fq0, a.fq1);
        relation[{i, j}] = -1;
        continue;
      }
      const int first = a.fq0 <= b.fq0 ? i : j;
      const int second = first == i ? j : i;
      const Unit& ua = units[static_cast<std::size_t>(owners[first].unit)];
      const Unit& ub = units[static_cast<std::size_t>(owners[second].unit)];
      const int p =
          boundary(ua.starts, ub.starts, ub.fq0, std::min(ua.fq1, ub.fq1));
      hi[first] = std::min(hi[first], p);
      lo[second] = std::max(lo[second], p);
      relation[{i, j}] = p;
    }
  }
  // A piece's score for `kept_anchors` of its chain's anchors, rounded as
  // minimap2 rounds a split chain's.
  const auto piece_score = [&](const Unit& unit, int kept_anchors) {
    if (kept_anchors == unit.n) return unit.score;
    const float share =
        static_cast<float>(kept_anchors) / static_cast<float>(unit.n);
    return static_cast<int>(unit.score * share + .499);
  };
  struct Held {
    int start = 0;
    int end = 0;
    int index = 0;  // in chain order
  };
  struct Piece {
    int owner = 0;
    std::vector<Held> held;  // forward order
  };
  std::vector<Piece> pieces;
  for (int o = 0; o < owner_count; ++o) {
    if (lo[o] >= hi[o]) continue;
    std::sort(holes[o].begin(), holes[o].end());
    std::vector<std::pair<int, int>> stretches;
    int cursor = lo[o];
    for (const auto& [h0, h1] : holes[o]) {
      if (h1 <= cursor) continue;
      if (h0 >= hi[o]) break;
      if (h0 > cursor) stretches.emplace_back(cursor, h0);
      cursor = std::max(cursor, h1);
    }
    if (cursor < hi[o]) stretches.emplace_back(cursor, hi[o]);
    const Unit& unit = units[static_cast<std::size_t>(owners[o].unit)];
    const DnaOwnershipPath& item =
        out.items[static_cast<std::size_t>(unit.item)];
    for (const auto& [plo, phi] : stretches) {
      Piece piece;
      for (int c = 0; c < unit.n; ++c) {
        const chaining::Anchor& anchor =
            (*item.anchors)[static_cast<std::size_t>(c)];
        const int start = forward_start(anchor, item.reverse, L, k);
        if (start < plo || start >= phi) continue;
        piece.held.push_back({start, start + anchor.span, c});
      }
      if (static_cast<int>(piece.held.size()) < kDnaOwnerMinAnchors) continue;
      if (item.reverse) std::reverse(piece.held.begin(), piece.held.end());
      piece.owner = o;
      pieces.push_back(std::move(piece));
    }
  }
  std::sort(pieces.begin(), pieces.end(),
            [](const Piece& left, const Piece& right) {
              return std::make_pair(left.held.front().start, left.owner) <
                     std::make_pair(right.held.front().start, right.owner);
            });
  if (pieces.empty()) return out;

  // Whole anchors. The bound between neighbours is their boundary, the
  // hole's edge when one owner contains the other, else the middle of the gap
  // between their anchors, moved between the end of the left piece's last
  // anchor and the start of the right piece's first. Where those two anchors
  // overlap, the piece of the later owner gives up the anchors that cross the
  // bound; a piece left under kDnaOwnerMinAnchors goes and the pass starts
  // again.
  const auto relation_bound = [&](const Piece& left, const Piece& right,
                                  int end_left, int start_right) {
    const int a = left.owner;
    const int b = right.owner;
    if (a != b) {
      const auto it = relation.find({std::min(a, b), std::max(a, b)});
      if (it != relation.end()) {
        if (it->second >= 0) return it->second;
        const int contained = std::min(a, b);
        return contained == b
                   ? units[static_cast<std::size_t>(owners[b].unit)].fq0
                   : units[static_cast<std::size_t>(owners[a].unit)].fq1;
      }
    }
    return static_cast<int>(
        (static_cast<std::int64_t>(end_left) + start_right) / 2);
  };
  const std::size_t piece_count = pieces.size();
  std::vector<int> first(piece_count, 0), last(piece_count);
  for (std::size_t p = 0; p < piece_count; ++p)
    last[p] = static_cast<int>(pieces[p].held.size());
  std::vector<int> live(piece_count);
  std::iota(live.begin(), live.end(), 0);
  std::vector<int> bounds;
  for (bool again = true; again;) {
    again = false;
    bounds.assign(live.size() + 1, 0);
    bounds.back() = L;
    for (std::size_t j = 0; j + 1 < live.size(); ++j) {
      const int l = live[j];
      const int r = live[j + 1];
      const std::vector<Held>& left = pieces[static_cast<std::size_t>(l)].held;
      const std::vector<Held>& right =
          pieces[static_cast<std::size_t>(r)].held;
      const int end_left = left[static_cast<std::size_t>(last[l] - 1)].end;
      const int start_right = right[static_cast<std::size_t>(first[r])].start;
      int at;
      if (end_left <= start_right) {
        at = std::min(
            std::max(relation_bound(pieces[static_cast<std::size_t>(l)],
                                    pieces[static_cast<std::size_t>(r)],
                                    end_left, start_right),
                     end_left),
            start_right);
      } else if (pieces[static_cast<std::size_t>(l)].owner <
                 pieces[static_cast<std::size_t>(r)].owner) {
        at = end_left;
        while (first[r] < last[r] &&
               right[static_cast<std::size_t>(first[r])].start < at)
          ++first[r];
      } else {
        at = start_right;
        while (last[l] > first[l] &&
               left[static_cast<std::size_t>(last[l] - 1)].end > at)
          --last[l];
      }
      const bool drop_left = last[l] - first[l] < kDnaOwnerMinAnchors;
      if (drop_left || last[r] - first[r] < kDnaOwnerMinAnchors) {
        live.erase(live.begin() +
                   static_cast<std::ptrdiff_t>(drop_left ? j : j + 1));
        again = true;
        break;
      }
      bounds[j + 1] = at;
    }
  }
  for (std::size_t j = 0; j < live.size(); ++j) {
    const int p = live[j];
    const Piece& piece = pieces[static_cast<std::size_t>(p)];
    const Owner& owner = owners[static_cast<std::size_t>(piece.owner)];
    const Unit& unit = units[static_cast<std::size_t>(owner.unit)];
    DnaOwnershipBlock block;
    block.item = unit.item;
    block.a0 = std::numeric_limits<int>::max();
    block.a1 = 0;
    for (int h = first[p]; h < last[p]; ++h) {
      const int index = piece.held[static_cast<std::size_t>(h)].index;
      block.a0 = std::min(block.a0, index);
      block.a1 = std::max(block.a1, index + 1);
    }
    block.forward_begin = bounds[j];
    block.forward_end = bounds[j + 1];
    block.kept = last[p] - first[p];
    block.score = piece_score(unit, block.kept);
    block.pool_subsc = owner.pool_subsc;
    block.pool_n_sub = owner.pool_n_sub;
    out.blocks.push_back(block);
  }
  return out;
}

}  // namespace fa::cpu::lr
