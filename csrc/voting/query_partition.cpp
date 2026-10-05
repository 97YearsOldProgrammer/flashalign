#include "query_partition.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fa::cpu::voting {
namespace {

// The support count saturates: once a block has two supporting tiles it may close, and
// more change nothing. This, not a bound on block length, keeps the state space finite.
constexpr int kSupportStates = 3; // 0, 1, and "at least the minimum"

// Every comparator key before the lexicographic tail is a sum over the same left-to-right
// steps, so a prefix is described by this key plus its end state.
struct Key {
  std::int64_t score = 0;
  std::int32_t non_null_blocks = 0;   // minimized
  std::int32_t candidate_changes = 0; // minimized
  std::int64_t chain_evidence = 0;    // maximized
  std::int64_t vote_evidence = 0;     // maximized
  std::int64_t vote_rank_sum = 0;     // minimized
};

int compare_key(const Key& left, const Key& right) noexcept {
  if (left.score != right.score)
    return left.score > right.score ? 1 : -1;
  if (left.non_null_blocks != right.non_null_blocks)
    return left.non_null_blocks < right.non_null_blocks ? 1 : -1;
  if (left.candidate_changes != right.candidate_changes)
    return left.candidate_changes < right.candidate_changes ? 1 : -1;
  if (left.chain_evidence != right.chain_evidence)
    return left.chain_evidence > right.chain_evidence ? 1 : -1;
  if (left.vote_evidence != right.vote_evidence)
    return left.vote_evidence > right.vote_evidence ? 1 : -1;
  if (left.vote_rank_sum != right.vote_rank_sum)
    return left.vote_rank_sum < right.vote_rank_sum ? 1 : -1;
  return 0;
}

// Deterministic tail: fewer/lexicographically smaller block identities first,
// then lexicographically smaller block end positions. Together those two
// sequences determine the assignment, so no further tie is reachable.
int compare_tail(const std::vector<CandidateId>& left,
                 const std::vector<CandidateId>& right) noexcept {
  const std::size_t size = left.size();
  // Block identities, in order, run-length compressed.
  std::size_t li = 0;
  std::size_t ri = 0;
  while (li < size && ri < size) {
    if (left[li] != right[ri])
      return left[li] < right[ri] ? 1 : -1;
    std::size_t ln = li + 1;
    while (ln < size && left[ln] == left[li])
      ++ln;
    std::size_t rn = ri + 1;
    while (rn < size && right[rn] == right[ri])
      ++rn;
    li = ln;
    ri = rn;
  }
  if (li != ri)
    return li < ri ? -1 : 1; // shorter id list sorts first
  // Identical identity sequences: compare block end positions.
  li = 0;
  ri = 0;
  while (li < size) {
    std::size_t ln = li + 1;
    while (ln < size && left[ln] == left[li])
      ++ln;
    std::size_t rn = ri + 1;
    while (rn < size && right[rn] == right[ri])
      ++rn;
    if (ln != rn)
      return ln < rn ? 1 : -1;
    li = ln;
    ri = rn;
  }
  return 0;
}

struct Provenance {
  std::int16_t state = -1;
  std::int8_t slot = -1;
};

// Key first, so the cell packs into 48 bytes. An absent cell is never read (push writes all
// members first, and an absent cell's backpointer is never followed), so the per-tile reset
// clears `present` alone.
struct Cell {
  Key key;
  Provenance from;
  bool present = false;
};

// Per-thread scratch for one solve, kept at capacity across solves. Every element is
// rewritten before it is read.
struct SolverScratch {
  std::vector<Cell> current;
  std::vector<Cell> next;
  std::vector<Provenance> back;
  std::vector<CandidateId> tail_left;
  std::vector<CandidateId> tail_right;
  std::vector<CandidateId> assignment;
  // The best assignment so far of the WinnerOnly finish.
  std::vector<CandidateId> best_assignment;
  // candidate_at(state) for every state of the current solve.
  std::vector<CandidateId> state_candidate;
};

SolverScratch& solver_scratch() {
  static thread_local SolverScratch scratch;
  return scratch;
}

// The relaxation helpers run on every edge; force inlining, which gcc otherwise skips for
// the larger ones.
#if defined(__GNUC__) || defined(__clang__)
#define FA_QP_INLINE inline __attribute__((always_inline))
#else
#define FA_QP_INLINE inline
#endif

class Solver {
public:
  explicit Solver(const QueryPartitionProblem& problem)
      : p_(problem), scratch_(solver_scratch()), current_(scratch_.current),
        next_(scratch_.next), back_(scratch_.back),
        tail_left_(scratch_.tail_left), tail_right_(scratch_.tail_right),
        state_candidate_(scratch_.state_candidate), tiles_(problem.tile_count),
        candidates_(static_cast<int>(problem.catalogue.candidates.size())),
        resume_cost_(problem.parameters.same_candidate_resume_cost),
        predecessors_(resume_cost_ >= 0 ? candidates_ + 1 : 1),
        owner_states_(1 + kSupportStates * candidates_),
        states_(owner_states_ * predecessors_),
        slots_(problem.rival == QueryPartitionRival::ExactSecondSolution ? 2
                                                                         : 1) {}

  QueryPartitionResult run();

private:
  // The WinnerOnly answer: the best final cell under compare_path. A cell's
  // key is its path's additive key, which evaluate() recomputes, so only cells
  // tied with the best on the key rebuild their assignments for the tail; the
  // winner alone is evaluated.
  QueryPartitionPath best_final_path();

  static constexpr int kNullOwnerState = 0;
  // "No preceding run", which is also what a preceding NULL run leaves behind:
  // neither can be resumed, so they need not be distinguished.
  static constexpr int kNoPredecessor = 0;

  FA_QP_INLINE bool resuming() const noexcept { return resume_cost_ >= 0; }

  // The predecessor slot a run owned by `candidate_index` leaves behind. With resumption
  // disabled there is one slot and the dimension collapses.
  FA_QP_INLINE int predecessor_of(int candidate_index) const noexcept {
    return resuming() ? candidate_index + 1 : kNoPredecessor;
  }
  FA_QP_INLINE int owner_state_of(int candidate_index,
                                  int support) const noexcept {
    return 1 + kSupportStates * candidate_index + support;
  }
  FA_QP_INLINE int state_of(int candidate_index, int support,
                            int predecessor) const noexcept {
    return owner_state_of(candidate_index, support) * predecessors_ +
           predecessor;
  }
  FA_QP_INLINE int null_state(int predecessor) const noexcept {
    return kNullOwnerState * predecessors_ + predecessor;
  }
  FA_QP_INLINE int predecessor_at(int state) const noexcept {
    return state % predecessors_;
  }
  // The candidate of a state; tabulated once per solve (candidate_at), since the tie-break
  // walk reads it on every step.
  CandidateId decode_candidate_at(int state) const noexcept {
    const int owner_state = state / predecessors_;
    if (owner_state == kNullOwnerState)
      return kNullCandidate;
    return p_.catalogue
        .candidates[static_cast<std::size_t>((owner_state - 1) /
                                             kSupportStates)]
        .id;
  }
  FA_QP_INLINE CandidateId candidate_at(int state) const noexcept {
    return state_candidate_[static_cast<std::size_t>(state)];
  }
  FA_QP_INLINE std::size_t cell_index(int state, int slot) const noexcept {
    return static_cast<std::size_t>(state) * static_cast<std::size_t>(slots_) +
           static_cast<std::size_t>(slot);
  }

  // Rebuild the tile assignment of the prefix that ends at layer `layer` in
  // (state, slot). Layer L means tiles [0, L) have been assigned.
  void reconstruct(int layer, int state, int slot,
                   std::vector<CandidateId>& out) const {
    out.assign(static_cast<std::size_t>(layer), kNullCandidate);
    for (int tile = layer - 1; tile >= 0; --tile) {
      out[static_cast<std::size_t>(tile)] = candidate_at(state);
      const Provenance from =
          back_[static_cast<std::size_t>(tile) *
                    static_cast<std::size_t>(states_ * slots_) +
                cell_index(state, slot)];
      state = from.state;
      slot = from.slot;
    }
  }

  FA_QP_INLINE Provenance predecessor(int tile,
                                      Provenance at) const noexcept {
    return back_[static_cast<std::size_t>(tile) *
                     static_cast<std::size_t>(states_ * slots_) +
                 cell_index(at.state, at.slot)];
  }

  // Order two prefixes that reach the same layer, optionally each extended by
  // one more tile owned by `append`. Only called when the additive key ties.
  //
  // Walks both backpointer chains in lockstep and stops where they meet: before that point
  // the prefixes are the same assignment and cannot decide the comparison. The shared tile
  // at the meeting point is kept, since the run straddling it may end at different
  // positions. Cost is the length of the divergence.
  int compare_prefix(int layer, Provenance left, Provenance right,
                     CandidateId append, bool extend) {
    ++tail_comparisons_;
    tail_left_.clear();
    tail_right_.clear();
    int at = layer;
    while (at > 0 && !(left.state == right.state && left.slot == right.slot)) {
      tail_left_.push_back(candidate_at(left.state));
      tail_right_.push_back(candidate_at(right.state));
      left = predecessor(at - 1, left);
      right = predecessor(at - 1, right);
      --at;
    }
    if (at > 0) {
      const CandidateId shared = candidate_at(left.state);
      tail_left_.push_back(shared);
      tail_right_.push_back(shared);
    }
    std::reverse(tail_left_.begin(), tail_left_.end());
    std::reverse(tail_right_.begin(), tail_right_.end());
    if (extend) {
      tail_left_.push_back(append);
      tail_right_.push_back(append);
    }
    tail_steps_ += tail_left_.size();
    return compare_tail(tail_left_, tail_right_);
  }

  FA_QP_INLINE void push(int layer, int to_state, const Key& key,
                         Provenance from, CandidateId to_candidate) {
    ++transitions_;
    const std::size_t base = cell_index(to_state, 0);
    for (int slot = 0; slot < slots_; ++slot) {
      Cell& held = next_[base + static_cast<std::size_t>(slot)];
      if (!held.present) {
        held = {key, from, true};
        return;
      }
      int order = compare_key(key, held.key);
      if (order == 0)
        order = compare_prefix(layer, from, held.from, to_candidate, true);
      if (order > 0) {
        // Shift the weaker entries down; distinct (source cell) provenances
        // are distinct assignments, so no duplicate can be introduced.
        Cell carry = {key, from, true};
        for (int at = slot; at < slots_; ++at)
          std::swap(carry, next_[base + static_cast<std::size_t>(at)]);
        return;
      }
    }
  }

  const QueryPartitionProblem& p_;
  SolverScratch& scratch_;
  std::vector<Cell>& current_;
  std::vector<Cell>& next_;
  std::vector<Provenance>& back_;
  std::vector<CandidateId>& tail_left_;
  std::vector<CandidateId>& tail_right_;
  std::vector<CandidateId>& state_candidate_;
  int tiles_;
  int candidates_;
  int resume_cost_;
  int predecessors_;
  int owner_states_;
  int states_;
  int slots_;
  std::uint64_t transitions_ = 0;
  std::uint64_t tail_comparisons_ = 0;
  std::uint64_t tail_steps_ = 0;
};

// Recomputes a finished assignment's score and counts directly from the problem.
QueryPartitionPath evaluate(const QueryPartitionProblem& p,
                            std::vector<CandidateId> assignment) {
  QueryPartitionPath path;
  path.assignment = std::move(assignment);
  CandidateId previous = kNullCandidate;
  // Candidate of the run before the current one. A null run leaves kNullCandidate
  // behind exactly as "no run yet" does, and neither can be resumed.
  CandidateId run_predecessor = kNullCandidate;
  const int resume_cost = p.parameters.same_candidate_resume_cost;
  for (int tile = 0; tile < p.tile_count; ++tile) {
    const CandidateId id = path.assignment[static_cast<std::size_t>(tile)];
    if (id == kNullCandidate) {
      if (p.valid_tiles.test(tile))
        path.score -= p.parameters.null_tile_cost;
      if (previous != id)
        run_predecessor = previous;
      previous = id;
      continue;
    }
    const QueryCandidate& value =
        p.catalogue.candidates[static_cast<std::size_t>(id)];
    if (id != previous) {
      // Same two clauses the DP charges: an A-B-A resumption of this candidate,
      // otherwise a full open.
      path.score -= (resume_cost >= 0 && run_predecessor == id)
                        ? resume_cost
                        : p.parameters.block_open_cost;
    }
    if (value.support.test(tile))
      path.score += p.parameters.supported_tile_reward;
    else
      path.score -= p.parameters.unsupported_ownership_cost;
    path.chain_evidence += value.chain_evidence;
    path.vote_evidence += value.vote_evidence;
    path.vote_rank_sum += value.catalogue_rank;
    if (previous != kNullCandidate && previous != id)
      ++path.candidate_changes;
    if (previous != id)
      run_predecessor = previous;
    previous = id;
  }
  for (int begin = 0; begin < p.tile_count;) {
    const CandidateId id = path.assignment[static_cast<std::size_t>(begin)];
    int end = begin + 1;
    while (end < p.tile_count &&
           path.assignment[static_cast<std::size_t>(end)] == id)
      ++end;
    QueryBlock block{id, begin, end, 0};
    if (id != kNullCandidate) {
      ++path.non_null_blocks;
      const QueryCandidate& value =
          p.catalogue.candidates[static_cast<std::size_t>(id)];
      for (int tile = begin; tile < end; ++tile)
        if (p.valid_tiles.test(tile))
          if (value.support.test(tile))
            ++block.supporting_tiles;
    }
    path.blocks.push_back(block);
    begin = end;
  }
  return path;
}

// Ranking used for the final answer and for the second solution. Mirrors the
// in-DP order, with the lexicographic tail resolved on the finished paths.
int compare_path(const QueryPartitionProblem&, const QueryPartitionPath& left,
                 const QueryPartitionPath& right) {
  Key lk{left.score,          left.non_null_blocks, left.candidate_changes,
         left.chain_evidence, left.vote_evidence,   left.vote_rank_sum};
  Key rk{right.score,          right.non_null_blocks, right.candidate_changes,
         right.chain_evidence, right.vote_evidence,   right.vote_rank_sum};
  const int order = compare_key(lk, rk);
  if (order != 0)
    return order;
  return compare_tail(left.assignment, right.assignment);
}

QueryPartitionResult Solver::run() {
  const std::size_t layer = static_cast<std::size_t>(states_ * slots_);
  current_.assign(layer, Cell{});
  next_.assign(layer, Cell{});
  // Every slot is written at the end of its tile before reconstruct() reads it, so the
  // table is sized, not cleared.
  back_.resize(layer * static_cast<std::size_t>(tiles_));
  state_candidate_.resize(static_cast<std::size_t>(states_));
  for (int state = 0; state < states_; ++state)
    state_candidate_[static_cast<std::size_t>(state)] =
        decode_candidate_at(state);
  current_[cell_index(null_state(kNoPredecessor), 0)].present = true;

  // Top entries among closable open blocks, used to relax every
  // "close one candidate, open another" edge without enumerating pairs. Four
  // entries are enough: at most `slots_` of them can belong to the excluded
  // candidate, so at least `slots_` survive the exclusion.
  struct Ranked {
    Key key;
    Provenance from;
    int candidate_index = -1;
    int predecessor = kNoPredecessor;
  };
  using RankedList = std::array<Ranked, 4>;
  RankedList closable{};
  int closable_size = 0;
  // With resumption on, a close/open edge lands in a state that remembers which candidate
  // closed, so the ranking is kept per closing candidate; empty otherwise.
  std::vector<RankedList> owner_closable;
  std::vector<int> owner_closable_size;
  if (resuming()) {
    owner_closable.assign(static_cast<std::size_t>(candidates_), RankedList{});
    owner_closable_size.assign(static_cast<std::size_t>(candidates_), 0);
  }

  for (int tile = 0; tile < tiles_; ++tile) {
    // An absent cell's key and provenance are never read (see Cell), so the
    // reset touches one byte per cell instead of rewriting the layer.
    for (Cell& cell : next_)
      cell.present = false;
    const bool valid = p_.valid_tiles.test(tile);

    closable_size = 0;
    std::fill(owner_closable_size.begin(), owner_closable_size.end(), 0);
    // Two closable blocks tied on the additive key always differ in candidate or end
    // positions, so the tail resolves here without an extension, in the order `push`
    // would apply after the shared next tile.
    const auto offer = [&](RankedList& list, int& size, const Key& key,
                           Provenance from, int candidate_index,
                           int predecessor) {
      const auto stronger = [&](const Key& lhs, Provenance lhs_from,
                                const Ranked& rhs) {
        const int order = compare_key(lhs, rhs.key);
        if (order != 0)
          return order > 0;
        return compare_prefix(tile, lhs_from, rhs.from, kNullCandidate, false) >
               0;
      };
      int at = size;
      while (at > 0 &&
             stronger(key, from, list[static_cast<std::size_t>(at) - 1]))
        --at;
      if (at >= static_cast<int>(list.size()))
        return;
      for (int move = std::min<int>(size, static_cast<int>(list.size()) - 1);
           move > at; --move)
        list[static_cast<std::size_t>(move)] =
            list[static_cast<std::size_t>(move) - 1];
      list[static_cast<std::size_t>(at)] = {key, from, candidate_index,
                                            predecessor};
      size = std::min<int>(size + 1, static_cast<int>(list.size()));
    };

    // Relax "stay null" and "close an open block" into the null state. A null
    // run keeps whatever run preceded it out of reach, so "stay null" carries
    // its predecessor slot forward unchanged.
    for (int predecessor = 0; predecessor < predecessors_; ++predecessor) {
      const int from_state = null_state(predecessor);
      for (int slot = 0; slot < slots_; ++slot) {
        const Cell& cell = current_[cell_index(from_state, slot)];
        if (!cell.present)
          continue;
        Key key = cell.key;
        if (valid)
          key.score -= p_.parameters.null_tile_cost;
        push(tile, from_state, key,
             Provenance{static_cast<std::int16_t>(from_state),
                        static_cast<std::int8_t>(slot)},
             kNullCandidate);
      }
    }
    for (int index = 0; index < candidates_; ++index) {
      // Closing this block makes it the predecessor of whatever comes next.
      const int to_null = null_state(predecessor_of(index));
      for (int predecessor = 0; predecessor < predecessors_; ++predecessor) {
        const int from_state = state_of(index, kSupportStates - 1, predecessor);
        for (int slot = 0; slot < slots_; ++slot) {
          const Cell& cell = current_[cell_index(from_state, slot)];
          if (!cell.present)
            continue;
          const Provenance from{static_cast<std::int16_t>(from_state),
                                static_cast<std::int8_t>(slot)};
          Key key = cell.key;
          if (valid)
            key.score -= p_.parameters.null_tile_cost;
          push(tile, to_null, key, from, kNullCandidate);
          if (!valid)
            continue;
          if (resuming())
            offer(owner_closable[static_cast<std::size_t>(index)],
                  owner_closable_size[static_cast<std::size_t>(index)],
                  cell.key, from, index, predecessor);
          else
            offer(closable, closable_size, cell.key, from, index, predecessor);
        }
      }
    }

    if (valid) {
      for (int index = 0; index < candidates_; ++index) {
        const QueryCandidate& value =
            p_.catalogue.candidates[static_cast<std::size_t>(index)];
        const bool supported = value.support.test(tile);
        const std::int64_t owned =
            supported
                ? static_cast<std::int64_t>(p_.parameters.supported_tile_reward)
                : -static_cast<std::int64_t>(
                      p_.parameters.unsupported_ownership_cost);
        const auto own = [&](Key key) {
          key.score += owned;
          key.chain_evidence += value.chain_evidence;
          key.vote_evidence += value.vote_evidence;
          key.vote_rank_sum += value.catalogue_rank;
          return key;
        };
        const int opened_support = supported ? 1 : 0;
        // The slot a resumption of THIS candidate must have come from.
        const int resume_predecessor = predecessor_of(index);

        // Continue the block that is already open on this candidate.
        for (int support = 0; support < kSupportStates; ++support) {
          const int to_support =
              std::min(kSupportStates - 1, support + (supported ? 1 : 0));
          for (int predecessor = 0; predecessor < predecessors_;
               ++predecessor) {
            const int from_state = state_of(index, support, predecessor);
            const int to_state = state_of(index, to_support, predecessor);
            for (int slot = 0; slot < slots_; ++slot) {
              const Cell& cell = current_[cell_index(from_state, slot)];
              if (!cell.present)
                continue;
              push(tile, to_state, own(cell.key),
                   Provenance{static_cast<std::int16_t>(from_state),
                              static_cast<std::int8_t>(slot)},
                   value.id);
            }
          }
        }
        // Open a block on an unowned tile. The run before the null run is what
        // this candidate would be resuming, and the block it opens is itself
        // preceded by that null run, so it leaves no resumable predecessor.
        const int opened_from_null =
            state_of(index, opened_support, kNoPredecessor);
        for (int predecessor = 0; predecessor < predecessors_; ++predecessor) {
          const int from_state = null_state(predecessor);
          const std::int64_t open_cost =
              resuming() && predecessor == resume_predecessor
                  ? resume_cost_
                  : static_cast<std::int64_t>(p_.parameters.block_open_cost);
          for (int slot = 0; slot < slots_; ++slot) {
            const Cell& cell = current_[cell_index(from_state, slot)];
            if (!cell.present)
              continue;
            Key key = cell.key;
            key.score -= open_cost;
            ++key.non_null_blocks;
            push(tile, opened_from_null, own(key),
                 Provenance{static_cast<std::int16_t>(from_state),
                            static_cast<std::int8_t>(slot)},
                 value.id);
          }
        }
        // Close a different candidate's block and open this one.
        const auto open_after = [&](const Ranked& ranked, int to_state,
                                    std::int64_t open_cost) {
          Key key = ranked.key;
          key.score -= open_cost;
          ++key.non_null_blocks;
          ++key.candidate_changes;
          push(tile, to_state, own(key), ranked.from, value.id);
        };
        if (!resuming()) {
          // Only the top few closable entries can win, so the pairwise edges
          // collapse into one ranked list.
          int taken = 0;
          for (int at = 0; at < closable_size && taken < slots_; ++at) {
            const Ranked& ranked = closable[static_cast<std::size_t>(at)];
            if (ranked.candidate_index == index)
              continue;
            ++taken;
            open_after(
                ranked, state_of(index, opened_support, kNoPredecessor),
                static_cast<std::int64_t>(p_.parameters.block_open_cost));
          }
          continue;
        }
        for (int closing = 0; closing < candidates_; ++closing) {
          if (closing == index)
            continue;
          const int to_state =
              state_of(index, opened_support, predecessor_of(closing));
          const RankedList& list =
              owner_closable[static_cast<std::size_t>(closing)];
          const int size =
              owner_closable_size[static_cast<std::size_t>(closing)];
          // A source whose own predecessor is this candidate is an A-B-A resumption and pays
          // `resume_cost_`; the others pay the full open and come from the ranked list.
          int taken = 0;
          for (int at = 0; at < size && taken < slots_; ++at) {
            const Ranked& ranked = list[static_cast<std::size_t>(at)];
            if (ranked.predecessor == resume_predecessor)
              continue;
            ++taken;
            open_after(
                ranked, to_state,
                static_cast<std::int64_t>(p_.parameters.block_open_cost));
          }
          const int resumed_state =
              state_of(closing, kSupportStates - 1, resume_predecessor);
          for (int slot = 0; slot < slots_; ++slot) {
            const Cell& cell = current_[cell_index(resumed_state, slot)];
            if (!cell.present)
              continue;
            open_after({cell.key,
                        Provenance{static_cast<std::int16_t>(resumed_state),
                                   static_cast<std::int8_t>(slot)},
                        closing, resume_predecessor},
                       to_state, resume_cost_);
          }
        }
      }
    }

    std::swap(current_, next_);
    for (int state = 0; state < states_; ++state) {
      for (int slot = 0; slot < slots_; ++slot) {
        back_[static_cast<std::size_t>(tile) * layer +
              cell_index(state, slot)] = current_[cell_index(state, slot)].from;
      }
    }
  }

  QueryPartitionResult result;
  result.dp_cells =
      static_cast<std::uint64_t>(tiles_) * static_cast<std::uint64_t>(states_);
  result.dp_transitions = transitions_;
  result.tail_comparisons = tail_comparisons_;
  result.tail_steps = tail_steps_;
  if (p_.rival == QueryPartitionRival::WinnerOnly) {
    result.selected = best_final_path();
    return result;
  }

  // A path may finish only in the null state or with a closable open block.
  std::vector<QueryPartitionPath> complete;
  std::vector<CandidateId>& assignment = scratch_.assignment;
  const auto collect = [&](int state) {
    for (int slot = 0; slot < slots_; ++slot) {
      if (!current_[cell_index(state, slot)].present)
        continue;
      reconstruct(tiles_, state, slot, assignment);
      complete.push_back(evaluate(p_, assignment));
    }
  };
  for (int predecessor = 0; predecessor < predecessors_; ++predecessor)
    collect(null_state(predecessor));
  for (int index = 0; index < candidates_; ++index)
    for (int predecessor = 0; predecessor < predecessors_; ++predecessor)
      collect(state_of(index, kSupportStates - 1, predecessor));
  std::sort(
      complete.begin(), complete.end(),
      [&](const QueryPartitionPath& left, const QueryPartitionPath& right) {
        return compare_path(p_, left, right) > 0;
      });

  if (!complete.empty()) {
    result.selected = complete.front();
  } else {
    result.selected.assignment.assign(static_cast<std::size_t>(tiles_),
                                      kNullCandidate);
  }

  if (p_.rival == QueryPartitionRival::ExactSecondSolution) {
    std::vector<const QueryPartitionPath*> rivals;
    for (const QueryPartitionPath& path : complete)
      if (path.assignment != result.selected.assignment)
        rivals.push_back(&path);
    std::sort(
        rivals.begin(), rivals.end(),
        [&](const QueryPartitionPath* left, const QueryPartitionPath* right) {
          return compare_path(p_, *left, *right) > 0;
        });
    if (!rivals.empty()) {
      result.runner_up = *rivals.front();
      result.has_runner_up = true;
      int begin = 0;
      while (begin < tiles_ &&
             result.selected.assignment[static_cast<std::size_t>(begin)] ==
                 result.runner_up.assignment[static_cast<std::size_t>(begin)])
        ++begin;
      if (begin < tiles_) {
        int end = begin + 1;
        while (end < tiles_ &&
               result.selected.assignment[static_cast<std::size_t>(end)] !=
                   result.runner_up.assignment[static_cast<std::size_t>(end)])
          ++end;
        result.first_difference_begin = begin;
        result.first_difference_end = end;
      }
    }
  }
  return result;
}

QueryPartitionPath Solver::best_final_path() {
  std::vector<CandidateId>& assignment = scratch_.assignment;
  std::vector<CandidateId>& best_assignment = scratch_.best_assignment;
  int best = -1;
  bool best_built = false;
  const auto consider = [&](int state) {
    const Cell& cell = current_[cell_index(state, 0)];
    if (!cell.present)
      return;
    if (best >= 0) {
      const int order = compare_key(cell.key, current_[cell_index(best, 0)].key);
      if (order < 0)
        return;
      if (order == 0) {
        if (!best_built) {
          reconstruct(tiles_, best, 0, best_assignment);
          best_built = true;
        }
        reconstruct(tiles_, state, 0, assignment);
        if (compare_tail(assignment, best_assignment) <= 0)
          return;
        best_assignment.swap(assignment);
        best = state;
        return;
      }
    }
    best = state;
    best_built = false;
  };
  // Distinct final cells hold distinct assignments, so the winner is unique
  // and the visiting order does not matter.
  for (int predecessor = 0; predecessor < predecessors_; ++predecessor)
    consider(null_state(predecessor));
  for (int index = 0; index < candidates_; ++index)
    for (int predecessor = 0; predecessor < predecessors_; ++predecessor)
      consider(state_of(index, kSupportStates - 1, predecessor));
  if (best < 0) {
    QueryPartitionPath path;
    path.assignment.assign(static_cast<std::size_t>(tiles_), kNullCandidate);
    return path;
  }
  if (!best_built)
    reconstruct(tiles_, best, 0, best_assignment);
  return evaluate(p_, best_assignment);
}

} // namespace

QueryPartitionResult
solve_query_partition(const QueryPartitionProblem& problem) {
  if (problem.tile_count < 1 || problem.tile_count > kMaxQueryTiles)
    throw std::invalid_argument(
        "query partition tile count must be within [1,4096]");
  if (problem.catalogue.candidates.empty() ||
      problem.catalogue.candidates.size() >
          static_cast<std::size_t>(2 * kMaxCatalogueLaneBound))
    throw std::invalid_argument(
        "query partition candidate count is out of range");
  if (problem.parameters.minimum_supported_tiles_per_non_null_block !=
      kSupportStates - 1)
    throw std::invalid_argument(
        "query partition currently requires two supporting tiles per block");
  for (std::size_t index = 0; index < problem.catalogue.candidates.size();
       ++index) {
    if (problem.catalogue.candidates[index].id !=
        static_cast<CandidateId>(index))
      throw std::invalid_argument(
          "query partition candidate IDs must be their catalogue positions");
  }
  return Solver(problem).run();
}

} // namespace fa::cpu::voting
