#include "backend.h"

#include "chain_mapq.h"
#include "context.h"
#include "family_projection.h"
#include "family_realization.h"
#include "mapq_routing.h"
#include "placement_chaining.h"
#include "placement.h"
#include "placement_family_adapter.h"
#include "postdp_scoring.h"
#include "query_seed_pool.h"
#include "record_divergence.h"
#include "record_family.h"
#include "residue_emission.h"
#include "residue_trigger.h"
#include "worker_scratch.h"
#include "../core/cigar.h"
#include "../core/sequence.h"
#include "../core/types.h"
#include "../index/index.h"
#include "../index/seed.h"
#include "../seeding/context.h"
#include "../seeding/scratch.h"
#include "../seeding/select.h"
#include "../seeding/syncmer.h"
#include "../seeding/tie_hash.h" // tie_read_seed (the vote's tie-break)
#include "../seeding/types.h"
#include "../voting/vote.h"
#include "../voting/vote_slope.h" // vote_slope_fit_winner (the winner's window)

#include <algorithm>
#include <tuple>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {
namespace lr {

namespace dna {

void zero_hypothesis_mapq(AlignResult& hypothesis) {
  hypothesis.mapq = 0;
  for (AlignResult& supplementary : hypothesis.supplementary)
    supplementary.mapq = 0;
}

namespace {

// Stands in for the reverse-complement query when no stage reads its bases
// (see reverse_query_bases_needed); consumers take the length from the read.
const std::vector<uint8_t> kNoQueryBases;

// The chain MAPQ rule bits in force for this run.
int chain_mapq_study_bits(const DnaContext& dctx) noexcept {
  return dna_chain_mapq_study_bits(dctx.opts.chain_mapq_hifi_margin);
}

// One catalogue rival's MAPQ evidence (see whole_query_chain). A rival with no
// accepted whole-query chain keeps only its vote and screening span.
DnaChainMapqRival
chain_mapq_rival_from_placement(const DnaPlacementFamily& catalogue,
                                const DnaPlacementChainingResult& placement,
                                ::fa::cpu::voting::CandidateId id) {
  DnaChainMapqRival rival;
  rival.candidate = id;
  if (const DnaPlacementCandidate* candidate = catalogue.find(id)) {
    rival.vote_evidence = candidate->vote_evidence;
    rival.contig = candidate->peak.chr;
    rival.reverse = candidate->peak.is_rc;
    rival.peak_diagonal = candidate->peak.raw_ref_start;
  }
  const DnaPlacementCandidateChain* whole = placement.whole_query_chain(id);
  const DnaPlacementCandidateChain* record = placement.find(id);
  if (whole != nullptr)
    rival.status = static_cast<int>(whole->status);
  else if (record != nullptr)
    rival.status = static_cast<int>(record->status);
  if (whole != nullptr && whole->status == DnaPlacementChainStatus::Accepted &&
      whole->chain_anchors > 0) {
    rival.chained = true;
    rival.chain_score = whole->chain_score;
    rival.chain_anchors = whole->chain_anchors;
    rival.q_begin = whole->forward_query_begin;
    rival.q_end = whole->forward_query_end;
    rival.ref_begin = whole->reference_begin;
    rival.ref_end = whole->reference_end;
  } else if (record != nullptr) {
    rival.q_begin = record->screening_forward_query_begin;
    rival.q_end = record->screening_forward_query_end;
  }
  return rival;
}

// The MAPQ evidence for one block owner of the committed family: the primary,
// or a supplementary scored on its own candidate's chain. Every owner sees the
// same rivals: candidates that own no selected block, are not
// residue-admitted and are not the winner. Only the primary has a dp2 owner.
DnaChainMapqEvidence build_chain_mapq_evidence(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement, bool placement_ran,
    const DnaPlacementCandidateChain* winner_chain,
    ::fa::cpu::voting::CandidateId winner_candidate,
    ::fa::cpu::voting::CandidateId dp2_owner, double dp1, double dp2,
    double dp1_raw, double dp2_raw, double identity, int read_len) {
  DnaChainMapqEvidence evidence;
  evidence.read_len = read_len;
  evidence.match_sc = dctx.opts.cigar_dp_match;
  // As minimap2: sub_diff = opt->a * 2 + opt->b.
  evidence.sub_diff =
      2 * dctx.opts.cigar_dp_match + dctx.opts.cigar_dp_mismatch;
  evidence.dp1 = dp1;
  evidence.dp2 = dp2;
  evidence.identity = identity;
  // HiFi margin rule: raw scores and the scoring row they are read against.
  evidence.hifi_margin_rule = dctx.opts.chain_mapq_hifi_margin;
  evidence.dp1_raw = dp1_raw;
  evidence.dp2_raw = dp2_raw;
  evidence.substitution_cost =
      dctx.opts.cigar_dp_match + dctx.opts.cigar_dp_mismatch;
  evidence.gap_open1 = dctx.opts.cigar_dp_gap_open1;
  evidence.gap_extend1 = dctx.opts.cigar_dp_gap_extend1;
  evidence.shadow_slack_bp = dctx.opts.cigar_local_interval_anchor_interval_pad;
  // The shadow refinements (chain_mapq.h), on by default.
  static_assert(ResolvedDnaOptions{}.dna_chain_mapq_shadow_vote_credit &&
                    ResolvedDnaOptions{}.dna_chain_mapq_chain_shadow_read_slack,
                "the dense-chain MAPQ locus rules are on by default");
  evidence.shadow_vote_credit = dctx.opts.dna_chain_mapq_shadow_vote_credit;
  evidence.chain_shadow_read_slack =
      dctx.opts.dna_chain_mapq_chain_shadow_read_slack;
  evidence.study_bits = chain_mapq_study_bits(dctx);
  if (winner_chain != nullptr) {
    evidence.f1 = winner_chain->chain_score;
    evidence.cnt = winner_chain->chain_anchors;
    evidence.sib_f2 = winner_chain->rival_chain_score;
    evidence.winner_q_begin = winner_chain->forward_query_begin;
    evidence.winner_q_end = winner_chain->forward_query_end;
    evidence.winner_ref_begin = winner_chain->reference_begin;
    evidence.winner_ref_end = winner_chain->reference_end;
  }
  if (const DnaPlacementCandidate* winner = catalogue.find(winner_candidate)) {
    evidence.winner_vote = winner->vote_evidence;
    evidence.winner_contig = winner->peak.chr;
    evidence.winner_reverse = winner->peak.is_rc;
    evidence.winner_peak_diagonal = winner->peak.raw_ref_start;
  }
  if (!placement_ran)
    return evidence;

  std::vector<::fa::cpu::voting::CandidateId> ids;
  ids.reserve(catalogue.candidates.size());
  for (const DnaPlacementCandidate& candidate : catalogue.candidates) {
    if (candidate.id == winner_candidate || candidate.residue_admitted)
      continue;
    bool owns_selected_block = false;
    for (const auto& block : catalogue.partition.selected.blocks) {
      if (block.candidate == candidate.id) {
        owns_selected_block = true;
        break;
      }
    }
    if (owns_selected_block)
      continue;
    ids.push_back(candidate.id);
  }
  // The realized rival always competes, even when it owns a block: after a
  // promotion it is the demoted incumbent, and its DP score is dp2.
  if (dp2_owner != ::fa::cpu::voting::kNullCandidate &&
      dp2_owner != winner_candidate &&
      std::find(ids.begin(), ids.end(), dp2_owner) == ids.end())
    ids.push_back(dp2_owner);

  evidence.rivals.reserve(ids.size());
  for (const ::fa::cpu::voting::CandidateId id : ids) {
    DnaChainMapqRival rival =
        chain_mapq_rival_from_placement(catalogue, placement, id);
    rival.owns_dp2 = id == dp2_owner;
    evidence.rivals.push_back(rival);
  }
  // Strongest vote first, then catalogue id; the formula does not depend on
  // the order.
  std::stable_sort(evidence.rivals.begin(), evidence.rivals.end(),
                   [](const DnaChainMapqRival& a, const DnaChainMapqRival& b) {
                     return std::tie(b.vote_evidence, a.candidate) <
                            std::tie(a.vote_evidence, b.candidate);
                   });
  return evidence;
}

// A record's gap-compressed divergence and its event denominator, as the PAF
// de:f tag computes them. False without base-level accounting.
bool record_event_divergence(const AlignResult& record, double& divergence,
                             double& denominator) {
  if (!record.alignment_accounting_valid || record.edit_distance < 0)
    return false;
  int gap_bases = 0;
  int gap_opens = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    if (op.second == 'I' || op.second == 'D') {
      gap_bases += op.first;
      ++gap_opens;
    }
  }
  const int events = record.block_len + std::max(0, record.ambiguities) -
                     gap_bases + gap_opens;
  if (events <= 0)
    return false;
  denominator = static_cast<double>(events);
  divergence = 1.0 - static_cast<double>(record.matches) / denominator;
  return true;
}

bool clears_dna_emission_floor(const AlignResult& record, int min_dp_max);

// Extend-then-demote changes records, not MAPQ: the MAPQ is scored and routed
// on the family as realized before the extension, since a merged read would
// otherwise escape R6's split-read cap. For the MAPQ stage the pre-extension
// family stands in: copies of the demoted records return to their slots and
// the primary's alignment fields are swapped with its pre-extension record
// (install_mapq_stand_in); after routing both are undone
// (remove_mapq_stand_in). The demoted records become secondaries with MAPQ 0.

// Swaps the alignment fields of two records (locus, CIGAR, spans, score,
// accounting), leaving the per-read fields (family vectors, MAPQ, s2, median
// occurrence, read length) in place.
void exchange_alignment_fields(AlignResult& left, AlignResult& right) {
  using std::swap;
  swap(left.score, right.score);
  swap(left.chromosome, right.chromosome);
  swap(left.pos, right.pos);
  swap(left.is_reverse, right.is_reverse);
  swap(left.origin, right.origin);
  swap(left.cigar, right.cigar);
  swap(left.target_regions, right.target_regions);
  swap(left.query_start, right.query_start);
  swap(left.query_end, right.query_end);
  swap(left.target_end, right.target_end);
  swap(left.matches, right.matches);
  swap(left.mismatches, right.mismatches);
  swap(left.insertions, right.insertions);
  swap(left.deletions, right.deletions);
  swap(left.ambiguities, right.ambiguities);
  swap(left.edit_distance, right.edit_distance);
  swap(left.block_len, right.block_len);
  swap(left.alignment_accounting_valid, right.alignment_accounting_valid);
  swap(left.cs, right.cs);
  swap(left.md, right.md);
}

// The demoted block records with their owners and former slots, and the
// primary's pre-extension record (see DnaFamilyRealizationOutcome).
struct DnaDemotedBlocks {
  std::vector<AlignResult> records;
  std::vector<::fa::cpu::voting::CandidateId> candidates;
  std::vector<int> positions;
  bool primary_extended = false;
  AlignResult pre_extension_primary;
};

// Reinserts copies of the demoted records into their slots, in ascending
// order, and swaps in the pre-extension primary. Returns the slots filled.
std::vector<std::size_t> install_mapq_stand_in(DnaDemotedBlocks& demoted,
                                               dna::Result& realized) {
  std::vector<std::size_t> slots;
  if (!demoted.records.empty() &&
      demoted.candidates.size() == demoted.records.size() &&
      demoted.positions.size() == demoted.records.size()) {
    std::vector<std::size_t> order(demoted.records.size());
    for (std::size_t k = 0; k < order.size(); ++k)
      order[k] = k;
    std::sort(order.begin(), order.end(),
              [&demoted](std::size_t left, std::size_t right) {
                return demoted.positions[left] < demoted.positions[right];
              });
    for (const std::size_t k : order) {
      if (demoted.positions[k] < 0)
        continue;
      const std::size_t slot = static_cast<std::size_t>(demoted.positions[k]);
      if (slot > realized.supplementary.size() ||
          slot > realized.supplementary_candidates.size())
        continue;
      realized.supplementary.insert(
          realized.supplementary.begin() + static_cast<std::ptrdiff_t>(slot),
          demoted.records[k]);
      realized.supplementary_candidates.insert(
          realized.supplementary_candidates.begin() +
              static_cast<std::ptrdiff_t>(slot),
          demoted.candidates[k]);
      slots.push_back(slot);
    }
  }
  if (demoted.primary_extended)
    exchange_alignment_fields(realized, demoted.pre_extension_primary);
  return slots;
}

// Undoes install_mapq_stand_in; the MAPQ routed to a copy leaves with it.
void remove_mapq_stand_in(DnaDemotedBlocks& demoted,
                          const std::vector<std::size_t>& slots,
                          dna::Result& realized) {
  std::vector<std::size_t> descending = slots;
  std::sort(descending.begin(), descending.end(),
            [](std::size_t left, std::size_t right) { return left > right; });
  for (const std::size_t slot : descending) {
    const auto at = static_cast<std::ptrdiff_t>(slot);
    realized.supplementary.erase(realized.supplementary.begin() + at);
    realized.supplementary_candidates.erase(
        realized.supplementary_candidates.begin() + at);
  }
  if (demoted.primary_extended)
    exchange_alignment_fields(realized, demoted.pre_extension_primary);
}

// Collects each record's divergence for the divergence contrast (index 0 the
// primary, i + 1 the i-th supplementary) and runs it. Inversion middles and
// records below the emission floor take no part.
DnaFamilyDivergence divergence_contrast(const dna::Result& realized,
                                        int min_dp_max) {
  DnaFamilyDivergence family;
  family.records.assign(1 + realized.supplementary.size(),
                        DnaRecordDivergence());
  std::vector<double> denominators(family.records.size(), 0.0);
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    const AlignResult& record = index == 0
                                    ? static_cast<const AlignResult&>(realized)
                                    : realized.supplementary[index - 1];
    if (record.origin == AlignmentOrigin::DnaLocalInversion)
      continue;
    if (!clears_dna_emission_floor(record, min_dp_max))
      continue;
    if (!record_event_divergence(record, family.records[index].divergence,
                                 denominators[index]))
      continue;
    family.records[index].accounted = true;
  }
  dna_divergence_contrast(family, denominators);
  return family;
}

// A record's divergence structure: mismatches and ambiguous bases from its
// replay accounting, indel events and bases and aligned columns from its
// CIGAR. All -1 without a CIGAR; the first two -1 without accounting.
struct DnaChainMapqRecordCensus {
  int mismatches = -1;
  int ambiguous = -1;
  int ins_events = -1;
  int ins_bases = -1;
  int del_events = -1;
  int del_bases = -1;
  int aligned_bases = -1;
};

DnaChainMapqRecordCensus chain_mapq_record_census(const AlignResult& record) {
  DnaChainMapqRecordCensus census;
  if (record.cigar.empty())
    return census;
  census.ins_events = 0;
  census.ins_bases = 0;
  census.del_events = 0;
  census.del_bases = 0;
  census.aligned_bases = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    if (op.first <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X':
      census.aligned_bases += op.first;
      break;
    case 'I':
      ++census.ins_events;
      census.ins_bases += op.first;
      break;
    case 'D':
      ++census.del_events;
      census.del_bases += op.first;
      break;
    default:
      break;
    }
  }
  if (record.alignment_accounting_valid) {
    census.mismatches = record.mismatches;
    census.ambiguous = record.ambiguities;
  }
  return census;
}

// The scoring the interval-matched scorer reads.
DpScoringParams chain_mapq_dp_params(const DnaContext& dctx) {
  DpScoringParams dp;
  dp.match = dctx.opts.cigar_dp_match;
  dp.mismatch = dctx.opts.cigar_dp_mismatch;
  dp.ambi = dctx.opts.cigar_dp_ambi;
  dp.gap_open1 = dctx.opts.cigar_dp_gap_open1;
  dp.gap_extend1 = dctx.opts.cigar_dp_gap_extend1;
  dp.gap_open2 = dctx.opts.cigar_dp_gap_open2;
  dp.gap_extend2 = dctx.opts.cigar_dp_gap_extend2;
  return dp;
}

// A record's dual-affine score split by the forward-query interval [lo, hi)
// (affine_score_over_query_interval); zero without a CIGAR, contig or
// encoded reference.
AffineIntervalScore chain_mapq_interval_score(
    const DnaContext& dctx, const AlignResult& record, int contig,
    const std::vector<uint8_t>& fwd, const std::vector<uint8_t>& rc,
    int read_len, int lo, int hi, const DpScoringParams& dp) {
  AffineIntervalScore score;
  if (dctx.ref.encoded == nullptr || contig < 0 ||
      static_cast<std::size_t>(contig) >= dctx.ref.encoded->size() ||
      record.pos < 0 || record.cigar.empty())
    return score;
  const std::vector<uint8_t>& reference =
      (*dctx.ref.encoded)[static_cast<std::size_t>(contig)];
  const std::vector<uint8_t>& query = record.is_reverse ? rc : fwd;
  return affine_score_over_query_interval(
      ::fa::cpu::output::parse_cigar_ops(record.cigar), query.data(),
      reference.data() + record.pos, read_len, record.is_reverse, lo, hi, dp);
}

// A record's start diagonal: reference start minus oriented query start.
std::int64_t chain_mapq_record_diagonal(const AlignResult& record,
                                        int read_len) {
  const std::int64_t oriented =
      record.is_reverse ? read_len - record.query_end : record.query_start;
  return static_cast<std::int64_t>(record.pos) - oriented;
}

// The contig index of a name, or -1 when unknown.
int chain_mapq_contig_index(const DnaContext& dctx, const std::string& name) {
  if (dctx.ref.names == nullptr)
    return -1;
  const std::vector<std::string>& names = *dctx.ref.names;
  for (std::size_t index = 0; index < names.size(); ++index)
    if (names[index] == name)
      return static_cast<int>(index);
  return -1;
}

// A committed record's alignment as a map from oriented query position to
// diagonal, one entry per M/=/X run of its CIGAR.
struct DnaChainMapqRecordDiagonals {
  struct Run {
    int q_begin = 0; // oriented query, half-open
    int q_end = 0;
    std::int64_t diagonal = 0;
  };
  std::vector<Run> runs;
};

DnaChainMapqRecordDiagonals
chain_mapq_record_diagonals(const AlignResult& record) {
  DnaChainMapqRecordDiagonals out;
  if (record.cigar.empty() || record.pos < 0)
    return out;
  int qoff = 0;
  std::int64_t toff = 0;
  for (const std::pair<int, char>& op :
       ::fa::cpu::output::parse_cigar_ops(record.cigar)) {
    const int len = op.first;
    if (len <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X':
      out.runs.push_back({qoff, qoff + len,
                          static_cast<std::int64_t>(record.pos) + toff - qoff});
      qoff += len;
      toff += len;
      break;
    case 'I':
    case 'S':
      qoff += len;
      break;
    case 'D':
    case 'N':
      toff += len;
      break;
    default: // H, P
      break;
    }
  }
  return out;
}

// The record's diagonal at oriented query position `oq`: the run holding it,
// else the nearest run boundary within `slack` (an insertion run between two
// aligned runs, or a clipped end). False when the record does not reach it.
bool chain_mapq_record_diagonal_at(const DnaChainMapqRecordDiagonals& map,
                                   std::int64_t oq, std::int64_t slack,
                                   std::int64_t& diagonal) {
  std::int64_t best_distance = -1;
  for (const DnaChainMapqRecordDiagonals::Run& run : map.runs) {
    std::int64_t distance = 0;
    if (oq < run.q_begin)
      distance = run.q_begin - oq;
    else if (oq >= run.q_end)
      distance = oq - (run.q_end - 1);
    if (best_distance < 0 || distance < best_distance) {
      best_distance = distance;
      diagonal = run.diagonal;
    }
    if (distance == 0)
      break;
  }
  return best_distance >= 0 && best_distance <= slack;
}

// A whole-query chain's contig, strand, forward-query span and reference
// span. Without a chain or either span it is invalid: it explains nothing and
// nothing explains it.
struct DnaChainMapqChainGeometry {
  int contig = -1;
  bool reverse = false;
  int q_begin = -1;
  int q_end = -1;
  int ref_begin = -1;
  int ref_end = -1;

  bool valid() const {
    return contig >= 0 && q_begin >= 0 && q_end > q_begin && ref_begin >= 0 &&
           ref_end > ref_begin;
  }
};

DnaChainMapqChainGeometry
chain_mapq_rival_geometry(const DnaChainMapqRival& rival) {
  DnaChainMapqChainGeometry chain;
  if (!rival.chained)
    return chain;
  chain.contig = rival.contig;
  chain.reverse = rival.reverse;
  chain.q_begin = rival.q_begin;
  chain.q_end = rival.q_end;
  chain.ref_begin = rival.ref_begin;
  chain.ref_end = rival.ref_end;
  return chain;
}

DnaChainMapqChainGeometry
chain_mapq_chain_geometry(const DnaPlacementCandidateChain& chain, int contig,
                          bool reverse) {
  DnaChainMapqChainGeometry out;
  out.contig = contig;
  out.reverse = reverse;
  out.q_begin = chain.forward_query_begin;
  out.q_end = chain.forward_query_end;
  out.ref_begin = chain.reference_begin;
  out.ref_end = chain.reference_end;
  return out;
}

// Whether one committed record explains a chain: both chain ends lie on the
// record's diagonal there, within `slack`, on the same contig and strand.
bool chain_mapq_record_explains(const DnaChainMapqRecordDiagonals& map,
                                int record_contig, bool record_reverse,
                                const DnaChainMapqChainGeometry& chain,
                                int read_len, std::int64_t slack) {
  if (!chain.valid() || chain.contig != record_contig ||
      chain.reverse != record_reverse)
    return false;
  const std::int64_t oriented_begin =
      chain.reverse ? read_len - chain.q_end : chain.q_begin;
  const std::int64_t oriented_end =
      chain.reverse ? read_len - chain.q_begin : chain.q_end;
  const std::int64_t start =
      static_cast<std::int64_t>(chain.ref_begin) - oriented_begin;
  const std::int64_t end =
      static_cast<std::int64_t>(chain.ref_end) - oriented_end;
  std::int64_t at_start = 0, at_end = 0;
  if (!chain_mapq_record_diagonal_at(map, oriented_begin, slack, at_start) ||
      !chain_mapq_record_diagonal_at(map, oriented_end - 1, slack, at_end))
    return false;
  return std::llabs(at_start - start) <= slack &&
         std::llabs(at_end - end) <= slack;
}

// The committed family's diagonal maps, built once per read (index 0 the
// primary, i + 1 the i-th supplementary). A record without CIGAR, position,
// known contig or aligned run explains nothing.
struct DnaChainMapqFamilyMap {
  struct Record {
    DnaChainMapqRecordDiagonals diagonals;
    int contig = -1;
    bool reverse = false;
    int q_begin = -1; // forward query, half-open
    int q_end = -1;
    bool inversion = false;
    // The candidate this record came from, or null for a terminal clip or
    // inversion middle. A record does not vouch for its own owner.
    ::fa::cpu::voting::CandidateId owner = ::fa::cpu::voting::kNullCandidate;
  };
  std::vector<Record> records;

  // Bit 8: the block's own record explains the chain.
  bool explains_by_record(std::size_t index,
                          const DnaChainMapqChainGeometry& chain, int read_len,
                          std::int64_t slack) const {
    if (index >= records.size())
      return false;
    const Record& record = records[index];
    return chain_mapq_record_explains(record.diagonals, record.contig,
                                      record.reverse, chain, read_len, slack);
  }

  // Bit 64: the chain's start lies on some record's diagonal and its end on
  // some record's, the same or another, as for a chain spanning a deletion
  // between two flank records. Inversion middles and the record `candidate`
  // owns (which trivially explains its own chain) take no part.
  bool explains_by_family(const DnaChainMapqChainGeometry& chain, int read_len,
                          std::int64_t slack,
                          ::fa::cpu::voting::CandidateId candidate) const {
    if (!chain.valid())
      return false;
    const std::int64_t oriented_begin =
        chain.reverse ? read_len - chain.q_end : chain.q_begin;
    const std::int64_t oriented_end =
        chain.reverse ? read_len - chain.q_begin : chain.q_end;
    const std::int64_t start =
        static_cast<std::int64_t>(chain.ref_begin) - oriented_begin;
    const std::int64_t end =
        static_cast<std::int64_t>(chain.ref_end) - oriented_end;
    bool start_held = false;
    bool end_held = false;
    for (const Record& record : records) {
      if (record.inversion || record.contig != chain.contig ||
          record.reverse != chain.reverse)
        continue;
      if (candidate != ::fa::cpu::voting::kNullCandidate &&
          record.owner == candidate)
        continue;
      std::int64_t at = 0;
      if (!start_held &&
          chain_mapq_record_diagonal_at(record.diagonals, oriented_begin, slack,
                                        at) &&
          std::llabs(at - start) <= slack)
        start_held = true;
      if (!end_held &&
          chain_mapq_record_diagonal_at(record.diagonals, oriented_end - 1,
                                        slack, at) &&
          std::llabs(at - end) <= slack)
        end_held = true;
      if (start_held && end_held)
        return true;
    }
    return false;
  }
};

DnaChainMapqFamilyMap chain_mapq_family_map(
    const DnaContext& dctx, const dna::Result& realized,
    ::fa::cpu::voting::CandidateId primary_candidate) {
  DnaChainMapqFamilyMap map;
  map.records.resize(1 + realized.supplementary.size());
  for (std::size_t index = 0; index < map.records.size(); ++index) {
    const AlignResult& member =
        index == 0 ? static_cast<const AlignResult&>(realized)
                   : realized.supplementary[index - 1];
    map.records[index].owner =
        index == 0 ? primary_candidate
        : index - 1 < realized.supplementary_candidates.size()
            ? realized.supplementary_candidates[index - 1]
            : ::fa::cpu::voting::kNullCandidate;
    if (member.cigar.empty() || member.pos < 0)
      continue;
    const int contig = chain_mapq_contig_index(dctx, member.chromosome);
    if (contig < 0)
      continue;
    DnaChainMapqFamilyMap::Record& record = map.records[index];
    record.diagonals = chain_mapq_record_diagonals(member);
    if (record.diagonals.runs.empty())
      continue;
    record.contig = contig;
    record.reverse = member.is_reverse;
    record.q_begin = member.query_start;
    record.q_end = member.query_end;
    record.inversion = member.origin == AlignmentOrigin::DnaLocalInversion;
  }
  return map;
}

// Bits 8 / 64: marks every chained rival the committed record (or family)
// explains as a record shadow, and adds the block owners it explains as
// rivals, so their vote is credited to the winner. `record_index` is the
// block's own record in the family map.
void mark_chain_mapq_record_shadows(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement,
    const DnaChainMapqFamilyMap& family_map, std::size_t record_index,
    ::fa::cpu::voting::CandidateId winner_candidate, int read_len,
    DnaChainMapqEvidence& evidence) {
  const int bits = chain_mapq_study_bits(dctx);
  if ((bits & (kDnaChainMapqStudyRecordShadow |
               kDnaChainMapqStudyFamilyVouch)) == 0)
    return;
  // The block's own record is asked first, then the family under bit 64.
  const bool family_arm = (bits & kDnaChainMapqStudyFamilyVouch) != 0;
  const std::int64_t slack = dna_chain_mapq_chain_shadow_slack(evidence);
  const auto vouch_for = [&](const DnaChainMapqRival& rival) {
    const DnaChainMapqChainGeometry chain = chain_mapq_rival_geometry(rival);
    if (family_map.explains_by_record(record_index, chain, read_len, slack))
      return 1;
    if (family_arm && family_map.explains_by_family(chain, read_len, slack,
                                                    rival.candidate))
      return 2;
    return 0;
  };
  for (DnaChainMapqRival& rival : evidence.rivals) {
    const int vouch = vouch_for(rival);
    rival.vouch = vouch;
    rival.record_shadow = vouch != 0;
  }
  bool appended = false;
  for (const DnaPlacementCandidate& candidate : catalogue.candidates) {
    if (candidate.id == winner_candidate || candidate.residue_admitted)
      continue;
    bool owns_selected_block = false;
    for (const auto& block : catalogue.partition.selected.blocks)
      if (block.candidate == candidate.id) {
        owns_selected_block = true;
        break;
      }
    if (!owns_selected_block)
      continue;
    bool present = false;
    for (const DnaChainMapqRival& rival : evidence.rivals)
      if (rival.candidate == candidate.id) {
        present = true;
        break;
      }
    if (present)
      continue;
    DnaChainMapqRival rival =
        chain_mapq_rival_from_placement(catalogue, placement, candidate.id);
    const int vouch = vouch_for(rival);
    if (vouch == 0)
      continue;
    rival.record_shadow = true;
    rival.vouch = vouch;
    evidence.rivals.push_back(rival);
    appended = true;
  }
  if (appended)
    std::stable_sort(
        evidence.rivals.begin(), evidence.rivals.end(),
        [](const DnaChainMapqRival& a, const DnaChainMapqRival& b) {
          return std::tie(b.vote_evidence, a.candidate) <
                 std::tie(a.vote_evidence, b.candidate);
        });
}

// The interval-matched scores of the incumbent and the alternative (R1), read
// before the commit and terminal-clip recovery change either family. [lo, hi)
// is the intersection of the two primaries' forward-query spans (hi < lo
// means disjoint). a_* sums the incumbent's primary and block-owning
// supplementaries inside / outside it; b_* is the alternative primary's.
struct DnaChainMapqAlternativeInterval {
  bool scored = false;
  int lo = -1;
  int hi = -1;
  long long a_inside = 0;
  long long a_outside = 0;
  long long b_inside = 0;
  long long b_outside = 0;
};

DnaChainMapqAlternativeInterval chain_mapq_alternative_interval(
    const DnaContext& dctx, const DnaFamilyRealizationOutcome& incumbent,
    const DnaFamilyRealizationOutcome& alternative,
    const std::vector<uint8_t>& fwd, const std::vector<uint8_t>& rc,
    int read_len) {
  DnaChainMapqAlternativeInterval interval;
  const AlignResult& win = incumbent.output;
  const AlignResult& alt = alternative.output;
  interval.scored = true;
  interval.lo = std::max(win.query_start, alt.query_start);
  interval.hi = std::min(win.query_end, alt.query_end);
  const DpScoringParams dp = chain_mapq_dp_params(dctx);
  const auto score_record = [&](const AlignResult& record) {
    return chain_mapq_interval_score(
        dctx, record, chain_mapq_contig_index(dctx, record.chromosome), fwd,
        rc, read_len, interval.lo, interval.hi, dp);
  };
  // Records that own no block (inversion middles) are left out.
  AffineIntervalScore a = score_record(win);
  for (std::size_t index = 0; index < win.supplementary.size(); ++index) {
    if (index < incumbent.output.supplementary_candidates.size() &&
        incumbent.output.supplementary_candidates[index] ==
            ::fa::cpu::voting::kNullCandidate)
      continue;
    const AffineIntervalScore part = score_record(win.supplementary[index]);
    a.inside += part.inside;
    a.outside += part.outside;
  }
  const AffineIntervalScore b = score_record(alt);
  interval.a_inside = a.inside;
  interval.a_outside = a.outside;
  interval.b_inside = b.inside;
  interval.b_outside = b.outside;
  return interval;
}

// The realized rivals (R1), HiFi presets only, run once placement,
// realization, the commit and terminal-clip recovery are final. Each rival
// goes through build_dna_rival_placement and realize_full_cigar_family into a
// local outcome; nothing from it is emitted. The results reach the formula as
// DnaChainMapqEvidence::sibling / block_rival.

// Realizes one rival chain of candidate `original` and measures it against
// `against`, the emitted record it competes with. The cost is recorded
// either way; returns whether the realization was accepted.
// [clip_begin, clip_end) confines it to a forward-query window (none when
// clip_end <= clip_begin); block rivals are clipped to their block's span, so
// their dp2_raw and record counts describe the clipped alignment.
bool realize_chain_mapq_rival(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    ::fa::cpu::voting::CandidateId original,
    const DnaPlacementCandidateChain& chain, const std::vector<uint8_t>& fwd,
    const std::vector<uint8_t>& rc, int read_len, const AlignResult& against,
    int against_contig, int rival_contig, DnaChainMapqRealizedRival& into,
    int clip_begin = 0, int clip_end = 0) {
  into.candidate = original;
  into.chain_score = chain.chain_score;
  into.chain_anchors = chain.chain_anchors;
  into.refusal = DnaChainMapqRealizeRefusal::RealizationFailed;
  const DnaPlacementChainingResult placement = build_dna_rival_placement(
      dctx, catalogue, original, chain, fwd, rc);
  if (!placement.accepted || !placement.family.valid ||
      !placement.family.original_candidate_id)
    return false;
  DnaFamilyRealizationRequest request;
  request.family = &placement.family;
  request.placement = &placement;
  request.forward_query = &fwd;
  request.reverse_query = &rc;
  request.clip_forward_begin = clip_begin;
  request.clip_forward_end = clip_end;
  const DnaFamilyRealizationOutcome outcome =
      realize_full_cigar_family(dctx, request);
  into.ksw2_attempts = outcome.ksw2_attempts;
  into.estimated_cells = outcome.estimated_cells;
  if (!outcome.accepted())
    return false;
  const AlignResult& rival = outcome.output;
  into.realized = true;
  into.refusal = DnaChainMapqRealizeRefusal::None;
  into.dp2_raw = static_cast<double>(outcome.decision_score);
  into.rec_score = rival.score;
  into.q_begin = rival.query_start;
  into.q_end = rival.query_end;
  into.contig = rival_contig;
  into.reverse = rival.is_reverse;
  into.ref_begin = rival.pos;
  into.ref_end = rival.target_end;
  if (rival_contig == against_contig && rival.is_reverse == against.is_reverse)
    into.diag_shift = chain_mapq_record_diagonal(rival, read_len) -
                      chain_mapq_record_diagonal(against, read_len);
  const DnaChainMapqRecordCensus census = chain_mapq_record_census(rival);
  into.mismatches = census.mismatches;
  into.ins_events = census.ins_events;
  into.del_events = census.del_events;
  into.aligned_bases = census.aligned_bases;
  // The intersection of the two forward spans; im_hi < im_lo means disjoint.
  into.im_lo = std::max(against.query_start, rival.query_start);
  into.im_hi = std::min(against.query_end, rival.query_end);
  const DpScoringParams dp = chain_mapq_dp_params(dctx);
  const AffineIntervalScore a = chain_mapq_interval_score(
      dctx, against, against_contig, fwd, rc, read_len, into.im_lo, into.im_hi,
      dp);
  const AffineIntervalScore b = chain_mapq_interval_score(
      dctx, rival, rival_contig, fwd, rc, read_len, into.im_lo, into.im_hi, dp);
  into.im_a_inside = a.inside;
  into.im_a_outside = a.outside;
  into.im_b_inside = b.inside;
  into.im_b_outside = b.outside;
  return true;
}

// What was realized for one read: the winner's sibling (single-block
// families) or one rival per block (index 0 the primary, i + 1 the i-th
// supplementary), and the total cost.
struct DnaChainMapqEvidenceRealizations {
  DnaChainMapqRealizedRival sibling;
  std::vector<DnaChainMapqRealizedRival> block_rivals;
  int ksw2_attempts = 0;
  std::int64_t estimated_cells = 0;
};

// Selects and realizes the rivals of one read.
//   (B) One block: the winner's rival sibling, when
//       max(sib_f2, min_chain_score) / f1 reaches
//       kDnaChainMapqSiblingRealizeMin, realized against the primary.
//   (A) Several blocks: for each block-owning record, its best rival
//       (select_dna_block_rival), realized clipped to the record's query
//       span, or unclipped if the clipped realization fails. At most
//       kDnaChainMapqBlockRivalMax per read, longest blocks first.
// `winner_chain` is the committed primary's whole-query chain; `family_map`
// feeds the own-locus verdict under kDnaChainMapqStudyLaneVerdict.
DnaChainMapqEvidenceRealizations realize_chain_mapq_evidence(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement,
    const DnaChainMapqFamilyMap& family_map,
    const DnaPlacementCandidateChain* winner_chain,
    ::fa::cpu::voting::CandidateId primary_candidate,
    const dna::Result& realized, const std::vector<uint8_t>& fwd,
    const std::vector<uint8_t>& rc, int read_len) {
  DnaChainMapqEvidenceRealizations out;
  out.block_rivals.assign(1 + realized.supplementary.size(),
                          DnaChainMapqRealizedRival());
  const auto charge = [&out](const DnaChainMapqRealizedRival& rival) {
    out.ksw2_attempts += rival.ksw2_attempts;
    out.estimated_cells += rival.estimated_cells;
  };
  struct Block {
    std::size_t record_index;
    ::fa::cpu::voting::CandidateId owner;
    const DnaPlacementCandidateChain* chain;
    const AlignResult* record;
  };
  std::vector<Block> blocks;
  blocks.push_back({0, primary_candidate, winner_chain,
                    &static_cast<const AlignResult&>(realized)});
  for (std::size_t index = 0; index < realized.supplementary.size(); ++index) {
    const ::fa::cpu::voting::CandidateId owner =
        index < realized.supplementary_candidates.size()
            ? realized.supplementary_candidates[index]
            : ::fa::cpu::voting::kNullCandidate;
    if (owner == ::fa::cpu::voting::kNullCandidate)
      continue;
    blocks.push_back({index + 1, owner, placement.find(owner),
                      &realized.supplementary[index]});
  }
  if (blocks.size() == 1) {
    // (B) the sibling path.
    DnaChainMapqRealizedRival& sibling = out.sibling;
    const DnaPlacementCandidate* winner = catalogue.find(primary_candidate);
    if (winner_chain == nullptr || winner == nullptr ||
        winner_chain->chain_score <= 0 || winner_chain->rival_sibling < 0) {
      sibling.refusal = DnaChainMapqRealizeRefusal::NoSibling;
      return out;
    }
    const double floor = static_cast<double>(
        std::max(winner_chain->rival_chain_score, kDnaMinChainScore));
    if (floor < kDnaChainMapqSiblingRealizeMin *
                    static_cast<double>(winner_chain->chain_score)) {
      sibling.refusal = DnaChainMapqRealizeRefusal::WeakSibling;
      return out;
    }
    const DnaPlacementCandidateChain sibling_chain = dna_sibling_rival_chain(
        *winner_chain, winner->peak.is_rc, read_len, catalogue.seed_length);
    if (sibling_chain.primary.empty()) {
      sibling.refusal = DnaChainMapqRealizeRefusal::NoSibling;
      return out;
    }
    realize_chain_mapq_rival(dctx, catalogue, primary_candidate, sibling_chain,
                             fwd, rc, read_len, realized, winner->peak.chr,
                             winner->peak.chr, sibling);
    charge(sibling);
    return out;
  }
  // (A) the per-block path.
  struct Planned {
    std::size_t block;
    DnaBlockRivalSelection selection;
    int query_len;
    ::fa::cpu::voting::CandidateId owner;
  };
  // Bit 128: a candidate that is this read's own placement is skipped before
  // it is realized. The slack is the evidence's chained shadow slack.
  const int study_bits = chain_mapq_study_bits(dctx);
  const bool lane_verdict =
      (study_bits & kDnaChainMapqStudyLaneVerdict) != 0;
  DnaChainMapqEvidence slack_frame;
  slack_frame.read_len = read_len;
  slack_frame.shadow_slack_bp =
      dctx.opts.cigar_local_interval_anchor_interval_pad;
  slack_frame.chain_shadow_read_slack =
      dctx.opts.dna_chain_mapq_chain_shadow_read_slack;
  const std::int64_t verdict_slack =
      dna_chain_mapq_chain_shadow_slack(slack_frame);
  std::vector<Planned> planned;
  for (std::size_t which = 0; which < blocks.size(); ++which) {
    const Block& block = blocks[which];
    DnaChainMapqRealizedRival& into = out.block_rivals[block.record_index];
    const DnaPlacementCandidate* owner = catalogue.find(block.owner);
    if (block.chain == nullptr || owner == nullptr) {
      into.refusal = DnaChainMapqRealizeRefusal::NotAttempted;
      continue;
    }
    const AlignResult& record = *block.record;
    // Own locus: a chain shadow of the block owner's chain on the same
    // contig and strand, or explained by the record / family (bits 8, 64).
    const std::size_t record_index = block.record_index;
    const int owner_contig = owner->peak.chr;
    const bool owner_reverse = owner->peak.is_rc;
    const DnaPlacementCandidateChain& owner_chain = *block.chain;
    std::function<bool(::fa::cpu::voting::CandidateId,
                       const DnaPlacementCandidateChain&)>
        own_locus;
    if (lane_verdict)
      own_locus = [&catalogue, &family_map, &owner_chain, owner_contig,
                   owner_reverse, record_index, read_len, verdict_slack,
                   study_bits](::fa::cpu::voting::CandidateId id,
                               const DnaPlacementCandidateChain& chain) {
        const DnaPlacementCandidate* rival = catalogue.find(id);
        if (rival == nullptr)
          return false;
        if (rival->peak.chr == owner_contig &&
            rival->peak.is_rc == owner_reverse &&
            dna_chain_mapq_chain_shadow(
                owner_chain.reference_begin, owner_chain.reference_end,
                owner_chain.forward_query_begin,
                owner_chain.forward_query_end, owner_reverse,
                chain.reference_begin, chain.reference_end,
                chain.forward_query_begin, chain.forward_query_end,
                rival->peak.is_rc, read_len, verdict_slack)
                .shadow)
          return true;
        const DnaChainMapqChainGeometry geometry = chain_mapq_chain_geometry(
            chain, rival->peak.chr, rival->peak.is_rc);
        if ((study_bits & kDnaChainMapqStudyRecordShadow) != 0 &&
            family_map.explains_by_record(record_index, geometry, read_len,
                                          verdict_slack))
          return true;
        if ((study_bits & kDnaChainMapqStudyFamilyVouch) == 0 ||
            !family_map.explains_by_family(geometry, read_len, verdict_slack,
                                           id))
          return false;
        return true;
      };
    Planned plan;
    plan.block = which;
    plan.selection = select_dna_block_rival(
        dctx, catalogue, placement, block.chain->chain_score,
        record.query_start, record.query_end, owner->peak.chr,
        owner->peak.is_rc, chain_mapq_record_diagonal(record, read_len),
        own_locus ? &own_locus : nullptr);
    into.refusal = plan.selection.refusal;
    if (plan.selection.candidate == ::fa::cpu::voting::kNullCandidate ||
        plan.selection.chain == nullptr)
      continue;
    plan.query_len = std::max(0, record.query_end - record.query_start);
    plan.owner = block.owner;
    planned.push_back(plan);
  }
  std::stable_sort(planned.begin(), planned.end(),
                   [](const Planned& a, const Planned& b) {
                     return std::tie(b.query_len, a.owner) <
                            std::tie(a.query_len, b.owner);
                   });
  for (std::size_t rank = 0; rank < planned.size(); ++rank) {
    const Planned& plan = planned[rank];
    const Block& block = blocks[plan.block];
    DnaChainMapqRealizedRival& into = out.block_rivals[block.record_index];
    if (rank >= static_cast<std::size_t>(kDnaChainMapqBlockRivalMax)) {
      into.refusal = DnaChainMapqRealizeRefusal::Capped;
      into.candidate = plan.selection.candidate;
      continue;
    }
    const DnaPlacementCandidate* owner = catalogue.find(block.owner);
    const DnaPlacementCandidate* rival =
        catalogue.find(plan.selection.candidate);
    if (owner == nullptr || rival == nullptr) {
      into.refusal = DnaChainMapqRealizeRefusal::RealizationFailed;
      continue;
    }
    // Realize over the block's query span, the only interval the rule
    // compares; if that fails (short blocks run out of anchors), realize the
    // whole chain instead. Both attempts are charged.
    if (!realize_chain_mapq_rival(
            dctx, catalogue, plan.selection.candidate, *plan.selection.chain,
            fwd, rc, read_len, *block.record, owner->peak.chr,
            rival->peak.chr, into, block.record->query_start,
            block.record->query_end)) {
      const int refused_attempts = into.ksw2_attempts;
      const long long refused_cells = into.estimated_cells;
      if (realize_chain_mapq_rival(dctx, catalogue, plan.selection.candidate,
                                   *plan.selection.chain, fwd, rc, read_len,
                                   *block.record, owner->peak.chr,
                                   rival->peak.chr, into))
        into.clip_arm = 2;
      into.ksw2_attempts += refused_attempts;
      into.estimated_cells += refused_cells;
    } else {
      into.clip_arm = 1;
    }
    charge(into);
  }
  return out;
}

// minimap2's emission floor (mm_filter_regs): no record with matched bases
// under min_chain_score or DP score under min_dp_max (-s). min_dp_max <= 0
// disables the DP half.
bool clears_dna_emission_floor(const AlignResult& record, int min_dp_max) {
  if (record.matches < kDnaMinChainScore)
    return false;
  return min_dp_max <= 0 || record.score >= min_dp_max;
}

// Applies the floor to one record family in place. Failing supplementaries are
// erased, since the writers emit every element regardless of mapped().
// `candidates`, when given, is kept parallel to the records. Returns false
// when nothing in the family clears the floor.
bool filter_dna_record_family(
    AlignResult& primary, std::vector<AlignResult>& supplementary,
    std::vector<::fa::cpu::voting::CandidateId>* candidates, int min_dp_max) {
  std::size_t kept = 0;
  for (std::size_t index = 0; index < supplementary.size(); ++index) {
    if (!clears_dna_emission_floor(supplementary[index], min_dp_max))
      continue;
    if (kept != index) {
      supplementary[kept] = std::move(supplementary[index]);
      if (candidates != nullptr && index < candidates->size())
        (*candidates)[kept] = (*candidates)[index];
    }
    ++kept;
  }
  supplementary.erase(supplementary.begin() + static_cast<std::ptrdiff_t>(kept),
                      supplementary.end());
  if (candidates != nullptr && candidates->size() > kept)
    candidates->resize(kept);
  if (clears_dna_emission_floor(primary, min_dp_max))
    return true;
  if (supplementary.empty())
    return false;
  // The primary failed but supplementaries remain: promote the one family
  // assembly would have named primary.
  const std::vector<AlignResult>::iterator promoted = std::max_element(
      supplementary.begin(), supplementary.end(), dna_primary_precedence_less);
  const std::size_t promoted_index =
      static_cast<std::size_t>(promoted - supplementary.begin());
  AlignResult carried = std::move(*promoted);
  supplementary.erase(promoted);
  if (candidates != nullptr && promoted_index < candidates->size())
    candidates->erase(candidates->begin() +
                      static_cast<std::ptrdiff_t>(promoted_index));
  // Read-level fields carry over; the promoted record keeps its own MAPQ.
  carried.read_len = primary.read_len;
  carried.median_occurrence = primary.median_occurrence;
  // Move the survivors and secondaries out first: `supplementary` aliases
  // primary.supplementary.
  std::vector<AlignResult> survivors = std::move(supplementary);
  std::vector<AlignResult> secondaries = std::move(primary.secondary);
  primary = std::move(carried);
  primary.supplementary = std::move(survivors);
  primary.secondary = std::move(secondaries);
  return true;
}

// Applies the emission floor to a read. Only a read left with no record at
// all is marked unmapped.
void apply_dna_emission_floor(int min_dp_max, dna::Result& realized) {
  for (std::vector<AlignResult>::iterator it = realized.secondary.begin();
       it != realized.secondary.end();) {
    if (filter_dna_record_family(*it, it->supplementary, nullptr, min_dp_max))
      ++it;
    else
      it = realized.secondary.erase(it);
  }
  if (filter_dna_record_family(static_cast<AlignResult&>(realized),
                               realized.supplementary,
                               &realized.supplementary_candidates, min_dp_max))
    return;
  demote_unmapped(realized);
  realized.supplementary.clear();
  realized.supplementary_candidates.clear();
  realized.secondary.clear();
}

// Map-only's floor. minimap2 keeps no chain under min_chain_score
// (mm_chain_dp) and does not test again what it splits off one
// (mm_split_reg), so a block record is judged by the whole-query chain its
// MAPQ is scored on, which a sibling slice prints only a share of. A record
// the MAPQ scores on no chain (a terminal clip, a secondary) is judged by its
// AS. The others print as before: a failing head gives its place to the next
// record. `owners`, when given, is parallel to `rest`. Returns false when no
// record is left.
bool keep_dna_min_chain_score(
    AlignResult& head, int head_chain_score, std::vector<AlignResult>& rest,
    std::vector<::fa::cpu::voting::CandidateId>* owners,
    const DnaPlacementChainingResult* placement) {
  std::size_t kept = 0;
  for (std::size_t index = 0; index < rest.size(); ++index) {
    const DnaPlacementCandidateChain* chain =
        owners != nullptr && placement != nullptr &&
                (*owners)[index] != ::fa::cpu::voting::kNullCandidate
            ? placement->find((*owners)[index])
            : nullptr;
    if ((chain != nullptr ? chain->chain_score : rest[index].score) <
        kDnaMinChainScore)
      continue;
    if (kept != index) {
      rest[kept] = std::move(rest[index]);
      if (owners != nullptr)
        (*owners)[kept] = (*owners)[index];
    }
    ++kept;
  }
  rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(kept), rest.end());
  if (owners != nullptr)
    owners->resize(kept);
  if (head_chain_score >= kDnaMinChainScore)
    return true;
  if (rest.empty())
    return false;
  // `rest` aliases head.supplementary: move everything out first.
  std::vector<AlignResult> survivors = std::move(rest);
  std::vector<AlignResult> secondaries = std::move(head.secondary);
  AlignResult next = std::move(survivors.front());
  survivors.erase(survivors.begin());
  if (owners != nullptr)
    owners->erase(owners->begin());
  // Read-level fields carry over; the record keeps its own MAPQ.
  next.read_len = head.read_len;
  next.median_occurrence = head.median_occurrence;
  head = std::move(next);
  head.supplementary = std::move(survivors);
  head.secondary = std::move(secondaries);
  return true;
}

// Applies it to a read, whose head the MAPQ scores on `winner_chain`. Only a
// read left with no record is unmapped.
void apply_dna_min_chain_score(const DnaPlacementCandidateChain* winner_chain,
                               const DnaPlacementChainingResult* placement,
                               dna::Result& realized) {
  for (std::vector<AlignResult>::iterator it = realized.secondary.begin();
       it != realized.secondary.end();) {
    if (keep_dna_min_chain_score(*it, it->score, it->supplementary, nullptr,
                                 nullptr))
      ++it;
    else
      it = realized.secondary.erase(it);
  }
  if (keep_dna_min_chain_score(
          static_cast<AlignResult&>(realized),
          winner_chain != nullptr ? winner_chain->chain_score : realized.score,
          realized.supplementary, &realized.supplementary_candidates,
          placement))
    return;
  demote_unmapped(realized);
  realized.supplementary.clear();
  realized.supplementary_candidates.clear();
  realized.secondary.clear();
}

// The whole-read winner: the best peak of fwd, then rc, by
// chain_peak_better_seeded under the read's tie seed, the first on an exact
// tie. It is catalogue candidate 0.
VotePeak* vote_slope_best_peak(std::vector<VotePeak>& fwd,
                               std::vector<VotePeak>& rc,
                               std::uint32_t tie_seed) {
  VotePeak* best = nullptr;
  for (std::vector<VotePeak>* lane : {&fwd, &rc})
    for (VotePeak& peak : *lane)
      if (best == nullptr || chain_peak_better_seeded(peak, *best, tie_seed))
        best = &peak;
  return best;
}

// Widens the winner's harvest window to its per-read line (vote_slope.h),
// fitted on the seeds its lane's exact refine walked. When the fit passes its
// gate, every chain pass harvests the winner over
// [min(raw, a) - pad, max(raw + L, a + L + stretch(L, b)) + pad]. Nothing else
// about the winner changes, and every other peak keeps the plain window.
void vote_slope_widen_winner(const LongReadSeedContext& seed_ctx,
                             std::uint32_t tie_seed, int span,
                             const ChainAnchorScratch& scratch,
                             std::vector<VotePeak>& fwd,
                             std::vector<VotePeak>& rc) {
  VotePeak* winner = vote_slope_best_peak(fwd, rc, tie_seed);
  if (winner == nullptr)
    return;
  const int lane = winner->is_rc ? 1 : 0;
  const std::uint64_t* chr_bounds = seed_ctx.index->chrom_offsets_data();
  std::vector<VoteSlopePair> pairs;
  std::vector<std::int64_t> residuals;
  const std::vector<DnaLongSeedView>& refine = scratch.slope_refine_views[lane];
  const VoteSlopeFit fit =
      !refine.empty()
          ? vote_slope_fit_winner(refine, chr_bounds, *winner, span, pairs,
                                  residuals)
          : vote_slope_fit_winner(
                (lane ? scratch.rc : scratch.fwd).retained_seeds, chr_bounds,
                *winner, span, pairs, residuals);
  if (!fit.gate)
    return;
  // b1_q20 is clamped to [-6 %, +10 %], well inside int32.
  const std::int64_t line_end =
      fit.a + vote_slope_stretch(span, static_cast<std::int32_t>(fit.b1_q20));
  winner->harvest_below = static_cast<std::int32_t>(
      std::max<std::int64_t>(0, winner->raw_ref_start - fit.a));
  winner->harvest_above = static_cast<std::int32_t>(
      std::max<std::int64_t>(0, line_end - winner->raw_ref_start));
}

// Map-only's alternative commit. The retained alternative, projected as the
// incumbent is, is promoted when its exact whole-query chain is strictly above
// the committed owner's, through the -c lane's
// commit_dna_alternative_hypothesis; otherwise it is the incumbent's MAPQ-0
// secondary.
struct DnaMapOnlyCommit {
  dna::Result primary;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  bool promoted = false;
  // The promoted record's pieces, for the map-only join.
  std::vector<DnaJoinPiece> pieces;
};

// `incumbent` is mapped and projected from `placement`.
DnaMapOnlyCommit commit_map_only_hypothesis(
    const DnaContext& dctx, const DnaPlacementChainingResult& placement,
    const DnaPlacementChainingResult& alternative_chaining,
    dna::Result incumbent,
    ::fa::cpu::voting::CandidateId incumbent_candidate) {
  DnaMapOnlyCommit out;
  out.primary_candidate = incumbent_candidate;
  if (alternative_chaining.accepted && alternative_chaining.family.valid &&
      alternative_chaining.family.original_candidate_id) {
    DnaFamilyProjectionResult projection =
        project_map_only_placement_family(dctx, alternative_chaining,
                                          incumbent);
    if (projection.committed && projection.output.mapped() &&
        projection.output.supplementary.empty()) {
      std::vector<DnaJoinPiece> pieces = std::move(projection.pieces);
      dna::Result alternative =
          settle_map_only_projection(std::move(projection));
      const DnaPlacementCandidateChain* incumbent_chain =
          placement.find(incumbent_candidate);
      const int incumbent_score =
          incumbent_chain != nullptr ? incumbent_chain->chain_score : 0;
      // The restricted family's only candidate has solver id 0.
      const DnaPlacementCandidateChain* alternative_chain =
          alternative_chaining.find(0);
      if (alternative_chain != nullptr &&
          alternative_chain->chain_score > incumbent_score) {
        DnaAlternativeCommit committed = commit_dna_alternative_hypothesis(
            std::move(incumbent), incumbent_candidate, incumbent_score,
            std::move(alternative),
            *alternative_chaining.family.original_candidate_id,
            alternative_chain->chain_score);
        out.primary = std::move(committed.primary);
        out.primary_candidate = committed.primary_candidate;
        out.promoted = committed.promoted;
        out.pieces = std::move(pieces);
        return out;
      }
      zero_hypothesis_mapq(alternative);
      incumbent.secondary.push_back(
          static_cast<AlignResult&&>(std::move(alternative)));
    }
  }
  out.primary = std::move(incumbent);
  return out;
}

} // namespace

AlignResult map_read(const DnaContext& base_dctx,
                     const LongReadSeedContext& seed_ctx,
                     DnaWorkerScratch& worker_scratch, const std::string& read,
                     const std::vector<uint8_t>& fwd_enc,
                     std::vector<uint8_t>& rc_enc,
                     const std::vector<QuerySeed>* shared_fwd_syncmer_seeds) {
  DnaContext dctx = base_dctx;
  auto effective_chain_segment_len = [&](int rl) {
    return ::fa::cpu::lr::effective_chain_segment_len(dctx, rl);
  };
  auto seed_window_estimate = [&]() {
    const int eff_s = effective_closed_syncmer_s(std::max(1, dctx.opts.k),
                                                 dctx.opts.chain_syncmer_s);
    return std::max(1, dctx.opts.k - eff_s + 1);
  };
  const int read_len = static_cast<int>(read.size());
  const int segment_len = effective_chain_segment_len(read_len);
  // The read's tie seed, from the name hash and the read length as minimap2;
  // every tie-break of the read uses it.
  dctx.vote_tie_seed = tie_read_seed(dctx.read_name_hash, read_len);
  dna::Result out{};
  auto finish = [](dna::Result result) {
    return static_cast<AlignResult&&>(std::move(result));
  };
  out.read_len = read_len;
  if (read_len < std::max(2 * dctx.opts.k, segment_len)) {
    return finish(std::move(out));
  }

  ChainAnchorScratch& scratch = worker_scratch.chain;
  worker_scratch.clear_for_read(static_cast<size_t>(read_len));
  auto ensure_rc_enc = [&]() -> const std::vector<uint8_t>& {
    if (rc_enc.size() != fwd_enc.size())
      rc_enc = reverse_complement_encoded_u8(fwd_enc);
    return rc_enc;
  };
  // Only CIGAR realization reads reverse-strand query bases; map-only never
  // builds the buffer and passes kNoQueryBases instead.
  const bool reverse_query_bases_needed = dctx.opts.enable_full_read_cigar;
  auto reverse_query_stream = [&]() -> const std::vector<uint8_t>& {
    return reverse_query_bases_needed ? ensure_rc_enc() : kNoQueryBases;
  };
  auto& raw = scratch.raw;
  raw.reserve(2);
  ChainSeedLookupCache* lookup_cache = &scratch.lookup_cache;
  const size_t expected_keys =
      static_cast<size_t>(
          std::max(1, read_len / std::max(1, seed_window_estimate()))) *
          2u +
      8u + 64u;
  lookup_cache->reset(expected_keys);
  if (shared_fwd_syncmer_seeds == nullptr ||
      shared_fwd_syncmer_seeds->empty()) {
    return finish(std::move(out));
  }
  DnaQuerySeedPool query_seed_pool(seed_ctx, read_len, *lookup_cache);
  // Per-seed captured views and cache slots, read by position later instead
  // of probing the cache again.
  query_seed_pool.bind_capture_stores(
      &scratch.captured_fwd_views, &scratch.captured_fwd_slots,
      &scratch.captured_rc_views, &scratch.captured_rc_slots);
  query_seed_pool.bind(DnaQuerySeedStrand::Forward, *shared_fwd_syncmer_seeds);
  if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Forward))
    return finish(std::move(out));
  const std::vector<QuerySeed>* shared_rc_syncmer_seeds = nullptr;
  auto ensure_shared_rc_stream = [&]() -> const std::vector<QuerySeed>* {
    if (!scratch.full_syncmer_rc_seeds.empty()) {
      query_seed_pool.bind(DnaQuerySeedStrand::Reverse,
                           scratch.full_syncmer_rc_seeds);
      if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Reverse))
        return nullptr;
      return &scratch.full_syncmer_rc_seeds;
    }
    project_query_seeds_to_rc_coordinates(*shared_fwd_syncmer_seeds, read_len,
                                          dctx.opts.k,
                                          scratch.full_syncmer_rc_seeds);
    query_seed_pool.bind(DnaQuerySeedStrand::Reverse,
                         scratch.full_syncmer_rc_seeds);
    if (!query_seed_pool.capture_posting_views(DnaQuerySeedStrand::Reverse))
      return nullptr;
    return &scratch.full_syncmer_rc_seeds;
  };
  auto prepare_shared_window_views =
      [&](DnaQuerySeedStrand strand, int seed_start, int seed_end,
          int downsample, std::vector<QuerySeed>& window_seeds,
          const std::vector<DnaLongSeedView>*& exact_refine_views)
      -> const std::vector<DnaLongSeedView>* {
    scratch.shadow_rc.clear_for_window();
    DnaLongSeedBundle& bundle = scratch.shadow_rc.seed_bundle;
    const bool selected = query_seed_pool.prepare_representatives(
        strand, seed_start, seed_end, downsample, window_seeds, bundle,
        scratch.shadow_rc);
    if (!selected)
      return nullptr;
    exact_refine_views = bundle.exact_refine_views.empty()
                             ? nullptr
                             : &bundle.exact_refine_views;
    return bundle.views_ready ? &bundle.selected_views : nullptr;
  };
  // The vote runs over one whole-read window [0, read_len).
  if (read_len >= dctx.opts.k) {
    const int start = 0;
    const int end = read_len;
    const int span = read_len;
    const int actual_downsample =
        std::max(1, dctx.opts.chain_syncmer_downsample);
    const std::vector<DnaLongSeedView>* fwd_seed_views = nullptr;
    const std::vector<DnaLongSeedView>* fwd_exact_refine_views = nullptr;
    fwd_seed_views = prepare_shared_window_views(
        DnaQuerySeedStrand::Forward, start, end, actual_downsample,
        scratch.shared_window_fwd_seeds, fwd_exact_refine_views);
    if (!fwd_seed_views)
      return finish(std::move(out));
    scratch.fwd.pending_exact_refine_views =
        fwd_exact_refine_views ? *fwd_exact_refine_views
                               : std::vector<DnaLongSeedView>{};
    // Kept for vote_slope_widen_winner: the reverse prepare below clears the
    // bundle these views live in.
    scratch.slope_refine_views[0] = scratch.fwd.pending_exact_refine_views;
    WindowAnchorPeakParams fwd_peak_params;
    fwd_peak_params.lookup_cache = lookup_cache;
    fwd_peak_params.seed_view_override = fwd_seed_views;
    fwd_peak_params.tie_seed = dctx.vote_tie_seed;
    auto fwd =
        window_anchor_peaks(seed_ctx, fwd_enc.data() + start, span,
                            /*is_rc=*/false, scratch.fwd, fwd_peak_params);

    const int rc_start = read_len - end;
    const std::vector<DnaLongSeedView>* rc_seed_views = nullptr;
    const std::vector<DnaLongSeedView>* rc_exact_refine_views = nullptr;
    shared_rc_syncmer_seeds = ensure_shared_rc_stream();
    if (!shared_rc_syncmer_seeds)
      return finish(std::move(out));
    rc_seed_views = prepare_shared_window_views(
        DnaQuerySeedStrand::Reverse, rc_start, rc_start + span,
        actual_downsample, scratch.shared_window_rc_seeds,
        rc_exact_refine_views);
    if (!rc_seed_views)
      return finish(std::move(out));
    scratch.rc.pending_exact_refine_views =
        rc_exact_refine_views ? *rc_exact_refine_views
                              : std::vector<DnaLongSeedView>{};
    scratch.slope_refine_views[1] = scratch.rc.pending_exact_refine_views;
    WindowAnchorPeakParams rc_peak_params;
    rc_peak_params.lookup_cache = lookup_cache;
    rc_peak_params.seed_view_override = rc_seed_views;
    rc_peak_params.tie_seed = dctx.vote_tie_seed;
    auto rc = window_anchor_peaks(seed_ctx, nullptr, span, /*is_rc=*/true,
                                  scratch.rc, rc_peak_params);
    // The winner's harvest window, widened to its line.
    vote_slope_widen_winner(seed_ctx, dctx.vote_tie_seed, span, scratch, fwd,
                            rc);
    raw.insert(raw.end(), fwd.begin(), fwd.end());
    raw.insert(raw.end(), rc.begin(), rc.end());
  }
  if (raw.empty())
    return finish(std::move(out));
  // The placement family comes from the whole-read vote alone.
  std::vector<VotePeak> selected_chain;
  // Only whether a best peak exists matters here; the family's own order
  // decides catalogue candidate 0.
  const VotePeak* best = best_peak(raw);
  if (best)
    selected_chain.push_back(*best);
  DnaPlacementFamily placement_family;
  DnaPlacementChainingResult placement_chaining;
  DnaPlacementChainingResult alternative_chaining;
  bool placement_chaining_ran = false;
  // The family before chaining, kept when the span rule widened an owner, for
  // the retry under the anchor-tile rule below. Invalid otherwise.
  DnaPlacementFamily unspanned_retry_family;
  const auto chain_placement_family = [&](DnaPlacementFamily family,
                                          DnaTileOwnership ownership) {
    return build_dna_placement_chains(
        dctx, std::move(family), fwd_enc, reverse_query_stream(),
        &scratch.fwd.retained_seeds, &scratch.rc.retained_seeds,
        shared_fwd_syncmer_seeds, shared_rc_syncmer_seeds, lookup_cache,
        query_seed_pool.captured_slots(DnaQuerySeedStrand::Forward),
        query_seed_pool.captured_slots(DnaQuerySeedStrand::Reverse),
        ownership);
  };
  if (best) {
    placement_family =
        build_dna_placement_family(dctx, fwd_enc, raw, scratch.fwd, scratch.rc);
    if (placement_family.valid) {
      const std::vector<uint8_t>& reverse_query = reverse_query_stream();
      placement_chaining =
          chain_placement_family(placement_family, DnaTileOwnership::Span);
      placement_chaining_ran = true;
      if (placement_chaining.span_widened)
        unspanned_retry_family = std::move(placement_family);
      placement_family = placement_chaining.family;
      alternative_chaining = build_dna_alternative_placement(
          dctx, placement_family, placement_chaining, fwd_enc, reverse_query);
    }
  }
  const auto& chain = selected_chain;
  if (chain.empty()) {
    return finish(std::move(out));
  }
  dna::Result realized;
  // Catalogue id 0 is the top-ranked vote candidate.
  ::fa::cpu::voting::CandidateId primary_candidate = 0;
  bool alternative_realization_accepted = false;
  bool alternative_promoted = false;
  // minimap2's dp_max / dp_max2 for this read, and the catalogue candidate
  // dp_max2 belongs to: the summed ksw2 scores of the primary and the retained
  // alternative.
  double chain_mapq_dp1 = 0.0;
  double chain_mapq_dp2 = 0.0;
  // The same pair as raw ksw2 scores, for the HiFi margin rule. Without
  // rescoring they equal the pair above; with it the pair above holds the
  // post-DP prices.
  double chain_mapq_dp1_raw = 0.0;
  double chain_mapq_dp2_raw = 0.0;
  // The primary/alternative interval-matched margin (chain_mapq.h R1), winner
  // minus rival; valid only when the two primaries' query spans intersect.
  // HiFi preset only.
  bool chain_mapq_alternative_im_valid = false;
  double chain_mapq_alternative_im_margin = 0.0;
  ::fa::cpu::voting::CandidateId chain_mapq_dp2_owner =
      ::fa::cpu::voting::kNullCandidate;
  // The post-DP prices. `ran` stays false without rescoring and when the
  // incumbent family was refused; consumers then read the raw decision scores.
  DnaPostDpOutcome postdp;
  // The block records demoted after the primary was extended, with their
  // owners, slots and the primary's pre-extension record. They join
  // `realized.secondary` after MAPQ is scored on the family they came from,
  // and leave with MAPQ 0.
  DnaDemotedBlocks demoted;
  if (dctx.opts.enable_full_read_cigar) {
    const std::vector<uint8_t>& reverse_query = ensure_rc_enc();
    DnaFamilyRealizationRequest request;
    request.family = &placement_family;
    request.placement = placement_chaining_ran ? &placement_chaining : nullptr;
    request.forward_query = &fwd_enc;
    request.reverse_query = &reverse_query;
    // Only the primary family's realization extends and demotes; the
    // alternative, the MAPQ rivals and terminal-clip recovery build their own
    // requests.
    request.primary_family = true;
    DnaFamilyRealizationOutcome family_outcome =
        realize_full_cigar_family(dctx, request);
    // Unspanned retry: the span rule can widen an owner into a block its
    // realization cannot fill (for example one whose interval holds none of
    // its chain's anchors). A widened family whose primary realization is not
    // accepted is chained once more under the anchor-tile rule and realized
    // again; the rest of the read (alternative, MAPQ evidence, terminal-clip
    // recovery) then uses the re-chained family. The request holds pointers,
    // so the second realization sees it. At most once per read.
    if (!family_outcome.accepted() && unspanned_retry_family.valid) {
      placement_chaining = chain_placement_family(
          std::move(unspanned_retry_family), DnaTileOwnership::AnchorTiles);
      placement_family = placement_chaining.family;
      alternative_chaining = build_dna_alternative_placement(
          dctx, placement_family, placement_chaining, fwd_enc, reverse_query);
      family_outcome = realize_full_cigar_family(dctx, request);
    }
    DnaFamilyRealizationOutcome alternative_outcome;
    if (family_outcome.accepted() && alternative_chaining.accepted &&
        alternative_chaining.family.valid &&
        alternative_chaining.family.original_candidate_id) {
      DnaFamilyRealizationRequest alternative_request;
      alternative_request.family = &alternative_chaining.family;
      alternative_request.placement = &alternative_chaining;
      alternative_request.forward_query = &fwd_enc;
      alternative_request.reverse_query = &reverse_query;
      alternative_outcome =
          realize_full_cigar_family(dctx, alternative_request);
      alternative_realization_accepted = alternative_outcome.accepted();
    }
    primary_candidate = family_outcome.primary_candidate;
    // Price both realized hypotheses before anything is moved out of them,
    // with the primary as it was before its extension. The stage declines
    // without rescoring or when the incumbent was refused.
    if (family_outcome.accepted()) {
      if (family_outcome.primary_extended)
        exchange_alignment_fields(family_outcome.output,
                                  family_outcome.pre_extension_primary);
      postdp = run_dna_postdp_scoring(
          dctx, family_outcome,
          alternative_realization_accepted ? &alternative_outcome : nullptr,
          fwd_enc, reverse_query, read_len);
      if (family_outcome.primary_extended)
        exchange_alignment_fields(family_outcome.output,
                                  family_outcome.pre_extension_primary);
    }
    // The interval-matched margin (chain_mapq.h R1): one CIGAR walk per
    // realized record while both outputs are whole. HiFi preset only.
    DnaChainMapqAlternativeInterval alternative_interval;
    const bool interval_wanted = dctx.opts.chain_mapq_hifi_margin;
    if (family_outcome.accepted() && alternative_realization_accepted &&
        interval_wanted)
      alternative_interval = chain_mapq_alternative_interval(
          dctx, family_outcome, alternative_outcome, fwd_enc, reverse_query,
          read_len);
    realized = std::move(family_outcome.output);
    if (alternative_outcome.accepted()) {
      // Promotion is decided on the raw decision score; post-DP prices feed
      // only the MAPQ.
      DnaAlternativeCommit committed = commit_dna_alternative_hypothesis(
          std::move(realized), primary_candidate, family_outcome.decision_score,
          std::move(alternative_outcome.output),
          *alternative_chaining.family.original_candidate_id,
          alternative_outcome.decision_score);
      alternative_promoted = committed.promoted;
      realized = std::move(committed.primary);
      primary_candidate = committed.primary_candidate;
    }
    // A family that demoted has more than one block, so it never retains an
    // alternative; the guard only states that a promoted alternative's family
    // is the one emitted.
    if (!alternative_promoted) {
      demoted.records = std::move(family_outcome.demoted);
      demoted.candidates = std::move(family_outcome.demoted_candidates);
      demoted.positions = std::move(family_outcome.demoted_positions);
      demoted.primary_extended = family_outcome.primary_extended;
      demoted.pre_extension_primary =
          std::move(family_outcome.pre_extension_primary);
    }
    realized.read_len = read_len;
    // Promotion swaps the roles: the promoted alternative owns dp_max and
    // the demoted incumbent owns dp_max2.
    chain_mapq_dp1 = static_cast<double>(
        alternative_promoted ? alternative_outcome.decision_score
                             : family_outcome.decision_score);
    if (alternative_realization_accepted) {
      chain_mapq_dp2 = static_cast<double>(
          alternative_promoted ? family_outcome.decision_score
                               : alternative_outcome.decision_score);
      chain_mapq_dp2_owner =
          alternative_promoted
              ? family_outcome.primary_candidate
              : *alternative_chaining.family.original_candidate_id;
    }
    chain_mapq_dp1_raw = chain_mapq_dp1;
    chain_mapq_dp2_raw = chain_mapq_dp2;
    // The incumbent's inside score minus the alternative's, negated after a
    // promotion so that it reads winner minus rival.
    if (dctx.opts.chain_mapq_hifi_margin && alternative_realization_accepted &&
        alternative_interval.scored &&
        alternative_interval.hi > alternative_interval.lo) {
      chain_mapq_alternative_im_valid = true;
      const double inside = static_cast<double>(
          alternative_interval.a_inside - alternative_interval.b_inside);
      chain_mapq_alternative_im_margin =
          alternative_promoted ? -inside : inside;
    }
    // The post-DP pair, when rescoring ran. It is priced in the incumbent's
    // frame (dp1 = the incumbent over the alternative's query span, dp2 = the
    // alternative primary's own price), so a promotion swaps it like the raw
    // pair. The dp2 owner does not depend on the scoring.
    if (postdp.ran) {
      chain_mapq_dp1 = static_cast<double>(
          alternative_promoted ? postdp.mapq_dp2 : postdp.mapq_dp1);
      if (alternative_realization_accepted)
        chain_mapq_dp2 = static_cast<double>(
            alternative_promoted ? postdp.mapq_dp1 : postdp.mapq_dp2);
    }
  } else {
    // Map-only: the placement family's projection replaces this carrier
    // wholesale, so no coarse vote geometry survives. A best peak without
    // vote support leaves the read unmapped without consulting the family.
    realized.read_len = read_len;
    realized.pos = 0;
    realized.score = std::max(0, chain.front().support);
  }
  std::vector<DnaJoinPiece> join_pieces;
  if (realized.mapped()) {
    if (!dctx.opts.enable_full_read_cigar) {
      // Placement chaining always ran here: a best peak exists past the
      // empty-chain return, so the catalogue is non-empty and the family is
      // valid by construction (placement_family_adapter.cpp).
      DnaFamilyProjectionResult projection =
          project_map_only_placement_family(dctx, placement_chaining, realized);
      primary_candidate = projection.primary_candidate;
      join_pieces = std::move(projection.pieces);
      realized = settle_map_only_projection(std::move(projection));
      if (realized.mapped()) {
        const ::fa::cpu::voting::CandidateId incumbent_candidate =
            primary_candidate;
        DnaMapOnlyCommit committed = commit_map_only_hypothesis(
            dctx, placement_chaining, alternative_chaining,
            std::move(realized), incumbent_candidate);
        realized = std::move(committed.primary);
        primary_candidate = committed.primary_candidate;
        if (committed.promoted) {
          // As after the -c lane's promotion: the promoted hypothesis owns
          // the winner chain and the demoted incumbent is the dp2 owner.
          alternative_promoted = true;
          join_pieces = std::move(committed.pieces);
          chain_mapq_dp2_owner = incumbent_candidate;
        }
      }
    }
    if (realized.mapped()) {
      // The selected record family is now committed. Terminal-clip records are
      // appended before MAPQ routing so they inherit read-level confidence
      // without entering selection.
      if (dna_residue_emission_boundary(
              realized.mapped(), placement_chaining_ran, alternative_promoted,
              dctx.opts.enable_full_read_cigar, !realized.cigar.empty())) {
        const std::vector<uint8_t>& reverse_query = reverse_query_stream();
        // Terminal-clip recovery: inspect only terminal query intervals the
        // committed records leave uncovered. It uses cached fine views, admits
        // at most kDnaResidueMaxAdmissionsPerRead records and never re-enters
        // family selection.
        dna_run_terminal_clip_recovery(dctx, placement_chaining, *lookup_cache,
                                       fwd_enc, reverse_query,
                                       dctx.opts.enable_full_read_cigar,
                                       realized);
      }
      // Map-only joins the block records the -c lane would bridge, so the
      // MAPQ below scores a joined record on its owner's chain.
      if (!dctx.opts.enable_full_read_cigar)
        dna_join_map_only_family(dctx, placement_chaining,
                                 std::move(join_pieces), realized,
                                 primary_candidate);
      // Confidence belongs to the committed hypothesis. Both output modes
      // reach this point, so MAPQ cannot affect placement, realization or
      // record-family decisions.
      const DnaPlacementCandidateChain* winner_chain =
          dna_committed_winner_chain(
              placement_chaining_ran ? &placement_chaining : nullptr,
              alternative_promoted ? &alternative_chaining : nullptr,
              primary_candidate, alternative_promoted);
      // The pre-extension family stands in for the MAPQ stage and its
      // routing (see exchange_alignment_fields). Terminal-clip recovery above
      // has already run on the extended family.
      const std::vector<std::size_t> stand_in_slots =
          install_mapq_stand_in(demoted, realized);
      // Confidence is minimap2's formula over the whole-query chains, scored
      // once per block-owning record (chain_mapq.h), after everything else is
      // committed.
      const double identity =
          realized.alignment_accounting_valid && realized.block_len > 0
              ? realized.identity()
              : 1.0;
      // The realized rivals (chain_mapq.h R1): the HiFi rule's extra DP, run
      // only under the HiFi preset with a committed primary CIGAR. The family
      // is final, so rivals are compared with the emitted records.
      DnaChainMapqEvidenceRealizations evidence_realizations;
      evidence_realizations.block_rivals.assign(
          1 + realized.supplementary.size(), DnaChainMapqRealizedRival());
      // The committed family's diagonal map, built once for the shadow rules
      // and the block rivals' own-locus test. HiFi preset only.
      const bool family_map_wanted = dctx.opts.chain_mapq_hifi_margin;
      DnaChainMapqFamilyMap family_map;
      if (family_map_wanted)
        family_map = chain_mapq_family_map(dctx, realized, primary_candidate);
      if (dctx.opts.chain_mapq_hifi_margin &&
          dctx.opts.enable_full_read_cigar && placement_chaining_ran &&
          !realized.cigar.empty()) {
        evidence_realizations = realize_chain_mapq_evidence(
            dctx, placement_family, placement_chaining, family_map,
            winner_chain, primary_candidate, realized, fwd_enc,
            ensure_rc_enc(), read_len);
      }
      // The committed family's block-owning records. The split-read caps
      // (chain_mapq.h R6) read this and each block record's divergence, so
      // the per-record counts are taken only where the caps can fire.
      int family_blocks = 1;
      for (const ::fa::cpu::voting::CandidateId owner :
           realized.supplementary_candidates)
        if (owner != ::fa::cpu::voting::kNullCandidate)
          ++family_blocks;
      const bool census_wanted =
          dctx.opts.chain_mapq_hifi_margin && family_blocks >= 2;
      const auto block_census = [census_wanted](const AlignResult& record) {
        return census_wanted ? chain_mapq_record_census(record)
                             : DnaChainMapqRecordCensus();
      };
      // For kDnaChainMapqStudyRelativeShareCap: each record's (mismatches,
      // aligned bases), over the records the divergence contrast collects, so
      // each block's evidence can carry the family's cleanest other record.
      // Index 0 is the primary.
      std::vector<std::pair<int, int>> family_mismatch_pairs;
      if ((chain_mapq_study_bits(dctx) & kDnaChainMapqStudyRelativeShareCap) !=
              0 &&
          census_wanted) {
        family_mismatch_pairs.assign(1 + realized.supplementary.size(),
                                     std::pair<int, int>(-1, -1));
        for (std::size_t index = 0; index < family_mismatch_pairs.size();
             ++index) {
          const AlignResult& member =
              index == 0 ? static_cast<const AlignResult&>(realized)
                         : realized.supplementary[index - 1];
          if (member.origin == AlignmentOrigin::DnaLocalInversion ||
              !member.alignment_accounting_valid ||
              !clears_dna_emission_floor(member, dctx.opts.cigar_dp_min_dp_max))
            continue;
          const DnaChainMapqRecordCensus member_census =
              chain_mapq_record_census(member);
          if (member_census.mismatches < 0 || member_census.aligned_bases <= 0)
            continue;
          family_mismatch_pairs[index] = {member_census.mismatches,
                                          member_census.aligned_bases};
        }
      }
      const auto rule_fields = [family_blocks, &family_mismatch_pairs](
                                   const DnaChainMapqRecordCensus& census,
                                   DnaChainMapqEvidence& into,
                                   std::size_t record_index) {
        into.family_blocks = family_blocks;
        into.record_mismatches = census.mismatches;
        into.record_ins_events = census.ins_events;
        into.record_del_events = census.del_events;
        into.record_aligned_bases = census.aligned_bases;
        // The cleanest other accounted record by mismatch rate.
        for (std::size_t other = 0; other < family_mismatch_pairs.size();
             ++other) {
          if (other == record_index || family_mismatch_pairs[other].second <= 0)
            continue;
          const bool cleaner =
              into.family_min_aligned <= 0 ||
              static_cast<double>(family_mismatch_pairs[other].first) *
                      static_cast<double>(into.family_min_aligned) <
                  static_cast<double>(into.family_min_mismatches) *
                      static_cast<double>(family_mismatch_pairs[other].second);
          if (cleaner) {
            into.family_min_mismatches = family_mismatch_pairs[other].first;
            into.family_min_aligned = family_mismatch_pairs[other].second;
          }
        }
      };
      DnaChainMapqEvidence evidence = build_chain_mapq_evidence(
          dctx, placement_family, placement_chaining, placement_chaining_ran,
          winner_chain, primary_candidate, chain_mapq_dp2_owner, chain_mapq_dp1,
          chain_mapq_dp2, chain_mapq_dp1_raw, chain_mapq_dp2_raw, identity,
          read_len);
      mark_chain_mapq_record_shadows(dctx, placement_family,
                                     placement_chaining, family_map, 0,
                                     primary_candidate, read_len, evidence);
      evidence.alternative_im_valid = chain_mapq_alternative_im_valid;
      evidence.alternative_im_margin = chain_mapq_alternative_im_margin;
      evidence.sibling = evidence_realizations.sibling;
      evidence.block_rival = evidence_realizations.block_rivals[0];
      const DnaChainMapqRecordCensus primary_census =
          block_census(static_cast<const AlignResult&>(realized));
      rule_fields(primary_census, evidence, 0);
      DnaChainMapqBreakdown breakdown;
      std::vector<DnaChainMapqRivalVerdict> verdicts;
      int mapq = dna_chain_mapq(evidence, &breakdown, &verdicts);
      // s2:i: the competing chain this record's MAPQ weighed. With f1 <= 0
      // the breakdown is default-constructed and "no evidence" is not "no
      // rival", so the record keeps -1 and writes no tag.
      if (evidence.f1 > 0)
        realized.secondary_chain_score = breakdown.f2;
      // Keep the retained alternative only if the MAPQ weighed it as an
      // admissible, non-shadow rival. A twin candidate of the committed locus
      // is a shadow; the writers derive md:i and XA:Z from `secondary` and
      // --secondary yes prints it, so such a twin must not appear as an
      // alternative placement. At most one alternative is retained, so one
      // verdict answers for the whole vector.
      if (!realized.secondary.empty()) {
        const ::fa::cpu::voting::CandidateId alternative_candidate =
            chain_mapq_dp2_owner != ::fa::cpu::voting::kNullCandidate
                ? chain_mapq_dp2_owner
                : alternative_chaining.family.original_candidate_id.value_or(
                      ::fa::cpu::voting::kNullCandidate);
        bool weighed = false;
        for (std::size_t index = 0;
             index < evidence.rivals.size() && index < verdicts.size();
             ++index) {
          if (evidence.rivals[index].candidate != alternative_candidate)
            continue;
          weighed = verdicts[index].admissible && !verdicts[index].shadow;
          break;
        }
        if (!weighed)
          realized.secondary.clear();
      }
      // Every supplementary that owns a selected block is scored on its own
      // candidate's whole-query chain, against the same rivals, with its own
      // realized score as dp1 (none in map-only) and its own identity. One
      // that owns no block (terminal-clip, inversion middle) keeps -1 and
      // inherits the primary's MAPQ during routing. A promoted alternative is
      // a single-record family, so the stable placement's chains suffice.
      std::vector<int> supplementary_mapq(realized.supplementary.size(), -1);
      for (std::size_t index = 0; index < realized.supplementary.size();
           ++index) {
        const ::fa::cpu::voting::CandidateId owner =
            index < realized.supplementary_candidates.size()
                ? realized.supplementary_candidates[index]
                : ::fa::cpu::voting::kNullCandidate;
        if (owner == ::fa::cpu::voting::kNullCandidate ||
            !placement_chaining_ran)
          continue;
        const DnaPlacementCandidateChain* block_chain =
            placement_chaining.find(owner);
        if (block_chain == nullptr)
          continue;
        const AlignResult& record = realized.supplementary[index];
        const double block_dp1 =
            dctx.opts.enable_full_read_cigar && record.score > 0
                ? static_cast<double>(record.score)
                : 0.0;
        const double block_identity =
            record.alignment_accounting_valid && record.block_len > 0
                ? record.identity()
                : 1.0;
        DnaChainMapqEvidence block_evidence = build_chain_mapq_evidence(
            dctx, placement_family, placement_chaining, placement_chaining_ran,
            block_chain, owner, ::fa::cpu::voting::kNullCandidate, block_dp1,
            0.0, block_dp1, 0.0, block_identity, read_len);
        mark_chain_mapq_record_shadows(dctx, placement_family,
                                       placement_chaining, family_map,
                                       index + 1, owner, read_len,
                                       block_evidence);
        // This block's realized rivals (the read-level sibling is never
        // realized on a multi-block family) and the R6 fields.
        block_evidence.sibling = evidence_realizations.sibling;
        block_evidence.block_rival =
            evidence_realizations.block_rivals[index + 1];
        const DnaChainMapqRecordCensus census = block_census(record);
        rule_fields(census, block_evidence, index + 1);
        DnaChainMapqBreakdown block_breakdown;
        std::vector<DnaChainMapqRivalVerdict> block_verdicts;
        supplementary_mapq[index] =
            dna_chain_mapq(block_evidence, &block_breakdown, &block_verdicts);
        // s2:i for this block's record, as for the primary. Written through
        // the vector because `record` is a const reference.
        if (block_evidence.f1 > 0)
          realized.supplementary[index].secondary_chain_score =
              block_breakdown.f2;
      }
      // Divergence contrast cap (chain_mapq.h): a record whose event
      // divergence sits far above the family's minimum is the wrong copy,
      // whatever its chain says. The cap applies to the MAPQ each record would
      // otherwise carry, inherited ones included.
      const DnaFamilyDivergence divergence =
          divergence_contrast(realized, dctx.opts.cigar_dp_min_dp_max);
      if (divergence.records[0].divergent)
        mapq = std::min(mapq, kDnaChainMapqDivergentRecordMapq);
      for (std::size_t index = 0; index < supplementary_mapq.size(); ++index) {
        if (!divergence.records[index + 1].divergent)
          continue;
        const int carried =
            supplementary_mapq[index] >= 0 ? supplementary_mapq[index] : mapq;
        supplementary_mapq[index] =
            std::min(carried, kDnaChainMapqDivergentRecordMapq);
      }
      // The demoted block records become secondaries before routing, which
      // gives every secondary 0. Their copies stay in their slots and the
      // primary stays pre-extension until routing is done, so the
      // inversion-middle rule sees that family too.
      for (AlignResult& record : demoted.records)
        realized.secondary.push_back(std::move(record));
      demoted.records.clear();
      route_dna_mapq(placement_family, primary_candidate, mapq,
                     supplementary_mapq, realized);
      // Restore the extended primary under the MAPQ computed above.
      remove_mapq_stand_in(demoted, stand_in_slots, realized);
    }
  }
  // The emitted family is final here, MAPQ included. Two demotions follow on
  // CIGAR output, and one on map-only, whose records carry chain scores and
  // approximate match counts that the CIGAR floor must not read.
  //
  // (1) A mapped primary with an empty CIGAR is a realization refusal. Emit
  // it unmapped, as minimap2 does for a read it cannot align, rather than as
  // a whole-read match or all-soft-clip. Such a primary has no
  // supplementaries.
  if (dctx.opts.enable_full_read_cigar && realized.mapped() &&
      realized.cigar.empty()) {
    demote_unmapped(realized);
  }
  // (2) minimap2's per-record emission floor (mm_filter_regs): every record
  // owns at least min_chain_score matched bases and min_dp_max (-s) DP score.
  // A failing supplementary or secondary is erased, a failing primary is
  // replaced by the widest surviving supplementary, and a family that keeps
  // nothing becomes unmapped.
  if (dctx.opts.enable_full_read_cigar && realized.mapped()) {
    apply_dna_emission_floor(dctx.opts.cigar_dp_min_dp_max, realized);
  }
  // (3) Map-only: no record under min_chain_score, judged on the chains the
  // MAPQ reads (apply_dna_min_chain_score).
  if (!dctx.opts.enable_full_read_cigar && realized.mapped()) {
    const DnaPlacementChainingResult* placement =
        placement_chaining_ran ? &placement_chaining : nullptr;
    apply_dna_min_chain_score(
        dna_committed_winner_chain(
            placement, alternative_promoted ? &alternative_chaining : nullptr,
            primary_candidate, alternative_promoted),
        placement, realized);
  }
  // ms:i: the max-scoring segment of each emitted record's CIGAR. Stamped
  // after the floor, which can erase records, promote a supplementary or
  // unmap the read (hence the second mapped() test). Map-only records have
  // no CIGAR and keep -1.
  if (dctx.opts.enable_full_read_cigar && realized.mapped()) {
    const std::vector<uint8_t>& reverse_query = ensure_rc_enc();
    auto stamp_dp_max_segment = [&](AlignResult& record) {
      record.dp_max_segment =
          dna_record_dp_max_segment(dctx, record, fwd_enc, reverse_query);
    };
    stamp_dp_max_segment(static_cast<AlignResult&>(realized));
    for (AlignResult& part : realized.supplementary)
      stamp_dp_max_segment(part);
    for (AlignResult& secondary : realized.secondary) {
      stamp_dp_max_segment(secondary);
      for (AlignResult& part : secondary.supplementary)
        stamp_dp_max_segment(part);
    }
  }
  return finish(std::move(realized));
}

} // namespace dna
} // namespace lr
} // namespace cpu
} // namespace fa
