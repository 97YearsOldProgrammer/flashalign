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
#include "region_realization.h"
#include "target_chaining.h"
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
// same rivals: candidates that own no selected block and are not the winner.
// Only the primary has a dp2 owner.
DnaChainMapqEvidence build_chain_mapq_evidence(
    const DnaContext& dctx, const DnaPlacementFamily& catalogue,
    const DnaPlacementChainingResult& placement, bool placement_ran,
    const DnaPlacementCandidateChain* winner_chain,
    ::fa::cpu::voting::CandidateId winner_candidate,
    ::fa::cpu::voting::CandidateId dp2_owner, double dp1, double dp2,
    double identity, int read_len) {
  DnaChainMapqEvidence evidence;
  evidence.min_chain_score = dctx.opts.min_chain_score;
  evidence.read_len = read_len;
  evidence.match_sc = dctx.opts.cigar_dp_match;
  // As minimap2: sub_diff = opt->a * 2 + opt->b.
  evidence.sub_diff =
      2 * dctx.opts.cigar_dp_match + dctx.opts.cigar_dp_mismatch;
  evidence.dp1 = dp1;
  evidence.dp2 = dp2;
  evidence.identity = identity;
  // HiFi margin rule: the scoring row its margin is read against.
  evidence.hifi_margin_rule = dctx.opts.chain_mapq_hifi_margin;
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
    if (candidate.id == winner_candidate)
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
  // The dp2 owner always competes, even when it owns a block: the candidate
  // of the record dp_max2 came from.
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

// A record printed from family.block_parts[part] is priced on that part
// (chain_mapq.h part_priced); its rivals stay.
void price_block_part(const DnaPlacementFamily& family, int part,
                      DnaChainMapqEvidence& evidence) {
  if (part < 0 || static_cast<std::size_t>(part) >= family.block_parts.size() ||
      evidence.f1 <= 0)
    return;
  const DnaBlockPart& block_part =
      family.block_parts[static_cast<std::size_t>(part)];
  evidence.f1 = block_part.item_score;
  evidence.cnt = block_part.item_anchors;
  // As minimap2 rounds a split piece's score.
  evidence.sib_f2 = static_cast<int>(block_part.pool_subsc + .499);
  evidence.part_priced = true;
  evidence.part_score = block_part.score;
  evidence.part_anchors = block_part.anchors;
  evidence.part_n_sub = block_part.pool_n_sub;
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

bool clears_dna_emission_floor(const AlignResult& record, int min_chain_score,
                               int min_dp_max);

// Collects each record's divergence for the divergence contrast (index 0 the
// primary, i + 1 the i-th supplementary) and runs it. Inversion middles and
// records below the emission floor take no part.
DnaFamilyDivergence divergence_contrast(const dna::Result& realized,
                                        int min_chain_score, int min_dp_max) {
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
    if (!clears_dna_emission_floor(record, min_chain_score, min_dp_max))
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
    // The candidate this record came from, or null for an inversion middle.
    // A record does not vouch for its own owner.
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
    if (candidate.id == winner_candidate)
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

// The MAPQ of the records region_realization.h emits, through the router:
// each primary record on its own chain with dp1 its dp_max and dp2 its
// dp_max2 (owned by the record that gave it, when that is another catalogue
// candidate's). The region stands in for the block part (price_block_part):
// f1 and cnt its whole chain's, score and cnt its own, sib_f2 and n_sub
// mm_set_parent's over the secondaries from its own candidate's chains. An
// inversion middle takes its flanks' MAPQ.
void price_dna_regions(const DnaContext& dctx,
                       const DnaPlacementFamily& family,
                       const DnaPlacementChainingResult& placement,
                       const std::vector<DnaRegionPrice>& prices, int read_len,
                       dna::Result& realized) {
  const std::size_t count = 1 + realized.supplementary.size();
  if (prices.size() != count)
    return;
  const auto member = [&realized](std::size_t index) -> AlignResult& {
    return index == 0 ? static_cast<AlignResult&>(realized)
                      : realized.supplementary[index - 1];
  };
  const bool family_map_wanted = dctx.opts.chain_mapq_hifi_margin;
  const ::fa::cpu::voting::CandidateId primary_candidate =
      prices[0].candidate >= 0
          ? static_cast<::fa::cpu::voting::CandidateId>(prices[0].candidate)
          : ::fa::cpu::voting::kNullCandidate;
  DnaChainMapqFamilyMap family_map;
  if (family_map_wanted)
    family_map = chain_mapq_family_map(dctx, realized, primary_candidate);
  int family_blocks = 0;
  for (const DnaRegionPrice& price : prices)
    if (!price.inversion)
      ++family_blocks;
  const bool census_wanted =
      dctx.opts.chain_mapq_hifi_margin && family_blocks >= 2;
  std::vector<std::pair<int, int>> family_mismatch_pairs;
  if ((chain_mapq_study_bits(dctx) & kDnaChainMapqStudyRelativeShareCap) !=
          0 &&
      census_wanted) {
    family_mismatch_pairs.assign(count, std::pair<int, int>(-1, -1));
    for (std::size_t index = 0; index < count; ++index) {
      const AlignResult& record = member(index);
      if (record.origin == AlignmentOrigin::DnaLocalInversion ||
          !record.alignment_accounting_valid ||
          !clears_dna_emission_floor(record, dctx.opts.min_chain_score,
                                     dctx.opts.cigar_dp_min_dp_max))
        continue;
      const DnaChainMapqRecordCensus census = chain_mapq_record_census(record);
      if (census.mismatches < 0 || census.aligned_bases <= 0)
        continue;
      family_mismatch_pairs[index] = {census.mismatches, census.aligned_bases};
    }
  }
  std::vector<int> mapqs(count, -1);
  for (std::size_t index = 0; index < count; ++index) {
    const DnaRegionPrice& price = prices[index];
    if (price.inversion)
      continue;
    AlignResult& record = member(index);
    const ::fa::cpu::voting::CandidateId owner =
        price.candidate >= 0
            ? static_cast<::fa::cpu::voting::CandidateId>(price.candidate)
            : ::fa::cpu::voting::kNullCandidate;
    const ::fa::cpu::voting::CandidateId dp2_owner =
        price.dp_max2_candidate >= 0 &&
                price.dp_max2_candidate != price.candidate
            ? static_cast<::fa::cpu::voting::CandidateId>(
                  price.dp_max2_candidate)
            : ::fa::cpu::voting::kNullCandidate;
    const double identity =
        record.alignment_accounting_valid && record.block_len > 0
            ? record.identity()
            : 1.0;
    const double dp1 = static_cast<double>(price.dp_max);
    const double dp2 = static_cast<double>(price.dp_max2);
    DnaChainMapqEvidence evidence = build_chain_mapq_evidence(
        dctx, family, placement, true,
        owner != ::fa::cpu::voting::kNullCandidate ? placement.find(owner)
                                                   : nullptr,
        owner, dp2_owner, dp1, dp2, identity, read_len);
    evidence.winner_q_begin = price.chain_q_begin;
    evidence.winner_q_end = price.chain_q_end;
    evidence.winner_ref_begin = price.chain_ref_begin;
    evidence.winner_ref_end = price.chain_ref_end;
    evidence.dp1_raw = static_cast<double>(price.dp_max0);
    evidence.dp2_raw = static_cast<double>(price.dp_max2_0);
    if (family_map_wanted)
      mark_chain_mapq_record_shadows(dctx, family, placement, family_map, index,
                                     owner, read_len, evidence);
    evidence.f1 = price.score0;
    evidence.cnt = price.anchors0;
    evidence.part_priced = true;
    evidence.part_score = price.score;
    evidence.part_anchors = price.cnt;
    evidence.sib_f2 = static_cast<int>(price.pool_subsc + .499);
    evidence.part_n_sub = price.pool_n_sub;
    const DnaChainMapqRecordCensus census =
        census_wanted ? chain_mapq_record_census(record)
                      : DnaChainMapqRecordCensus();
    evidence.family_blocks = family_blocks;
    evidence.record_mismatches = census.mismatches;
    evidence.record_ins_events = census.ins_events;
    evidence.record_del_events = census.del_events;
    evidence.record_aligned_bases = census.aligned_bases;
    for (std::size_t other = 0; other < family_mismatch_pairs.size(); ++other) {
      if (other == index || family_mismatch_pairs[other].second <= 0)
        continue;
      const bool cleaner =
          evidence.family_min_aligned <= 0 ||
          static_cast<double>(family_mismatch_pairs[other].first) *
                  static_cast<double>(evidence.family_min_aligned) <
              static_cast<double>(evidence.family_min_mismatches) *
                  static_cast<double>(family_mismatch_pairs[other].second);
      if (cleaner) {
        evidence.family_min_mismatches = family_mismatch_pairs[other].first;
        evidence.family_min_aligned = family_mismatch_pairs[other].second;
      }
    }
    DnaChainMapqBreakdown breakdown;
    mapqs[index] = dna_chain_mapq(evidence, &breakdown, nullptr);
    if (evidence.f1 > 0) {
      record.chain_anchors = price.cnt;
      record.chain_score = price.score;
      record.secondary_chain_score = breakdown.f2;
    }
  }
  const DnaFamilyDivergence divergence = divergence_contrast(
      realized, dctx.opts.min_chain_score, dctx.opts.cigar_dp_min_dp_max);
  for (std::size_t index = 0; index < count; ++index) {
    if (!divergence.records[index].divergent)
      continue;
    const int carried = mapqs[index] >= 0 ? mapqs[index] : mapqs[0];
    mapqs[index] = std::min(carried, kDnaChainMapqDivergentRecordMapq);
  }
  const int mapq = std::max(0, mapqs[0]);
  mapqs.erase(mapqs.begin());
  route_dna_mapq(family, primary_candidate, mapq, mapqs, realized);
}

// minimap2's emission floor (mm_filter_regs): no record with matched bases
// under min_chain_score or DP score under min_dp_max (-s). min_dp_max <= 0
// disables the DP half.
bool clears_dna_emission_floor(const AlignResult& record, int min_chain_score,
                               int min_dp_max) {
  if (record.matches < min_chain_score)
    return false;
  return min_dp_max <= 0 || record.score >= min_dp_max;
}

// Map-only's floor. minimap2 keeps no chain under
// min_chain_score (mm_chain_dp) and does not test again what it splits off
// one (mm_split_reg), so each record is judged by the score of the block part
// it prints from. A failing head gives its place to the next record. `owners`
// and `parts` are parallel to `rest`, `head_part` is the head's part. Returns
// false when no record is left.
bool keep_dna_min_chain_score(
    AlignResult& head, int head_chain_score, std::vector<AlignResult>& rest,
    std::vector<::fa::cpu::voting::CandidateId>& owners,
    const DnaPlacementFamily& family, int min_chain_score,
    std::vector<int>& parts, int& head_part) {
  std::size_t kept = 0;
  for (std::size_t index = 0; index < rest.size(); ++index) {
    if (family.block_parts[static_cast<std::size_t>(parts[index])].score <
        min_chain_score)
      continue;
    if (kept != index) {
      rest[kept] = std::move(rest[index]);
      owners[kept] = owners[index];
      parts[kept] = parts[index];
    }
    ++kept;
  }
  rest.erase(rest.begin() + static_cast<std::ptrdiff_t>(kept), rest.end());
  owners.resize(kept);
  parts.resize(kept);
  if (head_chain_score >= min_chain_score)
    return true;
  if (rest.empty())
    return false;
  // `rest` aliases head.supplementary: move everything out first.
  std::vector<AlignResult> survivors = std::move(rest);
  std::vector<AlignResult> secondaries = std::move(head.secondary);
  AlignResult next = std::move(survivors.front());
  survivors.erase(survivors.begin());
  owners.erase(owners.begin());
  head_part = parts.front();
  parts.erase(parts.begin());
  // Read-level fields carry over; the record keeps its own MAPQ.
  next.read_len = head.read_len;
  next.median_occurrence = head.median_occurrence;
  head = std::move(next);
  head.supplementary = std::move(survivors);
  head.secondary = std::move(secondaries);
  return true;
}

// Applies it to a map-only read. Only a read left with no record is unmapped.
void apply_dna_min_chain_score(const DnaPlacementFamily& family,
                               int min_chain_score, dna::Result& realized) {
  if (keep_dna_min_chain_score(
          static_cast<AlignResult&>(realized),
          family.block_parts[static_cast<std::size_t>(realized.primary_part)]
              .score,
          realized.supplementary, realized.supplementary_candidates, family,
          min_chain_score, realized.supplementary_parts,
          realized.primary_part))
    return;
  demote_unmapped(realized);
  realized.supplementary.clear();
  realized.supplementary_candidates.clear();
  realized.supplementary_parts.clear();
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
// [min(raw, a) - pad, max(raw + L, a + L + stretch(L, b)) + pad] and bands it
// on the line (VotePeak::line_a, line_b_q20). Nothing else about the winner
// changes, and every other peak keeps the plain window and band.
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
  winner->line_a = fit.a;
  winner->line_b_q20 = static_cast<std::int32_t>(fit.b1_q20);
  winner->line_gate = true;
}

// The all-chains lane (options/dna_profile.h all_chains). Every chain of
// `chains`, parallel to family.candidates, that reaches min_chain_score is
// projected whole as one record, best score first. A record on the contig and
// strand of a kept one, overlapping its reference span, is the same placement
// and is dropped. The best record is the read's head and the others its
// supplementary records, all MAPQ 0 and tp:A:S, as minimap2's -P prints them.
// A read left with none is unmapped.
void emit_all_chains(const DnaContext& dctx, const DnaPlacementFamily& family,
                     const std::vector<DnaPlacementCandidateChain>& chains,
                     dna::Result& out) {
  std::vector<std::size_t> order;
  for (std::size_t index = 0; index < chains.size(); ++index)
    if (chains[index].status == DnaPlacementChainStatus::Accepted &&
        chains[index].chain_score >= dctx.opts.min_chain_score)
      order.push_back(index);
  // Catalogue order breaks a tie on score.
  std::stable_sort(order.begin(), order.end(),
                   [&chains](std::size_t left, std::size_t right) {
                     return chains[left].chain_score >
                            chains[right].chain_score;
                   });
  // The whole read as one block, so the projection trims nothing.
  ::fa::cpu::voting::QueryBlock whole;
  whole.query_tile_begin = 0;
  whole.query_tile_end = family.tile_count;
  whole.supporting_tiles = family.tile_count;
  std::vector<AlignResult> records;
  std::vector<std::size_t> kept;
  for (const std::size_t index : order) {
    const DnaPlacementCandidate& candidate = family.candidates[index];
    whole.candidate = candidate.id;
    AlignResult record;
    if (!dna_project_chained_block(dctx, family, chains[index], whole, record))
      continue;
    bool duplicate = false;
    for (std::size_t slot = 0; slot < kept.size() && !duplicate; ++slot) {
      const VotePeak& other = family.candidates[kept[slot]].peak;
      duplicate = other.chr == candidate.peak.chr &&
                  other.is_rc == candidate.peak.is_rc &&
                  records[slot].pos < record.target_end &&
                  record.pos < records[slot].target_end;
    }
    if (duplicate)
      continue;
    record.origin = AlignmentOrigin::DnaAllChains;
    record.chain_anchors = chains[index].chain_anchors;
    record.chain_score = chains[index].chain_score;
    records.push_back(std::move(record));
    kept.push_back(index);
  }
  if (records.empty())
    return;
  static_cast<AlignResult&>(out) = std::move(records.front());
  for (std::size_t slot = 1; slot < records.size(); ++slot) {
    out.supplementary.push_back(std::move(records[slot]));
    out.supplementary_candidates.push_back(family.candidates[kept[slot]].id);
  }
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
  if (dctx.opts.all_chains)
    dctx.self_contig = seed_ctx.self_contig;
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
    // The winner's harvest window and band, on its line.
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
  bool placement_chaining_ran = false;
  if (best) {
    placement_family =
        build_dna_placement_family(dctx, fwd_enc, raw, scratch.fwd, scratch.rc);
    // The all-chains lane: the catalogue's chains are the records, and nothing
    // below runs.
    if (dctx.opts.all_chains) {
      if (placement_family.valid)
        emit_all_chains(
            dctx, placement_family,
            build_dna_target_chains(
                dctx, placement_family, fwd_enc, reverse_query_stream(),
                &scratch.fwd.retained_seeds, &scratch.rc.retained_seeds,
                shared_fwd_syncmer_seeds, shared_rc_syncmer_seeds,
                lookup_cache,
                query_seed_pool.captured_slots(DnaQuerySeedStrand::Forward),
                query_seed_pool.captured_slots(DnaQuerySeedStrand::Reverse)),
            out);
      return finish(std::move(out));
    }
    if (placement_family.valid) {
      const std::vector<uint8_t>& reverse_query = reverse_query_stream();
      placement_chaining = build_dna_placement_chains(
          dctx, placement_family, fwd_enc, reverse_query,
          &scratch.fwd.retained_seeds, &scratch.rc.retained_seeds,
          shared_fwd_syncmer_seeds, shared_rc_syncmer_seeds, lookup_cache,
          query_seed_pool.captured_slots(DnaQuerySeedStrand::Forward),
          query_seed_pool.captured_slots(DnaQuerySeedStrand::Reverse),
          worker_scratch.seed_density);
      placement_chaining_ran = true;
      dctx.inversion_gate_seeds = placement_chaining.inversion_gate_seeds;
      placement_family = placement_chaining.family;
    }
  }
  const auto& chain = selected_chain;
  if (chain.empty()) {
    return finish(std::move(out));
  }
  dna::Result realized;
  // Catalogue id 0 is the top-ranked vote candidate.
  ::fa::cpu::voting::CandidateId primary_candidate = 0;
  // The -c lane: the controller's records (region_realization.h) and what
  // the MAPQ reads for each primary.
  std::vector<DnaRegionPrice> region_prices;
  if (dctx.opts.enable_full_read_cigar) {
    const std::vector<uint8_t>& reverse_query = ensure_rc_enc();
    DnaFamilyRealizationRequest request;
    request.family = &placement_family;
    request.placement = placement_chaining_ran ? &placement_chaining : nullptr;
    request.forward_query = &fwd_enc;
    request.reverse_query = &reverse_query;
    DnaRegionOutcome regions = realize_dna_regions(dctx, request);
    DnaFamilyRealizationOutcome& outcome = regions.family;
    primary_candidate = outcome.primary_candidate;
    realized = std::move(outcome.output);
    realized.read_len = read_len;
    region_prices = std::move(regions.prices);
    if (!realized.mapped()) {
      if (regions.floored)
        demote_unmapped(realized);
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
    }
    if (dctx.opts.enable_full_read_cigar) {
      price_dna_regions(dctx, placement_family, placement_chaining,
                        region_prices, read_len, realized);
    } else if (realized.mapped()) {
      // Map-only joins the block records that pass the bridge test, so the
      // MAPQ below scores a joined record on its owner's chain.
      dna_join_map_only_family(dctx, placement_chaining,
                               std::move(join_pieces), realized,
                               primary_candidate);
      // Confidence belongs to the committed hypothesis, so MAPQ cannot affect
      // placement or record-family decisions.
      const DnaPlacementCandidateChain* winner_chain =
          placement_chaining_ran ? placement_chaining.find(primary_candidate)
                                 : nullptr;
      // Confidence is minimap2's formula over the whole-query chains, scored
      // once per block-owning record (chain_mapq.h), after everything else is
      // committed.
      const double identity =
          realized.alignment_accounting_valid && realized.block_len > 0
              ? realized.identity()
              : 1.0;
      // The committed family's diagonal map, built once for the shadow rules.
      // HiFi preset only.
      const bool family_map_wanted = dctx.opts.chain_mapq_hifi_margin;
      DnaChainMapqFamilyMap family_map;
      if (family_map_wanted)
        family_map = chain_mapq_family_map(dctx, realized, primary_candidate);
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
              !clears_dna_emission_floor(member, dctx.opts.min_chain_score,
                                         dctx.opts.cigar_dp_min_dp_max))
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
          winner_chain, primary_candidate, ::fa::cpu::voting::kNullCandidate,
          0.0, 0.0, identity, read_len);
      mark_chain_mapq_record_shadows(dctx, placement_family,
                                     placement_chaining, family_map, 0,
                                     primary_candidate, read_len, evidence);
      price_block_part(placement_family, realized.primary_part, evidence);
      const DnaChainMapqRecordCensus primary_census =
          block_census(static_cast<const AlignResult&>(realized));
      rule_fields(primary_census, evidence, 0);
      DnaChainMapqBreakdown breakdown;
      std::vector<DnaChainMapqRivalVerdict> verdicts;
      int mapq = dna_chain_mapq(evidence, &breakdown, &verdicts);
      // cm:i and s1:i: the block part the record prints from; s2:i: the
      // competing chain this record's MAPQ weighed. With f1 <= 0 there is no
      // chain, the breakdown is default-constructed and "no evidence" is not
      // "no rival", so the record keeps -1 and writes none of the three.
      if (evidence.f1 > 0) {
        const DnaBlockPart& part = placement_family.block_parts
            [static_cast<std::size_t>(realized.primary_part)];
        realized.chain_anchors = part.anchors;
        realized.chain_score = part.score;
        realized.secondary_chain_score = breakdown.f2;
      }
      // Every supplementary that owns a selected block is scored on its own
      // candidate's whole-query chain, against the same rivals. One that owns
      // no block keeps -1 and inherits the primary's MAPQ during routing.
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
        const double block_identity =
            record.alignment_accounting_valid && record.block_len > 0
                ? record.identity()
                : 1.0;
        DnaChainMapqEvidence block_evidence = build_chain_mapq_evidence(
            dctx, placement_family, placement_chaining, placement_chaining_ran,
            block_chain, owner, ::fa::cpu::voting::kNullCandidate, 0.0, 0.0,
            block_identity, read_len);
        mark_chain_mapq_record_shadows(dctx, placement_family,
                                       placement_chaining, family_map,
                                       index + 1, owner, read_len,
                                       block_evidence);
        price_block_part(placement_family,
                         index < realized.supplementary_parts.size()
                             ? realized.supplementary_parts[index]
                             : -1,
                         block_evidence);
        const DnaChainMapqRecordCensus census = block_census(record);
        rule_fields(census, block_evidence, index + 1);
        DnaChainMapqBreakdown block_breakdown;
        std::vector<DnaChainMapqRivalVerdict> block_verdicts;
        supplementary_mapq[index] =
            dna_chain_mapq(block_evidence, &block_breakdown, &block_verdicts);
        // cm:i, s1:i and s2:i for this block's record, as for the primary.
        // Written through the vector because `record` is a const reference.
        if (block_evidence.f1 > 0) {
          AlignResult& stamped = realized.supplementary[index];
          const DnaBlockPart& part = placement_family.block_parts
              [static_cast<std::size_t>(realized.supplementary_parts[index])];
          stamped.chain_anchors = part.anchors;
          stamped.chain_score = part.score;
          stamped.secondary_chain_score = block_breakdown.f2;
        }
      }
      // Divergence contrast cap (chain_mapq.h): a record whose event
      // divergence sits far above the family's minimum is the wrong copy,
      // whatever its chain says. The cap applies to the MAPQ each record would
      // otherwise carry, inherited ones included.
      const DnaFamilyDivergence divergence =
          divergence_contrast(realized, dctx.opts.min_chain_score,
                              dctx.opts.cigar_dp_min_dp_max);
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
      route_dna_mapq(placement_family, primary_candidate, mapq,
                     supplementary_mapq, realized);
    }
  }
  // The emitted family is final here, MAPQ included; the -c lane's emission
  // floor ran in its controller. One demotion follows on CIGAR output, and
  // one on map-only.
  //
  // (1) A mapped primary with an empty CIGAR is a realization refusal. Emit
  // it unmapped, as minimap2 does for a read it cannot align, rather than as
  // a whole-read match or all-soft-clip. Such a primary has no
  // supplementaries.
  if (dctx.opts.enable_full_read_cigar && realized.mapped() &&
      realized.cigar.empty()) {
    demote_unmapped(realized);
  }
  // (2) Map-only: no record under min_chain_score, judged on the block parts
  // the MAPQ reads (apply_dna_min_chain_score).
  if (!dctx.opts.enable_full_read_cigar && realized.mapped())
    apply_dna_min_chain_score(placement_chaining.family,
                              dctx.opts.min_chain_score, realized);
  // ms:i: the max-scoring segment of each emitted record's CIGAR, after the
  // demotion above (hence the second mapped() test). Map-only records have
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
