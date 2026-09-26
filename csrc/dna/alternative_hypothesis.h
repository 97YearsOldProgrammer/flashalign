#pragma once

#include "chain_mapq.h" // DnaChainMapqRealizeRefusal
#include "placement_family_adapter.h"
#include "result.h"

#include <cstdint>
#include <functional>

namespace fa::cpu::lr {

struct DnaContext;
struct DnaPlacementCandidateChain;
struct DnaPlacementChainingResult;

enum class DnaAlternativeRefusal : std::uint8_t {
  None,
  NoStableFamily,
  NotSingleBlock,
  NoCandidate,
  NoLocus,
  NoChain,
  ShadowOnly,
  DisjointOnly,
  IncredibleOnly,
  ExactRestoreFailed,
  StabilizationFailed,
  ProjectionFailed,
  RealizationFailed,
};

struct DnaAlternativeSelection {
  ::fa::cpu::voting::CandidateId incumbent =
      ::fa::cpu::voting::kNullCandidate;
  ::fa::cpu::voting::CandidateId candidate =
      ::fa::cpu::voting::kNullCandidate;
  DnaAlternativeRefusal refusal = DnaAlternativeRefusal::NoStableFamily;
};

// Picks at most one alternative to the owner of a single-block family, from
// the screening chains. A candidate qualifies when it is neither selected nor
// residue-admitted, is not a shadow of a selected locus (same contig and strand,
// vote start within max(2 kb, read length / 10)), overlaps the block by at
// least half the shorter query span, and scores at least 4/5 of the owner or
// within 2k of it. Best by (score, anchors, overlap, lowest id).
DnaAlternativeSelection select_dna_alternative_hypothesis(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable);

// The best rival of one block of a family, by the same rules applied to
// whole-query chains: the candidate owns no selected block and is not
// residue-admitted; its accepted whole-query chain competes with the block's
// forward query span [block_q_begin, block_q_end) under
// dna_chain_mapq_spans_compete; its start diagonal is not within the shadow
// window of `block_diagonal` on the same contig and strand; and it is credible
// against `owner_chain_score`. Best by (chain score, anchors, overlap, lowest
// id). `refusal` is None on success, else the furthest stage any candidate
// reached.
struct DnaBlockRivalSelection {
  ::fa::cpu::voting::CandidateId candidate =
      ::fa::cpu::voting::kNullCandidate;
  const DnaPlacementCandidateChain* chain = nullptr;
  int overlap = 0;
  DnaChainMapqRealizeRefusal refusal = DnaChainMapqRealizeRefusal::NoRival;
};

// `own_locus`, when set, is asked last: a candidate it accepts is the read's
// own placement rather than a rival and is skipped (refusal OwnLocus).
DnaBlockRivalSelection select_dna_block_rival(
    const DnaContext& context, const DnaPlacementFamily& family,
    const DnaPlacementChainingResult& stable, int owner_chain_score,
    int block_q_begin, int block_q_end, int block_contig, bool block_reverse,
    std::int64_t block_diagonal,
    const std::function<bool(::fa::cpu::voting::CandidateId,
                             const DnaPlacementCandidateChain&)>* own_locus =
        nullptr);

struct DnaAlternativeCommit {
  dna::Result primary;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  bool retained = false;
  bool promoted = false;
};

// CIGAR output: the higher decision score becomes primary (a tie keeps the
// incumbent) and the other is kept as a secondary with MAPQ 0. Each record's
// emitted score stays the score of its own CIGAR.
DnaAlternativeCommit commit_dna_alternative_hypothesis(
    dna::Result incumbent,
    ::fa::cpu::voting::CandidateId incumbent_candidate,
    int incumbent_decision_score,
    dna::Result alternative,
    ::fa::cpu::voting::CandidateId alternative_candidate,
    int alternative_decision_score);

}  // namespace fa::cpu::lr
