#pragma once

#include "../core/types.h"
#include "../voting/candidate_catalogue.h"

#include <vector>

namespace fa::cpu::lr {

struct DnaSegmentRecord {
  AlignResult alignment;
  ::fa::cpu::voting::CandidateId candidate = ::fa::cpu::voting::kNullCandidate;
  // The family.block_parts entry, -1 for none.
  int part = -1;
};

struct DnaRecordFamily {
  AlignResult primary;
  ::fa::cpu::voting::CandidateId primary_candidate =
      ::fa::cpu::voting::kNullCandidate;
  int primary_part = -1;
  // Each supplementary keeps the candidate whose block it came from, because
  // MAPQ is scored per block owner.
  std::vector<DnaSegmentRecord> supplementary;
  bool valid = false;
};

// Strict-weak ordering whose greatest element is the primary: widest query
// span, then higher score, then the smallest (query_start, chromosome, pos,
// is_reverse). Shared with the emission floor's promotion.
bool dna_primary_precedence_less(const AlignResult& left,
                                 const AlignResult& right);

// Orders records by query, rejects any query overlap, picks the primary and
// orders the supplementaries. Changes no geometry and assigns no MAPQ.
DnaRecordFamily
assemble_dna_record_family(std::vector<DnaSegmentRecord> records);

} // namespace fa::cpu::lr
