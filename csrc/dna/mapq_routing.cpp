#include "mapq_routing.h"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

namespace fa::cpu::lr {
namespace {

bool inversion_record(const AlignResult& record) noexcept {
  return record.origin == AlignmentOrigin::DnaLocalInversion;
}

// As minimap2's mm_set_inv_mapq: records are ordered by (contig, position)
// and an inversion record between two others takes the smaller of their
// MAPQs. Secondaries take no part.
void apply_inversion_mapq(dna::Result& realized) {
  std::vector<AlignResult*> ordered;
  ordered.reserve(1 + realized.supplementary.size());
  bool has_inversion = false;
  ordered.push_back(static_cast<AlignResult*>(&realized));
  for (AlignResult& supplementary : realized.supplementary) {
    has_inversion = has_inversion || inversion_record(supplementary);
    ordered.push_back(&supplementary);
  }
  if (!has_inversion || ordered.size() < 3)
    return;
  std::stable_sort(ordered.begin(), ordered.end(),
                   [](const AlignResult* left, const AlignResult* right) {
                     return std::tie(left->chromosome, left->pos) <
                            std::tie(right->chromosome, right->pos);
                   });
  for (std::size_t index = 1; index + 1 < ordered.size(); ++index) {
    if (!inversion_record(*ordered[index]))
      continue;
    ordered[index]->mapq =
        std::min(ordered[index - 1]->mapq, ordered[index + 1]->mapq);
  }
}

} // namespace

void route_dna_mapq(const DnaPlacementFamily& catalogue,
                    ::fa::cpu::voting::CandidateId primary_candidate, int mapq,
                    const std::vector<int>& supplementary_mapq,
                    dna::Result& realized) {
  if (!realized.mapped())
    return;

  realized.mapq = mapq;

  const DnaPlacementCandidate* candidate = catalogue.find(primary_candidate);
  realized.median_occurrence =
      candidate == nullptr ? 0 : candidate->peak.anchor.median_occurrence;
  for (std::size_t index = 0; index < realized.supplementary.size(); ++index) {
    AlignResult& supplementary = realized.supplementary[index];
    const int own =
        index < supplementary_mapq.size() ? supplementary_mapq[index] : -1;
    supplementary.mapq =
        inversion_record(supplementary) ? 0 : own >= 0 ? own : mapq;
  }
  for (AlignResult& secondary : realized.secondary) {
    secondary.mapq = 0;
    for (AlignResult& supplementary : secondary.supplementary)
      supplementary.mapq = 0;
  }
  // Last, as in minimap2's mm_set_mapq2.
  apply_inversion_mapq(realized);
}

} // namespace fa::cpu::lr
