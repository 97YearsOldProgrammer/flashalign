// The divergence contrast between the records of one read. All records of a
// read share one error rate, so a record whose divergence sits far above the
// family's minimum is likely aligned to the wrong copy of a repeat.
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include the DNA divergence contrast"
#endif

#include "chain_mapq.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace fa::cpu::lr {
namespace dna {

// `accounted` and `divergence` (gap-compressed, as PAF de:f) are inputs: the
// caller sets `accounted` for a record with base-level accounting that clears
// the emission floor and is not an inversion middle. `z` and `divergent` are
// outputs.
struct DnaRecordDivergence {
  bool accounted = false;
  bool divergent = false;
  double divergence = 0.0;
  double z = 0.0;
};

// The family: index 0 is the primary, index i + 1 the i-th supplementary.
struct DnaFamilyDivergence {
  int accounted = 0;
  std::vector<DnaRecordDivergence> records;
};

// Fills `family.accounted` and each accounted record's `z` and `divergent`.
// `denominators` holds each record's event denominator L, in record order.
// Against the accounted record with the lowest divergence, record i scores the
// pooled two-sample rate difference
//   z = (de_i - de_min) / sqrt(p (1 - p) (1 / L_i + 1 / L_min)),
//   p = (de_i L_i + de_min L_min) / (L_i + L_min),
// and is divergent above kDnaChainMapqDivergentRecordZ. With fewer than two
// accounted records nothing is flagged.
inline void dna_divergence_contrast(DnaFamilyDivergence& family,
                                    const std::vector<double>& denominators) {
  family.accounted = 0;
  std::size_t lowest = family.records.size();
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    if (!family.records[index].accounted)
      continue;
    if (lowest == family.records.size() ||
        family.records[index].divergence < family.records[lowest].divergence)
      lowest = index;
    ++family.accounted;
  }
  if (family.accounted < 2)
    return;
  const double de_min = family.records[lowest].divergence;
  const double length_min = denominators[lowest];
  const double events_min = de_min * length_min;
  for (std::size_t index = 0; index < family.records.size(); ++index) {
    DnaRecordDivergence& record = family.records[index];
    if (!record.accounted || index == lowest)
      continue;
    const double length = denominators[index];
    const double events = record.divergence * length;
    const double pooled = (events + events_min) / (length + length_min);
    const double variance = std::max(pooled * (1.0 - pooled), 1e-12) *
                            (1.0 / length + 1.0 / length_min);
    record.z = (record.divergence - de_min) / std::sqrt(variance);
    record.divergent = record.z > kDnaChainMapqDivergentRecordZ;
  }
}

} // namespace dna
} // namespace fa::cpu::lr
