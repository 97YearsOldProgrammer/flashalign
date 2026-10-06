// minimap2's DP terms of one emitted record, read from its CIGAR: the
// maximum-scoring segment under the log-gap cost (mm_update_extra), which is
// the ms:i tag and the -c lane's dp_max, the event identity
// (mm_event_identity) and dp_max recomputed under b2 (mm_update_dp_max).
#pragma once

#ifdef FLASHALIGN_BUILDING_RNA
#error "flashalign_rna may not include the DNA record DP terms"
#endif

#include "context.h"
#include "result.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

// The ms:i value of one record: the maximum-scoring segment of its CIGAR under
// the log-gap cost, as minimap2 computes ms. Returns -1 (no tag) for an
// unmapped record, an empty CIGAR or an unknown contig.
int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc);

// The terms of minimap2's mm_update_dp_max (align.c) for one record, under
// the same pricing as dna_record_dp_max_segment: its event identity
// (mm_event_identity) and its dp_max recomputed under b2 (mm_recal_max_dp, a
// gap costing b2 + log2(1 + len), floored at 0). False, and -1, for a record
// that pricing cannot read.
bool dna_record_event_identity(const DnaContext& context,
                               const AlignResult& record,
                               const std::vector<std::uint8_t>& fwd,
                               const std::vector<std::uint8_t>& rc,
                               double& identity);
int dna_record_recal_dp_max(const DnaContext& context,
                            const AlignResult& record,
                            const std::vector<std::uint8_t>& fwd,
                            const std::vector<std::uint8_t>& rc, double b2);
// mm_update_dp_max's b2 from the best record's event identity.
double dna_rank_b2(double identity, int match, int mismatch) noexcept;

} // namespace fa::cpu::lr
