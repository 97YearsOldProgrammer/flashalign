// Counters for the steps of verified-geometry realization. DP cells are not
// counted here; they go to ksw2_attempts / estimated_cells.
#pragma once

#include <cstdint>

namespace fa::cpu::lr {

struct DnaGeometryWork {
  int verified_regions = 0;        // verified regions emitted, no DP
  std::int64_t verified_bases = 0; // their bases
  int seams = 0;                   // fills between two verified neighbours
  int stretches = 0;               // stretches whose first piece fired
  int stretch_pieces = 0;          // pieces fired, one gap fill each
  std::int64_t stretch_bases = 0;  // those stretches' query bases

  // Regions of two or more anchors that went to a stretch because a gap failed
  // the certificate: their anchors, query bases and failing gaps.
  int failed_regions = 0;
  int failed_region_anchors = 0;
  std::int64_t failed_region_bases = 0;
  int failing_gaps = 0;

  // Seams and pieces whose right corner is a long join, filled at
  // minimap2's max(q_span, r_span) band.
  int join_calls = 0;

  void add(const DnaGeometryWork& other) {
    verified_regions += other.verified_regions;
    verified_bases += other.verified_bases;
    seams += other.seams;
    stretches += other.stretches;
    stretch_pieces += other.stretch_pieces;
    stretch_bases += other.stretch_bases;
    failed_regions += other.failed_regions;
    failed_region_anchors += other.failed_region_anchors;
    failed_region_bases += other.failed_region_bases;
    failing_gaps += other.failing_gaps;
    join_calls += other.join_calls;
  }
};

} // namespace fa::cpu::lr
