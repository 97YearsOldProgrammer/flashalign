#pragma once

#include "../core/cigar.h"
#include "../core/types.h"

#include <cstdint>
#include <vector>

namespace fa::cpu::lr {

// Replays a fixed CIGAR and fills the record's coordinates and accounting;
// false if the CIGAR does not replay. `request` asks for cs:Z / MD:Z / =X
// output (minimap2 --cs / --MD / --eqx), built during the same replay; with
// --eqx the CIGAR is replaced by its =/X spelling.
inline bool commit_fixed_cigar_geometry(
    AlignResult& result, const std::vector<std::uint8_t>& query,
    const std::vector<std::uint8_t>& reference,
    const ::fa::cpu::output::CigarReplayRequest& request = {}) {
  if (result.cigar.empty() || result.read_len < 0 ||
      static_cast<int>(query.size()) < result.read_len || result.pos < 0) {
    return false;
  }
  const ::fa::cpu::output::CigarReplay replay = ::fa::cpu::output::replay_cigar(
      result.cigar, query.data(), result.read_len, reference.data(),
      static_cast<int>(reference.size()), result.pos, request);
  const ::fa::cpu::output::ForwardQueryInterval query_interval =
      ::fa::cpu::output::forward_query_interval_from_cigar(
          result.cigar, result.read_len, result.is_reverse);
  if (!replay.valid || replay.query_consumed != result.read_len ||
      replay.target_start != result.pos || !query_interval.valid) {
    return false;
  }

  result.query_start = query_interval.q_lo;
  result.query_end = query_interval.q_hi;
  result.target_end = replay.target_end;
  result.matches = replay.matches;
  result.mismatches = replay.mismatches;
  result.insertions = replay.insertions;
  result.deletions = replay.deletions;
  result.ambiguities = replay.ambiguities;
  result.edit_distance = replay.edit_distance;
  result.block_len = replay.block_len;
  result.cs = replay.cs;
  result.md = replay.md;
  if (request.eqx)
    result.cigar = replay.eqx_cigar;
  result.alignment_accounting_valid = true;
  return true;
}

}  // namespace fa::cpu::lr
