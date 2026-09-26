// Aligned-block accounting shared by the PAF, SAM and BAM writers: the query
// interval, block length, matches, edit distance and gap counts behind PAF
// columns 3, 4, 10 and 11 and the NM:i and de:f tags.
#pragma once

#include "../core/cigar.h"
#include "../core/types.h"

#include <algorithm>
#include <cstdio>
#include <string>

namespace fa {
namespace cpu {
namespace output {

struct BlockAccounting {
  int qstart = 0;
  int qend = 0;
  int block_len = 0;
  int matches = 0;
  int nm = -1; // -1 => not computable (the record carries no CIGAR)
  int ambiguities = 0; // N bases the replay excluded from block_len/matches
  int gap_bases = 0;
  int gap_opens = 0;
};

// Uses the engine's exact accounting when alignment_accounting_valid is set.
// Otherwise the query interval comes from the CIGAR's clips (mirrored for a
// reverse alignment, whose CIGAR is reference-forward), and NM is
// block_len - matches: mismatches plus inserted and deleted bases.
inline BlockAccounting block_accounting(const AlignResult& result, int qlen) {
  BlockAccounting acc;
  acc.qstart = std::max(0, result.query_start);
  acc.qend = result.alignment_accounting_valid
                 ? std::max(acc.qstart, std::min(qlen, result.query_end))
             : result.query_end > result.query_start
                 ? std::min(qlen, result.query_end)
                 : qlen;
  if (!result.alignment_accounting_valid && !result.cigar.empty()) {
    const auto qi = forward_query_interval_from_cigar(result.cigar, qlen,
                                                      result.is_reverse);
    if (qi.valid) {
      acc.qstart = qi.q_lo;
      acc.qend = qi.q_hi;
    }
  }
  acc.block_len = result.alignment_accounting_valid ? result.block_len
                  : result.block_len > 0 ? result.block_len
                                         : std::max(0, acc.qend - acc.qstart);
  acc.matches = result.alignment_accounting_valid ? result.matches
                : result.matches > 0
                    ? result.matches
                    : std::max(0, std::min(acc.block_len, result.score));
  acc.ambiguities = std::max(0, result.ambiguities);
  if (result.cigar.empty())
    return acc;

  int cigar_block = 0;
  int cigar_matches = 0;
  const auto ops = parse_cigar_ops(result.cigar);
  for (const auto& op : ops) {
    const int len = op.first;
    const char code = op.second;
    if (code == 'M' || code == '=' || code == 'X' || code == 'I' ||
        code == 'D') {
      cigar_block += len;
    }
    if (code == 'M' || code == '=') {
      cigar_matches += len;
    }
    if (code == 'I' || code == 'D') {
      acc.gap_bases += len;
      ++acc.gap_opens;
    }
  }
  if (!result.alignment_accounting_valid && cigar_block > 0)
    acc.block_len = cigar_block;
  // Counting every M as a match overstates identity, so this is only a
  // fallback when the producer gave no match count.
  if (!result.alignment_accounting_valid && cigar_matches > 0 &&
      result.matches <= 0)
    acc.matches = cigar_matches;
  acc.nm = result.alignment_accounting_valid
               ? result.edit_distance
               : std::max(0, acc.block_len - acc.matches);
  return acc;
}

// minimap2's de:f, gap-compressed divergence: each indel counts once. As
// mm_event_identity, identity is mlen / (blen + n_ambi - n_gap + n_gapo), with
// ambiguous bases added back to the denominator. False when there is no CIGAR
// or the denominator is empty; then no tag is written.
inline bool event_divergence(const BlockAccounting& acc, double& divergence) {
  if (acc.nm < 0)
    return false;
  const int denominator =
      acc.block_len + acc.ambiguities - acc.gap_bases + acc.gap_opens;
  if (denominator <= 0)
    return false;
  divergence =
      1.0 - static_cast<double>(acc.matches) / static_cast<double>(denominator);
  return true;
}

inline std::string format_paf_de(double de) {
  de = std::max(0.0, std::min(1.0, de));
  if (de == 0.0)
    return "0";
  char buf[16];
  std::snprintf(buf, sizeof(buf), "%.4f", de);
  return std::string(buf);
}

} // namespace output
} // namespace cpu
} // namespace fa
