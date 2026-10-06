#include "postdp_scoring.h"

#include "../core/cigar.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <utility>

namespace fa::cpu::lr {

namespace {

// C's (int) conversion truncates toward zero; saturate instead of UB.
int truncate_to_int(double value) noexcept {
  if (!(value == value))
    return 0;
  if (value >= 2147483647.0)
    return 2147483647;
  if (value <= -2147483648.0)
    return -2147483648;
  return static_cast<int>(value);
}

// The terms minimap2's mm_update_extra derives from one alignment. Bases are
// encoded 0..3; above 3 is ambiguous.
struct DpMaxSweep {
  int dp_max = 0; // max-scoring sub-segment, log-gap cost, (int)(max + .499)
  int mlen = 0;   // matched bases, ambiguous excluded
  int blen = 0;   // aligned block, ambiguous excluded
  int n_ambi = 0;
  int n_gap = 0;  // gap bases
  int n_gapo = 0; // gap opens
  // Sum over gaps of log2(1 + len): minibwa's gap cost is b2 * this.
  double log_gap_sum = 0.0;
};

// mm_update_extra's sweep with log gaps: add the substitution score per
// aligned base, subtract q + e * log2(1 + len) per gap, reset at 0 and keep
// the maximum. `query` is the oriented query from base 0 and `reference`
// starts at the record's pos.
DpMaxSweep dp_max_sweep(const std::vector<std::pair<int, char>>& ops,
                        const std::uint8_t* query,
                        const std::uint8_t* reference, int match, int mismatch,
                        int gap_open, int gap_extend) noexcept {
  DpMaxSweep out;
  if (query == nullptr || reference == nullptr)
    return out;
  // Signs as in ksw_gen_simple_mat.
  const double a = std::abs(match);
  const double b = -std::abs(mismatch);
  const double ambi = -1.0;
  const double q = static_cast<double>(gap_open);
  const double e = static_cast<double>(gap_extend);
  double s = 0.0, max = 0.0;
  std::size_t qoff = 0, toff = 0;
  for (const auto& op : ops) {
    const int len = op.first;
    if (len <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X': {
      int n_ambi = 0, n_diff = 0;
      for (int l = 0; l < len; ++l) {
        const int cq = query[qoff + static_cast<std::size_t>(l)];
        const int ct = reference[toff + static_cast<std::size_t>(l)];
        if (ct > 3 || cq > 3) {
          ++n_ambi;
          s += ambi;
        } else if (ct != cq) {
          ++n_diff;
          s += b;
        } else {
          s += a;
        }
        if (s < 0)
          s = 0;
        else
          max = max > s ? max : s;
      }
      out.blen += len - n_ambi;
      out.mlen += len - (n_ambi + n_diff);
      out.n_ambi += n_ambi;
      qoff += static_cast<std::size_t>(len);
      toff += static_cast<std::size_t>(len);
      break;
    }
    case 'I': {
      int n_ambi = 0;
      for (int l = 0; l < len; ++l)
        if (query[qoff + static_cast<std::size_t>(l)] > 3)
          ++n_ambi;
      out.blen += len - n_ambi;
      out.n_ambi += n_ambi;
      const double lg = std::log2(1.0 + static_cast<double>(len));
      s -= q + e * lg;
      if (s < 0)
        s = 0;
      out.n_gap += len;
      out.n_gapo += 1;
      out.log_gap_sum += lg;
      qoff += static_cast<std::size_t>(len);
      break;
    }
    case 'D': {
      int n_ambi = 0;
      for (int l = 0; l < len; ++l)
        if (reference[toff + static_cast<std::size_t>(l)] > 3)
          ++n_ambi;
      out.blen += len - n_ambi;
      out.n_ambi += n_ambi;
      const double lg = std::log2(1.0 + static_cast<double>(len));
      s -= q + e * lg;
      if (s < 0)
        s = 0;
      out.n_gap += len;
      out.n_gapo += 1;
      out.log_gap_sum += lg;
      toff += static_cast<std::size_t>(len);
      break;
    }
    case 'N':
      toff += static_cast<std::size_t>(len);
      break;
    case 'S':
      qoff += static_cast<std::size_t>(len);
      break;
    default: // H, P: consume nothing
      break;
    }
  }
  out.dp_max = truncate_to_int(max + 0.499);
  return out;
}

// minimap2's mm_event_identity: mlen / (blen + n_ambi - n_gap + n_gapo), or
// 0 on an empty denominator.
double event_identity(const DpMaxSweep& sweep) noexcept {
  const int denominator =
      sweep.blen + sweep.n_ambi - sweep.n_gap + sweep.n_gapo;
  if (denominator <= 0)
    return 0.0;
  return static_cast<double>(sweep.mlen) / static_cast<double>(denominator);
}

// b2 as in minimap2's mm_update_dp_max: 0.5 / max(0.02, 1 - identity), with
// its b2 * a < b guard.
double rescore_b2(double best_identity, int match, int mismatch) noexcept {
  double div = 1.0 - best_identity;
  if (div < 0.02)
    div = 0.02;
  double b2 = 0.5 / div;
  if (b2 * static_cast<double>(match) < static_cast<double>(mismatch))
    b2 = static_cast<double>(match) / static_cast<double>(mismatch);
  return b2;
}

// The record's contig index by name, or -1 when unknown.
int record_contig_index(const DnaContext& context, const std::string& name) {
  if (context.ref.names == nullptr)
    return -1;
  const auto& names = *context.ref.names;
  for (std::size_t index = 0; index < names.size(); ++index)
    if (names[index] == name)
      return static_cast<int>(index);
  return -1;
}

struct RecordPrice {
  DpMaxSweep sweep;
  // False when the record has no CIGAR or its contig is unknown, so that
  // dna_record_dp_max_segment can report -1 instead of 0.
  bool swept = false;
};

RecordPrice price_record(const DnaContext& context, const AlignResult& record,
                         const std::vector<std::uint8_t>& fwd,
                         const std::vector<std::uint8_t>& rc) {
  RecordPrice price;
  const int contig = record_contig_index(context, record.chromosome);
  if (context.ref.encoded != nullptr && contig >= 0 &&
      static_cast<std::size_t>(contig) < context.ref.encoded->size() &&
      record.pos >= 0 && !record.cigar.empty()) {
    const std::vector<std::uint8_t>& reference =
        (*context.ref.encoded)[static_cast<std::size_t>(contig)];
    const std::vector<std::uint8_t>& query = record.is_reverse ? rc : fwd;
    price.sweep = dp_max_sweep(::fa::cpu::output::parse_cigar_ops(record.cigar),
                               query.data(), reference.data() + record.pos,
                               context.opts.cigar_dp_match,
                               context.opts.cigar_dp_mismatch,
                               context.opts.cigar_dp_gap_open1,
                               context.opts.cigar_dp_gap_extend1);
    price.swept = true;
  }
  return price;
}

} // namespace

int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc) {
  if (!record.mapped())
    return -1;
  const RecordPrice price = price_record(context, record, fwd, rc);
  return price.swept ? price.sweep.dp_max : -1;
}

bool dna_record_event_identity(const DnaContext& context,
                               const AlignResult& record,
                               const std::vector<std::uint8_t>& fwd,
                               const std::vector<std::uint8_t>& rc,
                               double& identity) {
  if (!record.mapped())
    return false;
  const RecordPrice price = price_record(context, record, fwd, rc);
  if (!price.swept)
    return false;
  identity = event_identity(price.sweep);
  return true;
}

int dna_record_recal_dp_max(const DnaContext& context,
                            const AlignResult& record,
                            const std::vector<std::uint8_t>& fwd,
                            const std::vector<std::uint8_t>& rc, double b2) {
  if (!record.mapped())
    return -1;
  const RecordPrice price = price_record(context, record, fwd, rc);
  if (!price.swept)
    return -1;
  const DpMaxSweep& sweep = price.sweep;
  const double gap_cost =
      b2 * static_cast<double>(sweep.n_gapo) + sweep.log_gap_sum;
  const int n_mis = sweep.blen + sweep.n_ambi - sweep.mlen - sweep.n_gap;
  const int recal = truncate_to_int(
      static_cast<double>(std::abs(context.opts.cigar_dp_match)) *
          (static_cast<double>(sweep.mlen) - b2 * n_mis - gap_cost) +
      .499);
  return std::max(0, recal);
}

double dna_rank_b2(double identity, int match, int mismatch) noexcept {
  return rescore_b2(identity, match, mismatch);
}

} // namespace fa::cpu::lr
