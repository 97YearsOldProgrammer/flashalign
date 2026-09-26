#include "postdp_scoring.h"

#include "../core/cigar.h"
#include "chain_mapq.h" // kDnaChainMapqMaskLevel

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <utility>

namespace fa::cpu::lr {

namespace {

// minimap2's rank_frac. minibwa's trigger: the best hit covers rank_frac of
// the read and the runner-up's span reaches sqrt(rank_frac) of the best's.
constexpr double kRankFrac = 0.9;

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

// minibwa's mb_recal_max_dp: a gap costs b2 * log2(1 + len) and each clipped
// query base adds 1 / b2 mismatches, a net -a per base. Floored at 0.
int rescored_dp_max(const DpMaxSweep& sweep, double b2, int match,
                    int read_len, int aligned_query_span) noexcept {
  const double a = static_cast<double>(match);
  int n_mis = sweep.blen + sweep.n_ambi - sweep.mlen - sweep.n_gap;
  n_mis += truncate_to_int(
      static_cast<double>(std::max(0, read_len - aligned_query_span)) / b2 +
      0.499);
  const int rescored = truncate_to_int(
      a * (static_cast<double>(sweep.mlen) - b2 * n_mis -
           b2 * sweep.log_gap_sum) +
      0.499);
  return std::max(0, rescored);
}

// mm_set_parent's overlap test without the uncovered term: two records
// compete when overlap / shorter span > mask_level. A record without a valid
// span competes with everything.
bool spans_compete(int a_begin, int a_end, int b_begin, int b_end) noexcept {
  const bool a_span = a_begin >= 0 && a_end > a_begin;
  const bool b_span = b_begin >= 0 && b_end > b_begin;
  if (!a_span || !b_span)
    return true;
  const int overlap = std::min(a_end, b_end) - std::max(a_begin, b_begin);
  if (overlap <= 0)
    return false;
  const int shorter = std::min(a_end - a_begin, b_end - b_begin);
  return static_cast<double>(overlap) / static_cast<double>(shorter) >
         kDnaChainMapqMaskLevel;
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

// One realized record's prices and its forward-query span.
struct RecordPrice {
  int dp_sweep = 0;
  int dp_rescored = 0; // 0 when the rescoring did not run
  int dp_used = 0;     // rescored if rescoring ran, else the sweep
  int q_begin = -1;
  int q_end = -1;
  DpMaxSweep sweep;
  // False when the record has no CIGAR or its contig is unknown, so that
  // dna_record_dp_max_segment can report -1 instead of 0.
  bool swept = false;
};

int record_span(const RecordPrice& price) noexcept {
  return std::max(0, price.q_end - price.q_begin);
}

RecordPrice price_record(const DnaContext& context, const AlignResult& record,
                         const std::vector<std::uint8_t>& fwd,
                         const std::vector<std::uint8_t>& rc) {
  RecordPrice price;
  price.q_begin = record.query_start;
  price.q_end = record.query_end;
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
  price.dp_sweep = price.sweep.dp_max;
  price.dp_used = price.dp_sweep;
  return price;
}

// Prices a family: the primary first, then the supplementaries in order.
void price_family(const DnaContext& context,
                  const DnaFamilyRealizationOutcome& outcome,
                  const std::vector<std::uint8_t>& fwd,
                  const std::vector<std::uint8_t>& rc,
                  std::vector<RecordPrice>& into) {
  into.push_back(price_record(context, outcome.output, fwd, rc));
  for (const AlignResult& part : outcome.output.supplementary)
    into.push_back(price_record(context, part, fwd, rc));
}

} // namespace

int dna_record_dp_max_segment(const DnaContext& context,
                              const AlignResult& record,
                              const std::vector<std::uint8_t>& fwd,
                              const std::vector<std::uint8_t>& rc) {
  // Same pricing as the family sums, so ms and dp1 / dp2 agree.
  if (!record.mapped())
    return -1;
  const RecordPrice price = price_record(context, record, fwd, rc);
  return price.swept ? price.dp_sweep : -1;
}

AffineIntervalScore affine_score_over_query_interval(
    const std::vector<std::pair<int, char>>& ops, const std::uint8_t* query,
    const std::uint8_t* reference, int read_len, bool is_reverse, int lo,
    int hi, const DpScoringParams& dp) noexcept {
  AffineIntervalScore out;
  if (query == nullptr || reference == nullptr || read_len <= 0)
    return out;
  // Magnitudes; the signs are applied below.
  const long long a = std::abs(dp.match);
  const long long b = std::abs(dp.mismatch);
  const long long ambi = std::abs(dp.ambi);
  const auto forward = [&](int oriented) {
    return is_reverse ? read_len - 1 - oriented : oriented;
  };
  const auto inside = [&](int oriented) {
    const int q = forward(oriented);
    return q >= lo && q < hi;
  };
  const auto gap_cost = [&](int len) {
    const long long one =
        static_cast<long long>(dp.gap_open1) +
        static_cast<long long>(dp.gap_extend1) * static_cast<long long>(len);
    const long long two =
        static_cast<long long>(dp.gap_open2) +
        static_cast<long long>(dp.gap_extend2) * static_cast<long long>(len);
    return std::min(one, two);
  };
  const auto credit = [&](bool in, long long delta) {
    (in ? out.inside : out.outside) += delta;
  };
  int qoff = 0;
  std::size_t toff = 0;
  for (const auto& op : ops) {
    const int len = op.first;
    if (len <= 0)
      continue;
    switch (op.second) {
    case 'M':
    case '=':
    case 'X': {
      // Never read past the query, even on a malformed CIGAR.
      if (qoff + len > read_len)
        return out;
      for (int l = 0; l < len; ++l) {
        const int cq = query[static_cast<std::size_t>(qoff + l)];
        const int ct = reference[toff + static_cast<std::size_t>(l)];
        long long s;
        if (cq > 3 || ct > 3)
          s = -ambi;
        else if (cq != ct)
          s = -b;
        else
          s = a;
        const bool in = inside(qoff + l);
        credit(in, s);
        ++(in ? out.aligned_inside : out.aligned_outside);
      }
      qoff += len;
      toff += static_cast<std::size_t>(len);
      break;
    }
    case 'I':
      // One event at the run's midpoint.
      credit(inside(qoff + len / 2), -gap_cost(len));
      qoff += len;
      break;
    case 'D':
      // One event at the query base the deletion precedes in the walk.
      credit(inside(std::min(qoff, read_len - 1)), -gap_cost(len));
      toff += static_cast<std::size_t>(len);
      break;
    case 'N':
      toff += static_cast<std::size_t>(len);
      break;
    case 'S':
      qoff += len;
      break;
    default: // H, P: consume nothing
      break;
    }
  }
  return out;
}

DnaPostDpOutcome
run_dna_postdp_scoring(const DnaContext& context,
                       const DnaFamilyRealizationOutcome& incumbent,
                       const DnaFamilyRealizationOutcome* alternative,
                       const std::vector<std::uint8_t>& fwd,
                       const std::vector<std::uint8_t>& rc, int read_len) {
  DnaPostDpOutcome out;
  if (!context.opts.postdp_rescoring || !incumbent.accepted())
    return out;
  out.ran = true;

  // Sweep every record of both families.
  std::vector<RecordPrice> incumbent_prices;
  std::vector<RecordPrice> alternative_prices;
  price_family(context, incumbent, fwd, rc, incumbent_prices);
  const bool has_alternative =
      alternative != nullptr && alternative->accepted();
  if (has_alternative)
    price_family(context, *alternative, fwd, rc, alternative_prices);
  out.incumbent_records = static_cast<int>(incumbent_prices.size());
  out.alternative_records = static_cast<int>(alternative_prices.size());
  out.incumbent_dp_raw = incumbent.decision_score;
  if (has_alternative)
    out.alternative_dp_raw = alternative->decision_score;
  for (const RecordPrice& price : incumbent_prices)
    out.incumbent_dp_sweep += price.dp_sweep;
  for (const RecordPrice& price : alternative_prices)
    out.alternative_dp_sweep += price.dp_sweep;

  // The rescoring trigger compares the two primaries, as minibwa's max /
  // max2. Without an alternative there is no runner-up and no rescoring.
  if (has_alternative) {
    const RecordPrice& incumbent_primary = incumbent_prices.front();
    const RecordPrice& alternative_primary = alternative_prices.front();
    // A tie keeps the incumbent as the best hit.
    const bool incumbent_best =
        incumbent_primary.dp_sweep >= alternative_primary.dp_sweep;
    const RecordPrice& best =
        incumbent_best ? incumbent_primary : alternative_primary;
    const RecordPrice& runner =
        incumbent_best ? alternative_primary : incumbent_primary;
    out.best_span = record_span(best);
    out.runner_span = record_span(runner);
    // As minibwa: no length floor and no score clause.
    out.triggered =
        static_cast<double>(out.best_span) >=
            static_cast<double>(read_len) * kRankFrac &&
        !(static_cast<double>(out.runner_span) <
          static_cast<double>(out.best_span) * std::sqrt(kRankFrac));
    if (out.triggered) {
      out.best_identity = event_identity(best.sweep);
      out.b2 = rescore_b2(out.best_identity, context.opts.cigar_dp_match,
                          context.opts.cigar_dp_mismatch);
      out.clip_bp = std::max(0, read_len - out.best_span);
      // One b2, from the best hit, for every record of both families.
      const auto rescore_all = [&](std::vector<RecordPrice>& family,
                                   int& total) {
        for (RecordPrice& price : family) {
          price.dp_rescored =
              rescored_dp_max(price.sweep, out.b2, context.opts.cigar_dp_match,
                              read_len, record_span(price));
          price.dp_used = price.dp_rescored;
          total += price.dp_rescored;
        }
      };
      rescore_all(incumbent_prices, out.incumbent_dp_rescored);
      rescore_all(alternative_prices, out.alternative_dp_rescored);
    }
  }

  // The MAPQ's dp1 / dp2.
  if (has_alternative) {
    // minimap2 judges a hit against rivals over its own query span, so dp1
    // sums only the incumbent records that compete with the alternative's
    // span, not supplementaries the rival never contested.
    const RecordPrice& alternative_primary = alternative_prices.front();
    int priced = 0;
    for (const RecordPrice& price : incumbent_prices)
      if (spans_compete(alternative_primary.q_begin, alternative_primary.q_end,
                        price.q_begin, price.q_end))
        priced += price.dp_used;
    // An alternative that competes with no record at all is priced against
    // the primary record, the hit the MAPQ is about.
    out.mapq_dp1 = priced > 0 ? priced : incumbent_prices.front().dp_used;
    out.mapq_dp2 = alternative_primary.dp_used;
  } else {
    out.mapq_dp1 = incumbent_prices.front().dp_used;
  }
  return out;
}

} // namespace fa::cpu::lr
