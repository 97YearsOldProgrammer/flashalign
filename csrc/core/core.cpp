// CIGAR building, parsing, clipping and replay; reverse complement of a 0-4 encoded sequence.
#include "core.h"

#include <algorithm>
#include <array>

/**** cigar ****/

namespace fa { namespace cpu { namespace output {

void append_cigar_run(std::string& cigar, int count, char op) {
    if (count <= 0) return;
    if (!cigar.empty() && cigar.back() == op) {
        size_t begin = cigar.size() - 1;
        while (begin > 0 && cigar[begin - 1] >= '0' && cigar[begin - 1] <= '9') {
            --begin;
        }
        const int prev = std::stoi(cigar.substr(begin, cigar.size() - 1 - begin));
        cigar.erase(begin);
        cigar += std::to_string(prev + count);
        cigar += op;
        return;
    }
    cigar += std::to_string(count);
    cigar += op;
}

void append_cigar_string(std::string& cigar, const std::string& suffix) {
    int count = 0;
    for (char c : suffix) {
        if (c >= '0' && c <= '9') {
            count = count * 10 + (c - '0');
        } else {
            append_cigar_run(cigar, count, c);
            count = 0;
        }
    }
}

namespace {

// Packed append_cigar_run.
void append_packed_run(PackedCigar& cigar, int count, int op) {
    if (count <= 0) return;
    const std::uint32_t code = static_cast<std::uint32_t>(op) & 0xfu;
    if (!cigar.empty() && (cigar.back() & 0xfu) == code) {
        cigar.back() += static_cast<std::uint32_t>(count) << 4;
        return;
    }
    cigar.push_back((static_cast<std::uint32_t>(count) << 4) | code);
}

} // namespace

void append_packed_cigar(PackedCigar& cigar, const PackedCigar& suffix) {
    for (const std::uint32_t run : suffix)
        append_packed_run(cigar, packed_run_count(run), packed_run_op(run));
}

void prepend_packed_cigar(PackedCigar& cigar, PackedCigar prefix) {
    if (cigar.empty() || packed_run_count(cigar.front()) <= 0) {
        append_packed_cigar(prefix, cigar);
        cigar = std::move(prefix);
        return;
    }
    append_packed_run(prefix, packed_run_count(cigar.front()),
                      packed_run_op(cigar.front()));
    prefix.insert(prefix.end(), cigar.begin() + 1, cigar.end());
    cigar = std::move(prefix);
}

namespace {

constexpr int kCigarMatch = 0;
constexpr int kCigarInsertion = 1;
constexpr int kCigarDeletion = 2;
constexpr int kCigarReferenceSkip = 3;
constexpr std::uint32_t kMaximumCigarRun = UINT32_MAX >> 4;

bool valid_span(int begin, int count, std::size_t length) noexcept {
  return count >= 0 && mapping::checked_count_range(
      begin, static_cast<std::uint64_t>(count),
      static_cast<std::uint64_t>(length));
}

} // namespace

bool normalize_cigar(std::vector<std::uint32_t>& cigar,
                     const std::uint8_t* query,
                     const std::uint8_t* reference,
                     int query_span, int reference_span,
                     bool mapping_reverse,
                     int& projected_query_begin,
                     int& projected_query_end,
                     int& reference_begin,
                     int& query_shift,
                     int& reference_shift) {
  query_shift = reference_shift = 0;
  if (cigar.empty()) return false;
  // As in minimap2, a single-operation CIGAR skips normalization, including the
  // leading-gap rule.
  if (cigar.size() <= 1) return true;
  {
    int reference_offset = 0;
    int query_offset = 0;
    bool shrink = false;
    for (std::size_t index = 0; index < cigar.size(); ++index) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int length = static_cast<int>(cigar[index] >> 4);
      if (length == 0) shrink = true;
      if (operation == kCigarMatch) {
        if (!valid_span(query_offset, length,
                        static_cast<std::size_t>(query_span)) ||
            !valid_span(reference_offset, length,
                        static_cast<std::size_t>(reference_span)))
          return false;
        reference_offset += length;
        query_offset += length;
      } else if (operation == kCigarInsertion ||
                 operation == kCigarDeletion) {
        if (index > 0 && index + 1 < cigar.size() &&
            (cigar[index - 1] & 0xfu) == kCigarMatch &&
            (cigar[index + 1] & 0xfu) == kCigarMatch) {
          const int previous_length =
              static_cast<int>(cigar[index - 1] >> 4);
          int shift = 0;
          if (operation == kCigarInsertion) {
            if (!valid_span(query_offset, length,
                            static_cast<std::size_t>(query_span)))
              return false;
            for (; shift < previous_length; ++shift) {
              if (query[query_offset - 1 - shift] !=
                  query[query_offset + length - 1 - shift])
                break;
            }
          } else {
            if (!valid_span(reference_offset, length,
                            static_cast<std::size_t>(reference_span)))
              return false;
            for (; shift < previous_length; ++shift) {
              if (reference[reference_offset - 1 - shift] !=
                  reference[reference_offset + length - 1 - shift])
                break;
            }
          }
          const std::uint32_t next_length = cigar[index + 1] >> 4;
          if (shift > 0 && static_cast<std::uint32_t>(shift) <=
                  kMaximumCigarRun - next_length) {
            cigar[index - 1] -= static_cast<std::uint32_t>(shift) << 4;
            cigar[index + 1] += static_cast<std::uint32_t>(shift) << 4;
            query_offset -= shift;
            reference_offset -= shift;
          }
          if (shift == previous_length && static_cast<std::uint32_t>(shift) <=
                  kMaximumCigarRun - next_length) shrink = true;
        }
        if (operation == kCigarInsertion) {
          if (!valid_span(query_offset, length,
                          static_cast<std::size_t>(query_span))) return false;
          query_offset += length;
        } else {
          if (!valid_span(reference_offset, length,
                          static_cast<std::size_t>(reference_span))) return false;
          reference_offset += length;
        }
      } else if (operation == kCigarReferenceSkip) {
        if (!valid_span(reference_offset, length,
                        static_cast<std::size_t>(reference_span))) return false;
        reference_offset += length;
      } else {
        return false;
      }
      if (query_offset > query_span || reference_offset > reference_span)
        return false;
    }
    if (query_offset != query_span || reference_offset != reference_span)
      return false;

    for (std::size_t index = 0; index + 2 < cigar.size();) {
      const int operation = static_cast<int>(cigar[index] & 0xfu);
      const int next_operation = static_cast<int>(cigar[index + 1] & 0xfu);
      if (operation > 0 && operation + next_operation == 3) {
        std::array<std::uint32_t, 3> totals{};
        std::size_t limit = index;
        for (; limit < cigar.size(); ++limit) {
          const int gap_operation = static_cast<int>(cigar[limit] & 0xfu);
          const std::uint32_t length = cigar[limit] >> 4;
          if (gap_operation == kCigarInsertion ||
              gap_operation == kCigarDeletion || length == 0) {
            if (gap_operation > kCigarDeletion) return false;
            totals[static_cast<std::size_t>(gap_operation)] += length;
          } else {
            break;
          }
        }
        if (totals[kCigarInsertion] > 0 && totals[kCigarDeletion] > 0 &&
            totals[kCigarInsertion] <= kMaximumCigarRun &&
            totals[kCigarDeletion] <= kMaximumCigarRun &&
            limit - index > 2) {
          cigar[index] = totals[kCigarInsertion] << 4 | kCigarInsertion;
          cigar[index + 1] = totals[kCigarDeletion] << 4 | kCigarDeletion;
          for (std::size_t zero = index + 2; zero < limit; ++zero)
            cigar[zero] &= 0xfu;
          shrink = true;
        }
        index = limit + 1;
      } else {
        ++index;
      }
    }
    if (shrink) {
      cigar.erase(std::remove_if(cigar.begin(), cigar.end(),
                                 [](std::uint32_t encoded) {
                                   return (encoded >> 4) == 0;
                                 }),
                  cigar.end());
      std::vector<std::uint32_t> merged;
      merged.reserve(cigar.size());
      for (std::uint32_t encoded : cigar) {
        const std::uint32_t encoded_length = encoded >> 4;
        const std::uint32_t merged_length = merged.empty()
            ? 0 : merged.back() >> 4;
        if (!merged.empty() &&
            (merged.back() & 0xfu) == (encoded & 0xfu) &&
            encoded_length <= kMaximumCigarRun - merged_length)
          merged.back() += encoded_length << 4;
        else
          merged.push_back(encoded);
      }
      cigar.swap(merged);
    }
  }
  if (cigar.empty()) return false;
  const int leading_operation = static_cast<int>(cigar.front() & 0xfu);
  if (leading_operation == kCigarInsertion ||
      leading_operation == kCigarDeletion) {
    const int length = static_cast<int>(cigar.front() >> 4);
    if (leading_operation == kCigarInsertion) {
      if (mapping_reverse) projected_query_end -= length;
      else projected_query_begin += length;
      query_shift = length;
    } else {
      reference_begin += length;
      reference_shift = length;
    }
    cigar.erase(cigar.begin());
  }
  return !cigar.empty();
}

int packed_cigar_query_consumed(const PackedCigar& cigar) {
    int query = 0;
    for (const std::uint32_t run : cigar) {
        const int count = packed_run_count(run);
        if (count <= 0) return -1;
        switch (packed_run_op(run)) {
        case 0:  // M
        case 1:  // I
        case 4:  // S
        case 7:  // =
        case 8:  // X
            query += count;
            break;
        default:
            break;
        }
    }
    return query;
}

namespace {

// Text rendering of a packed CIGAR.
void append_packed_cigar_text(std::string& out, const PackedCigar& cigar) {
    char digits[12];
    for (const std::uint32_t run : cigar) {
        unsigned count = run >> 4;
        int n = 0;
        do {
            digits[n++] = static_cast<char>('0' + count % 10u);
            count /= 10u;
        } while (count != 0u);
        while (n > 0) out += digits[--n];
        const int op = packed_run_op(run);
        out += kPackedCigarOps[op > 9 ? 0 : op];
    }
}

} // namespace

std::string packed_cigar_text(const PackedCigar& cigar) {
    std::string out;
    out.reserve(cigar.size() * 4);
    append_packed_cigar_text(out, cigar);
    return out;
}

std::vector<std::pair<int, char>> parse_cigar_ops(const std::string& cigar) {
    std::vector<std::pair<int, char>> ops;
    int count = 0;
    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            count = count * 10 + (c - '0');
            continue;
        }
        if (count <= 0) return {};
        ops.push_back({count, c});
        count = 0;
    }
    if (count != 0) return {};
    return ops;
}

TerminalHardClip hard_clip_terminal(const std::string& cigar) {
    TerminalHardClip out;
    auto ops = parse_cigar_ops(cigar);
    if (ops.empty()) { out.cigar = cigar; return out; }
    if (ops.front().second == 'S') { out.lead = ops.front().first; ops.front().second = 'H'; }
    if (ops.size() > 1 && ops.back().second == 'S') {
        out.trail = ops.back().first;
        ops.back().second = 'H';
    }
    for (const auto& [len, op] : ops) append_cigar_run(out.cigar, len, op);
    return out;
}

ForwardQueryInterval forward_query_interval_from_cigar(
    const std::string& cigar, int read_len, bool is_reverse) {
    ForwardQueryInterval out;
    const auto ops = parse_cigar_ops(cigar);
    if (ops.empty() || read_len < 0) return out;
    int lead = 0, trail = 0;
    for (size_t i = 0; i < ops.size(); ++i) {
        const char op = ops[i].second;
        if (op != 'S' && op != 'H') break;
        lead += ops[i].first;
    }
    for (size_t i = ops.size(); i > 0; --i) {
        const char op = ops[i - 1].second;
        if (op != 'S' && op != 'H') break;
        trail += ops[i - 1].first;
    }
    const int q0 = is_reverse ? trail : lead;
    const int q1 = read_len - (is_reverse ? lead : trail);
    if (q0 < 0 || q1 < q0 || q1 > read_len) return out;
    out.q_lo = q0;
    out.q_hi = q1;
    out.valid = true;
    return out;
}

std::string clipped_reference_oriented_cigar(
    const std::string& reference_cigar, int read_len, int forward_q_lo,
    int forward_q_hi, bool is_reverse) {
    if (read_len < 0 || forward_q_lo < 0 || forward_q_hi < forward_q_lo ||
        forward_q_hi > read_len) return {};
    auto ops = parse_cigar_ops(reference_cigar);
    if (ops.empty()) return {};
    while (!ops.empty() && (ops.front().second == 'S' || ops.front().second == 'H'))
        ops.erase(ops.begin());
    while (!ops.empty() && (ops.back().second == 'S' || ops.back().second == 'H'))
        ops.pop_back();
    if (ops.empty()) return {};
    int body_query = 0;
    for (const auto& op : ops) {
        if (op.second == 'M' || op.second == 'I' || op.second == '=' ||
            op.second == 'X') body_query += op.first;
    }
    if (body_query != forward_q_hi - forward_q_lo) return {};
    const auto clips = reference_oriented_terminal_clips(
        read_len, forward_q_lo, forward_q_hi, is_reverse);
    std::string out;
    append_cigar_run(out, clips.lead, 'S');
    for (const auto& op : ops) append_cigar_run(out, op.first, op.second);
    append_cigar_run(out, clips.trail, 'S');
    return out;
}

namespace {

// minimap2's cs/MD base letters.
constexpr char kCsLower[] = "acgtn";
constexpr char kCsUpper[] = "ACGTN";

inline char cs_lower(uint8_t code) { return kCsLower[code > 4 ? 4 : code]; }
inline char cs_upper(uint8_t code) { return kCsUpper[code > 4 ? 4 : code]; }

} // namespace

CigarReplay replay_cigar(const std::string& cigar,
                                const uint8_t* qseq, int qlen,
                                const uint8_t* tseq, int tlen,
                                int ref_start) {
  return replay_cigar(cigar, qseq, qlen, tseq, tlen, ref_start,
                      CigarReplayRequest{});
}

CigarReplay replay_cigar(const std::string& cigar, const uint8_t* qseq,
                         int qlen, const uint8_t* tseq, int tlen, int ref_start,
                         const CigarReplayRequest& request) {
  CigarReplay out;
  const bool want_cs = request.cs != CigarReplayRequest::Cs::None;
  const bool want_cs_long = request.cs == CigarReplayRequest::Cs::Long;
  const bool want_md = request.md;
  const bool want_eqx = request.eqx;
  // Pending cs identical run and the MD match count.
  std::string cs_run;
  int cs_run_len = 0;
  int md_run = 0;
  // No merging: minimap2 rewrites each op separately, so "5M5M" becomes "5=5=", not "10=".
  const auto emit_eqx_run = [&](int count, char op) {
    out.eqx_cigar += std::to_string(count);
    out.eqx_cigar += op;
  };
  const auto flush_cs_run = [&]() {
    if (cs_run_len == 0)
      return;
    if (want_cs_long) {
      out.cs += '=';
      out.cs += cs_run;
    } else {
      out.cs += ':';
      out.cs += std::to_string(cs_run_len);
    }
    cs_run.clear();
    cs_run_len = 0;
  };
  const auto ops = parse_cigar_ops(cigar);
  if (ops.empty() || qseq == nullptr || tseq == nullptr || qlen < 0 ||
      tlen < 0 || ref_start < 0 || ref_start > tlen) {
    return out;
  }

  size_t core_begin = 0;
  size_t core_end = ops.size();
  while (core_begin < core_end &&
         (ops[core_begin].second == 'S' || ops[core_begin].second == 'H')) {
    out.lead_clip += ops[core_begin].first;
    ++core_begin;
  }
  while (core_end > core_begin &&
         (ops[core_end - 1].second == 'S' || ops[core_end - 1].second == 'H')) {
    out.trail_clip += ops[core_end - 1].first;
    --core_end;
  }
  if (core_begin == core_end || out.lead_clip > qlen ||
      out.trail_clip > qlen - out.lead_clip) {
    return out;
  }
  if (want_eqx) {
    for (size_t k = 0; k < core_begin; ++k)
      emit_eqx_run(ops[k].first, ops[k].second);
  }

  int qi = out.lead_clip;
  int ti = ref_start;
  out.target_start = ref_start;
  for (size_t k = core_begin; k < core_end; ++k) {
    const int len = ops[k].first;
    const char op = ops[k].second;
    if (len <= 0)
      return CigarReplay{};
    if (op == 'M' || op == '=' || op == 'X') {
      if (len > qlen - qi || len > tlen - ti)
        return CigarReplay{};
      for (int x = 0; x < len; ++x) {
        const uint8_t q = qseq[qi + x];
        const uint8_t t = tseq[ti + x];
        if (q > 3 || t > 3) {
          ++out.ambiguities;
        } else {
          ++out.block_len;
          if (q == t)
            ++out.matches;
          else
            ++out.mismatches;
        }
      }
      if (want_cs || want_md) {
        // Code comparison, as minimap2: N against N extends an identical run.
        for (int x = 0; x < len; ++x) {
          const uint8_t q = qseq[qi + x];
          const uint8_t t = tseq[ti + x];
          if (q == t) {
            if (want_cs) {
              if (want_cs_long)
                cs_run += cs_upper(q);
              ++cs_run_len;
            }
            if (want_md)
              ++md_run;
            continue;
          }
          if (want_cs) {
            flush_cs_run();
            out.cs += '*';
            out.cs += cs_lower(t);
            out.cs += cs_lower(q);
          }
          if (want_md) {
            out.md += std::to_string(md_run);
            out.md += cs_upper(t);
            md_run = 0;
          }
        }
        // minimap2 ends the run at every op, so "5M5M" gives ":5:5".
        if (want_cs)
          flush_cs_run();
      }
      if (want_eqx) {
        int x = 0;
        while (x < len) {
          const bool equal = qseq[qi + x] == tseq[ti + x];
          int run = 1;
          while (x + run < len &&
                 (qseq[qi + x + run] == tseq[ti + x + run]) == equal)
            ++run;
          emit_eqx_run(run, equal ? '=' : 'X');
          x += run;
        }
      }
      qi += len;
      ti += len;
      out.aligned_query_consumed += len;
    } else if (op == 'I') {
      if (len > qlen - qi)
        return CigarReplay{};
      out.insertions += len;
      out.aligned_query_consumed += len;
      for (int x = 0; x < len; ++x) {
        if (qseq[qi + x] > 3)
          ++out.ambiguities;
        else
          ++out.block_len;
      }
      if (want_eqx)
        emit_eqx_run(len, 'I');
      if (want_cs) {
        // MD ignores I and keeps its count.
        out.cs += '+';
        for (int x = 0; x < len; ++x)
          out.cs += cs_lower(qseq[qi + x]);
      }
      qi += len;
    } else if (op == 'D') {
      if (len > tlen - ti)
        return CigarReplay{};
      out.deletions += len;
      for (int x = 0; x < len; ++x) {
        if (tseq[ti + x] > 3)
          ++out.ambiguities;
        else
          ++out.block_len;
      }
      if (want_eqx)
        emit_eqx_run(len, 'D');
      if (want_cs) {
        out.cs += '-';
        for (int x = 0; x < len; ++x)
          out.cs += cs_lower(tseq[ti + x]);
      }
      if (want_md) {
        out.md += std::to_string(md_run);
        out.md += '^';
        for (int x = 0; x < len; ++x)
          out.md += cs_upper(tseq[ti + x]);
        md_run = 0;
      }
      ti += len;
    } else if (op == 'N') {
      if (len > tlen - ti)
        return CigarReplay{};
      if (want_eqx)
        emit_eqx_run(len, 'N');
      if (want_cs) {
        // Forward-strand bases even for a minus-strand transcript. A 1 bp N repeats its
        // base rather than reading past the op.
        const int step = len >= 2 ? 1 : 0;
        out.cs += '~';
        out.cs += cs_lower(tseq[ti]);
        out.cs += cs_lower(tseq[ti + step]);
        out.cs += std::to_string(len);
        out.cs += cs_lower(tseq[ti + len - 1 - step]);
        out.cs += cs_lower(tseq[ti + len - 1]);
      }
      // MD ignores N and keeps its count.
      ti += len;
    } else {
      // Interior S/H, P and unknown operations.
      return CigarReplay{};
    }
  }

  if (out.trail_clip > qlen - qi || qi + out.trail_clip != qlen) {
    return CigarReplay{};
  }
  if (want_eqx) {
    for (size_t k = core_end; k < ops.size(); ++k)
      emit_eqx_run(ops[k].first, ops[k].second);
  }
  qi += out.trail_clip;
  out.query_consumed = qi;
  out.target_consumed = ti - ref_start;
  out.target_end = ti;
  out.edit_distance = out.block_len - out.matches + out.ambiguities;
  // Always written, even "0", as the SAM MD grammar requires; minimap2 omits a zero.
  if (want_md)
    out.md += std::to_string(md_run);
  out.valid = true;
  return out;
}

}}}  // namespace fa::cpu::output

/**** sequence ****/

namespace fa {
namespace cpu {
namespace lr {

std::vector<uint8_t>
reverse_complement_encoded_u8(const std::vector<uint8_t> &fwd_enc) {
  std::vector<uint8_t> rc_enc(fwd_enc.size());
  for (size_t i = 0; i < fwd_enc.size(); i++) {
    const uint8_t b = fwd_enc[fwd_enc.size() - 1 - i];
    rc_enc[i] = (b > 3) ? 4 : static_cast<uint8_t>(3 - b);
  }
  return rc_enc;
}

} // namespace lr
} // namespace cpu
} // namespace fa
