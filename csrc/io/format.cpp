// SAM and PAF output: records, the header, aux tags, cs and MD, block divergence.
#include "io.h"

/**** emitted_role ****/

namespace fa::cpu::output {

static bool emitted_role_is_secondary(EmittedRole role) noexcept {
  return role == EmittedRole::Secondary ||
         role == EmittedRole::SecondarySupplementary;
}

static bool emitted_role_is_supplementary(EmittedRole role) noexcept {
  return role == EmittedRole::Supplementary ||
         role == EmittedRole::SecondarySupplementary;
}

static int emitted_role_sam_flag(EmittedRole role) noexcept {
  return (emitted_role_is_secondary(role) ? 0x100 : 0) |
         (emitted_role_is_supplementary(role) ? 0x800 : 0);
}

static char emitted_role_paf_type(EmittedRole role) noexcept {
  return emitted_role_is_secondary(role) ? 'S' : 'P';
}

// The tp:A type as minimap2's write_tags sets it:
//   type = r->id == r->parent ? (r->inv? 'I' : 'P') : (r->inv? 'i' : 'S');
// Only DNA local-inversion records (the mm_align1_inv port) get 'I' or 'i'.
// An all-chains record is 'S' in any role: minimap2's -P sets no parent.
static char emitted_role_paf_type(EmittedRole role,
                                  const AlignResult& result) noexcept {
  if (result.origin == AlignmentOrigin::DnaLocalInversion)
    return emitted_role_is_secondary(role) ? 'i' : 'I';
  if (result.origin == AlignmentOrigin::DnaAllChains)
    return 'S';
  return emitted_role_paf_type(role);
}

void validate_one_level_hypotheses(const AlignResult& result) {
  for (const AlignResult& secondary : result.secondary) {
    if (!secondary.secondary.empty()) {
      throw std::invalid_argument(
          "alignment secondary hypotheses may not be nested");
    }
    for (const AlignResult& supplementary : secondary.supplementary) {
      if (!supplementary.secondary.empty() ||
          !supplementary.supplementary.empty()) {
        throw std::invalid_argument(
            "alignment secondary supplementary segments must be leaves");
      }
    }
  }
  for (const AlignResult& supplementary : result.supplementary) {
    if (!supplementary.secondary.empty() ||
        !supplementary.supplementary.empty()) {
      throw std::invalid_argument(
          "alignment supplementary segments must be leaves");
    }
  }
}

}  // namespace fa::cpu::output

/**** alignment_aux_tags ****/

namespace fa::cpu::output {

static AlignmentAuxiliaryTags alignment_auxiliary_tags(
    const AlignResult& result) {
  AlignmentAuxiliaryTags tags;
  if (!result.mapped()) return tags;

  tags.alignment_score = result.score;
  if (result.alignment_accounting_valid) {
    tags.edit_distance = result.edit_distance;
    tags.ambiguities = result.ambiguities;
  }
  // md is minibwa's margin: this record's ms minus the best ms among its
  // mapped alternatives (0 when there are none). It may be negative.
  if (result.dp_max_segment >= 0) {
    tags.max_segment_score = result.dp_max_segment;
    int best_alternative = 0;
    for (const AlignResult& alternative : result.secondary) {
      if (!alternative.mapped() || alternative.dp_max_segment < 0)
        continue;
      best_alternative = std::max(best_alternative, alternative.dp_max_segment);
    }
    tags.max_score_margin = result.dp_max_segment - best_alternative;
  }
  if (result.chain_anchors >= 0)
    tags.chain_anchors = result.chain_anchors;
  if (result.chain_score >= 0)
    tags.chain_score = result.chain_score;
  if (result.secondary_chain_score >= 0)
    tags.secondary_chain_score = result.secondary_chain_score;
  return tags;
}

static int factual_sa_edit_distance(const AlignResult& result) {
  const AlignmentAuxiliaryTags tags = alignment_auxiliary_tags(result);
  if (!tags.edit_distance) {
    throw std::logic_error(
        "mapped split-family member lacks exact alignment accounting");
  }
  return *tags.edit_distance;
}

}  // namespace fa::cpu::output

/**** block_divergence ****/

namespace fa {
namespace cpu {
namespace output {

// Uses the engine's exact accounting when alignment_accounting_valid is set.
// Otherwise the query interval comes from the CIGAR's clips (mirrored for a
// reverse alignment, whose CIGAR is reference-forward), and NM is
// block_len - matches: mismatches plus inserted and deleted bases.
static BlockAccounting block_accounting(const AlignResult& result, int qlen) {
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
static bool event_divergence(const BlockAccounting& acc, double& divergence) {
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

static std::string format_paf_de(double de) {
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

/**** sam ****/

namespace fa { namespace cpu { namespace output {

static int clamp_sam_mapq(int mapq) {
    return std::clamp(mapq, 0, 255);
}

static std::string clean_sam_field(const std::string& field) {
    if (field.empty()) return "*";
    std::string out = field;
    for (char& ch : out) {
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = '_';
    }
    return out;
}

static std::string cigar_or_fallback(const AlignResult& result) {
    if (!result.mapped()) return "*";
    return result.cigar;
}

static int reference_span_from_cigar(const std::string& cigar, int fallback_len) {
    if (cigar.empty() || cigar == "*") return std::max(0, fallback_len);

    int64_t span = 0;
    int64_t count = 0;
    bool saw_op = false;
    bool malformed = false;

    for (char ch : cigar) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            count = count * 10 + static_cast<int64_t>(ch - '0');
            if (count > INT_MAX) malformed = true;
            continue;
        }

        if (count <= 0) malformed = true;
        switch (ch) {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                span += count;
                break;
            case 'I':
            case 'S':
            case 'H':
            case 'P':
                break;
            default:
                malformed = true;
                break;
        }
        count = 0;
        saw_op = true;
    }

    if (malformed || count != 0 || !saw_op) return std::max(0, fallback_len);
    if (span > INT_MAX) return INT_MAX;
    return static_cast<int>(span);
}

int reference_span(const AlignResult& result) {
    return reference_span_from_cigar(cigar_or_fallback(result), result.read_len);
}

static int single_segment_flag(const AlignResult& result) {
    int flag = 0;
    if (!result.mapped()) flag |= 0x4;
    if (result.mapped() && result.is_reverse) flag |= 0x10;
    return flag;
}

static void append_sam_int(std::string& out, int64_t value) {
    char buf[32];
    auto result = std::to_chars(buf, buf + sizeof(buf), value);
    if (result.ec == std::errc()) {
        out.append(buf, result.ptr);
    } else {
        out += std::to_string(value);
    }
}

static void append_clean_sam_field(std::string& out, const std::string& field) {
    if (field.empty()) {
        out.push_back('*');
        return;
    }
    for (char ch : field) {
        out.push_back((ch == '\t' || ch == '\n' || ch == '\r') ? '_' : ch);
    }
}

static void append_sam_cigar(std::string& out, const AlignResult& result) {
    if (!result.mapped()) {
        out.push_back('*');
        return;
    }
    append_clean_sam_field(out, result.cigar);
}

// Complement of each byte: A, C, G, T in either case; any other byte gives 'N'.
const char* complement_table() {
    static const std::array<char, 256> table = [] {
        std::array<char, 256> t{};
        t.fill('N');
        t[static_cast<unsigned char>('A')] = 'T';
        t[static_cast<unsigned char>('T')] = 'A';
        t[static_cast<unsigned char>('C')] = 'G';
        t[static_cast<unsigned char>('G')] = 'C';
        t[static_cast<unsigned char>('a')] = 't';
        t[static_cast<unsigned char>('t')] = 'a';
        t[static_cast<unsigned char>('c')] = 'g';
        t[static_cast<unsigned char>('g')] = 'c';
        return t;
    }();
    return table.data();
}

void append_reverse_complement(std::string& out, const std::string& seq) {
    const std::size_t base = out.size();
    out.resize(base + seq.size());
    const char* table = complement_table();
    char* dst = out.data() + base;
    const std::size_t n = seq.size();
    for (std::size_t i = 0; i < n; ++i) {
        dst[i] = table[static_cast<unsigned char>(seq[n - 1 - i])];
    }
}

// Appends `text` reversed, into pre-sized storage.
void append_reversed(std::string& out, const std::string& text) {
    const std::size_t base = out.size();
    out.resize(base + text.size());
    char* dst = out.data() + base;
    const std::size_t n = text.size();
    for (std::size_t i = 0; i < n; ++i) dst[i] = text[n - 1 - i];
}

void append_sam_record_line(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::string& extra_tags,
    int flag_extra,
    const std::string* cigar_override,
    bool seq_preoriented,
    const std::string& read_group_id,
    EmittedRole role,
    const std::string& record_tags
) {
    const bool mapped = result.mapped();
    append_clean_sam_field(out, read_name);
    out.push_back('\t');
    append_sam_int(out, single_segment_flag(result) | (mapped ? flag_extra : 0));
    out.push_back('\t');
    if (mapped) append_clean_sam_field(out, result.chromosome);
    else out.push_back('*');
    out.push_back('\t');
    append_sam_int(out, mapped ? result.pos + 1 : 0);
    out.push_back('\t');
    append_sam_int(out, mapped ? clamp_sam_mapq(result.mapq) : 0);
    out.push_back('\t');
    if (cigar_override) append_clean_sam_field(out, *cigar_override);
    else append_sam_cigar(out, result);
    out += "\t*\t0\t0\t";
    // SEQ/QUAL in reference orientation, unless the caller already did it.
    const bool rc = mapped && result.is_reverse && !seq_preoriented;
    if (seq.empty()) out.push_back('*');
    else if (rc) append_reverse_complement(out, seq);
    else out += seq;
    out.push_back('\t');
    if (qual.empty()) out.push_back('*');
    else if (rc) append_reversed(out, qual);
    else out += qual;
    // RG:Z comes first, as in minimap2, on mapped and unmapped records alike.
    if (!read_group_id.empty()) {
        out += "\tRG:Z:";
        out += read_group_id;
    }
    const AlignmentAuxiliaryTags auxiliary =
        alignment_auxiliary_tags(result);
    if (auxiliary.alignment_score) {
        out += "\tAS:i:";
        append_sam_int(out, *auxiliary.alignment_score);
    }
    // ms:i and md:i follow AS:i, as in minimap2 and minibwa. md compares a
    // record with its alternatives, so secondaries do not carry it.
    if (auxiliary.max_segment_score) {
        out += "\tms:i:";
        append_sam_int(out, *auxiliary.max_segment_score);
    }
    if (auxiliary.max_score_margin && !emitted_role_is_secondary(role)) {
        out += "\tmd:i:";
        append_sam_int(out, *auxiliary.max_score_margin);
    }
    if (auxiliary.edit_distance) {
        out += "\tNM:i:";
        append_sam_int(out, *auxiliary.edit_distance);
    }
    if (auxiliary.ambiguities) {
        out += "\tnn:i:";
        append_sam_int(out, *auxiliary.ambiguities);
    }
    // tp:A and de:f, as minimap2 writes them in SAM too, computed as for PAF.
    if (mapped) {
        out += "\ttp:A:";
        out.push_back(emitted_role_paf_type(role, result));
        // cm:i, s1:i and s2:i between tp:A and de:f, as in minimap2; s2 not
        // on secondaries.
        if (auxiliary.chain_anchors) {
            out += "\tcm:i:";
            append_sam_int(out, *auxiliary.chain_anchors);
        }
        if (auxiliary.chain_score) {
            out += "\ts1:i:";
            append_sam_int(out, *auxiliary.chain_score);
        }
        if (auxiliary.secondary_chain_score &&
            !emitted_role_is_secondary(role)) {
            out += "\ts2:i:";
            append_sam_int(out, *auxiliary.secondary_chain_score);
        }
        const BlockAccounting acc =
            block_accounting(result, result.read_len);
        double divergence = 0.0;
        if (event_divergence(acc, divergence)) {
            out += "\tde:f:";
            out += format_paf_de(divergence);
        }
    }
    // RNA transcript strand (minimap2 ts:A); '\0' on DNA records.
    if (mapped && (result.transcript_strand == '+' ||
                   result.transcript_strand == '-')) {
        out += "\tts:A:";
        out.push_back(result.transcript_strand);
    }
    if (!extra_tags.empty()) {
        out.push_back('\t');
        out += extra_tags;
    }
    // cs:Z and MD:Z after extra_tags (which hold SA:Z), in minimap2's order.
    if (mapped && result.alignment_accounting_valid) {
      if (!result.cs.empty()) {
        out += "\tcs:Z:";
        out += result.cs;
      }
      if (!result.md.empty()) {
        out += "\tMD:Z:";
        out += result.md;
      }
    }
    // Last, where minimap2 puts the -y comment.
    if (!record_tags.empty()) {
        out.push_back('\t');
        out += record_tags;
    }
    out.push_back('\n');
}

// SA:Z for segment `self`, one "rname,pos,strand,CIGAR,mapQ,NM;" entry per
// other mapped segment. Throws if a segment lacks exact accounting for NM.
// Entries carry the full CIGAR, not minimap2's condensed form.
static std::string build_sa_tag(const std::vector<const AlignResult*>& segs, size_t self) {
    std::string sa = "SA:Z:";
    bool any = false;
    for (size_t i = 0; i < segs.size(); ++i) {
        if (i == self) continue;
        const AlignResult& s = *segs[i];
        if (!s.mapped()) continue;
        const int nm = factual_sa_edit_distance(s);
        sa += clean_sam_field(s.chromosome);
        sa += ',';
        sa += std::to_string(s.pos + 1);
        sa += s.is_reverse ? ",-," : ",+,";
        sa += clean_sam_field(cigar_or_fallback(s));
        sa += ',';
        sa += std::to_string(clamp_sam_mapq(s.mapq));
        sa += ',';
        sa += std::to_string(nm);
        sa += ';';
        any = true;
    }
    return any ? sa : std::string();
}

// XA:Z from the read's alternatives, as bwa and minibwa write it:
// "chr,{+|-}pos,CIGAR,NM;" per mapped secondary head, POS 1-based. Entries
// without a CIGAR or exact accounting are skipped.
// Empty when there are none.
static std::string build_xa_tag(const AlignResult& result) {
    std::string xa = "XA:Z:";
    bool any = false;
    for (const AlignResult& alternative : result.secondary) {
        if (!alternative.mapped()) continue;
        const std::string& cigar = alternative.cigar;
        if (cigar.empty() || cigar == "*") continue;
        // An alternative without exact accounting has no NM and is skipped.
        const AlignmentAuxiliaryTags auxiliary =
            alignment_auxiliary_tags(alternative);
        if (!auxiliary.edit_distance) continue;
        xa += clean_sam_field(alternative.chromosome);
        xa += ',';
        xa += alternative.is_reverse ? '-' : '+';
        xa += std::to_string(alternative.pos + 1);
        xa += ',';
        xa += cigar;
        xa += ',';
        xa += std::to_string(*auxiliary.edit_distance);
        xa += ';';
        any = true;
    }
    return any ? xa : std::string();
}

// Turns terminal soft clips into hard clips and drops the clipped bases from
// SEQ/QUAL. The returned seq/qual are reference-oriented. An absent or "*"
// read field gives an empty slice.
static HardClippedRecord hard_clip_supplementary(
    const std::string& soft_cigar, const std::string& fwd_seq,
    const std::string& fwd_qual, bool is_reverse) {
    HardClippedRecord r;
    const auto hc = hard_clip_terminal(soft_cigar);
    r.cigar = hc.cigar;
    std::string oseq, oqual;
    if (!fwd_seq.empty() && fwd_seq != "*") {
        if (is_reverse) append_reverse_complement(oseq, fwd_seq);
        else oseq = fwd_seq;
    }
    if (!fwd_qual.empty() && fwd_qual != "*") {
        if (is_reverse) append_reversed(oqual, fwd_qual);
        else oqual = fwd_qual;
    }
    const auto slice = [](const std::string& s, int lead, int trail) -> std::string {
        if (s.empty()) return s;
        const int n = static_cast<int>(s.size());
        const int a = std::min(std::max(lead, 0), n);
        const int b = std::max(a, n - std::max(trail, 0));
        return s.substr(static_cast<size_t>(a), static_cast<size_t>(b - a));
    };
    r.seq = slice(oseq, hc.lead, hc.trail);
    r.qual = slice(oqual, hc.lead, hc.trail);
    return r;
}

void append_sam_records(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const PassthroughTags& passthrough,
    SamEmitOptions opts,
    const std::string& record_tags
) {
    const auto emit_hypothesis = [&](const AlignResult& hypothesis,
                                     bool secondary_hypothesis) {
        std::vector<const AlignResult*> segs;
        segs.reserve(1 + hypothesis.supplementary.size());
        segs.push_back(&hypothesis);
        for (const auto& s : hypothesis.supplementary) segs.push_back(&s);
        for (size_t i = 0; i < segs.size(); ++i) {
            const EmittedRole role = secondary_hypothesis
                ? (i == 0 ? EmittedRole::Secondary
                          : EmittedRole::SecondarySupplementary)
                : (i == 0 ? EmittedRole::Primary
                          : EmittedRole::Supplementary);
            AlignResult emitted = *segs[i];
            if (secondary_hypothesis) emitted.mapq = 0;
            const bool hard_clip = !secondary_hypothesis && i > 0 &&
                                   opts.hard_clip_supp && emitted.mapped() &&
                                   !emitted.cigar.empty();
            // Per-base tags only where SEQ is the whole read.
            std::string tags = (secondary_hypothesis || hard_clip)
                                   ? passthrough.clipped
                                   : passthrough.whole;
            // SA:Z links the primary's segments only.
            if (!secondary_hypothesis) {
                const std::string sa = build_sa_tag(segs, i);
                if (!sa.empty()) {
                    if (!tags.empty()) tags += '\t';
                    tags += sa;
                }
            }
            // XA:Z after SA:Z, on the primary only, whether or not the
            // alternatives are also written as records.
            if (!secondary_hypothesis && i == 0 && segs[i]->mapped()) {
                const std::string xa = build_xa_tag(*segs[i]);
                if (!xa.empty()) {
                    if (!tags.empty()) tags += '\t';
                    tags += xa;
                }
            }
            const int flag_extra = emitted_role_sam_flag(role);
            if (hard_clip) {
                const HardClippedRecord hc = hard_clip_supplementary(
                    emitted.cigar, seq, qual, emitted.is_reverse);
                append_sam_record_line(out, read_name, emitted, hc.seq, hc.qual,
                                       tags, flag_extra, &hc.cigar,
                                       /*seq_preoriented=*/true,
                                       opts.read_group_id, role, record_tags);
            } else {
                // As in minimap2, secondaries keep the soft-clipped CIGAR but
                // omit SEQ/QUAL.
                append_sam_record_line(
                    out, read_name, emitted,
                    secondary_hypothesis ? std::string() : seq,
                    secondary_hypothesis ? std::string() : qual,
                    tags, flag_extra, /*cigar_override=*/nullptr,
                    /*seq_preoriented=*/false, opts.read_group_id, role,
                    record_tags);
            }
        }
    };
    emit_hypothesis(result, /*secondary_hypothesis=*/false);
    // An unmapped read has no alternatives to print.
    if (opts.emit_secondary && result.mapped()) {
        for (const AlignResult& secondary : result.secondary) {
            if (secondary.mapped())
                emit_hypothesis(secondary, /*secondary_hypothesis=*/true);
        }
    }
}

std::vector<std::string> sam_header(
    const std::vector<SamReference>& references,
    const std::string& program_name,
    const std::string& version,
    const std::string& command_line,
    const std::string& program_id,
    const std::string& previous_program
) {
    std::vector<std::string> lines;
    lines.push_back("@HD\tVN:1.6\tSO:unsorted");
    for (const SamReference& ref : references) {
        lines.push_back("@SQ\tSN:" + clean_sam_field(ref.name) + "\tLN:" + std::to_string(ref.length));
    }
    std::string program = "@PG\tID:" + clean_sam_field(program_id.empty() ? program_name : program_id)
                          + "\tPN:" + clean_sam_field(program_name);
    if (!previous_program.empty())
        program += "\tPP:" + clean_sam_field(previous_program);
    program += "\tVN:" + clean_sam_field(version);
    if (!command_line.empty())
        program += "\tCL:" + clean_sam_field(command_line);
    lines.push_back(std::move(program));
    return lines;
}

}}}  // namespace fa::cpu::output

/**** paf ****/

namespace fa { namespace cpu { namespace output {

static int target_end_from_result(const AlignResult& result) {
    if (!result.mapped()) return 0;
    if (result.alignment_accounting_valid)
        return std::max(result.pos, result.target_end);
    if (result.target_end > result.pos) return result.target_end;
    if (!result.target_regions.empty()) {
        int end = result.pos;
        for (const auto& region : result.target_regions) {
            end = std::max(end, region.second);
        }
        return end;
    }
    return result.pos + reference_span(result);
}

static std::string paf_cigar_from_sam_cigar(const std::string& cigar) {
    const auto ops = parse_cigar_ops(cigar);
    if (ops.empty()) return cigar;
    std::size_t begin = 0;
    std::size_t end = ops.size();
    if (ops[begin].second == 'S') ++begin;
    if (end > begin && ops[end - 1].second == 'S') --end;
    std::string out;
    for (std::size_t i = begin; i < end; ++i) {
        append_cigar_run(out, ops[i].first, ops[i].second);
    }
    return out;
}

void write_paf_record(
    std::ostream& out,
    const ::fa::cpu::io::FastxRecord& read,
    const AlignResult& result,
    const std::unordered_map<std::string, int64_t>& ref_lengths,
    bool with_cigar,
    EmittedRole role,
    bool copy_comment
) {
    if (!result.mapped()) return;
    const int qlen = static_cast<int>(read.seq.size());
    // Shared with SAM (io/io.h), so NM:i and de:f agree.
    const BlockAccounting acc = block_accounting(result, qlen);
    const char strand = result.is_reverse ? '-' : '+';
    const auto found = ref_lengths.find(result.chromosome);
    const int64_t tlen = found == ref_lengths.end() ? 0 : found->second;
    const int tstart = std::max(0, result.pos);
    const int tend = std::max(tstart, target_end_from_result(result));
    out << read.name << '\t'
        << qlen << '\t'
        << acc.qstart << '\t'
        << acc.qend << '\t'
        << strand << '\t'
        << result.chromosome << '\t'
        << tlen << '\t'
        << tstart << '\t'
        << tend << '\t'
        << acc.matches << '\t'
        << acc.block_len << '\t'
        << (emitted_role_is_secondary(role) ? 0 : clamp_sam_mapq(result.mapq));
    // An all-chains record carries minimap2's overlap tags, which have no AS.
    if (result.origin != AlignmentOrigin::DnaAllChains)
        out << "\tAS:i:" << result.score;
    // md and s2 compare a record with its alternatives, so secondaries do not
    // carry them, as in minimap2 and minibwa.
    const AlignmentAuxiliaryTags auxiliary = alignment_auxiliary_tags(result);
    if (auxiliary.max_segment_score)
        out << "\tms:i:" << *auxiliary.max_segment_score;
    if (auxiliary.max_score_margin && !emitted_role_is_secondary(role))
        out << "\tmd:i:" << *auxiliary.max_score_margin;
    if (acc.nm >= 0) out << "\tNM:i:" << acc.nm;
    if (result.alignment_accounting_valid)
        out << "\tnn:i:" << result.ambiguities;
    out << "\ttp:A:" << emitted_role_paf_type(role, result);
    // cm:i and s1:i, the record's own chain, on any record that holds one.
    if (auxiliary.chain_anchors)
        out << "\tcm:i:" << *auxiliary.chain_anchors;
    if (auxiliary.chain_score)
        out << "\ts1:i:" << *auxiliary.chain_score;
    if (auxiliary.secondary_chain_score && !emitted_role_is_secondary(role))
        out << "\ts2:i:" << *auxiliary.secondary_chain_score;
    double divergence = 0.0;
    if (with_cigar && event_divergence(acc, divergence))
        out << "\tde:f:" << format_paf_de(divergence);
    if (with_cigar && !result.cigar.empty()) {
        const std::string paf_cigar = paf_cigar_from_sam_cigar(result.cigar);
        if (!paf_cigar.empty()) out << "\tcg:Z:" << paf_cigar;
    }
    // cs:Z and MD:Z follow cg:Z, in minimap2's order.
    if (result.alignment_accounting_valid) {
      if (!result.cs.empty())
        out << "\tcs:Z:" << result.cs;
      if (!result.md.empty())
        out << "\tMD:Z:" << result.md;
    }
    // RNA transcript strand (minimap2 ts:A); '\0' on DNA records.
    if (result.transcript_strand == '+' || result.transcript_strand == '-') {
        out << "\tts:A:" << result.transcript_strand;
    }
    // -y: the FASTA/Q comment, verbatim and last, as minimap2 writes it.
    if (copy_comment && !read.comment.empty()) out << '\t' << read.comment;
    out << '\n';
}

void write_paf_no_hit_record(
    std::ostream& out, const ::fa::cpu::io::FastxRecord& read) {
    out << read.name << '\t' << read.seq.size()
        << "\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n";
}

}}}  // namespace fa::cpu::output
