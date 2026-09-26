#include "sam.h"

namespace fa { namespace cpu { namespace output {

int clamp_sam_mapq(int mapq) {
    return std::clamp(mapq, 0, 255);
}

std::string clean_sam_field(const std::string& field) {
    if (field.empty()) return "*";
    std::string out = field;
    for (char& ch : out) {
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = '_';
    }
    return out;
}

std::string cigar_or_fallback(const AlignResult& result) {
    if (!result.mapped()) return "*";
    return result.cigar;
}

int reference_span_from_cigar(const std::string& cigar, int fallback_len) {
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

void add_int_tag(std::vector<std::string>& tags, const std::string& tag, int value) {
    tags.push_back(tag + ":i:" + std::to_string(value));
}

void add_string_tag(std::vector<std::string>& tags, const std::string& tag, const std::string& value) {
    tags.push_back(tag + ":Z:" + clean_sam_field(value));
}

int single_segment_flag(const AlignResult& result) {
    int flag = 0;
    if (!result.mapped()) flag |= 0x4;
    if (result.mapped() && result.is_reverse) flag |= 0x10;
    return flag;
}

void append_sam_int(std::string& out, int64_t value) {
    char buf[32];
    auto result = std::to_chars(buf, buf + sizeof(buf), value);
    if (result.ec == std::errc()) {
        out.append(buf, result.ptr);
    } else {
        out += std::to_string(value);
    }
}

void append_clean_sam_field(std::string& out, const std::string& field) {
    if (field.empty()) {
        out.push_back('*');
        return;
    }
    for (char ch : field) {
        out.push_back((ch == '\t' || ch == '\n' || ch == '\r') ? '_' : ch);
    }
}

void append_sam_cigar(std::string& out, const AlignResult& result) {
    if (!result.mapped()) {
        out.push_back('*');
        return;
    }
    append_clean_sam_field(out, result.cigar);
}

char complement_base(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default:  return 'N';
    }
}

// Table form of complement_base(), for the single-threaded writer's hottest
// loop.
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
        // s2:i between tp:A and de:f, as in minimap2; not on secondaries.
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

std::string build_sa_tag(const std::vector<const AlignResult*>& segs, size_t self) {
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

std::string build_xa_tag(const AlignResult& result) {
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

HardClippedRecord hard_clip_supplementary(
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
    const std::string& extra_tags,
    SamEmitOptions opts,
    const std::string& record_tags
) {
    const auto emit_hypothesis = [&](const AlignResult& hypothesis,
                                     bool secondary_hypothesis,
                                     const std::string& root_tags) {
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
            std::string tags = (i == 0 && !secondary_hypothesis)
                                   ? root_tags
                                   : std::string();
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
            AlignResult emitted = *segs[i];
            if (secondary_hypothesis) emitted.mapq = 0;
            const int flag_extra = emitted_role_sam_flag(role);
            if (!secondary_hypothesis && i > 0 && opts.hard_clip_supp &&
                emitted.mapped() && !emitted.cigar.empty()) {
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
    emit_hypothesis(result, /*secondary_hypothesis=*/false, extra_tags);
    // An unmapped read has no alternatives to print.
    if (opts.emit_secondary && result.mapped()) {
        for (const AlignResult& secondary : result.secondary) {
            if (secondary.mapped())
                emit_hypothesis(secondary, /*secondary_hypothesis=*/true, {});
        }
    }
}

std::vector<std::string> sam_header(
    const std::vector<SamReference>& references,
    const std::string& program_name,
    const std::string& version,
    const std::string& command_line
) {
    std::vector<std::string> lines;
    lines.push_back("@HD\tVN:1.6\tSO:unsorted");
    for (const SamReference& ref : references) {
        lines.push_back("@SQ\tSN:" + clean_sam_field(ref.name) + "\tLN:" + std::to_string(ref.length));
    }
    std::string program = "@PG\tID:" + clean_sam_field(program_name)
                          + "\tPN:" + clean_sam_field(program_name)
                          + "\tVN:" + clean_sam_field(version);
    if (!command_line.empty())
        program += "\tCL:" + clean_sam_field(command_line);
    lines.push_back(std::move(program));
    return lines;
}

}}}  // namespace fa::cpu::output
