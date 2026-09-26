#include "bam.h"

#include <stdexcept>

namespace fa { namespace cpu { namespace output {

void bam_put_i32(std::string& b, int32_t v) {
    const uint32_t u = static_cast<uint32_t>(v);
    for (int i = 0; i < 4; ++i) b.push_back(static_cast<char>((u >> (8 * i)) & 0xff));
}

void bam_put_u32(std::string& b, uint32_t u) {
    for (int i = 0; i < 4; ++i) b.push_back(static_cast<char>((u >> (8 * i)) & 0xff));
}

void bam_put_u16(std::string& b, uint16_t u) {
    b.push_back(static_cast<char>(u & 0xff));
    b.push_back(static_cast<char>((u >> 8) & 0xff));
}

void bam_tag_i(std::string& tags, char a, char b, int32_t v) {
    tags.push_back(a); tags.push_back(b); tags.push_back('i'); bam_put_i32(tags, v);
}

void bam_tag_A(std::string& tags, char a, char b, char v) {
    tags.push_back(a); tags.push_back(b); tags.push_back('A'); tags.push_back(v);
}

void bam_tag_f(std::string& tags, char a, char b, float v) {
    tags.push_back(a); tags.push_back(b); tags.push_back('f');
    char raw[sizeof(float)];
    std::memcpy(raw, &v, sizeof(float));
    tags.append(raw, sizeof(float));
}

void bam_tag_Z(std::string& tags, char a, char b, const std::string& v) {
    tags.push_back(a); tags.push_back(b); tags.push_back('Z');
    tags += v; tags.push_back('\0');
}

const std::array<uint8_t, 256>& seq_nt16_table() {
    static const std::array<uint8_t, 256> table = [] {
        std::array<uint8_t, 256> a{};
        a.fill(15);  // unknown -> N
        const char* s = "=ACMGRSVTWYHKDBN";
        for (int i = 0; i < 16; ++i) {
            a[static_cast<unsigned char>(s[i])] = static_cast<uint8_t>(i);
            a[static_cast<unsigned char>(std::tolower(static_cast<unsigned char>(s[i])))] =
                static_cast<uint8_t>(i);
        }
        return a;
    }();
    return table;
}

int bam_cigar_op_code(char op) {
    switch (op) {
        case 'M': return 0; case 'I': return 1; case 'D': return 2; case 'N': return 3;
        case 'S': return 4; case 'H': return 5; case 'P': return 6; case '=': return 7;
        case 'X': return 8; default: return -1;
    }
}

int bam_reg2bin(int64_t beg, int64_t end) {
    // As the SAMv1 section 5.3 listing, with 64-bit operands.
    --end;
    if ((beg >> 14) == (end >> 14))
        return static_cast<int>(((1 << 15) - 1) / 7 + (beg >> 14));
    if ((beg >> 17) == (end >> 17))
        return static_cast<int>(((1 << 12) - 1) / 7 + (beg >> 17));
    if ((beg >> 20) == (end >> 20))
        return static_cast<int>(((1 << 9) - 1) / 7 + (beg >> 20));
    if ((beg >> 23) == (end >> 23))
        return static_cast<int>(((1 << 6) - 1) / 7 + (beg >> 23));
    if ((beg >> 26) == (end >> 26))
        return static_cast<int>(((1 << 3) - 1) / 7 + (beg >> 26));
    return 0;
}

std::string encode_bam_header(
    const std::vector<SamReference>& refs,
    const std::string& sam_header_text) {
    std::string h;
    h += "BAM\1";
    bam_put_i32(h, static_cast<int32_t>(sam_header_text.size()));
    h += sam_header_text;
    bam_put_i32(h, static_cast<int32_t>(refs.size()));
    for (const SamReference& r : refs) {
        // SAMv1 limits l_ref to < 2^31. The index already refuses longer
        // contigs; refuse here too rather than wrap to a negative length.
        if (r.length < 0 || r.length > 2147483647LL) {
            throw std::runtime_error(
                "BAM header: reference " + r.name + " is " +
                std::to_string(r.length) +
                " bp, outside the [0, 2147483647] range the BAM l_ref field "
                "and the signed BAM POS can represent");
        }
        bam_put_i32(h, static_cast<int32_t>(r.name.size() + 1));
        h += r.name; h.push_back('\0');
        bam_put_i32(h, static_cast<int32_t>(r.length));
    }
    return h;
}

void encode_bam_record(
    std::string& out,
    const std::string& qname,
    int32_t refID,
    int32_t pos0,          // 0-based; -1 if unmapped
    int mapq,
    int flag,
    const std::string& cigar_text,   // AlignResult.cigar ("" or "*" = none)
    const std::string& seq_in,
    const std::string& qual_in,
    bool is_reverse,
    const std::string& extra_tags) {
    // SEQ and QUAL in reference orientation, as in SAM.
    std::string oseq, oqual;
    const bool have_seq = !(seq_in.empty() || seq_in == "*");
    if (have_seq) {
        if (is_reverse) append_reverse_complement(oseq, seq_in);
        else oseq = seq_in;
    }
    const bool have_qual = !(qual_in.empty() || qual_in == "*");
    if (have_qual) {
        if (is_reverse) oqual.assign(qual_in.rbegin(), qual_in.rend());
        else oqual = qual_in;
    }
    const int32_t l_seq = have_seq ? static_cast<int32_t>(oseq.size()) : 0;

    // CIGAR ops + reference span (for bin).
    std::vector<std::pair<int, char>> ops;
    if (!cigar_text.empty() && cigar_text != "*") ops = parse_cigar_ops(cigar_text);
    int ref_span = 0;
    for (const auto& [len, op] : ops) {
        if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') ref_span += len;
    }

    // CG:B,I escape when the op count overflows the 16-bit n_cigar_op field.
    std::string cg_tag;
    std::vector<std::pair<int, char>> emit_ops = ops;
    if (ops.size() > 65535) {
        cg_tag.push_back('C'); cg_tag.push_back('G'); cg_tag.push_back('B'); cg_tag.push_back('I');
        bam_put_i32(cg_tag, static_cast<int32_t>(ops.size()));
        for (const auto& [len, op] : ops) {
            bam_put_u32(cg_tag, (static_cast<uint32_t>(len) << 4) |
                                 static_cast<uint32_t>(bam_cigar_op_code(op)));
        }
        // Placeholder CIGAR: <l_seq>S<ref_span>N (htslib convention).
        emit_ops.clear();
        emit_ops.push_back({l_seq, 'S'});
        emit_ops.push_back({ref_span, 'N'});
    }

    // SAMv1 section 4.2.1: BIN is reg2bin() of the start and the CIGAR end; a
    // record with no reference span counts as length one, and an unmapped one
    // gets 4680. Beyond 512 Mbp the bin does not fit 16 bits and is truncated,
    // as htslib does; readers recompute it, and such contigs need a CSI index.
    const int64_t begin = pos0;
    const int64_t end = begin + (ref_span > 0 ? ref_span : 1);
    const int bin =
        (refID < 0 || pos0 < 0) ? 4680 : bam_reg2bin(begin, end);

    std::string rec;  // everything after block_size
    bam_put_i32(rec, refID);
    bam_put_i32(rec, pos0);
    rec.push_back(static_cast<char>(static_cast<uint8_t>(qname.size() + 1)));  // l_read_name
    rec.push_back(static_cast<char>(static_cast<uint8_t>(mapq < 0 ? 0 : (mapq > 255 ? 255 : mapq))));
    bam_put_u16(rec, static_cast<uint16_t>(bin));
    bam_put_u16(rec, static_cast<uint16_t>(emit_ops.size()));
    bam_put_u16(rec, static_cast<uint16_t>(flag));
    bam_put_i32(rec, l_seq);
    bam_put_i32(rec, -1);  // next_refID
    bam_put_i32(rec, -1);  // next_pos
    bam_put_i32(rec, 0);   // tlen

    rec += qname; rec.push_back('\0');

    for (const auto& [len, op] : emit_ops) {
        bam_put_u32(rec, (static_cast<uint32_t>(len) << 4) |
                          static_cast<uint32_t>(bam_cigar_op_code(op)));
    }

    // 4-bit packed SEQ (high nibble first).
    const auto& nt16 = seq_nt16_table();
    for (int32_t i = 0; i < l_seq; i += 2) {
        uint8_t hi = nt16[static_cast<unsigned char>(oseq[i])];
        uint8_t lo = (i + 1 < l_seq) ? nt16[static_cast<unsigned char>(oseq[i + 1])] : 0;
        rec.push_back(static_cast<char>((hi << 4) | lo));
    }

    // QUAL: raw Phred (= ASCII - 33). 0xff*l_seq when absent.
    if (have_qual) {
        for (char c : oqual) rec.push_back(static_cast<char>(static_cast<uint8_t>(c) - 33));
    } else {
        rec.append(static_cast<size_t>(l_seq), static_cast<char>(0xff));
    }

    rec += extra_tags;
    rec += cg_tag;

    bam_put_i32(out, static_cast<int32_t>(rec.size()));
    out += rec;
}

std::vector<std::string> encode_bam_records(
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::unordered_map<std::string, int>& reference_ids,
    const std::string& primary_extra_tags,
    SamEmitOptions opts,
    const std::string& record_extra_tags) {
    std::vector<std::string> records;
    std::size_t record_count = 1 + result.supplementary.size();
    for (const AlignResult& secondary : result.secondary)
        record_count += 1 + secondary.supplementary.size();
    records.reserve(record_count);
    const auto emit_hypothesis = [&](const AlignResult& hypothesis,
                                     bool secondary_hypothesis) {
      std::vector<const AlignResult*> segments;
      segments.reserve(1 + hypothesis.supplementary.size());
      segments.push_back(&hypothesis);
      for (const AlignResult& supplementary : hypothesis.supplementary)
          segments.push_back(&supplementary);
      for (std::size_t index = 0; index < segments.size(); ++index) {
        AlignResult segment = *segments[index];
        if (secondary_hypothesis) segment.mapq = 0;
        const EmittedRole role = secondary_hypothesis
            ? (index == 0 ? EmittedRole::Secondary
                          : EmittedRole::SecondarySupplementary)
            : (index == 0 ? EmittedRole::Primary
                          : EmittedRole::Supplementary);
        int reference_id = -1;
        if (segment.mapped()) {
            const auto found = reference_ids.find(segment.chromosome);
            if (found != reference_ids.end()) reference_id = found->second;
        }

        std::string tags;
        // RG:Z first, mirroring the SAM writer and minimap2.
        if (!opts.read_group_id.empty())
            bam_tag_Z(tags, 'R', 'G', opts.read_group_id);
        const AlignmentAuxiliaryTags auxiliary =
            alignment_auxiliary_tags(segment);
        if (auxiliary.alignment_score)
            bam_tag_i(tags, 'A', 'S', *auxiliary.alignment_score);
        // Tags in the SAM writer's order. md and s2 are not written on
        // secondaries.
        if (auxiliary.max_segment_score)
            bam_tag_i(tags, 'm', 's', *auxiliary.max_segment_score);
        if (auxiliary.max_score_margin && !emitted_role_is_secondary(role))
            bam_tag_i(tags, 'm', 'd', *auxiliary.max_score_margin);
        if (auxiliary.edit_distance)
            bam_tag_i(tags, 'N', 'M', *auxiliary.edit_distance);
        if (auxiliary.ambiguities)
            bam_tag_i(tags, 'n', 'n', *auxiliary.ambiguities);
        if (segment.mapped()) {
            bam_tag_A(tags, 't', 'p', emitted_role_paf_type(role, segment));
            if (auxiliary.secondary_chain_score &&
                !emitted_role_is_secondary(role))
                bam_tag_i(tags, 's', '2', *auxiliary.secondary_chain_score);
            const BlockAccounting accounting =
                block_accounting(segment, segment.read_len);
            double divergence = 0.0;
            if (event_divergence(accounting, divergence)) {
                // The rounded value SAM prints, so the two agree.
                bam_tag_f(tags, 'd', 'e',
                          std::strtof(format_paf_de(divergence).c_str(),
                                      nullptr));
            }
        }
        // RNA transcript strand (minimap2 ts:A); none when unknown.
        if (segment.mapped() && (segment.transcript_strand == '+' ||
                                 segment.transcript_strand == '-'))
          bam_tag_A(tags, 't', 's', segment.transcript_strand);
        if (!secondary_hypothesis && segments.size() > 1) {
            const std::string sa = build_sa_tag(segments, index);
            if (sa.size() > 5)
                bam_tag_Z(tags, 'S', 'A', sa.substr(5));
        }
        // XA:Z on the primary only. build_xa_tag's text carries the "XA:Z:"
        // prefix, which is stripped.
        if (!secondary_hypothesis && index == 0 && segment.mapped()) {
            const std::string xa = build_xa_tag(segment);
            if (xa.size() > 5)
                bam_tag_Z(tags, 'X', 'A', xa.substr(5));
        }
        // cs:Z and MD:Z after SA:Z, in minimap2's SAM order.
        if (segment.mapped() && segment.alignment_accounting_valid) {
          if (!segment.cs.empty())
            bam_tag_Z(tags, 'c', 's', segment.cs);
          if (!segment.md.empty())
            bam_tag_Z(tags, 'M', 'D', segment.md);
        }
        if (!secondary_hypothesis && index == 0) tags += primary_extra_tags;
        // Tags for every record of the read (-y), last.
        tags += record_extra_tags;

        const int flag =
            single_segment_flag(segment) |
            (segment.mapped() ? emitted_role_sam_flag(role) : 0);
        std::string cigar = cigar_or_fallback(segment);
        const std::string empty;
        const std::string* emitted_seq = secondary_hypothesis ? &empty : &seq;
        const std::string* emitted_qual = secondary_hypothesis ? &empty : &qual;
        bool reverse = segment.is_reverse;
        HardClippedRecord hard_clipped;
        if (!secondary_hypothesis && index > 0 && opts.hard_clip_supp && segment.mapped() &&
            !segment.cigar.empty()) {
            hard_clipped = hard_clip_supplementary(
                segment.cigar, seq, qual, segment.is_reverse);
            cigar = hard_clipped.cigar;
            emitted_seq = &hard_clipped.seq;
            emitted_qual = &hard_clipped.qual;
            reverse = false; // hard_clipped sequence is already oriented
        }

        records.emplace_back();
        encode_bam_record(
            records.back(), read_name, reference_id,
            segment.mapped() ? segment.pos : -1,
            segment.mapped() ? segment.mapq : 0, flag, cigar, *emitted_seq,
            *emitted_qual, reverse, tags);
      }
    };
    emit_hypothesis(result, /*secondary_hypothesis=*/false);
    // As in the SAM writer: no alternative record under an unmapped primary.
    if (opts.emit_secondary && result.mapped()) {
      for (const AlignResult& secondary : result.secondary) {
        if (secondary.mapped())
          emit_hypothesis(secondary, /*secondary_hypothesis=*/true);
      }
    }
    return records;
}

}}}  // namespace fa::cpu::output
