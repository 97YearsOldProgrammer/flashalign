#include "paf.h"

namespace fa { namespace cpu { namespace output {

int target_end_from_result(const AlignResult& result) {
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

std::string paf_cigar_from_sam_cigar(const std::string& cigar) {
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
    // Shared with SAM (io/block_divergence.h), so NM:i and de:f agree.
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
