#include "output.h"

#include <algorithm>

namespace fa { namespace cpu { namespace lr { namespace rna {

namespace {

void stamp_chain_tags(AlignResult& record, const RnaChainTags& tags) {
    record.chain_anchors = tags.chain_anchors;
    record.chain_score = tags.chain_score;
    record.secondary_chain_score = tags.secondary_chain_score;
}

}  // namespace

bool commit_realized_alignment(
    const RnaSpliceRealizationResult& realized,
    int read_len,
    const int* chain_mapq,
    const RnaChainTags& chain_tags,
    const std::vector<RnaChimericFamily>& chimeric,
    Result& out) {
    if (realized.refused || realized.segments.empty()) return false;

    AlignResult best = realized.segments.front();
    for (size_t i = 1; i < realized.segments.size(); ++i)
        best.supplementary.push_back(realized.segments[i]);

    static_cast<AlignResult&>(out) = std::move(best);
    out.read_len = read_len;
    if (chain_mapq) {
    // Chain-formula MAPQ, computed upstream in the committed frame.
        out.mapq = *chain_mapq;
    } else {
        // Only the MAPQ input accounting failed: emit a mapped Q0 row, as the formula
        // does for an undecidable placement. The caller counts this.
        out.mapq = 0;
    }
    // Segments of one spliced family share the read-level MAPQ and the chain tags. A
    // chimeric family explains different query bases on its own chain against
    // different rivals, so it carries its own.
    stamp_chain_tags(out, chain_tags);
    for (auto& supplementary : out.supplementary) {
        supplementary.mapq = out.mapq;
        stamp_chain_tags(supplementary, chain_tags);
    }
    // Chimeric families are appended flat (io/emitted_role.h requires every
    // supplementary to be a leaf), in the caller's order, each segment stamped with
    // its family's MAPQ.
    for (const RnaChimericFamily& family : chimeric) {
        if (family.result == nullptr) continue;
        if (family.result->refused || family.result->segments.empty()) continue;
        for (const AlignResult& segment : family.result->segments) {
            out.supplementary.push_back(segment);
            out.supplementary.back().mapq = family.mapq;
            stamp_chain_tags(out.supplementary.back(), family.chain_tags);
        }
    }
    if (realized.transcript_strand == 1)
        out.transcript_strand = '+';
    else if (realized.transcript_strand == 2)
        out.transcript_strand = '-';
    return true;
}

void attach_runner_up_family(RnaSpliceRealizationResult& runner_up,
                             const RnaChainTags& chain_tags,
                             AlignResult& out) {
    if (runner_up.refused || runner_up.segments.empty()) return;
    // Flat, as above: segment 0 is the head and its continuations are leaf
    // supplementaries. MAPQ 0 on every record, as DNA does for a retained alternative:
    // the confidence belongs to the committed record. The chain tags are the family's
    // own, as for a chimeric family.
    AlignResult head = std::move(runner_up.segments.front());
    head.mapq = 0;
    stamp_chain_tags(head, chain_tags);
    head.supplementary.reserve(runner_up.segments.size() - 1);
    for (size_t i = 1; i < runner_up.segments.size(); ++i) {
        head.supplementary.push_back(std::move(runner_up.segments[i]));
        head.supplementary.back().mapq = 0;
        stamp_chain_tags(head.supplementary.back(), chain_tags);
    }
    out.secondary.push_back(std::move(head));
}

}}}}  // namespace fa::cpu::lr::rna
