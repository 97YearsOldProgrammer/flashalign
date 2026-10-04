// PAF record output.
#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <ostream>
#include <string>
#include <unordered_map>

#include "../core/types.h"   // AlignResult
#include "../core/cigar.h"   // parse_cigar_ops
#include "block_divergence.h"
#include "fastx.h"           // FastxRecord
#include "sam.h"             // reference_span, clamp_sam_mapq
#include "emitted_role.h"

namespace fa { namespace cpu { namespace output {

int target_end_from_result(const AlignResult& result);

std::string paf_cigar_from_sam_cigar(const std::string& cigar);

// `copy_comment` (-y) appends the read's FASTA/Q comment verbatim as the last
// field, as minimap2 does.
void write_paf_record(
    std::ostream& out,
    const ::fa::cpu::io::FastxRecord& read,
    const AlignResult& result,
    const std::unordered_map<std::string, int64_t>& ref_lengths,
    bool with_cigar,
    EmittedRole role = EmittedRole::Primary,
    bool copy_comment = false
);

// minimap2's --paf-no-hit row for a read with no alignment.
void write_paf_no_hit_record(
    std::ostream& out, const ::fa::cpu::io::FastxRecord& read);

}}}  // namespace fa::cpu::output
