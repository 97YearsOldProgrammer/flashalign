// SAM record and header output.
#pragma once

#include "../core/cigar.h"
#include "emitted_role.h"
#include "../core/types.h"
#include "alignment_aux_tags.h"
#include "block_divergence.h"
#include "fa_version.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <climits>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace output {

struct SamReference {
    std::string name;
    int64_t length = 0;
};

int clamp_sam_mapq(int mapq);

std::string clean_sam_field(const std::string& field);

std::string cigar_or_fallback(const AlignResult& result);

int reference_span_from_cigar(const std::string& cigar, int fallback_len);

int reference_span(const AlignResult& result);

void add_int_tag(std::vector<std::string>& tags, const std::string& tag, int value);

void add_string_tag(std::vector<std::string>& tags, const std::string& tag, const std::string& value);

int single_segment_flag(const AlignResult& result);

void append_sam_int(std::string& out, int64_t value);

void append_clean_sam_field(std::string& out, const std::string& field);

void append_sam_cigar(std::string& out, const AlignResult& result);

char complement_base(char c);

// Appends the reverse complement of `seq`, for reverse-strand SEQ.
void append_reverse_complement(std::string& out, const std::string& seq);

void append_sam_record_line(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::string& extra_tags = {},
    int flag_extra = 0,  // OR-ed into FLAG when mapped (e.g. 0x800 supplementary)
    // For hard-clipped supplementaries: `cigar_override` replaces
    // result.cigar, and `seq_preoriented` means seq/qual are already
    // reference-oriented and clipped.
    const std::string* cigar_override = nullptr,
    bool seq_preoriented = false,
    // RG:Z value; empty emits no tag.
    const std::string& read_group_id = {},
    // For tp:A; consistent with flag_extra.
    EmittedRole role = EmittedRole::Primary,
    // Tag text written last, after cs/MD, on every record of the read (the -y
    // comment). extra_tags, by contrast, belongs to this record only.
    const std::string& record_tags = {}
);

// SA:Z for segment `self`, one "rname,pos,strand,CIGAR,mapQ,NM;" entry per
// other mapped segment. Throws if a segment lacks exact accounting for NM.
// Entries carry the full CIGAR, not minimap2's condensed form.
std::string build_sa_tag(const std::vector<const AlignResult*>& segs, size_t self);

// XA:Z from the read's alternatives, as bwa and minibwa write it:
// "chr,{+|-}pos,CIGAR,NM;" per mapped secondary head, POS 1-based. Entries
// without a CIGAR or exact accounting are skipped. Empty when there are none.
std::string build_xa_tag(const AlignResult& result);

// A supplementary segment hard-clipped as minimap2 does.
struct HardClippedRecord {
    std::string cigar;  // terminal soft clips rewritten to hard clips ('S' -> 'H')
    std::string seq;    // reference-oriented read, aligned portion only (clips dropped)
    std::string qual;   // matching QUAL slice; empty when the read has no QUAL
};

// Turns terminal soft clips into hard clips and drops the clipped bases from
// SEQ/QUAL. The returned seq/qual are reference-oriented. An absent or "*"
// read field gives an empty slice.
HardClippedRecord hard_clip_supplementary(
    const std::string& soft_cigar, const std::string& fwd_seq,
    const std::string& fwd_qual, bool is_reverse);

struct SamEmitOptions {
    bool hard_clip_supp = false;  // the CLI sets this unless -Y (minimap2's default)
    // One record per mapped alternative (FLAG 0x100), as --secondary yes.
    bool emit_secondary = false;
    // RG:Z for every record of the read; empty for none, including a read
    // that already carries its own RG:Z.
    std::string read_group_id;
};

// All SAM records of one read: the primary, then its supplementaries (FLAG
// 0x800), each with an SA:Z listing the others, then secondaries if enabled.
// Only supplementaries are hard-clipped; SA:Z lists the soft CIGARs, as
// minimap2 does. extra_tags (the uBAM passthrough) goes on the primary only,
// record_tags on every record.
void append_sam_records(
    std::string& out,
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::string& extra_tags = {},
    SamEmitOptions opts = {},
    const std::string& record_tags = {}
);

// @HD, @SQ and @PG lines. A non-empty `command_line` becomes @PG CL:, with
// tabs and newlines neutralized so an argument cannot forge a header line.
std::vector<std::string> sam_header(
    const std::vector<SamReference>& references,
    const std::string& program_name = "flashalign",
    const std::string& version = FA_VERSION,
    const std::string& command_line = {}
);

}}} // namespace fa::cpu::output
