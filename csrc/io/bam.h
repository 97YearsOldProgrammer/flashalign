// BAM header and record encoding, per the SAMv1 spec. A CIGAR of more than
// 65535 operations is stored in a CG:B,I tag, as the spec requires.
#pragma once

#include "../core/cigar.h"
#include "sam.h"
#include "emitted_role.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace output {

// Little-endian field appenders.
void bam_put_i32(std::string& b, int32_t v);
void bam_put_u32(std::string& b, uint32_t u);
void bam_put_u16(std::string& b, uint16_t u);

// Tag-block appenders: two tag characters, a type byte, then the value.
void bam_tag_i(std::string& tags, char a, char b, int32_t v);
void bam_tag_A(std::string& tags, char a, char b, char v);
void bam_tag_f(std::string& tags, char a, char b, float v);
void bam_tag_Z(std::string& tags, char a, char b, const std::string& v);

// Base character to 4-bit code, per "=ACMGRSVTWYHKDBN".
const std::array<uint8_t, 256>& seq_nt16_table();

int bam_cigar_op_code(char op);

// reg2bin from SAMv1 section 5.3, with int64 arguments as htslib's hts_reg2bin
// so beg + span cannot overflow. Unmapped records use reg2bin(-1, 0) = 4680.
// The scheme addresses only 2^29 bp (see encode_bam_record).
int bam_reg2bin(int64_t beg, int64_t end);

// BAM header: magic, SAM text, then n_ref x (l_name, name\0, l_ref).
std::string encode_bam_header(
    const std::vector<SamReference>& refs,
    const std::string& sam_header_text);

// Appends one block_size-prefixed BAM record to `out`. seq_in/qual_in are the
// read as sequenced; a reverse-strand record stores the reverse complement and
// reversed qualities, as in SAM. `extra_tags` is a raw tag block.
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
    const std::string& extra_tags = {});

// BAM counterpart of sam.h's append_sam_records, with the same records and
// tags. `primary_extra_tags` (the uBAM passthrough) goes on the primary only;
// `record_extra_tags` (the -y comment) goes last on every record. Both are raw
// BAM tag blocks.
std::vector<std::string> encode_bam_records(
    const std::string& read_name,
    const AlignResult& result,
    const std::string& seq,
    const std::string& qual,
    const std::unordered_map<std::string, int>& reference_ids,
    const std::string& primary_extra_tags = {},
    SamEmitOptions opts = {},
    const std::string& record_extra_tags = {});

}}}  // namespace fa::cpu::output
