// SAM tag text to raw BAM aux bytes, the inverse of bam_reader.h's
// bam_tags_to_sam_text. Used by -y with BAM output, where the FASTA/Q comment
// must already be tag text (as `samtools fastq -T` writes it). Types: A, i
// (stored in the smallest integer code that fits, as htslib does), f, Z, H,
// and B with subtypes c C s S i I f.
#pragma once

#include <string>
#include <string_view>

namespace fa { namespace cpu { namespace io {

// True when `text` is tab-separated SAM tag fields (TAG:TYPE:VALUE with
// TAG = [A-Za-z][A-Za-z0-9] and TYPE one of A i f Z H B). Empty text is valid
// (no fields). On false, a non-null *offending_field receives the first
// field that is not tag text.
bool sam_tag_text_valid(std::string_view text,
                        std::string* offending_field = nullptr);

// Encode tag text as a raw BAM tag block. `text` must already have passed
// sam_tag_text_valid; a field that is not tag text throws std::runtime_error.
std::string sam_tag_text_to_bam(std::string_view text);

}}}  // namespace fa::cpu::io
