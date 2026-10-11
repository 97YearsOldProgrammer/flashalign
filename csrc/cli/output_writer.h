#pragma once

#include "../core/types.h"
#include "../io/fastx.h"
#include "../io/reads.h"
#include "../io/sam.h"

#include <atomic>
#include <cstdint>
#include <iosfwd>
#include <streambuf>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace fa::cpu::cli {

// The tag keys FlashAlign writes on SAM records. An input tag with one of
// them is dropped from the passthrough in favour of ours, with one warning
// per key.
class OwnTagKeys {
 public:
    explicit OwnTagKeys(std::vector<std::string> keys);
    // True when the two characters at `key` are one of ours.
    bool drop(const char* key) const;

 private:
    std::vector<std::string> keys_;
    std::unique_ptr<std::atomic<bool>[]> warned_;
};

struct AlignmentOutputContext {
    std::ostream& text;
    const std::unordered_map<std::string, int64_t>* reference_lengths = nullptr;
    std::string format;
    bool include_unmapped = false;
    bool paf_cigar = false;
    // -y: copy the read's FASTA/Q comment onto the record (minimap2 -y).
    bool copy_comment = false;
    output::SamEmitOptions sam;
    const OwnTagKeys* own_tag_keys = nullptr;
};

// The header lines a run carries over from its BAM inputs, in input order:
// every @PG line (the first of each ID), the union of the @RG lines by ID,
// and every @CO line.
struct InputHeader {
    std::vector<std::string> program_lines;
    // @RG and @CO lines, in the order read.
    std::vector<std::string> passthrough;
    // The ID of the last @PG line of the first input that has one: our @PG's PP.
    std::string previous_program;
};

// Reads the header of every input that is a BAM file. Two inputs declaring
// one read group ID with different @RG lines is an error naming both.
InputHeader read_input_headers(const std::vector<std::string>& paths);

class AlignmentOutputWriter {
 public:
    AlignmentOutputWriter(
        const std::string& path,
        const std::string& format,
        const io::GenomeLoad& genome,
        // @PG CL: value; empty omits the field.
        const std::string& command_line = {},
        // @RG line to declare, already validated and unescaped; empty means
        // no read group.
        const std::string& read_group_line = {});
    ~AlignmentOutputWriter();

    AlignmentOutputWriter(const AlignmentOutputWriter&) = delete;
    AlignmentOutputWriter& operator=(const AlignmentOutputWriter&) = delete;

    void write_header(
        const InputHeader& input_header,
        bool suppress_sam_header);
    // `own_tag_keys`: the SAM tag keys this run writes (see OwnTagKeys).
    AlignmentOutputContext context(
        bool include_unmapped,
        bool paf_cigar,
        bool copy_comment,
        output::SamEmitOptions sam,
        std::vector<std::string> own_tag_keys);
    void close();
    void flush();

 private:
    std::string format_;
    std::string command_line_;
    std::string read_group_line_;
    std::string read_group_id_;
    // Declared first, so the stream below is destroyed (and flushed) before it.
    std::unique_ptr<std::streambuf> stdout_sink_;
    std::unique_ptr<std::ostream> owner_;
    std::ostream* text_ = nullptr;
    std::vector<output::SamReference> references_;
    std::unordered_map<std::string, int64_t> reference_lengths_;
    std::unique_ptr<OwnTagKeys> own_tag_keys_;
};

// Writes one read's records in the requested format. Returns the number of
// records written.
int64_t write_alignment_records(
    const AlignmentOutputContext& context,
    const io::FastxRecord& read,
    const AlignResult& result);

}  // namespace fa::cpu::cli
