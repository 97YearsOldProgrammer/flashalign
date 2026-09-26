#pragma once

#include "../core/types.h"
#include "../io/fastx.h"
#include "../io/reads.h"
#include "../io/sam.h"

#include <cstdint>
#include <iosfwd>
#include <streambuf>
#include <memory>
#include <string>
#include <unordered_map>

namespace fa::cpu::io {
class BgzfWriter;
}

namespace fa::cpu::cli {

struct AlignmentOutputContext {
    std::ostream& text;
    io::BgzfWriter* bam_writer = nullptr;
    const std::unordered_map<std::string, int>* bam_reference_ids = nullptr;
    const std::unordered_map<std::string, int64_t>* reference_lengths = nullptr;
    std::string format;
    bool include_unmapped = false;
    bool paf_cigar = false;
    // -y: copy the read's FASTA/Q comment onto the record (minimap2 -y).
    bool copy_comment = false;
    output::SamEmitOptions sam;
};

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
        const std::string& source_header_text,
        bool suppress_sam_header);
    AlignmentOutputContext context(
        bool include_unmapped,
        bool paf_cigar,
        bool copy_comment,
        output::SamEmitOptions sam);
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
    std::unordered_map<std::string, int> bam_reference_ids_;
    std::unique_ptr<io::BgzfWriter> bam_writer_;
};

// Writes one read's records in the requested format. Returns the number of
// records written.
int64_t write_alignment_records(
    const AlignmentOutputContext& context,
    const io::FastxRecord& read,
    const AlignResult& result);

}  // namespace fa::cpu::cli
