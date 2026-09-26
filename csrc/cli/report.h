#pragma once

#include "cli/parse.h"       // AlignOptions
#include "api/aligner.h"     // BackendMode, ResolvedOptions

#include <cstddef>
#include <cstdint>
#include <optional>
#include <ostream>
#include <string>

namespace fa::cpu::cli {

// What --stats prints. Over a multi-part index the read counts sum over the
// passes; bases are counted once.
struct RunSummary {
  std::string preset;
  std::string output_format;
  bool rna = false;
  int64_t reads = 0;
  int64_t bases = 0;
  int64_t mapped_reads = 0;
  int64_t records = 0; // every record written, unmapped SAM/BAM rows included
  int64_t supplementary_records = 0;
  int64_t secondary_records = 0; // written only with --secondary yes
  // RNA: reads whose primary alignment has an intron, and the introns (N
  // operators) of all primary alignments.
  int64_t spliced_reads = 0;
  int64_t junctions = 0;
  double wall_seconds = 0.0; // from before the index is loaded
  double map_seconds = 0.0;  // the mapping passes alone, for throughput
  int64_t peak_rss_bytes = 0;
};

// Peak resident set size of this process, in bytes.
int64_t peak_rss_bytes();

// --stats: the run summary block, written to stderr.
void write_stats(const RunSummary& summary, std::ostream& out);

// The closing stderr lines of every align and index run: version, the command
// line as typed, and real time, CPU time and peak RSS.
void write_run_footer(int argc, char** argv, double real_seconds);

// --show-config: print the resolved configuration, each value tagged with its
// source; the caller then exits 0. `reference_sequences` is the sequence count
// when known (from an index header); nullopt omits the row.
void show_config(const AlignOptions& opt,
                 const fa::cpu::options::ResolvedOptions& cfg,
                 fa::cpu::api::BackendMode mode,
                 std::optional<std::size_t> reference_sequences,
                 std::ostream& out);

}  // namespace fa::cpu::cli
