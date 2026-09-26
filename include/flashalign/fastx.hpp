#pragma once

#include <memory>
#include <string>

namespace flashalign {

/**
 * One FASTA/FASTQ record. `qual` is empty for FASTA. `comment` is the rest of the header
 * line after the name and the whitespace that follows it, as kseq.h splits it.
 */
struct FastxRecord {
  std::string name;
  std::string comment;
  std::string seq;
  std::string qual;
};

/**
 * Streaming FASTA/FASTQ reader, the one the CLI uses. Reads plain, gzip and bgzf files
 * (multi-member gzip included) and "-" for stdin; not unaligned BAM. Sequences are returned
 * uppercased; names, comments and qualities as they are.
 */
class FastxReader {
 public:
  /** @throws std::runtime_error if the path cannot be opened */
  explicit FastxReader(const std::string& path);
  ~FastxReader();
  FastxReader(FastxReader&&) noexcept;
  FastxReader& operator=(FastxReader&&) noexcept;
  FastxReader(const FastxReader&) = delete;
  FastxReader& operator=(const FastxReader&) = delete;

  /**
   * Reads the next record into `out`.
   * @return true on success, false at end of file
   * @throws std::runtime_error on a read error, a truncated FASTQ or a quality/sequence
   *         length mismatch
   */
  bool next(FastxRecord& out);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace flashalign
