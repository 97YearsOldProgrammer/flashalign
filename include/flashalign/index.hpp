#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace flashalign {

class Aligner;

/**
 * A closed-syncmer seed index with its reference embedded. Copies share one index, so an
 * Index can be handed to several Aligners at no cost. Thread-safe for concurrent reads.
 */
class Index {
 public:
  /** An empty index; an Aligner rejects it. */
  Index();
  ~Index();
  Index(const Index&);
  Index& operator=(const Index&);
  Index(Index&&) noexcept;
  Index& operator=(Index&&) noexcept;

  /**
   * Loads a .faix written by `flashalign index` or save().
   * @throws std::runtime_error when the file cannot be loaded
   */
  static Index load(const std::string& path);
  /**
   * Builds an index over {name, sequence} pairs, embedding the reference. Contigs are sorted
   * by name and empty sequences skipped. Use preset_seeding() for a preset's k and s;
   * `threads` 0 uses every available core. The index records no preset.
   * @throws std::runtime_error when the build fails
   */
  static Index build(
      const std::vector<std::pair<std::string, std::string>>& references,
      int k = 21, int syncmer_s = 9, int threads = 0);
  /**
   * Builds an index from a FASTA file (plain, gzip or bgzf, or "-" for stdin), as build().
   * It maps the same as a `flashalign index` build with the same k and s.
   * @throws std::runtime_error when the file cannot be read, is not FASTA or FASTQ, or holds
   * no sequence
   */
  static Index build_from_fasta(
      const std::string& path, int k = 21, int syncmer_s = 9,
      int threads = 0);
  /**
   * Whether `path` is a .faix, judged from its header, not its name. Never throws; false for
   * "-" and for a path that cannot be opened.
   */
  static bool is_index_file(const std::string& path);

  /** Writes the index as a .faix. @throws std::runtime_error on failure */
  void save(const std::string& path) const;
  /** Seed k-mer length. */
  int k() const;
  /** Closed-syncmer s. */
  int syncmer_s() const;
  /** Whether the reference sequence is embedded (sequence() needs it). */
  bool has_reference() const;
  /** Total reference bases. */
  std::int64_t total_bases() const;
  /** Memory held by the index, in megabytes. */
  double memory_megabytes() const;
  /**
   * The preset recorded in the .faix header by `flashalign index -x`, or "" when none is
   * recorded (always for build() and build_from_fasta()). `flashalign align` without -x maps
   * under it.
   */
  std::string preset() const;
  /** Reference sequence names, in index order. */
  std::vector<std::string> reference_names() const;
  /** Reference sequence lengths, in reference_names() order. */
  std::vector<std::int64_t> reference_lengths() const;
  /**
   * The reference slice [start, end) of `name`, uppercase, with every non-ACGT base as 'N'.
   * `end == -1` means the contig end. Clamped as mappy's Aligner.seq(): a negative start is
   * 0, an end past the contig is its length, and start >= end gives "".
   *
   * The first call unpacks the whole reference at one byte per base (about 3 GB for a human
   * genome) and caches it for every copy of this Index; later calls only copy.
   * @throws std::out_of_range for an unknown name
   * @throws std::logic_error when the index embeds no reference
   */
  std::string sequence(
      const std::string& name, std::int64_t start = 0,
      std::int64_t end = -1) const;

 private:
  struct Impl;
  explicit Index(std::shared_ptr<Impl>);
  std::shared_ptr<Impl> impl_;
  friend class Aligner;
};

}  // namespace flashalign
