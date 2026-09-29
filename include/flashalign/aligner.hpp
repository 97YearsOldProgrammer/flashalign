#pragma once

#include "alignment.hpp"
#include "config.hpp"
#include "index.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace flashalign {

/**
 * Per-call cs/MD request, like mappy's map(cs=True, MD=True). An empty field keeps the
 * configured value (Config::cs, Config::emit_md). The request applies to one call only and
 * config() never reports it.
 *
 * A request that differs from the currently installed one reconfigures the aligner before
 * mapping: cheap (no index reload), but paid on every change, so a caller that wants cs or MD
 * on every read should set them in the Config instead. Every preset accepts both strings;
 * spliced records carry minimap2's `~` intron token in cs.
 */
struct MapRequest {
  std::optional<CsMode> cs;  ///< empty: Config::cs
  std::optional<bool> md;    ///< empty: Config::emit_md
};

/**
 * Maps reads against one reference, in the manner of mappy's Aligner.
 *
 * Thread safety: map(), map_batch(), paf(), config(), reference_names() and
 * reference_lengths() may be called concurrently on one Aligner. reconfigure(), and a map
 * call whose MapRequest differs from the installed one, change the configuration under an
 * exclusive lock; a concurrent map may then run under either configuration. Threads that
 * need different cs/MD settings at the same time should each use their own Aligner over a
 * shared Index.
 */
class Aligner {
 public:
  /**
   * Aligner over a loaded or built index. k and s always come from the index.
   * @throws std::invalid_argument for an empty index, an unknown preset or a Config the
   *         option resolver rejects
   */
  Aligner(Index index, Config config = {});
  /**
   * Builds an index over `references` ({name, sequence} pairs; empty sequences are skipped)
   * with the Config's k and s (the preset's by default), then an aligner over it.
   * @throws std::invalid_argument as the constructor above
   */
  Aligner(
      const std::vector<std::pair<std::string, std::string>>& references,
      Config config = {});
  ~Aligner();
  Aligner(Aligner&&) noexcept;
  Aligner& operator=(Aligner&&) noexcept;
  Aligner(const Aligner&) = delete;
  Aligner& operator=(const Aligner&) = delete;

  /**
   * Maps one read. The result is the primary alignment, with its supplementary segments
   * and secondary hypotheses nested inside; check Alignment::mapped(). Without a read name,
   * ties between equally good loci are broken from the read length alone, as minimap2 does
   * without a query name.
   */
  Alignment map(const std::string& read) const;
  /** Maps reads in parallel on Config::threads threads; one result per read, in order. */
  std::vector<Alignment> map_batch(
      const std::vector<std::string>& reads) const;
  /** map() with a per-call cs/MD request (see MapRequest). */
  Alignment map(const std::string& read, const MapRequest& request) const;
  /** map_batch() with one cs/MD request for the whole batch. */
  std::vector<Alignment> map_batch(
      const std::vector<std::string>& reads,
      const MapRequest& request) const;
  /**
   * The PAF lines the CLI writes for one read: the primary and its supplementary segments,
   * newline-terminated, without secondary rows. `with_cigar` adds the cg:Z tag (the CLI's
   * -c). The name also seeds the tie-break, as in the CLI.
   * @return the lines, or "" for an unmapped read
   */
  std::string paf(
      const std::string& name, const std::string& read,
      bool with_cigar = false) const;
  /**
   * The resolved configuration. Passed back to reconfigure() it maps as before, with the
   * resolved values now explicit settings.
   */
  Config config() const;
  /**
   * Installs a new configuration; k and s stay the index's.
   * @throws std::invalid_argument when the option resolver rejects `config`
   */
  void reconfigure(Config config);
  /** Reference sequence names, in index order. */
  std::vector<std::string> reference_names() const;
  /** Reference sequence lengths, in reference_names() order. */
  std::vector<std::int64_t> reference_lengths() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace flashalign
