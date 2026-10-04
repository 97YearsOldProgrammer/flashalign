#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace flashalign {

/** How an alignment record was produced. */
enum class AlignmentOrigin : std::uint8_t {
  Ordinary,
  DnaLocalInversion,  ///< the inverted middle of a DNA alignment (PAF tp:A:I / i)
  DnaQueryPartition,  ///< DNA map-only placement, without a CIGAR
  DnaCigarFamily,     ///< DNA placement with a base-level CIGAR
};

/**
 * One alignment record. The object map() returns is the primary; its supplementary
 * segments and secondary hypotheses are nested inside it. Coordinates are 0-based and
 * half-open; query coordinates are on the forward read.
 */
struct Alignment {
  int score = 0;                ///< alignment score (AS:i), or the chain score without a CIGAR
  std::string chromosome;       ///< reference sequence name
  int pos = -1;                 ///< reference start; -1 when unmapped
  int read_len = 0;             ///< query length
  int mapq = 255;               ///< mapping quality, 0-60
  std::uint32_t median_occurrence = 0;  ///< median index occurrence of the placing seeds
  bool is_reverse = false;      ///< the read maps to the reverse strand
  AlignmentOrigin origin = AlignmentOrigin::Ordinary;
  std::string cigar;            ///< SAM CIGAR with soft clips; empty without base alignment
  std::vector<std::pair<int, int>> target_regions;  ///< reference intervals of a CIGAR-free record
  int query_start = 0;          ///< query start on the forward read
  int query_end = 0;            ///< query end on the forward read
  int target_end = 0;           ///< reference end
  int matches = 0;              ///< matching bases (PAF column 10)
  int mismatches = 0;
  int insertions = 0;           ///< inserted query bases
  int deletions = 0;            ///< deleted reference bases
  int ambiguities = 0;          ///< bases aligned to or from N, left out of matches and block_len
  int edit_distance = -1;       ///< NM:i; -1 without a CIGAR
  int block_len = 0;            ///< alignment block length (PAF column 11)
  bool alignment_accounting_valid = false;  ///< the counts above come from a base-level CIGAR
  /// cs and MD strings without their "cs:Z:" / "MD:Z:" prefixes; empty unless requested
  /// (Config::cs, Config::emit_md or a MapRequest) and the record has a CIGAR.
  std::string cs;
  std::string md;
  char transcript_strand = '\0';  ///< RNA: '+' or '-', '\0' when unknown
  std::vector<Alignment> supplementary;  ///< supplementary segments of this record
  /// Alternative placements weighed by MAPQ. Each may have supplementary segments but no
  /// secondary hypotheses of its own.
  std::vector<Alignment> secondary;
  /// Length of `chromosome` (mappy's ctg_len). Set by Aligner on every record it returns;
  /// 0 on an unmapped record.
  std::int64_t target_len = 0;
  /// True for a secondary hypothesis and for its supplementary segments; false for the
  /// primary and the primary's supplementary segments. Set by Aligner.
  bool is_secondary = false;
  /// cm:i, the number of anchors on the record's chain; -1 when the record has none.
  int chain_anchors = -1;
  /// s1:i, the score of the record's chain; -1 when the record has none.
  int chain_score = -1;
  /// s2:i, the chain score of the best competing chain; -1 when not computed.
  int secondary_chain_score = -1;
  /// ms:i, the score of the maximum-scoring segment of the CIGAR; -1 without a CIGAR.
  int dp_max_segment = -1;
  /// True for a DNA secondary beyond the first alternative (-N 2 and above): output
  /// only, weighed by no MAPQ, and left out of the primary's md:i and XA:Z.
  bool output_only = false;

  /** True when the record has a placement. */
  bool mapped() const {
    return pos >= 0 && (score > 0 || alignment_accounting_valid);
  }

  /** Not secondary, as mappy's is_primary; the primary's supplementary segments count. */
  bool is_primary() const { return !is_secondary; }

  /** matches / block_len, or 0 without base-level accounting. */
  double identity() const {
    return alignment_accounting_valid && block_len > 0
               ? static_cast<double>(matches) /
                     static_cast<double>(block_len)
               : 0.0;
  }
};

}  // namespace flashalign
