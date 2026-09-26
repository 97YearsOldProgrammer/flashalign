#include "bindings.h"

#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

// mappy's operation alphabet; Hit.cigar holds indices into it.
constexpr const char* kCigarOps = "MIDNSHP=XB";

std::vector<std::pair<int, char>> split_cigar(const std::string& cigar) {
  std::vector<std::pair<int, char>> ops;
  int length = 0;
  bool digits = false;
  for (const char letter : cigar) {
    if (letter >= '0' && letter <= '9') {
      length = length * 10 + (letter - '0');
      digits = true;
      continue;
    }
    if (!digits)
      return {};
    ops.emplace_back(length, letter);
    length = 0;
    digits = false;
  }
  return digits ? std::vector<std::pair<int, char>>{} : ops;
}

// The SAM CIGAR without its terminal soft clips, as in PAF cg:Z and mappy's Alignment.cigar.
std::string clipped_cigar(const std::string& sam_cigar) {
  const auto ops = split_cigar(sam_cigar);
  if (ops.empty())
    return sam_cigar;
  std::size_t begin = 0;
  std::size_t end = ops.size();
  if (ops[begin].second == 'S')
    ++begin;
  if (end > begin && ops[end - 1].second == 'S')
    --end;
  std::string out;
  for (std::size_t i = begin; i < end; ++i) {
    out += std::to_string(ops[i].first);
    out += ops[i].second;
  }
  return out;
}

std::vector<std::pair<int, int>> cigar_tuples(const std::string& sam_cigar) {
  std::vector<std::pair<int, int>> tuples;
  for (const auto& op : split_cigar(clipped_cigar(sam_cigar))) {
    const char* found = nullptr;
    for (const char* scan = kCigarOps; *scan != '\0'; ++scan) {
      if (*scan == op.second) {
        found = scan;
        break;
      }
    }
    if (found == nullptr) {
      throw std::runtime_error(
          std::string("flashalign: CIGAR operation not in \"MIDNSHP=XB\": ") +
          op.second);
    }
    tuples.emplace_back(op.first, static_cast<int>(found - kCigarOps));
  }
  return tuples;
}

int trans_strand_of(const flashalign::Alignment& record) {
  if (record.transcript_strand == '+')
    return 1;
  if (record.transcript_strand == '-')
    return -1;
  return 0;
}

} // namespace

std::vector<Hit> flatten_hits(flashalign::Alignment record) {
  std::vector<Hit> hits;
  if (!record.mapped())
    return hits;
  std::vector<flashalign::Alignment> segments;
  segments.swap(record.supplementary);
  hits.push_back(Hit{std::move(record), false});
  for (flashalign::Alignment& segment : segments) {
    if (!segment.mapped())
      continue;
    hits.push_back(Hit{std::move(segment), true});
  }
  return hits;
}

void bind_alignment(nb::module_& module) {
  nb::class_<Hit>(
      module, "Hit",
      "One alignment record, with mappy's field names.\n"
      "\n"
      "Hits are read-only and come only from Aligner.map() and map_batch(),\n"
      "in the CLI's record order: the primary, then its supplementary\n"
      "segments. Coordinates are 0-based, half-open and on the forward\n"
      "strand of both sequences, as in PAF and mappy.")
      .def_prop_ro(
          "ctg", [](const Hit& hit) { return hit.record.chromosome; },
          "Target (reference) contig name. PAF column 6.")
      .def_prop_ro(
          "ctg_len", [](const Hit& hit) { return hit.record.target_len; },
          "Length of `ctg`, in bases. PAF column 7.")
      .def_prop_ro(
          "r_st", [](const Hit& hit) { return hit.record.pos; },
          "Target start, 0-based. PAF column 8.")
      .def_prop_ro(
          "r_en", [](const Hit& hit) { return hit.record.target_end; },
          "Target end, exclusive. PAF column 9.")
      .def_prop_ro(
          "q_st", [](const Hit& hit) { return hit.record.query_start; },
          "Query start on the forward strand, 0-based. PAF column 3.")
      .def_prop_ro(
          "q_en", [](const Hit& hit) { return hit.record.query_end; },
          "Query end on the forward strand, exclusive. PAF column 4.")
      .def_prop_ro(
          "strand",
          [](const Hit& hit) { return hit.record.is_reverse ? -1 : 1; },
          "+1 if the read maps forward, -1 if reverse. PAF column 5.")
      .def_prop_ro(
          "mapq", [](const Hit& hit) { return hit.record.mapq; },
          "Mapping quality, 0-60. PAF column 12.")
      .def_prop_ro(
          "NM", [](const Hit& hit) { return hit.record.edit_distance; },
          "Edit distance over the aligned block (NM:i); -1 without a CIGAR.")
      .def_prop_ro(
          "blen", [](const Hit& hit) { return hit.record.block_len; },
          "Aligned block length including gaps. PAF column 11.")
      .def_prop_ro(
          "mlen", [](const Hit& hit) { return hit.record.matches; },
          "Residue matches in the aligned block. PAF column 10.")
      .def_prop_ro(
          "is_primary", [](const Hit& hit) { return hit.record.is_primary(); },
          "True unless this record is a rival hypothesis (tp:A:S).")
      .def_prop_ro(
          "trans_strand", &trans_strand_of,
          "Transcript strand: +1 / -1 / 0 for unknown (ts:A, RNA only).")
      .def_prop_ro(
          "cigar",
          [](const Hit& hit) { return cigar_tuples(hit.record.cigar); },
          "Aligned block as [(length, op), ...]; op indexes \"MIDNSHP=XB\".")
      .def_prop_ro(
          "cigar_str",
          [](const Hit& hit) { return clipped_cigar(hit.record.cigar); },
          "The same CIGAR as a string; the PAF cg:Z value (clips dropped).")
      .def_prop_ro(
          "cs", [](const Hit& hit) { return hit.record.cs; },
          "cs:Z difference string; \"\" unless Config.cs or map(cs=...) asks\n"
          "for it.")
      .def_prop_ro(
          "MD", [](const Hit& hit) { return hit.record.md; },
          "MD:Z difference string; \"\" unless Config.emit_md or map(MD=True)\n"
          "asks for it.")
      .def_prop_ro(
          "score", [](const Hit& hit) { return hit.record.score; },
          "Alignment score (the AS:i tag), or the chain score without a CIGAR.")
      .def_prop_ro(
          "read_len", [](const Hit& hit) { return hit.record.read_len; },
          "Length of the query this record was produced from.")
      .def_prop_ro(
          "is_supplementary", [](const Hit& hit) { return hit.supplementary; },
          "True for a supplementary segment of the record before it.")
      .def_prop_ro(
          "identity", [](const Hit& hit) { return hit.record.identity(); },
          "mlen / blen, or 0.0 without a base-level alignment.")
      .def_prop_ro(
          "mismatches", [](const Hit& hit) { return hit.record.mismatches; },
          "Substituted bases in the aligned block.")
      .def_prop_ro(
          "insertions", [](const Hit& hit) { return hit.record.insertions; },
          "Inserted query bases in the aligned block.")
      .def_prop_ro(
          "deletions", [](const Hit& hit) { return hit.record.deletions; },
          "Deleted reference bases in the aligned block.")
      .def_prop_ro(
          "ambiguities", [](const Hit& hit) { return hit.record.ambiguities; },
          "Ambiguous (N) bases excluded from mlen and blen (the nn:i tag).")
      .def_prop_ro(
          "is_mapped", [](const Hit& hit) { return hit.record.mapped(); },
          "True for every Hit map() returns; unmapped reads yield no Hit.")
      .def_prop_ro(
          "secondary",
          [](const Hit& hit) {
            std::vector<Hit> alternatives;
            if (hit.supplementary || hit.record.is_secondary)
              return alternatives;
            for (const flashalign::Alignment& alternative :
                 hit.record.secondary) {
              for (Hit& nested : flatten_hits(alternative))
                alternatives.push_back(std::move(nested));
            }
            return alternatives;
          },
          "Rival placement hypotheses, on the primary only.\n"
          "\n"
          "These are what MAPQ was computed against, not records of their\n"
          "own: the CLI writes secondary rows only under --secondary yes, and\n"
          "map() never returns them. Each rival is flattened like the primary\n"
          "(the rival, then its supplementary segments). Empty on a\n"
          "supplementary Hit and on a rival.")
      .def(
          "__str__",
          [](const Hit& hit) {
            const flashalign::Alignment& record = hit.record;
            std::ostringstream out;
            out << record.query_start << '\t' << record.query_end << '\t'
                << (record.is_reverse ? '-' : '+') << '\t' << record.chromosome
                << '\t' << record.target_len << '\t' << record.pos << '\t'
                << record.target_end << '\t' << record.matches << '\t'
                << record.block_len << '\t' << record.mapq << '\t'
                << (record.is_primary() ? "tp:A:P" : "tp:A:S") << '\t';
            const int strand = trans_strand_of(record);
            out << (strand > 0   ? "ts:A:+"
                    : strand < 0 ? "ts:A:-"
                                 : "ts:A:.")
                << "\tcg:Z:" << clipped_cigar(record.cigar);
            if (!record.cs.empty())
              out << "\tcs:Z:" << record.cs;
            if (!record.md.empty())
              out << "\tMD:Z:" << record.md;
            return out.str();
          },
          "mappy's line: the PAF row without its first two columns.\n"
          "\n"
          "q_st q_en strand ctg ctg_len r_st r_en mlen blen mapq tp:A ts:A\n"
          "cg:Z [cs:Z] [MD:Z], tab-separated. Aligner.paf() gives the full\n"
          "PAF row, query name and length included.")
      .def(
          "__repr__",
          [](const Hit& hit) {
            std::ostringstream out;
            out << "<flashalign.Hit " << hit.record.chromosome << ':'
                << hit.record.pos << '-' << hit.record.target_end << ' '
                << (hit.record.is_reverse ? '-' : '+')
                << " mapq=" << hit.record.mapq;
            if (hit.supplementary)
              out << " supplementary";
            if (hit.record.is_secondary)
              out << " secondary";
            out << '>';
            return out.str();
          },
          "A short one-line summary; str(hit) gives the PAF-style line.");
}
