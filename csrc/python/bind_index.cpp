#include "bindings.h"

#include <flashalign/index.hpp>

#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

std::optional<std::string> index_sequence(const flashalign::Index& index,
                                          const std::string& name,
                                          std::int64_t start,
                                          std::int64_t end) {
  try {
    return index.sequence(name, start, end);
  } catch (const std::out_of_range&) {
    // As mappy's Aligner.seq(): None for an unknown contig.
    return std::nullopt;
  }
}

} // namespace

void bind_index(nb::module_& module) {
  nb::class_<flashalign::Index>(
      module, "Index",
      "A closed-syncmer seed index, with its reference embedded.\n"
      "\n"
      "An Index is reference-counted, so several Aligners can share one\n"
      "copy; Aligner.index returns the one an aligner uses.\n"
      "\n"
      "    idx = fa.Index.load('ref.faix')       # a .faix written by the CLI\n"
      "    idx = fa.Index.build_from_fasta('ref.fa', "
      "*fa.preset_seeding('lr'))\n"
      "\n"
      "k and syncmer_s are fixed when the index is built, and mapping always\n"
      "uses the index's pair, never a preset's.")
      .def_static(
          "load",
          [](const std::filesystem::path& path) {
            return flashalign::Index::load(path.string());
          },
          nb::arg("path"), nb::call_guard<nb::gil_scoped_release>(),
          "Load a .faix seed index written by `flashalign index`.")
      .def_static(
          "build",
          [](const std::vector<std::pair<std::string, std::string>>& sequences,
             int k, int syncmer_s, int threads) {
            return flashalign::Index::build(sequences, k, syncmer_s, threads);
          },
          nb::arg("sequences"), nb::arg("k") = 21, nb::arg("syncmer_s") = 9,
          nb::arg("threads") = 0,
          nb::call_guard<nb::gil_scoped_release>(),
          "Build an index from [(name, seq), ...], embedding the reference.\n"
          "\n"
          "Contigs are sorted by name, as the CLI sorts them, so the\n"
          "alignments match the CLI's; empty sequences are skipped. The index\n"
          "records no preset name.")
      .def_static(
          "build_from_fasta",
          [](const std::filesystem::path& path, int k, int syncmer_s,
             int threads) {
            return flashalign::Index::build_from_fasta(path.string(), k,
                                                       syncmer_s, threads);
          },
          nb::arg("path"), nb::arg("k") = 21, nb::arg("syncmer_s") = 9,
          nb::arg("threads") = 0,
          nb::call_guard<nb::gil_scoped_release>(),
          "Build an index from a reference FASTA on disk.\n"
          "\n"
          "Plain, gzip or bgzf, or '-' for stdin, read with the CLI's reader.\n"
          "Pass *preset_seeding(name) for k and syncmer_s if the aligner will\n"
          "run under a preset other than 'lr'.")
      .def_static(
          "is_index_file",
          [](const std::filesystem::path& path) {
            return flashalign::Index::is_index_file(path.string());
          },
          nb::arg("path"), nb::call_guard<nb::gil_scoped_release>(),
          "Is this path a .faix seed index?\n"
          "\n"
          "Decided from the file header, not the name, as `flashalign align`\n"
          "decides. False, never an exception, for '-', for a path that\n"
          "cannot be opened and for a file that is not a .faix.")
      .def(
          "save",
          [](const flashalign::Index& self, const std::filesystem::path& path) {
            self.save(path.string());
          },
          nb::arg("path"), nb::call_guard<nb::gil_scoped_release>(),
          "Write this index to `path` as a .faix.")
      .def("seq", &index_sequence, nb::arg("name"), nb::arg("start") = 0,
           nb::arg("end") = -1, nb::call_guard<nb::gil_scoped_release>(),
           "One slice of the embedded reference, or None for an unknown "
           "contig.\n"
           "\n"
           "Uppercase, with every non-ACGT base as 'N'. end=-1 means the\n"
           "contig end, and the range is clamped like mappy's Aligner.seq():\n"
           "start<0 becomes 0, end past the contig becomes its length, and\n"
           "start>=end gives ''.\n"
           "\n"
           "The first call unpacks the whole reference at one byte per base\n"
           "(about 3 GB for a human genome) and caches it on the index; later\n"
           "calls only copy.")
      .def_prop_ro("k", &flashalign::Index::k,
                   "Seed k-mer length this index was built with.")
      .def_prop_ro("syncmer_s", &flashalign::Index::syncmer_s,
                   "Closed-syncmer s this index was built with.")
      .def_prop_ro("has_reference", &flashalign::Index::has_reference,
                   "True when the index embeds its reference (seq() needs it).")
      .def_prop_ro("memory_mb", &flashalign::Index::memory_megabytes,
                   "Memory held by the index, in megabytes.")
      .def_prop_ro("total_bp", &flashalign::Index::total_bases,
                   "Total reference bases indexed.")
      .def_prop_ro("preset", &flashalign::Index::preset,
                   "The preset name the .faix header records, or ''.\n"
                   "\n"
                   "`flashalign index -x NAME` records it and\n"
                   "`flashalign align ref.faix reads.fq` maps under it. An\n"
                   "index built here records none, and Aligner then uses\n"
                   "'lr' unless a preset is given.")
      .def_prop_ro("seq_names", &flashalign::Index::reference_names,
                   "Contig names, in index order (sorted by name).")
      .def_prop_ro("seq_lengths", &flashalign::Index::reference_lengths,
                   "Contig lengths in bases, in seq_names order.")
      .def_prop_ro(
          "n_seq",
          [](const flashalign::Index& self) {
            return self.reference_names().size();
          },
          "Number of contigs in the index.")
      .def(
          "__repr__",
          [](const flashalign::Index& self) {
        std::string out = "<flashalign.Index k=" + std::to_string(self.k()) +
                          " s=" + std::to_string(self.syncmer_s()) + " n_seq=" +
                          std::to_string(self.reference_names().size()) +
                          " total_bp=" + std::to_string(self.total_bases());
        const std::string preset = self.preset();
        if (!preset.empty())
          out += " preset=" + preset;
        return out + ">";
          },
          "The seeding, the contig count, the size and any recorded preset.");
}
