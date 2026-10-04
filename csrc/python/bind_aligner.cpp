#include "bindings.h"

#include <flashalign/aligner.hpp>
#include <flashalign/index.hpp>

#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

// The public Aligner plus its Index, so Aligner.index can share the index with another
// Aligner.
struct PyAligner {
  flashalign::Index index;
  std::unique_ptr<flashalign::Aligner> aligner;
};

namespace {

[[noreturn]] void raise_not_found(const std::string& path) {
  const std::string message = "flashalign: no such file: " + path;
  PyErr_SetString(PyExc_FileNotFoundError, message.c_str());
  throw nb::python_error();
}

std::optional<std::string> optional_string(nb::handle value, const char* what) {
  if (value.is_none())
    return std::nullopt;
  if (!nb::isinstance<nb::str>(value)) {
    throw nb::type_error(
        (std::string("flashalign: ") + what + " must be a string or None")
            .c_str());
  }
  return nb::cast<std::string>(value);
}

// None keeps the configured value; False turns it off and True/'short'/'long' on, for this
// call.
flashalign::MapRequest make_request(nb::handle cs, nb::handle md) {
  flashalign::MapRequest request;
  if (!cs.is_none()) request.cs = parse_cs_mode(cs);
  if (!md.is_none()) {
    if (!nb::isinstance<nb::bool_>(md))
      throw nb::type_error("flashalign: MD must be a bool or None");
    request.md = nb::cast<bool>(md);
  }
  return request;
}

} // namespace

void bind_aligner(nb::module_& module) {
  nb::class_<PyAligner>(
      module, "Aligner",
      "An aligner over one reference, in the mappy manner.\n"
      "\n"
      "    a = fa.Aligner('ref.fa', preset='lr:hq')   # or 'ref.faix', or an "
      "Index\n"
      "    for h in a.map(read_seq):\n"
      "        print(h.ctg, h.r_st, h.r_en, h.strand, h.mapq, h.cigar_str)\n"
      "\n"
      "`target` is a path (str or os.PathLike) or an Index. A path is judged\n"
      "by its content, not its name: a .faix is loaded, anything else is\n"
      "indexed as a reference FASTA with the preset's (k, s).\n"
      "\n"
      "`preset` defaults to the preset the .faix records, or 'lr' when it\n"
      "records none, as `flashalign align ref.faix reads.fq` does. An\n"
      "explicit preset wins over the index's.\n"
      "\n"
      "`config` replaces the preset-derived configuration and is used as\n"
      "given; `threads` then only sizes an index build from a FASTA, since\n"
      "the Config has its own. A `config` and a `preset` that disagree are\n"
      "a ValueError.\n"
      "\n"
      "Construction raises FileNotFoundError for a missing path, ValueError\n"
      "for an unknown preset or a rejected configuration, and RuntimeError\n"
      "otherwise. It never returns a broken object, so mappy's\n"
      "`if not a: raise` idiom is unnecessary: __bool__ is always True.\n"
      "\n"
      "map(), map_batch() and paf() may be called concurrently on one\n"
      "Aligner from several Python threads, and release the GIL for the\n"
      "whole C++ call. reconfigure() may also run meanwhile; a map() that\n"
      "overlaps it uses either the old or the new configuration.")
      .def(
          "__init__",
          [](PyAligner* self, nb::object target, nb::object preset_obj,
             int threads, nb::object config_obj) {
            const std::optional<std::string> preset =
                optional_string(preset_obj, "preset");
            std::optional<flashalign::Config> config;
            if (!config_obj.is_none()) {
              if (!nb::isinstance<flashalign::Config>(config_obj)) {
                throw nb::type_error(
                    "flashalign: config must be a flashalign.Config or None");
              }
              config = nb::cast<flashalign::Config>(config_obj);
            }
            if (preset && config && config->preset != *preset) {
              throw std::invalid_argument(
                  "flashalign: preset='" + *preset + "' and config.preset='" +
                  config->preset +
                  "' disagree; pass one or the other, not both");
            }
            // Seeding for a FASTA target: the explicit preset, else the Config's, else lr.
            const std::string build_preset = preset   ? *preset
                                             : config ? config->preset
                                                      : std::string("lr");

            std::optional<flashalign::Index> index;
            if (nb::isinstance<flashalign::Index>(target)) {
              index = nb::cast<flashalign::Index>(target);
            } else {
              std::filesystem::path path;
              try {
                path = nb::cast<std::filesystem::path>(target);
              } catch (const nb::cast_error&) {
                throw nb::type_error(
                    "flashalign: target must be a path (str or os.PathLike) "
                    "or a flashalign.Index");
              }
              const std::string spelling = path.string();
              if (!std::filesystem::exists(path)) raise_not_found(spelling);
              // By content, as the CLI: load a .faix, index anything else as FASTA.
              const bool is_faix = flashalign::Index::is_index_file(spelling);
              const std::pair<int, int> seeding =
                  is_faix ? std::pair<int, int>{0, 0}
                          : flashalign::preset_seeding(build_preset);
              nb::gil_scoped_release release;
              index = is_faix ? flashalign::Index::load(spelling)
                              : flashalign::Index::build_from_fasta(
                                    spelling, seeding.first, seeding.second,
                                    threads);
            }

            if (!config) {
              // As the CLI: the explicit preset, else the index's, else lr.
              const std::string recorded = index->preset();
              flashalign::Config derived;
              derived.preset = preset            ? *preset
                               : recorded.empty() ? std::string("lr")
                                                  : recorded;
              derived.threads = threads;
              config = derived;
            }
            std::unique_ptr<flashalign::Aligner> aligner;
            {
              nb::gil_scoped_release release;
              aligner =
                  std::make_unique<flashalign::Aligner>(*index, *config);
            }
            // Construct last, so a refused configuration leaves no half-built object.
            new (self) PyAligner{std::move(*index), std::move(aligner)};
          },
          nb::arg("target"), nb::arg("preset") = nb::none(), nb::kw_only(),
          nb::arg("threads") = 0, nb::arg("config") = nb::none(),
          nb::sig("def __init__(self, target: str | os.PathLike | Index, "
                  "preset: str | None = None, *, threads: int = 0, "
                  "config: Config | None = None) -> None"),
          "Aligner(target, preset=None, *, threads=0, config=None)")
      .def(
          "map",
          [](const PyAligner& self, const std::string& seq, nb::handle cs,
             nb::handle md) {
            const flashalign::MapRequest request = make_request(cs, md);
            flashalign::Alignment record;
            {
              nb::gil_scoped_release release;
              record = self.aligner->map(seq, request);
            }
            return flatten_hits(std::move(record));
          },
          nb::arg("seq"), nb::kw_only(), nb::arg("cs") = nb::none(),
          nb::arg("MD") = nb::none(),
          nb::sig("def map(self, seq: str, *, cs: bool | str | None = None, "
                  "MD: bool | None = None) -> list[Hit]"),
          "Align one read; the CLI's record set, in the CLI's order.\n"
          "\n"
          "Returns [primary, *supplementary segments] as separate Hits, or\n"
          "[] for an unmapped read. Rival hypotheses stay on the primary's\n"
          "`.secondary`; the CLI writes them as rows only under\n"
          "--secondary yes.\n"
          "\n"
          "cs and MD request minimap2's difference strings for this call.\n"
          "None (the default) keeps Config.cs and Config.emit_md, both off\n"
          "unless set; True (or 'short' / 'long') and MD=True turn them on,\n"
          "False turns them off. A call whose request differs from the one\n"
          "in force reconfigures the aligner first, so a loop that\n"
          "alternates cs=True with plain calls reconfigures every time; to\n"
          "get cs or MD on every read, set them in the Config.\n"
          "\n"
          "The GIL is released for the whole C++ call.")
      .def(
          "map_batch",
          [](const PyAligner& self, const std::vector<std::string>& seqs,
             nb::handle cs, nb::handle md) {
            const flashalign::MapRequest request = make_request(cs, md);
            std::vector<flashalign::Alignment> records;
            {
              nb::gil_scoped_release release;
              records = self.aligner->map_batch(seqs, request);
            }
            std::vector<std::vector<Hit>> out;
            out.reserve(records.size());
            for (flashalign::Alignment& record : records)
              out.push_back(flatten_hits(std::move(record)));
            return out;
          },
          nb::arg("seqs"), nb::kw_only(), nb::arg("cs") = nb::none(),
          nb::arg("MD") = nb::none(),
          nb::sig("def map_batch(self, seqs: Sequence[str], *, "
                  "cs: bool | str | None = None, MD: bool | None = None) "
                  "-> list[list[Hit]]"),
          "Align many reads; one inner list per read, in input order.\n"
          "\n"
          "Each inner list is what map() would return for that read; cs and\n"
          "MD mean the same as for map() and apply to the whole batch. The\n"
          "batch runs on `threads` workers and releases the GIL once, so\n"
          "this is the fast path for many reads.")
      .def(
          "paf",
          [](const PyAligner& self, const std::string& name,
             const std::string& seq, bool cigar) {
            nb::gil_scoped_release release;
            return self.aligner->paf(name, seq, cigar);
          },
          nb::arg("name"), nb::arg("seq"), nb::kw_only(),
          nb::arg("cigar") = false,
          "The CLI's PAF record(s) for one read, as one string.\n"
          "\n"
          "Complete PAF lines, newline-terminated, with the query name and\n"
          "length; cigar=True adds the cg:Z tag (the CLI's -c). Under the\n"
          "default configuration, cigar=True gives what\n"
          "`flashalign align -c` writes for the read. The CLI's plain\n"
          "PAF maps without base-level alignment;\n"
          "Config(full_read_cigar=False) with cigar=False gives that output.\n"
          "No secondary rows are written, as under the CLI's default. An\n"
          "unmapped read gives ''.\n"
          "\n"
          "The name seeds the tie-break between equally good loci, as a\n"
          "query name does in minimap2; map() and map_batch() have no name\n"
          "and seed it from the read length alone, as mappy's map() does.")
      .def(
          "seq",
          [](const PyAligner& self, const std::string& name, std::int64_t start,
             std::int64_t end) -> std::optional<std::string> {
            nb::gil_scoped_release release;
            try {
              return self.index.sequence(name, start, end);
            } catch (const std::out_of_range&) {
              return std::nullopt;
            }
          },
          nb::arg("name"), nb::arg("start") = 0, nb::arg("end") = -1,
          "A slice of the reference, or None for a contig the index lacks.\n"
          "\n"
          "Index.seq() by another name; the same clamping, the same uppercase\n"
          "ACGT/N alphabet, and the same first-call cost (the whole embedded\n"
          "reference is unpacked once and cached).")
      .def(
          "reconfigure",
          [](PyAligner& self, const flashalign::Config& config) {
            nb::gil_scoped_release release;
            self.aligner->reconfigure(config);
          },
          nb::arg("config"),
          "Install a new configuration on this aligner.\n"
          "\n"
          "k and s stay the index's. Safe to call while other threads map;\n"
          "a map() that overlaps it uses either the old or the new\n"
          "configuration. ValueError for a rejected configuration.")
      .def_prop_ro(
          "seq_names",
          [](const PyAligner& self) { return self.aligner->reference_names(); },
          "Reference contig names, in index order.")
      .def_prop_ro(
          "seq_lengths",
          [](const PyAligner& self) {
            return self.aligner->reference_lengths();
          },
          "Reference contig lengths, in seq_names order.")
      .def_prop_ro(
          "n_seq",
          [](const PyAligner& self) {
            return self.aligner->reference_names().size();
          },
          "Number of reference contigs.")
      .def_prop_ro(
          "k", [](const PyAligner& self) { return self.index.k(); },
          "Seed k-mer length, which is always the index's.")
      .def_prop_ro(
          "s", [](const PyAligner& self) { return self.index.syncmer_s(); },
          "Closed-syncmer s, which is always the index's.")
      .def_prop_ro(
          "preset",
          [](const PyAligner& self) { return self.aligner->config().preset; },
          "The preset this aligner resolved its configuration from.")
      .def_prop_ro(
          "config",
          [](const PyAligner& self) { return self.aligner->config(); },
          "A copy of the resolved configuration from the constructor or\n"
          "reconfigure(); a per-call cs/MD request never changes it.")
      .def_prop_ro(
          "index", [](const PyAligner& self) { return self.index; },
          "The Index this aligner maps against.\n"
          "\n"
          "Reference-counted: a second Aligner built on it shares the same\n"
          "copy.")
      .def(
          "__bool__", [](const PyAligner&) { return true; },
          "Always True: a failed construction raises instead of returning a\n"
          "falsy object.")
      .def(
          "__repr__",
          [](const PyAligner& self) {
        const flashalign::Config config = self.aligner->config();
        return "<flashalign.Aligner preset=" + config.preset +
               " k=" + std::to_string(self.index.k()) +
               " s=" + std::to_string(self.index.syncmer_s()) + " n_seq=" +
               std::to_string(self.aligner->reference_names().size()) +
               " threads=" + std::to_string(config.threads) + ">";
          },
          "The preset, the seeding in force, the contig count and the thread "
          "count.");
}
