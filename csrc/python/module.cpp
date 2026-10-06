#include "bindings.h"

#include <flashalign/index.hpp>
#include <flashalign/sequence.hpp>
#include <flashalign/version.hpp>

#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>

#include <exception>
#include <stdexcept>
#include <string>

namespace nb = nanobind;

namespace {

// Maps std::out_of_range to KeyError and file-open failures to FileNotFoundError. Registered
// last, so nanobind tries it first; anything rethrown goes to the default translator.
void translate_exception(const std::exception_ptr& error, void*) {
  try {
    std::rethrow_exception(error);
  } catch (const std::out_of_range& e) {
    // An unknown name: a lookup by key, so not IndexError.
    PyErr_SetString(PyExc_KeyError, e.what());
  } catch (const std::runtime_error& e) {
    const std::string message = e.what();
    if (message.rfind("failed to open", 0) == 0 ||
        message.find("No such file") != std::string::npos) {
      PyErr_SetString(PyExc_FileNotFoundError, message.c_str());
      return;
    }
    throw;
  }
}

} // namespace

NB_MODULE(_flashalign, module) {
  module.doc() =
      "FlashAlign's compiled extension. Import `flashalign` instead: the "
      "package re-exports every name defined here.";
  module.attr("__version__") = std::string(flashalign::version);

  bind_config(module);
  bind_alignment(module);
  bind_index(module);
  bind_aligner(module);
  bind_fastx(module);

  module.def("revcomp", &flashalign::revcomp, nb::arg("seq"),
             "Reverse complement, as mappy's revcomp().\n"
             "\n"
             "IUPAC codes are complemented (R<->Y, K<->M, B<->V, D<->H; S, W\n"
             "and N map to themselves), case is kept, U becomes A, and any\n"
             "other character is kept as is, in reversed position.");

  module.def(
      "preset_seeding",
      [](const std::string& preset) {
        return flashalign::preset_seeding(preset);
      },
      nb::arg("preset"),
      "The (k, syncmer_s) a preset builds its index with.\n"
      "\n"
      "    fa.Index.build_from_fasta('ref.fa', *fa.preset_seeding('lr:hq'))\n"
      "\n"
      "Mapping always uses the index's own seeding, so an index built with\n"
      "Index.build_from_fasta for a preset other than 'lr' needs this pair.\n"
      "ValueError for an unknown preset, and for 'ava-ont' and 'ava-hifi',\n"
      "which the binding refuses.");

  nb::register_exception_translator(translate_exception);
}
