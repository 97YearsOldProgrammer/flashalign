#include "bindings.h"

#include <flashalign/fastx.hpp>

#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>

#include <filesystem>
#include <memory>
#include <string>
#include <utility>

namespace nb = nanobind;

namespace {

// The iterator fastx_read() returns; it keeps the file open while it lives.
struct PyFastxReader {
  std::unique_ptr<flashalign::FastxReader> reader;
  bool read_comment = false;
};

} // namespace

void bind_fastx(nb::module_& module) {
  nb::class_<PyFastxReader>(
      module, "FastxReader",
      "An iterator over the records of one FASTA/FASTQ file.\n"
      "\n"
      "Made by fastx_read(); yields (name, seq, qual) tuples, or\n"
      "(name, seq, qual, comment) with read_comment=True. Reads plain, gzip\n"
      "and bgzf files (multi-member gzip included) with the CLI's reader,\n"
      "and '-' for stdin; not unaligned BAM. Sequences are uppercased.")
      .def(
          "__iter__", [](nb::object self) { return self; },
          nb::sig("def __iter__(self) -> FastxReader"),
          "The reader is its own iterator.")
      .def(
          "__next__",
          [](PyFastxReader& self) {
            flashalign::FastxRecord record;
            bool more = false;
            {
              nb::gil_scoped_release release;
              more = self.reader->next(record);
            }
            if (!more)
              throw nb::stop_iteration();
            nb::object qual =
                record.qual.empty() ? nb::none() : nb::cast(record.qual);
            if (!self.read_comment)
              return nb::make_tuple(record.name, record.seq, qual);
            nb::object comment =
                record.comment.empty() ? nb::none() : nb::cast(record.comment);
            return nb::make_tuple(record.name, record.seq, qual, comment);
          },
          "The next record, or StopIteration at end of file.\n"
          "\n"
          "The GIL is released around each read. `qual` is None for FASTA,\n"
          "and `comment` is None when the header carries none.");

  module.def(
      "fastx_read",
      [](const std::filesystem::path& path, bool read_comment) {
        PyFastxReader reader;
        reader.read_comment = read_comment;
        {
          nb::gil_scoped_release release;
          reader.reader =
              std::make_unique<flashalign::FastxReader>(path.string());
        }
        return reader;
      },
      nb::arg("path"), nb::arg("read_comment") = false,
      "Iterate one FASTA/FASTQ file: (name, seq, qual) per record.\n"
      "\n"
      "    for name, seq, qual in fa.fastx_read('reads.fq.gz'):\n"
      "        ...\n"
      "\n"
      "`qual` is None for FASTA. read_comment=True yields four-tuples\n"
      "(name, seq, qual, comment), with comment None when the header has\n"
      "none. Plain, gzip and bgzf are all read (multi-member gzip included),\n"
      "and path='-' reads stdin; unaligned BAM is not read here. Sequences\n"
      "are uppercased.\n"
      "\n"
      "Unlike mappy, a path that cannot be opened raises FileNotFoundError\n"
      "rather than returning None.");
}
