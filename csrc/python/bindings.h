// Shared declarations for the nanobind module. The bindings use only the public API
// (<flashalign/*.hpp>) and only change representation: CIGARs become mappy's (length, op)
// tuples and each record tree is flattened into the CLI's record order.
#pragma once

#include <nanobind/nanobind.h>

#include <flashalign/alignment.hpp>
#include <flashalign/config.hpp>

#include <string>
#include <vector>

// One record as Aligner.map() returns it: a node of the Alignment tree, plus whether it is a
// supplementary segment of the record before it.
struct Hit {
  flashalign::Alignment record;
  bool supplementary = false;
};

// The primary, then its supplementary segments; empty for an unmapped read. Secondary
// hypotheses stay on the primary Hit.
std::vector<Hit> flatten_hits(flashalign::Alignment record);

// False/"none" -> None, True/"short" -> Short, "long" -> Long; anything else throws
// (ValueError or TypeError).
flashalign::CsMode parse_cs_mode(nanobind::handle value);

// "none", "short" or "long".
const char* cs_mode_name(flashalign::CsMode mode);

void bind_alignment(nanobind::module_& module);
void bind_aligner(nanobind::module_& module);
void bind_config(nanobind::module_& module);
void bind_fastx(nanobind::module_& module);
void bind_index(nanobind::module_& module);
