#include "bindings.h"

#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <optional>
#include <string_view>
#include <sstream>
#include <stdexcept>
#include <string>

namespace nb = nanobind;

const char* cs_mode_name(flashalign::CsMode mode) {
  switch (mode) {
  case flashalign::CsMode::Short:
    return "short";
  case flashalign::CsMode::Long:
    return "long";
  case flashalign::CsMode::None:
    break;
  }
  return "none";
}

flashalign::CsMode parse_cs_mode(nb::handle value) {
  if (value.is_none())
    return flashalign::CsMode::None;
  if (nb::isinstance<nb::bool_>(value)) {
    return nb::cast<bool>(value) ? flashalign::CsMode::Short
                                 : flashalign::CsMode::None;
  }
  if (nb::isinstance<nb::str>(value)) {
    const std::string form = nb::cast<std::string>(value);
    if (form == "none")
      return flashalign::CsMode::None;
    if (form == "short")
      return flashalign::CsMode::Short;
    if (form == "long")
      return flashalign::CsMode::Long;
    throw std::invalid_argument(
        "flashalign: cs expects False, True, 'none', 'short' or 'long' (got '" +
        form + "')");
  }
  throw nb::type_error(
      "flashalign: cs expects False, True, 'none', 'short' or 'long'");
}

namespace {

void append(std::ostringstream& out, const char* name, int value) {
  out << name << '=' << value;
}
void append(std::ostringstream& out, const char* name, double value) {
  out << name << '=' << value;
}
void append(std::ostringstream& out, const char* name, bool value) {
  out << name << '=' << (value ? "True" : "False");
}
void append(std::ostringstream& out, const char* name,
            const std::string& value) {
  out << name << "='" << value << '\'';
}
void append(std::ostringstream& out, const char* name,
            const std::optional<int>& value) {
  out << name << '=';
  if (value)
    out << *value;
  else
    out << "None";
}

} // namespace

// Every Config field except cs, with its docstring. Expanded for the properties, the keyword
// constructor and __repr__.
#define FA_CONFIG_FIELDS(X)                                                    \
  X(preset,                                                                    \
    "Preset name: 'lr', 'lr:hq', 'splice', 'splice:hq' [-x / --preset].")      \
  X(k, "Seed k-mer length [-k]; mapping always uses the index's.")             \
  X(min_support, "Minimum anchor support [--min-support]. -1 = preset-owned.") \
  X(max_query_seeds, "Query seeds kept per strand. -1 = preset-owned.")        \
  X(full_read_cigar,                                                           \
    "Realize base-level CIGARs; False maps only (the CLI's plain -f paf).")    \
  X(syncmer_s, "Closed-syncmer s [-s]; mapping always uses the index's.")      \
  X(syncmer_downsample, "Query seed downsampling stride; 1 = none.")           \
  X(vote_diag_bin_width, "Vote diagonal bin width [--dw]. -1 = preset-owned.") \
  X(long_occ_cap,                                                              \
    "Seed occurrence cap floor; the index may raise it. -1 = preset-owned.")   \
  X(long_primary_occ_cap,                                                      \
    "Cap of the fixed occurrence policy that --max-vote-occ selects;\n"        \
    "unused under the presets' policy. -1 = preset-owned.")                    \
  X(threads, "Worker threads [-t]; 0 = one per available CPU.")                \
  X(dp_match, "DP matching score [-A]. -1 = preset-owned.")                    \
  X(dp_mismatch, "DP mismatch penalty [-B]. -1 = preset-owned.")               \
  X(dp_score_n, "DP ambiguous-base score [--score-N]. -1 = preset-owned.")     \
  X(dp_gap_open1, "First gap-open penalty [-O]. -1 = preset-owned.")           \
  X(dp_gap_extend1, "First gap-extension penalty [-E]. -1 = preset-owned.")    \
  X(dp_gap_open2, "Second gap-open penalty [-O]. -1 = preset-owned.")          \
  X(dp_gap_extend2, "Second gap-extension penalty [-E]. -1 = preset-owned.")   \
  X(dp_zdrop, "Z-drop score [-z]. -1 = preset-owned.")                         \
  X(dp_zdrop_inv,                                                              \
    "Inversion Z-drop, the second -z value (row dp_inversion_zdrop).\n"        \
    "-1 = dp_zdrop; on a copy of Aligner.config, set it with dp_zdrop.")       \
  X(dp_end_bonus, "Alignment end bonus [--end-bonus]. -2 = preset-owned.")     \
  X(rna_min_intron, "Minimum intron length [--min-intron]; splice only.")      \
  X(rna_max_intron, "Maximum intron length [-G]; splice only.")                \
  X(rna_strand_mode,                                                           \
    "Transcript strand [-u]: -1 unset, 0 auto, 1 forward, 2 reverse, 3 none.") \
  X(tile_score_hit,                                                            \
    "Query-tile supported reward [--tile-score]. None = preset-owned.")        \
  X(tile_score_block,                                                          \
    "Query-tile block-open cost [--tile-score]. None = preset-owned.")         \
  X(tile_score_null,                                                           \
    "Query-tile null cost [--tile-score]. None = preset-owned.")               \
  X(tile_score_miss,                                                           \
    "Query-tile unsupported-ownership cost [--tile-score]. None = "            \
    "preset-owned.")                                                           \
  X(rna_junction_bed,                                                          \
    "Known-junction BED6/BED12 path [--junc-bed]; splice only.")               \
  X(rna_junction_bonus,                                                        \
    "Known-junction endpoint bonus [--junc-bonus]. -1 = preset-owned.")        \
  X(rna_rival_pri_ratio,                                                       \
    "Rival-to-primary chain score ratio [-p]. -1 = preset-owned.")             \
  X(rna_max_loci,                                                              \
    "Candidate loci retained per read [-N + 1]; at most one fewer are\n"       \
    "realized. -1 = preset-owned.")                                            \
  X(dp_min_score, "Minimum DP alignment score [-S]. -1 = preset-owned.")       \
  X(dp_bw,                                                                     \
    "Chaining and alignment bandwidth [-r]. -1 = preset-owned; a splice\n"     \
    "preset rejects it.")                                                      \
  X(dp_bw_long,                                                                \
    "Long-join bandwidth, the second -r value. -1 = preset-owned; a\n"         \
    "splice preset rejects it.")                                               \
  X(emit_md, "Emit the MD:Z difference string [--MD].")                        \
  X(emit_eqx, "Write =/X CIGAR operators instead of M [--eqx].")

void bind_config(nb::module_& module) {
  nb::class_<flashalign::Config> config(
      module, "Config",
      "Mapping options, named as `flashalign align --show-config` prints\n"
      "them; a dot there (tile_score.hit) is an underscore here.\n"
      "\n"
      "-1 (or None, or \"\") means \"take the preset's value\"; anything else\n"
      "is an explicit setting, which may be rejected with a ValueError.\n"
      "\n"
      "    fa.Config(preset='lr:hq', threads=8, dp_min_score=50)\n"
      "\n"
      "Aligner.config returns the resolved configuration, with every field\n"
      "set to the value in use, or -1 / None where the preset decides it.");

  config.def(
      "__init__",
      [](flashalign::Config* self, const std::string& preset,
         nb::kwargs fields) {
        new (self) flashalign::Config();
        self->preset = preset;
        for (auto item : fields) {
          const std::string name = nb::cast<std::string>(item.first);
          nb::handle value = item.second;
          if (name == "cs") {
            self->cs = parse_cs_mode(value);
            continue;
          }
#define FA_CONFIG_KWARG(field, doc)                                            \
  if (name == #field) {                                                        \
    self->field = nb::cast<decltype(self->field)>(value);                      \
    continue;                                                                  \
  }
          FA_CONFIG_FIELDS(FA_CONFIG_KWARG)
#undef FA_CONFIG_KWARG
          throw std::invalid_argument(
              "flashalign: Config has no field '" + name +
              "' (see help(flashalign.Config) for the field names)");
        }
      },
      nb::arg("preset") = "lr", nb::arg("fields"),
      "Config(preset='lr', **fields): a Config with those fields set.\n"
      "\n"
      "Every keyword must be one of the fields below; an unknown name is a\n"
      "ValueError.");

#define FA_CONFIG_PROP(field, doc)                                             \
  config.def_rw(#field, &flashalign::Config::field, doc);
  FA_CONFIG_FIELDS(FA_CONFIG_PROP)
#undef FA_CONFIG_PROP

  config.def_prop_rw(
      "cs",
      [](const flashalign::Config& self) {
        return std::string(cs_mode_name(self.cs));
      },
      [](flashalign::Config& self, nb::handle value) {
        self.cs = parse_cs_mode(value);
      },
      "cs:Z difference string form [--cs]: 'none', 'short' or 'long'.\n"
      "Accepts False/True for 'none'/'short'.");

  config.def(
      "__repr__",
      [](const flashalign::Config& self) {
        const flashalign::Config base;
        std::ostringstream out;
        out << "flashalign.Config(";
        append(out, "preset", self.preset);
#define FA_CONFIG_REPR(field, doc)                                             \
  if (std::string_view(#field) != "preset" && !(self.field == base.field)) {   \
    out << ", ";                                                               \
    append(out, #field, self.field);                                           \
  }
        FA_CONFIG_FIELDS(FA_CONFIG_REPR)
#undef FA_CONFIG_REPR
        if (self.cs != base.cs) {
          out << ", cs='" << cs_mode_name(self.cs) << '\'';
        }
        out << ')';
        return out.str();
      },
      "The preset, plus every field that differs from the struct default.");
}
