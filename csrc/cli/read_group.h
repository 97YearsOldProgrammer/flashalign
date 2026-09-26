// -R/--rg validation and unescaping, as minimap2's sam_write_rg_line and
// mm_escape, except that an unrecognized escape is an error instead of being
// silently dropped.
#pragma once

#include "cli/errors.h" // UsageError (a malformed -R/--rg is a usage error)

#include <stdexcept>
#include <string>

namespace fa::cpu::cli {

struct ReadGroup {
  std::string line; // unescaped @RG header line (real tab characters)
  std::string id;   // the ID: field value, for RG:Z and duplicate detection
};

// Throws UsageError with a diagnostic naming the defect. An empty
// spec is rejected like any other non-@RG line.
inline ReadGroup parse_read_group(const std::string& spec) {
  if (spec.rfind("@RG", 0) != 0) {
    throw UsageError(
        "-R/--rg: the read group line is not started with @RG");
  }
  if (spec.find('\t') != std::string::npos) {
    throw UsageError(
        "-R/--rg: the read group line contained literal <tab> characters -- "
        "replace literal tabs with escaped tabs: \\t");
  }
  ReadGroup group;
  for (std::size_t i = 0; i < spec.size(); ++i) {
    if (spec[i] != '\\') {
      group.line.push_back(spec[i]);
      continue;
    }
    if (i + 1 == spec.size()) {
      throw UsageError(
          "-R/--rg: the read group line ends with a dangling backslash");
    }
    const char escaped = spec[i + 1];
    if (escaped == 't')
      group.line.push_back('\t');
    else if (escaped == '\\')
      group.line.push_back('\\');
    else
      throw UsageError(
          std::string("-R/--rg: unsupported escape '\\") + escaped +
          "' in the read group line; only \\t and \\\\ are recognized");
    ++i;
  }
  const std::size_t id_field = group.line.find("\tID:");
  if (id_field == std::string::npos) {
    throw UsageError("-R/--rg: no ID within the read group line");
  }
  const std::size_t begin = id_field + 4;
  std::size_t end = begin;
  while (end < group.line.size() && group.line[end] != '\t' &&
         group.line[end] != '\n') {
    ++end;
  }
  group.id = group.line.substr(begin, end - begin);
  if (group.id.size() > 255) {
    throw UsageError(
        "-R/--rg: @RG:ID is longer than 255 characters");
  }
  if (group.id.empty()) {
    throw UsageError("-R/--rg: @RG:ID is empty");
  }
  return group;
}

// The ID declared by an @RG header line (e.g. a uBAM passthrough line), or an
// empty string when the line declares none.
inline std::string read_group_line_id(const std::string& line) {
  const std::size_t id_field = line.find("\tID:");
  if (id_field == std::string::npos)
    return std::string();
  const std::size_t begin = id_field + 4;
  std::size_t end = begin;
  while (end < line.size() && line[end] != '\t' && line[end] != '\n')
    ++end;
  return line.substr(begin, end - begin);
}

} // namespace fa::cpu::cli
