#pragma once

// UsageError means the command line is wrong (main() exits 2 with a hint);
// any other exception means the run failed (main() exits 1).

#include <stdexcept>
#include <string>

namespace fa::cpu::cli {

class UsageError : public std::invalid_argument {
public:
  explicit UsageError(const std::string& what) : std::invalid_argument(what) {}
  explicit UsageError(const char* what) : std::invalid_argument(what) {}
};

} // namespace fa::cpu::cli
