#include "config.h"

#include <cstring>

namespace fa {
namespace cpu {
namespace lr {
namespace rna {

StrandMode parse_strand_mode(const char* text, StrandMode fallback) {
  if (!text || !*text) return fallback;
  if (std::strcmp(text, "unknown") == 0 || std::strcmp(text, "auto") == 0)
    return StrandMode::Unknown;
  if (std::strcmp(text, "forward") == 0 ||
      std::strcmp(text, "fwd") == 0 || std::strcmp(text, "f") == 0 ||
      std::strcmp(text, "+") == 0)
    return StrandMode::Forward;
  if (std::strcmp(text, "reverse") == 0 ||
      std::strcmp(text, "rev") == 0 || std::strcmp(text, "r") == 0 ||
      std::strcmp(text, "-") == 0)
    return StrandMode::Reverse;
  // minimap2 -u n: score no splice motif at all.
  if (std::strcmp(text, "none") == 0 || std::strcmp(text, "n") == 0)
    return StrandMode::None;
  return fallback;
}

}  // namespace rna
}  // namespace lr
}  // namespace cpu
}  // namespace fa
