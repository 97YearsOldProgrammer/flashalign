#include "aligner.h"

#include "../index/format.h"
#include "../io/byte_source.h"

#include <zlib.h>

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fa {
namespace cpu {
namespace api {

ReferenceKind reference_kind(const std::string &path) {
  if (path == "-")
    return ReferenceKind::Sequences;
  const auto refuse = [&](const std::string &why) {
    return std::runtime_error(path + ": " + why);
  };
  // stat() first: opening a named pipe would block, and closing it could end its writer.
  struct stat st;
  if (::stat(path.c_str(), &st) != 0)
    throw refuse(std::strerror(errno));
  if (S_ISDIR(st.st_mode))
    throw refuse("is a directory");
  if (!S_ISREG(st.st_mode))
    return ReferenceKind::Sequences;
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw refuse(std::strerror(errno));

  char magic[8] = {};
  if (::pread(fd, magic, sizeof(magic), 0) == static_cast<ssize_t>(sizeof(magic)) &&
      (std::memcmp(magic, kFaixMagic, sizeof(magic)) == 0 ||
       std::memcmp(magic, kFaixDevelopmentMagic, sizeof(magic)) == 0)) {
    ::close(fd);
    return ReferenceKind::Index;
  }

  // zlib reads plain input as is. gzclose() closes the descriptor.
  gzFile in = gzdopen(fd, "r");
  if (in == nullptr) {
    ::close(fd);
    throw refuse("cannot read the file");
  }
  char buffer[4096];
  char first = 0;
  while (first == 0) {
    const int n = gzread(in, buffer, sizeof(buffer));
    if (n < 0) {
      int code = 0;
      const std::string why = gzerror(in, &code);
      gzclose(in);
      throw refuse("cannot decompress: " + why);
    }
    if (n == 0) {
      // A gzip file cut before any content: say so before the refusal below.
      int code = Z_OK;
      (void)gzerror(in, &code);
      if (code == Z_BUF_ERROR)
        ::fa::cpu::io::warn_truncated_gzip(path);
      break;
    }
    for (int i = 0; i < n; ++i) {
      const char c = buffer[i];
      if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
        first = c;
        break;
      }
    }
  }
  gzclose(in);
  if (first == 0)
    throw refuse("the file is empty");
  if (first != '>' && first != '@')
    throw refuse("not a FASTA or FASTQ file, nor a FlashAlign index");
  return ReferenceKind::Sequences;
}

} // namespace api
} // namespace cpu
} // namespace fa
