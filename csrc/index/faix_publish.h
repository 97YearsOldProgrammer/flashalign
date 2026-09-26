#pragma once

#include "format.h"

#include <atomic>
#include <cerrno>
#include <cstdio>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fa {
namespace cpu {

// Sequential transactional file: bytes go to a private partial in the same directory,
// renamed into place only after an exact-size check and fsync. An existing destination
// survives any failure.
class FaixTransactionalFile final : public FaixByteSink {
public:
  explicit FaixTransactionalFile(const std::string& final_path)
      : final_path_(final_path) {
    static std::atomic<uint64_t> sequence{0};
    for (int attempt = 0; attempt < 1024 && fd_ < 0; ++attempt) {
      partial_path_ =
          final_path + ".partial." +
          std::to_string(static_cast<long long>(getpid())) + "." +
          std::to_string(sequence.fetch_add(1, std::memory_order_relaxed));
      fd_ = open(partial_path_.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_TRUNC,
                 0644);
      if (fd_ < 0 && errno != EEXIST)
        break;
    }
    if (fd_ < 0)
      partial_path_.clear();
  }

  ~FaixTransactionalFile() override {
    if (fd_ >= 0)
      close(fd_);
    if (!published_ && !partial_path_.empty())
      unlink(partial_path_.c_str());
  }

  FaixTransactionalFile(const FaixTransactionalFile&) = delete;
  FaixTransactionalFile& operator=(const FaixTransactionalFile&) = delete;

  bool opened() const noexcept { return fd_ >= 0; }

  bool write_bytes(const void* data, uint64_t size) override {
    if (fd_ < 0 || (size != 0 && data == nullptr))
      return false;
    const char* cursor = static_cast<const char*>(data);
    while (size > 0) {
      const size_t chunk =
          static_cast<size_t>(std::min<uint64_t>(size, 1ULL << 20));
      const ssize_t n = write(fd_, cursor, chunk);
      if (n < 0) {
        if (errno == EINTR)
          continue;
        return false;
      }
      if (n == 0)
        return false;
      cursor += n;
      size -= static_cast<uint64_t>(n);
      written_ += static_cast<uint64_t>(n);
    }
    return true;
  }

  // Overwrites `size` bytes at `offset` within what has already been written; it never
  // extends the file.
  bool patch_bytes(uint64_t offset, const void* data, uint64_t size) {
    if (fd_ < 0 || (size != 0 && data == nullptr) || offset > written_ ||
        size > written_ - offset)
      return false;
    const char* cursor = static_cast<const char*>(data);
    while (size > 0) {
      const size_t chunk =
          static_cast<size_t>(std::min<uint64_t>(size, 1ULL << 20));
      const ssize_t n = pwrite(fd_, cursor, chunk, static_cast<off_t>(offset));
      if (n < 0) {
        if (errno == EINTR)
          continue;
        return false;
      }
      if (n == 0)
        return false;
      cursor += n;
      offset += static_cast<uint64_t>(n);
      size -= static_cast<uint64_t>(n);
    }
    return true;
  }

  uint64_t written() const noexcept { return written_; }

  bool publish(uint64_t expected_bytes) {
    if (fd_ < 0 || written_ != expected_bytes)
      return false;
    struct stat st;
    if (fstat(fd_, &st) != 0 ||
        static_cast<uint64_t>(st.st_size) != expected_bytes ||
        fsync(fd_) != 0) {
      return false;
    }
    if (close(fd_) != 0) {
      fd_ = -1;
      return false;
    }
    fd_ = -1;
    if (rename(partial_path_.c_str(), final_path_.c_str()) != 0)
      return false;
    published_ = true;

    const size_t slash = final_path_.find_last_of('/');
    const std::string parent =
        slash == std::string::npos
            ? std::string(".")
            : final_path_.substr(0, slash == 0 ? 1 : slash);
    const int parent_fd = open(parent.c_str(), O_RDONLY);
    if (parent_fd >= 0) {
      (void)fsync(parent_fd);
      close(parent_fd);
    }
    return true;
  }

private:
  std::string final_path_;
  std::string partial_path_;
  int fd_ = -1;
  uint64_t written_ = 0;
  bool published_ = false;
};

} // namespace cpu
} // namespace fa
