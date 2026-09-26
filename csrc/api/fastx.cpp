// The public FASTA/FASTQ reader and revcomp(). The reader drives the io layer's kseq
// directly, rather than io::FastxReader, because it also returns the comment.
#include <flashalign/fastx.hpp>
#include <flashalign/sequence.hpp>

#include "../io/fastx.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace flashalign {

struct FastxReader::Impl {
  std::string path;
  std::unique_ptr<fa::cpu::io::ByteSource> source;
  fa::cpu::io::kseq_t* ks = nullptr;

  explicit Impl(const std::string& value) : path(value) {
    // Throws std::runtime_error when the path cannot be opened.
    source = fa::cpu::io::make_byte_source(path);
    ks = fa::cpu::io::kseq_init(source.get());
    if (ks == nullptr)
      throw std::runtime_error("failed to open FASTA/FASTQ: " + path);
  }
  ~Impl() {
    if (ks != nullptr) fa::cpu::io::kseq_destroy(ks);
    // kseq does not own the ByteSource.
    source.reset();
  }
  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;
};

FastxReader::FastxReader(const std::string& path)
    : impl_(std::make_unique<Impl>(path)) {}
FastxReader::~FastxReader() = default;
FastxReader::FastxReader(FastxReader&&) noexcept = default;
FastxReader& FastxReader::operator=(FastxReader&&) noexcept = default;

bool FastxReader::next(FastxRecord& out) {
  if (impl_ == nullptr || impl_->ks == nullptr) return false;
  const std::int64_t length =
      fa::cpu::io::read_whole_record(impl_->ks, *impl_->source);
  if (length >= 0) {
    const fa::cpu::io::kseq_t* record = impl_->ks;
    out.name.assign(record->name.s, record->name.l);
    if (record->comment.l != 0)
      out.comment.assign(record->comment.s, record->comment.l);
    else
      out.comment.clear();
    out.seq.assign(record->seq.s, record->seq.l);
    if (record->qual.l != 0)
      out.qual.assign(record->qual.s, record->qual.l);
    else
      out.qual.clear();
    // Uppercase, as the io reader does, so soft-masked bases match the index.
    fa::cpu::io::uppercase_inplace(out.seq);
    return true;
  }
  if (length == -3)
    throw std::runtime_error("error reading FASTA/FASTQ: " + impl_->path);
  if (length == -2) {
    throw std::runtime_error(
        "truncated FASTQ or qual/seq length mismatch in " + impl_->path);
  }
  return false;  // length == -1: end of file.
}

namespace {

// IUPAC complement with case kept and U -> A, as mappy; other characters are unchanged.
constexpr char complement_of(char base) {
  switch (base) {
    case 'A': return 'T';
    case 'B': return 'V';
    case 'C': return 'G';
    case 'D': return 'H';
    case 'G': return 'C';
    case 'H': return 'D';
    case 'K': return 'M';
    case 'M': return 'K';
    case 'N': return 'N';
    case 'R': return 'Y';
    case 'S': return 'S';
    case 'T': return 'A';
    case 'U': return 'A';
    case 'V': return 'B';
    case 'W': return 'W';
    case 'Y': return 'R';
    case 'a': return 't';
    case 'b': return 'v';
    case 'c': return 'g';
    case 'd': return 'h';
    case 'g': return 'c';
    case 'h': return 'd';
    case 'k': return 'm';
    case 'm': return 'k';
    case 'n': return 'n';
    case 'r': return 'y';
    case 's': return 's';
    case 't': return 'a';
    case 'u': return 'a';
    case 'v': return 'b';
    case 'w': return 'w';
    case 'y': return 'r';
    default: return base;
  }
}

}  // namespace

std::string revcomp(std::string_view seq) {
  std::string out;
  out.resize(seq.size());
  for (std::size_t i = 0; i < seq.size(); ++i)
    out[seq.size() - 1 - i] = complement_of(seq[i]);
  return out;
}

}  // namespace flashalign
