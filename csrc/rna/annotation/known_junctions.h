#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace fa::cpu::lr::rna {

struct KnownIntron {
  std::uint32_t start = 0;
  std::uint32_t end = 0;
  std::uint8_t strand_mask = 0; // 1 plus, 2 minus
};

struct KnownJunctionLoadCounters {
  std::uint64_t rows = 0;
  std::uint64_t introns = 0;
  std::uint64_t duplicates_merged = 0;
  std::uint64_t single_block_rows = 0;
};

struct KnownJunctionContigView {
  const KnownIntron* begin = nullptr;
  const KnownIntron* end = nullptr;
  bool empty() const noexcept { return begin == end; }
};

class KnownJunctionStore {
public:
  static KnownJunctionStore
  load_bed(const std::string& path, const std::vector<std::string>& names,
           const std::vector<std::vector<std::uint8_t>>& references);
  KnownJunctionContigView contig(std::size_t id) const noexcept;
  const KnownJunctionLoadCounters& counters() const noexcept {
    return counters_;
  }

private:
  std::vector<std::vector<KnownIntron>> contigs_;
  KnownJunctionLoadCounters counters_;
};

std::uint64_t
materialize_known_junction_mask(KnownJunctionContigView junctions,
                                int reference_begin, int reference_end,
                                std::vector<std::uint8_t>& output);

} // namespace fa::cpu::lr::rna
