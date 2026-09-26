#include "known_junctions.h"

#include <zlib.h>
#include <algorithm>
#include <charconv>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace fa::cpu::lr::rna {
namespace {
void strip_carriage_return(std::string& line) {
  if (!line.empty() && line.back() == '\r')
    line.pop_back();
}

std::vector<std::string_view> fields(std::string_view line,
                                     char delimiter = '\t') {
  std::vector<std::string_view> out;
  while (true) {
    const auto pos = line.find(delimiter);
    out.push_back(line.substr(0, pos));
    if (pos == std::string_view::npos)
      break;
    line.remove_prefix(pos + 1);
  }
  return out;
}
std::uint32_t number(std::string_view text, const std::string& path,
                     std::uint64_t line) {
  std::uint64_t value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (text.empty() || result.ec != std::errc{} ||
      result.ptr != text.data() + text.size() ||
      value > std::numeric_limits<std::uint32_t>::max())
    throw std::invalid_argument(path + ":" + std::to_string(line) +
                                ": invalid BED coordinate");
  return static_cast<std::uint32_t>(value);
}
std::vector<std::string> read_lines(const std::string& path) {
  std::vector<std::string> out;
  if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") {
    gzFile file = gzopen(path.c_str(), "rb");
    if (!file)
      throw std::invalid_argument("cannot open RNA junction BED: " + path);
    std::string line;
    char buffer[65536];
    while (gzgets(file, buffer, sizeof buffer)) {
      line += buffer;
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
        strip_carriage_return(line);
        out.push_back(line);
        line.clear();
      }
    }
    int gz_error = Z_OK;
    (void)gzerror(file, &gz_error);
    const int error = gzclose(file);
    if (!line.empty()) {
      strip_carriage_return(line);
      out.push_back(line);
    }
    if ((gz_error != Z_OK && gz_error != Z_STREAM_END) || error != Z_OK)
      throw std::invalid_argument("cannot decompress RNA junction BED: " +
                                  path);
  } else {
    std::ifstream file(path);
    if (!file)
      throw std::invalid_argument("cannot open RNA junction BED: " + path);
    std::string line;
    while (std::getline(file, line)) {
      strip_carriage_return(line);
      out.push_back(line);
    }
    if (file.bad())
      throw std::invalid_argument("cannot read RNA junction BED: " + path);
  }
  return out;
}
} // namespace

KnownJunctionStore KnownJunctionStore::load_bed(
    const std::string& path, const std::vector<std::string>& names,
    const std::vector<std::vector<std::uint8_t>>& references) {
  if (names.size() != references.size())
    throw std::invalid_argument(
        "RNA junction BED: invalid reference dictionary");
  std::unordered_map<std::string, std::size_t> ids;
  for (std::size_t i = 0; i < names.size(); ++i)
    if (!ids.emplace(names[i], i).second)
      throw std::invalid_argument(
          "RNA junction BED: duplicate reference contig: " + names[i]);
  KnownJunctionStore store;
  store.contigs_.resize(names.size());
  std::uint64_t line_number = 0;
  for (const std::string& line : read_lines(path)) {
    ++line_number;
    if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 ||
        line.rfind("browser", 0) == 0)
      continue;
    ++store.counters_.rows;
    const auto f = fields(line);
    if (f.size() != 6 && f.size() != 12)
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": expected BED6 or BED12");
    const auto found = ids.find(std::string(f[0]));
    if (found == ids.end())
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": unknown contig: " + std::string(f[0]));
    if (f[5] != "+" && f[5] != "-")
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": strand must be + or -");
    const std::uint8_t strand = f[5] == "+" ? 1 : 2;
    const auto chrom_start = number(f[1], path, line_number);
    const auto chrom_end = number(f[2], path, line_number);
    if (chrom_end <= chrom_start ||
        chrom_end > references[found->second].size())
      throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                  ": BED interval outside reference");
    auto add = [&](std::uint64_t start, std::uint64_t end) {
      if (end <= start || end > references[found->second].size() ||
          end > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                    ": intron outside reference");
      store.contigs_[found->second].push_back(
          {static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(end),
           strand});
      ++store.counters_.introns;
    };
    if (f.size() == 6)
      add(chrom_start, chrom_end);
    else {
      const auto count = number(f[9], path, line_number);
      const auto sizes = fields(f[10], ',');
      const auto starts = fields(f[11], ',');
      auto usable = [](const std::vector<std::string_view>& v) {
        return !v.empty() && v.back().empty() ? v.size() - 1 : v.size();
      };
      if (count == 0 || usable(sizes) != count || usable(starts) != count)
        throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                    ": malformed BED12 blocks");
      if (count == 1)
        ++store.counters_.single_block_rows;
      std::uint64_t previous_end = 0;
      for (std::uint32_t i = 0; i < count; ++i) {
        const auto block_start = number(starts[i], path, line_number);
        const auto block_size = number(sizes[i], path, line_number);
        const std::uint64_t exon_start =
            std::uint64_t(chrom_start) + block_start;
        const std::uint64_t exon_end = exon_start + block_size;
        if (block_size == 0 || exon_end > chrom_end ||
            (i && exon_start < previous_end))
          throw std::invalid_argument(path + ":" + std::to_string(line_number) +
                                      ": invalid BED12 block geometry");
        if (i)
          add(previous_end, exon_start);
        previous_end = exon_end;
      }
    }
  }
  for (auto& contig : store.contigs_) {
    std::sort(contig.begin(), contig.end(), [](const auto& a, const auto& b) {
      return a.start < b.start || (a.start == b.start && a.end < b.end);
    });
    std::vector<KnownIntron> merged;
    for (const auto& intron : contig) {
      if (!merged.empty() && merged.back().start == intron.start &&
          merged.back().end == intron.end) {
        merged.back().strand_mask |= intron.strand_mask;
        ++store.counters_.duplicates_merged;
      } else
        merged.push_back(intron);
    }
    contig = std::move(merged);
  }
  return store;
}
KnownJunctionContigView
KnownJunctionStore::contig(std::size_t id) const noexcept {
  if (id >= contigs_.size())
    return {};
  return {contigs_[id].data(), contigs_[id].data() + contigs_[id].size()};
}
std::uint64_t materialize_known_junction_mask(KnownJunctionContigView view,
                                              int begin, int end,
                                              std::vector<std::uint8_t>& out) {
  if (begin < 0 || end <= begin) {
    out.clear();
    return 0;
  }
  out.assign(static_cast<std::size_t>(end - begin), 0);
  if (view.empty())
    return 0;
  auto it = std::lower_bound(view.begin, view.end, begin,
                             [](const KnownIntron& i, int p) {
                               return i.start < static_cast<std::uint32_t>(p);
                             });
  for (; it != view.end && it->start < static_cast<std::uint32_t>(end); ++it) {
    if (it->end > static_cast<std::uint32_t>(end))
      continue;
    auto& donor = out[it->start - begin];
    auto& acceptor = out[it->end - 1 - begin];
    if (it->strand_mask & 1) {
      donor |= 1;
      acceptor |= 2;
    }
    if (it->strand_mask & 2) {
      donor |= 8;
      acceptor |= 4;
    }
  }
  return static_cast<std::uint64_t>(std::count_if(
      out.begin(), out.end(), [](std::uint8_t value) { return value != 0; }));
}
} // namespace fa::cpu::lr::rna
