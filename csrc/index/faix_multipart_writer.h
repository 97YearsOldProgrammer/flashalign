// Writer for the multi-part container of `flashalign index -I` (format.h). One
// transactional file in layout order: the header and global tables, then each part's image
// as the caller builds it (so only one part's index is held at a time), and finally the
// part table patched in at the front. Nothing appears at `path` until finish() publishes; a
// failure leaves the destination as it was.
#pragma once

#include "faix.h"
#include "faix_publish.h"

#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace fa {
namespace cpu {

class FaixMultipartWriter {
public:
  // What the container declares up front. The parts' image fields are filled in as images
  // are appended; everything else is checked against each part as it arrives.
  struct Plan {
    int k = 0;
    int syncmer_s = 0;
    bool has_reference = true; // every part embeds its 4-bit reference
    std::string preset;        // recorded on the container as on the parts
    std::vector<uint64_t> chr_offsets; // global prefix sums, chrom_count + 1
    std::vector<std::string> names;    // global, in contig-id order
    std::vector<FaixPartEntry> parts;  // first_contig / contig_count set
  };

  FaixMultipartWriter(const std::string& path, Plan plan)
      : plan_(std::move(plan)), out_(path) {}

  bool opened() const noexcept { return out_.opened(); }

  // The header, the global offsets, a zeroed part table and the global names. finish()
  // patches the part table.
  bool begin() {
    if (!out_.opened() || begun_ || plan_.parts.empty() ||
        plan_.names.empty() ||
        plan_.chr_offsets.size() != plan_.names.size() + 1 ||
        plan_.chr_offsets.front() != 0) {
      return false;
    }
    // The parts must tile the contig range in order.
    uint64_t next = 0;
    for (FaixPartEntry& e : plan_.parts) {
      if (e.contig_count == 0 || e.first_contig != next ||
          e.contig_count > plan_.names.size() - next) {
        return false;
      }
      e.image_offset = 0;
      e.image_bytes = 0;
      next += e.contig_count;
    }
    if (next != plan_.names.size())
      return false;

    FaixHeader header;
    header.flags = kFaixFlagMultipart |
                    (plan_.has_reference ? kFaixFlagHasReference : 0u);
    header.k = static_cast<uint32_t>(plan_.k);
    header.syncmer_s = static_cast<uint32_t>(plan_.syncmer_s);
    header.contig_count = plan_.names.size();
    header.total_bp = plan_.chr_offsets.back();
    header.part_count = plan_.parts.size();
    for (const std::string& name : plan_.names) {
      if (name.size() >
          static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
        return false;
      }
      header.name_bytes +=
          sizeof(uint32_t) + static_cast<uint64_t>(name.size());
    }
    if (!faix_header_set_preset(header, plan_.preset))
      return false;

    char head[kFaixHeaderBytes] = {};
    faix_header_to_bytes(header, head);
    if (!out_.write_bytes(head, kFaixHeaderBytes))
      return false;
    if (!out_.write_bytes(plan_.chr_offsets.data(),
                          plan_.chr_offsets.size() * sizeof(uint64_t))) {
      return false;
    }
    part_table_offset_ = out_.written();
    const std::vector<FaixPartEntry> zeroed(plan_.parts.size());
    if (!out_.write_bytes(zeroed.data(),
                          zeroed.size() * sizeof(FaixPartEntry))) {
      return false;
    }
    for (const std::string& name : plan_.names) {
      const uint32_t len = static_cast<uint32_t>(name.size());
      if (!out_.write_bytes(&len, sizeof(len)) ||
          !out_.write_bytes(name.data(), name.size())) {
        return false;
      }
    }
    begun_ = true;
    return true;
  }

  // The next part in order: the index of exactly the contigs the plan gave it, with the
  // plan's names, lengths, seeding, reference policy and preset.
  bool append_part(const FaixIndex& part) {
    if (!begun_ || next_part_ >= plan_.parts.size())
      return false;
    FaixPartEntry& entry = plan_.parts[next_part_];
    const size_t first = static_cast<size_t>(entry.first_contig);
    const size_t count = static_cast<size_t>(entry.contig_count);
    const uint64_t part_bp =
        plan_.chr_offsets[first + count] - plan_.chr_offsets[first];
    if (part.empty() || part.k() != plan_.k ||
        part.build_syncmer_s() != plan_.syncmer_s ||
        part.has_reference_payload() != plan_.has_reference ||
        part.chrom_count() != entry.contig_count ||
        part.total_bp() != part_bp || part.build_preset() != plan_.preset ||
        part.chromosome_names().size() != count) {
      return false;
    }
    const uint64_t* offsets = part.chrom_offsets_data();
    if (offsets == nullptr)
      return false;
    for (size_t j = 0; j < count; ++j) {
      if (part.chromosome_names()[j] != plan_.names[first + j] ||
          offsets[j + 1] - offsets[j] !=
              plan_.chr_offsets[first + j + 1] - plan_.chr_offsets[first + j]) {
        return false;
      }
    }
    const uint64_t at = out_.written();
    uint64_t expected = 0;
    if (!part.write_image(out_, expected) || out_.written() - at != expected) {
      return false;
    }
    entry.image_offset = at;
    entry.image_bytes = expected;
    ++next_part_;
    return true;
  }

  // Every part appended: patch the part table, publish.
  bool finish() {
    if (!begun_ || next_part_ != plan_.parts.size())
      return false;
    if (!out_.patch_bytes(part_table_offset_, plan_.parts.data(),
                          plan_.parts.size() * sizeof(FaixPartEntry))) {
      return false;
    }
    return out_.publish(out_.written());
  }

private:
  Plan plan_;
  FaixTransactionalFile out_;
  uint64_t part_table_offset_ = 0;
  size_t next_part_ = 0;
  bool begun_ = false;
};

} // namespace cpu
} // namespace fa
