#include "record_family.h"

#include <algorithm>
#include <tuple>
#include <utility>

namespace fa::cpu::lr {
namespace {

using RecordIterator = std::vector<DnaSegmentRecord>::iterator;

RecordIterator widest_primary(std::vector<DnaSegmentRecord>& records) {
  return std::max_element(
      records.begin(), records.end(),
      [](const DnaSegmentRecord& left, const DnaSegmentRecord& right) {
        return dna_primary_precedence_less(left.alignment, right.alignment);
      });
}

} // namespace

bool dna_primary_precedence_less(const AlignResult& left,
                                 const AlignResult& right) {
  const int left_span = left.query_end - left.query_start;
  const int right_span = right.query_end - right.query_start;
  if (left_span != right_span)
    return left_span < right_span;
  if (left.score != right.score)
    return left.score < right.score;
  return std::tie(left.query_start, left.chromosome, left.pos,
                  left.is_reverse) > std::tie(right.query_start,
                                              right.chromosome, right.pos,
                                              right.is_reverse);
}

DnaRecordFamily
assemble_dna_record_family(std::vector<DnaSegmentRecord> records) {
  DnaRecordFamily family;
  if (records.empty())
    return family;

  std::sort(
      records.begin(), records.end(),
      [](const DnaSegmentRecord& left, const DnaSegmentRecord& right) {
        return std::tie(left.alignment.query_start, left.alignment.query_end,
                        left.candidate, left.alignment.chromosome,
                        left.alignment.pos, left.alignment.is_reverse) <
               std::tie(right.alignment.query_start, right.alignment.query_end,
                        right.candidate, right.alignment.chromosome,
                        right.alignment.pos, right.alignment.is_reverse);
      });
  for (std::size_t index = 1; index < records.size(); ++index) {
    if (records[index].alignment.query_start <
        records[index - 1].alignment.query_end)
      return family;
  }

  const RecordIterator primary = widest_primary(records);

  family.primary = primary->alignment;
  family.primary_candidate = primary->candidate;
  family.supplementary.reserve(records.size() - 1);
  for (const DnaSegmentRecord& record : records) {
    if (&record != &*primary)
      family.supplementary.push_back(record);
  }
  std::sort(family.supplementary.begin(), family.supplementary.end(),
            [](const DnaSegmentRecord& left, const DnaSegmentRecord& right) {
              return std::tie(left.alignment.query_start,
                              left.alignment.query_end,
                              left.alignment.chromosome, left.alignment.pos) <
                     std::tie(right.alignment.query_start,
                              right.alignment.query_end,
                              right.alignment.chromosome,
                              right.alignment.pos);
            });
  family.valid = true;
  return family;
}

} // namespace fa::cpu::lr
