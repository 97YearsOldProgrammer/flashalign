// The DNA lane's catalogue step: each vote peak's factual tile support, then the catalogue over them.
#include "dna.h"

#include "../index/format.h"
#include "../vote/vote.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

/**** placement_family_adapter ****/

namespace fa::cpu::lr {
namespace {

using ::fa::cpu::voting::CandidateCatalogue;
using ::fa::cpu::voting::CandidateId;
using ::fa::cpu::voting::CandidateInput;
using ::fa::cpu::voting::QueryTileMask;

// `packed` is the posting word `global` was decoded from.
bool posting_supports_peak(
    const DnaContext& context,
    const std::vector<std::uint8_t>& forward_query,
    const DnaLongSeedView& seed_view, const VotePeak& peak,
    std::uint64_t global, PackedRefPos packed) {
  // Uses only the index offsets, never reference bases.
  if (context.ref.index == nullptr) return false;
  const std::uint64_t* offsets =
      context.ref.index->chrom_offsets_data();
  const int chromosome_count = context.ref.contig_count();
  if (peak.chr < 0 || peak.chr >= chromosome_count)
    return false;
  const std::uint64_t chromosome_begin =
      offsets[static_cast<std::size_t>(peak.chr)];
  const std::uint64_t chromosome_end =
      offsets[static_cast<std::size_t>(peak.chr + 1)];
  if (global < chromosome_begin || global >= chromosome_end)
    return false;

  const std::int64_t local =
      static_cast<std::int64_t>(global - chromosome_begin);
  const int seed_length = std::max(1, context.opts.k);
  const int oriented_query = seed_view.seed.read_pos;
  if (oriented_query < 0 ||
      oriented_query + seed_length >
          static_cast<int>(forward_query.size()) ||
      local < 0 ||
      local + seed_length > context.ref.contig_length(peak.chr))
    return false;

  const int width = std::max(1, peak.anchor.ref_start_bin_width);
  const int peak_bin =
      peak.anchor.ref_start_bin_width > 0
          ? peak.anchor.ref_start_bin
          : vote_floor_div(static_cast<int>(peak.raw_ref_start), width);
  const int reference_start =
      static_cast<int>(local) - oriented_query;
  if (std::abs(vote_floor_div(reference_start, width) - peak_bin) > 1)
    return false;

  // The posting was found under this seed's canonical key, so the k-mers
  // match exactly when their orientation bits agree with the strand.
  const bool supported = packed_ref_orientation_compatible(
      seed_view.seed.z, packed, peak.is_rc);
  return supported;
}

// First index of `view` whose global position is at least `global`; the
// positions ascend.
std::uint32_t posting_lower_bound(const KmerPostingView& view,
                                  std::uint32_t from, std::uint64_t global) {
  std::uint32_t low = from;
  std::uint32_t high = view.count;
  while (low < high) {
    const std::uint32_t middle = low + (high - low) / 2;
    if (static_cast<std::uint64_t>(view.positions[middle]) < global)
      low = middle + 1;
    else
      high = middle;
  }
  return low;
}

// The postings of a seed at oriented position `read_pos` that
// posting_supports_peak can accept form one run [first, last) of the list:
// those on the peak's contig at a local position below both the contig's end
// and its length minus the seed length, and at a diagonal bin within one of
// the peak's, which is a local interval since the bin is non-decreasing in
// the local position. False when the bin test's int cast of the local position
// could wrap (a contig past INT_MAX): the whole list must then be walked.
bool peak_posting_run(const DnaContext& context, const VotePeak& peak,
                      int seed_length, int read_length, int read_pos,
                      const KmerPostingView& view, std::uint32_t& first,
                      std::uint32_t& last) {
  first = 0;
  last = 0;
  if (context.ref.index == nullptr || peak.chr < 0 ||
      peak.chr >= context.ref.contig_count() || read_pos < 0 ||
      read_pos + seed_length > read_length)
    return true;
  const std::uint64_t* offsets = context.ref.index->chrom_offsets_data();
  const std::uint64_t chromosome_begin =
      offsets[static_cast<std::size_t>(peak.chr)];
  const std::uint64_t chromosome_end =
      offsets[static_cast<std::size_t>(peak.chr + 1)];
  if (chromosome_end <= chromosome_begin) return true;
  const std::int64_t local_end = std::min<std::int64_t>(
      static_cast<std::int64_t>(chromosome_end - chromosome_begin),
      context.ref.contig_length(peak.chr) - seed_length + 1);
  if (local_end <= 0) return true;
  if (local_end - 1 > std::numeric_limits<int>::max()) return false;
  const int width = std::max(1, peak.anchor.ref_start_bin_width);
  const int peak_bin =
      peak.anchor.ref_start_bin_width > 0
          ? peak.anchor.ref_start_bin
          : vote_floor_div(static_cast<int>(peak.raw_ref_start), width);
  const std::int64_t low = std::max<std::int64_t>(
      0, (static_cast<std::int64_t>(peak_bin) - 1) * width + read_pos);
  const std::int64_t high = std::min<std::int64_t>(
      local_end, (static_cast<std::int64_t>(peak_bin) + 2) * width + read_pos);
  if (low >= high) return true;
  first = posting_lower_bound(
      view, 0, chromosome_begin + static_cast<std::uint64_t>(low));
  last = posting_lower_bound(
      view, first, chromosome_begin + static_cast<std::uint64_t>(high));
  return true;
}

// The forward stretch [begin, begin + length) of the read whose tile grid a
// mask is on: the whole read, or under a vote per stretch one stretch's
// window.
struct TileFrame {
  int begin = 0;
  int length = 0;
};

// dna_forward_query_tile on `frame`'s grid; on the whole read, the same.
int frame_tile(const TileFrame& frame, int oriented_seed_position,
               bool reverse, int read_length, int seed_length,
               int tile_count) noexcept {
  const int forward_position =
      reverse ? read_length - seed_length - oriented_seed_position
              : oriented_seed_position;
  return ::fa::cpu::voting::query_tile_for_position(
      forward_position - frame.begin, frame.length, seed_length, tile_count);
}

QueryTileMask factual_forward_support(
    const DnaContext& context,
    const std::vector<std::uint8_t>& forward_query,
    const VotePeak& peak, const ChainWindowPeakScratch& scratch,
    std::uint64_t& posting_tests, int tile_count, const TileFrame& frame) {
  QueryTileMask support;
  const int read_length = static_cast<int>(forward_query.size());
  const int seed_length = std::max(1, context.opts.k);
  const auto test_view = [&](const DnaLongSeedView& view) {
    // Only the run can hold a supporting posting. posting_tests counts what
    // the walk from the list's start would have tested.
    std::uint32_t posting = 0;
    std::uint32_t last = view.view.count;
    if (view.view.found() &&
        !peak_posting_run(context, peak, seed_length, read_length,
                          view.seed.read_pos, view.view, posting, last)) {
      posting = 0;
      last = view.view.count;
    }
    bool contributed = false;
    for (; posting < last; ++posting) {
      if (posting_supports_peak(
              context, forward_query, view, peak,
              view.view.positions[posting],
              view.view.positions.packed_at(posting))) {
        contributed = true;
        break;
      }
    }
    posting_tests += contributed ? posting + 1 : view.view.count;
    if (contributed) {
      support.set(frame_tile(frame, view.seed.read_pos, peak.is_rc,
                             read_length, seed_length, tile_count));
    }
  };
  if (!scratch.retained_seeds.empty()) {
    for (const ChainWindowRetainedSeed& retained :
         scratch.retained_seeds) {
      test_view({retained.seed, retained.view});
    }
  } else {
    const std::size_t count = std::min(
        scratch.seeds.size(),
        std::min(scratch.seed_views.size(), scratch.seed_used.size()));
    for (std::size_t index = 0; index < count; ++index) {
      if (scratch.seed_used[index])
        test_view({scratch.seeds[index], scratch.seed_views[index]});
    }
  }
  return support;
}

// Computes a ranked peak's tile mask on demand, so peaks that ratio admission
// rejects from their coarse range never walk their postings.
class RankedPeakMaskSource final
    : public ::fa::cpu::voting::CandidateMaskSource {
 public:
  RankedPeakMaskSource(const DnaContext& context,
                       const std::vector<std::uint8_t>& forward_query,
                       const std::vector<VotePeak>& ranked,
                       const ChainWindowPeakScratch& forward_scratch,
                       const ChainWindowPeakScratch& reverse_scratch,
                       std::uint64_t& posting_tests, int tile_count,
                       const TileFrame& frame)
      : context_(context),
        forward_query_(forward_query),
        ranked_(ranked),
        forward_scratch_(forward_scratch),
        reverse_scratch_(reverse_scratch),
        posting_tests_(posting_tests),
        tile_count_(tile_count),
        frame_(frame),
        masks_(ranked.size()),
        ready_(ranked.size(), 0) {}

  const QueryTileMask& mask(int slot) override {
    const std::size_t index = static_cast<std::size_t>(slot);
    if (!ready_[index]) {
      const VotePeak& peak = ranked_[index];
      masks_[index] = factual_forward_support(
          context_, forward_query_, peak,
          peak.is_rc ? reverse_scratch_ : forward_scratch_, posting_tests_,
          tile_count_, frame_);
      ready_[index] = 1;
    }
    return masks_[index];
  }

 private:
  const DnaContext& context_;
  const std::vector<std::uint8_t>& forward_query_;
  const std::vector<VotePeak>& ranked_;
  const ChainWindowPeakScratch& forward_scratch_;
  const ChainWindowPeakScratch& reverse_scratch_;
  std::uint64_t& posting_tests_;
  int tile_count_;
  TileFrame frame_;
  std::vector<QueryTileMask> masks_;
  std::vector<char> ready_;
};

struct GeometryKey {
  int chromosome = -1;
  bool reverse = false;
  int diagonal_bin = 0;

  friend bool operator==(const GeometryKey& left,
                         const GeometryKey& right) {
    return std::tie(left.chromosome, left.reverse, left.diagonal_bin) ==
           std::tie(right.chromosome, right.reverse, right.diagonal_bin);
  }
};

GeometryKey geometry_key(const VotePeak& peak) {
  const int width = std::max(1, peak.anchor.ref_start_bin_width);
  const int bin = peak.anchor.ref_start_bin_width > 0
                      ? peak.anchor.ref_start_bin
                      : vote_floor_div(
                            static_cast<int>(peak.raw_ref_start), width);
  return {peak.chr, peak.is_rc, bin};
}

}  // namespace

int dna_forward_query_tile(int oriented_seed_position, bool reverse,
                           int read_length, int seed_length,
                           int tile_count) noexcept {
  const int forward_position =
      reverse ? read_length - seed_length - oriented_seed_position
              : oriented_seed_position;
  return ::fa::cpu::voting::query_tile_for_position(
      forward_position, read_length, seed_length, tile_count);
}

const DnaPlacementCandidate* DnaPlacementFamily::find(
    CandidateId id) const noexcept {
  const auto found = std::find_if(
      candidates.begin(), candidates.end(),
      [id](const DnaPlacementCandidate& candidate) {
        return candidate.id == id;
      });
  return found == candidates.end() ? nullptr : &*found;
}

DnaPlacementFamily build_dna_placement_family(
    const DnaContext& context, const std::vector<std::uint8_t>& forward_query,
    const std::vector<VotePeak>& whole_read_peaks,
    const ChainWindowPeakScratch& forward_scratch,
    const ChainWindowPeakScratch& reverse_scratch, int frame_begin,
    int frame_length) {
  DnaPlacementFamily family;
  family.read_length = static_cast<int>(forward_query.size());
  family.seed_length = std::max(1, context.opts.k);
  family.tile_count = context.opts.query_tiles;
  // The grid the admission compares masks on; candidates keep the family's.
  const TileFrame whole{0, family.read_length};
  const bool framed = frame_length > 0;
  const TileFrame frame = framed ? TileFrame{frame_begin, frame_length} : whole;
  std::vector<VotePeak> ranked = whole_read_peaks;
  // Peaks tied on every evidence key are ordered by a read-seeded hash of
  // their locus, as minimap2 breaks ties; under the vote per stretch by their
  // read_lo first (map.cpp stretch_to_contig). Catalogue ranks and ids follow
  // this order, so every later tie-break on rank or id inherits it.
  const std::uint32_t seed = context.vote_tie_seed;
  std::stable_sort(ranked.begin(), ranked.end(),
                   [seed](const VotePeak& left, const VotePeak& right) {
                     return chain_peak_better_seeded(left, right, seed);
                   });

  // Ratio admission (the default): masks are computed lazily and the
  // admission gets each peak's coarse tile range, which contains its mask.
  // With --vote-ratio 0 the count rule applies and masks are computed eagerly,
  // except in the all-chains lane, which reads none.
  // Ratio admission compares masks on a grid of at most
  // kMaxAdmissionQueryTiles tiles; the partition keeps the family's.
  const bool ratio_admission = context.opts.vote_admission_ratio > 0.0;
  const int admission_tiles = std::min(
      family.tile_count, ::fa::cpu::voting::kMaxAdmissionQueryTiles);
  std::vector<GeometryKey> keys;
  std::vector<VotePeak> key_peaks;
  std::vector<CandidateInput> inputs;
  int rank = 0;
  for (std::size_t peak_index = 0; peak_index < ranked.size(); ++peak_index) {
    const VotePeak& peak = ranked[peak_index];
    const GeometryKey key = geometry_key(peak);
    auto found = std::find(keys.begin(), keys.end(), key);
    std::size_t key_index = static_cast<std::size_t>(
        std::distance(keys.begin(), found));
    if (found == keys.end()) {
      keys.push_back(key);
      key_peaks.push_back(peak);
      key_index = keys.size() - 1;
    }
    CandidateInput input;
    input.equivalence_key = key_index + 1;
    input.lane = peak.is_rc ? 1 : 0;
    input.catalogue_rank = rank++;
    input.vote_evidence = std::max(peak.vote_score, peak.support);
    if (ratio_admission) {
      input.mask_slot = static_cast<int>(peak_index);
      if (peak.evidence_read_lo <= peak.evidence_read_hi) {
        const int tile_a = frame_tile(
            frame, peak.evidence_read_lo, peak.is_rc, family.read_length,
            family.seed_length, admission_tiles);
        const int tile_b = frame_tile(
            frame, peak.evidence_read_hi, peak.is_rc, family.read_length,
            family.seed_length, admission_tiles);
        input.coarse_tile_lo = std::min(tile_a, tile_b);
        input.coarse_tile_hi = std::max(tile_a, tile_b);
      }
    } else if (!context.opts.all_chains) {
      input.support = factual_forward_support(
          context, forward_query, peak,
          peak.is_rc ? reverse_scratch : forward_scratch,
          family.exact_posting_tests, family.tile_count, frame);
    }
    inputs.push_back(input);
  }

  // Under ratio admission the per-lane count is only a cost ceiling; the
  // ratio does the selecting.
  RankedPeakMaskSource mask_source(context, forward_query, ranked,
                                   forward_scratch, reverse_scratch,
                                   family.exact_posting_tests,
                                   admission_tiles, frame);
  CandidateCatalogue catalogue = ::fa::cpu::voting::build_candidate_catalogue(
      std::move(inputs), context.opts.catalogue_lane_bound,
      ratio_admission ? context.opts.vote_admission_ratio : 0.0,
      ratio_admission ? &mask_source : nullptr, admission_tiles,
      context.opts.all_chains);
  // Admitted on the coarser grid, or on a stretch's, a candidate's support is
  // rebuilt on the family's.
  if (framed || (ratio_admission && admission_tiles != family.tile_count)) {
    for (auto& generic : catalogue.candidates) {
      const VotePeak& peak =
          key_peaks[static_cast<std::size_t>(generic.equivalence_key - 1)];
      generic.support = factual_forward_support(
          context, forward_query, peak,
          peak.is_rc ? reverse_scratch : forward_scratch,
          family.exact_posting_tests, family.tile_count, whole);
    }
  }
  for (const auto& generic : catalogue.candidates) {
    const std::size_t key_index =
        static_cast<std::size_t>(generic.equivalence_key - 1);
    if (key_index >= key_peaks.size()) return family;
    DnaPlacementCandidate candidate;
    candidate.id = generic.id;
    candidate.peak = key_peaks[key_index];
    candidate.support = generic.support;
    candidate.equivalence_key = generic.equivalence_key;
    candidate.lane = generic.lane;
    candidate.catalogue_rank = generic.catalogue_rank;
    candidate.vote_evidence = generic.vote_evidence;
    family.candidates.push_back(std::move(candidate));
    if (key_peaks[key_index].is_rc)
      ++family.reverse_candidates;
    else
      ++family.forward_candidates;
  }
  family.valid = !family.candidates.empty();
  return family;
}

// Each stretch's catalogue is its own: its ratio admission, lane bound and
// key dedup see only its peaks, and a place named in two stretches stays two
// candidates. Ids, catalogue ranks and keys continue the contig's, so a tie
// on rank breaks in stretch order.
void append_stretch_family(DnaPlacementFamily& family,
                           DnaPlacementFamily stretch) {
  if (family.read_length == 0) {
    family.read_length = stretch.read_length;
    family.seed_length = stretch.seed_length;
    family.tile_count = stretch.tile_count;
  }
  int rank_base = 0;
  std::uint64_t key_base = 0;
  for (const DnaPlacementCandidate& candidate : family.candidates) {
    rank_base = std::max(rank_base, candidate.catalogue_rank + 1);
    key_base = std::max(key_base, candidate.equivalence_key);
  }
  for (DnaPlacementCandidate& candidate : stretch.candidates) {
    candidate.id = static_cast<CandidateId>(family.candidates.size());
    candidate.catalogue_rank += rank_base;
    candidate.equivalence_key += key_base;
    family.candidates.push_back(std::move(candidate));
  }
  family.forward_candidates += stretch.forward_candidates;
  family.reverse_candidates += stretch.reverse_candidates;
  family.exact_posting_tests += stretch.exact_posting_tests;
  family.valid = !family.candidates.empty();
}

std::string CandidateIdOverflow::refusal(const std::string* name, int length,
                                         std::size_t candidates,
                                         int lane_bound) {
  return (name ? "query contig '" + *name + "'"
               : std::string("a query contig")) +
         " (" + std::to_string(length) + " bp) has " +
         std::to_string(candidates) +
         " candidates in its catalogue per stretch, more than a candidate id "
         "holds (" +
         std::to_string(std::numeric_limits<CandidateId>::max()) +
         "); lower --max-cands (" + std::to_string(lane_bound) +
         " per strand here) to map it";
}

}  // namespace fa::cpu::lr
