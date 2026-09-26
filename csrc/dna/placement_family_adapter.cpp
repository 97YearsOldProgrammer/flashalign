#include "placement_family_adapter.h"

#include "../index/format.h"
#include "../voting/state.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <tuple>
#include <vector>

namespace fa::cpu::lr {
namespace {

using ::fa::cpu::voting::CandidateCatalogue;
using ::fa::cpu::voting::CandidateId;
using ::fa::cpu::voting::CandidateInput;
using ::fa::cpu::voting::QueryPartitionProblem;
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

QueryTileMask factual_forward_support(
    const DnaContext& context,
    const std::vector<std::uint8_t>& forward_query,
    const VotePeak& peak, const ChainWindowPeakScratch& scratch,
    std::uint64_t& posting_tests) {
  QueryTileMask support;
  const int read_length = static_cast<int>(forward_query.size());
  const int seed_length = std::max(1, context.opts.k);
  const auto test_view = [&](const DnaLongSeedView& view) {
    bool contributed = false;
    for (std::uint32_t posting = 0; posting < view.view.count; ++posting) {
      ++posting_tests;
      if (posting_supports_peak(
              context, forward_query, view, peak,
              view.view.positions[posting],
              view.view.positions.packed_at(posting))) {
        contributed = true;
        break;
      }
    }
    if (contributed) {
      support.set(dna_forward_query_tile(
          view.seed.read_pos, peak.is_rc, read_length, seed_length));
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
                       std::uint64_t& posting_tests)
      : context_(context),
        forward_query_(forward_query),
        ranked_(ranked),
        forward_scratch_(forward_scratch),
        reverse_scratch_(reverse_scratch),
        posting_tests_(posting_tests),
        masks_(ranked.size()),
        ready_(ranked.size(), 0) {}

  const QueryTileMask& mask(int slot) override {
    const std::size_t index = static_cast<std::size_t>(slot);
    if (!ready_[index]) {
      const VotePeak& peak = ranked_[index];
      masks_[index] = factual_forward_support(
          context_, forward_query_, peak,
          peak.is_rc ? reverse_scratch_ : forward_scratch_, posting_tests_);
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
                           int read_length, int seed_length) noexcept {
  const int forward_position =
      reverse ? read_length - seed_length - oriented_seed_position
              : oriented_seed_position;
  return ::fa::cpu::voting::query_tile_for_position(
      forward_position, read_length, seed_length);
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
    const ChainWindowPeakScratch& reverse_scratch) {
  DnaPlacementFamily family;
  family.read_length = static_cast<int>(forward_query.size());
  family.seed_length = std::max(1, context.opts.k);
  std::vector<VotePeak> ranked = whole_read_peaks;
  // Peaks tied on every evidence key are ordered by a read-seeded hash of
  // their locus, as minimap2 breaks ties. Catalogue ranks and ids follow this
  // order, so every later tie-break on rank or id inherits it.
  const std::uint32_t seed = context.vote_tie_seed;
  std::stable_sort(ranked.begin(), ranked.end(),
                   [seed](const VotePeak& left, const VotePeak& right) {
                     return chain_peak_better_seeded(left, right, seed);
                   });

  // Ratio admission (the default): masks are computed lazily and the
  // admission gets each peak's coarse tile range, which contains its mask.
  // With --vote-ratio 0 the count rule applies and masks are computed eagerly.
  const bool ratio_admission = context.opts.vote_admission_ratio > 0.0;
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
        const int tile_a = dna_forward_query_tile(
            peak.evidence_read_lo, peak.is_rc, family.read_length,
            family.seed_length);
        const int tile_b = dna_forward_query_tile(
            peak.evidence_read_hi, peak.is_rc, family.read_length,
            family.seed_length);
        input.coarse_tile_lo = std::min(tile_a, tile_b);
        input.coarse_tile_hi = std::max(tile_a, tile_b);
      }
    } else {
      input.support = factual_forward_support(
          context, forward_query, peak,
          peak.is_rc ? reverse_scratch : forward_scratch,
          family.exact_posting_tests);
    }
    inputs.push_back(input);
  }

  QueryPartitionProblem problem;
  problem.tile_count = ::fa::cpu::voting::kQueryTileCount;
  for (int tile = 0; tile < problem.tile_count; ++tile) {
    if (::fa::cpu::voting::query_tile_end(
            tile + 1, family.read_length, family.seed_length) >
        ::fa::cpu::voting::query_tile_begin(
            tile, family.read_length, family.seed_length)) {
      problem.valid_tiles.set(tile);
    }
  }
  problem.parameters = context.opts.query_partition;
  // Under ratio admission the per-lane count is only a cost ceiling; the
  // ratio does the selecting.
  RankedPeakMaskSource mask_source(context, forward_query, ranked,
                                   forward_scratch, reverse_scratch,
                                   family.exact_posting_tests);
  problem.catalogue = ::fa::cpu::voting::build_candidate_catalogue(
      std::move(inputs),
      ratio_admission ? ::fa::cpu::voting::kCatalogueLaneBound
                      : ::fa::cpu::voting::kCountAdmissionLaneBound,
      ratio_admission ? context.opts.vote_admission_ratio : 0.0,
      ratio_admission ? &mask_source : nullptr);
  for (const auto& generic : problem.catalogue.candidates) {
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
  family.partition = ::fa::cpu::voting::solve_query_partition(problem);
  family.valid = !family.candidates.empty();
  return family;
}

}  // namespace fa::cpu::lr
