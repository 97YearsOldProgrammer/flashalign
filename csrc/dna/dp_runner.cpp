#include "dp_runner.h"

#include "../dp/ksw2_align.h" // ::fa::cpu::ksw2_simple_mat + KSW flags
#include "../dp/params.h"     // ::fa::cpu::DpMapOpt defaults/effective bands

#include <algorithm>
#include <utility>

namespace fa {
namespace cpu {
namespace lr {

// As in minimap2: flank extensions use the band bw_eff() and gap fills
// bw_long_eff() (max(q_span, r_span) across a long join); a gap fill runs the
// KSW_EZ_APPROX_MAX pass first; a fill whose matrix exceeds opt.max_sw_mat is
// not run but returned Z-dropped.
realization::RealizationOutcome run_dna_long_realization(
    realization::RealizationRole role, const DpScoringParams &dp,
    realization::SequenceSlice query, realization::SequenceSlice target,
    int caller_requested_band, int interior_zdrop, bool long_join, realization::ExecutionContext *execution_context,
    bool inversion_probe_enabled, const VerifiedRegionInputs *verified_region,
    const realization::InversionProbeGate *inversion_probe_gate) {
  realization::RealizationRequest request;
  request.role = role;
  request.objective = realization::objective_for_role(role);
  request.query = query;
  request.target = target;
  const bool left =
      request.objective == realization::EndpointObjective::OpenEndedLeft;
  request.query.reversed_for_executor = left;
  request.target.reversed_for_executor = left;

  request.scoring.match = dp.match;
  request.scoring.mismatch = dp.mismatch;
  request.scoring.ambiguity = dp.ambi;
  ::fa::cpu::ksw2_simple_mat(request.scoring.substitution_matrix.data(),
                             dp.match, dp.mismatch, dp.ambi);
  request.scoring.gap1 = {dp.gap_open1, dp.gap_extend1};
  request.scoring.gap2 = {dp.gap_open2, dp.gap_extend2};
  request.scoring.model =
      dp.gap_open1 == dp.gap_open2 && dp.gap_extend1 == dp.gap_extend2
          ? realization::ScoreModel::SingleAffine
          : realization::ScoreModel::DualAffine;

  ::fa::cpu::DpMapOpt opt;
  opt.bw = dp.bw;
  opt.bw_long = dp.bw_long;
  request.caller_requested_band = caller_requested_band;
  request.configured_band = dp.bw;
  request.configured_long_band = dp.bw_long;
  request.selected_band =
      role == realization::RealizationRole::DnaLocalInversionMiddle
          ? static_cast<int>(dp.bw * 1.5)
          : request.objective == realization::EndpointObjective::PinnedGlobal
          ? (long_join ? std::max(query.length, target.length)
                       : opt.bw_long_eff())
          : opt.bw_eff();
  request.zdrop =
      request.objective == realization::EndpointObjective::PinnedGlobal
          ? interior_zdrop
          : dp.zdrop;
  request.inversion_zdrop = dp.inversion_zdrop;
  const bool inversion_role =
      role == realization::RealizationRole::DnaInternalFill ||
      role == realization::RealizationRole::DnaSupplementaryFill;
  request.inversion_probe_enabled =
      inversion_role && interior_zdrop >= 0 && inversion_probe_enabled;
  request.inversion_max_gap = dp.inversion_max_gap;
  request.inversion_min_chain_score = dp.inversion_min_chain_score;
  request.inversion_min_dp_max = dp.inversion_min_dp_max;
  if (inversion_probe_gate != nullptr)
    request.inversion_probe_gate = *inversion_probe_gate;
  request.end_bonus =
      request.objective == realization::EndpointObjective::PinnedGlobal ||
              role == realization::RealizationRole::DnaLocalInversionMiddle
          ? -1
          : dp.tail_end_bonus;
  request.matrix_cell_cap = opt.max_sw_mat;
  request.long_join = long_join;
  request.ksw_flags =
      request.objective == realization::EndpointObjective::PinnedGlobal
          ? KSW_EZ_APPROX_MAX
          : KSW_EZ_EXTZ_ONLY;
  if (left)
    request.ksw_flags |= KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR;
  if (verified_region != nullptr) {
    request.verified = true;
    request.gaps = verified_region->gaps;
    request.gap_count = verified_region->gap_count;
  }

  realization::ExecutionContext local_context;
  realization::ExecutionContext &context =
      execution_context == nullptr ? local_context : *execution_context;
  return realization::run(request, realization::production_policy(),
                          realization::production_executors(), context);
}

} // namespace lr
} // namespace cpu
} // namespace fa
