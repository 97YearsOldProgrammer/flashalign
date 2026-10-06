// Resolved long-read seeding parameters (LongReadSeedContext) and the vote occurrence-cap
// policy.
#pragma once

#include "../index/index.h"       // SeedIndex
#include "../index/format.h"      // KmerPostingView

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <string>
#include <vector>

namespace fa { namespace cpu { namespace lr {

enum class LongOccPolicy {
    Off,
    Fixed,
    Platform,
};

struct LongOccPolicyConfig {
    LongOccPolicy policy = LongOccPolicy::Off;
    int primary_occ_cap = 0;
    int platform_occ_cap = 12;
};

// The DNA vote occurrence cap. Under LongOccPolicy::Platform, the presets' policy, it is
//     max(platform_occ_cap, mm_idx_cal_max_occ(attached index, kLongOccQuantileF))
// limited by the preset's ceiling, if any. The quantile only raises the preset floor, so a
// repeat-rich reference gets a higher cap; f is chosen so a human reference resolves to the
// floor. --max-vote-occ selects LongOccPolicy::Fixed (or Off at 0), which skips the rule
// and its index scan.
inline constexpr double kLongOccQuantileF = 1.81e-4;

// How the cap was resolved, filled at index-attach time for reporting only. `resolved` is
// false until an index is attached.
struct LongOccIndexResolution {
    bool resolved = false;
    double quantile_f = 0.0;
    std::uint64_t distinct_keys = 0;
    std::uint64_t total_postings = 0;
    std::int64_t raw_cap = 0;   // faix_cal_max_occ, before the floor
    int floor_cap = 0;          // the preset/caller cap the quantile floors to
    int ceiling_cap = 0;        // the ceiling in force, 0 = none
    int resolved_cap = 0;       // min(ceiling_cap, max(floor_cap, raw_cap))
};

inline std::string long_occ_policy_token(std::string value) {
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](unsigned char c) {
            const char lower = static_cast<char>(std::tolower(c));
            return lower == '-' ? '_' : lower;
        });
    return value;
}

inline std::string long_occ_policy_name(LongOccPolicy policy) {
    switch (policy) {
        case LongOccPolicy::Fixed:
            return "fixed";
        case LongOccPolicy::Platform:
            return "platform";
        case LongOccPolicy::Off:
        default:
            return "off";
    }
}

inline LongOccPolicy long_occ_policy_from_string(
    const std::string& value,
    LongOccPolicy fallback
) {
    const std::string mode = long_occ_policy_token(value);
    if (mode == "off" || mode == "none" || mode == "0") {
        return LongOccPolicy::Off;
    }
    if (mode == "fixed" || mode == "cap" || mode == "cap_fixed") {
        return LongOccPolicy::Fixed;
    }
    if (mode == "platform" || mode == "platform_fixed") {
        return LongOccPolicy::Platform;
    }
    return fallback;
}

inline bool long_occ_policy_string_valid(const std::string& value) {
    const std::string mode = long_occ_policy_token(value);
    return mode == "off" || mode == "none" || mode == "0" ||
           mode == "fixed" || mode == "cap" || mode == "cap_fixed" ||
           mode == "platform" || mode == "platform_fixed";
}

inline int effective_long_primary_occ_cap(
    const LongOccPolicyConfig& cfg
) {
    switch (cfg.policy) {
        case LongOccPolicy::Fixed:
            return std::max(0, cfg.primary_occ_cap);
        case LongOccPolicy::Platform:
            // The rule resolves into platform_occ_cap; before an index is attached this
            // is the preset floor, never 0 (no cap).
            return std::max(0, cfg.platform_occ_cap);
        case LongOccPolicy::Off:
        default:
            return 0;
    }
}

inline bool seed_allowed_by_long_occ_policy(
    const KmerPostingView& view,
    const LongOccPolicyConfig& cfg
) {
    const int cap = effective_long_primary_occ_cap(cfg);
    if (cap <= 0) return true;
    if (!view.found() || view.count == 0) return false;
    return view.occurrence > 0 &&
           view.occurrence <= static_cast<uint32_t>(cap);
}

// The same test for a found seed, on a given occurrence.
inline bool occurrence_allowed_by_long_occ_policy(
    uint32_t occurrence,
    const LongOccPolicyConfig& cfg
) {
    const int cap = effective_long_primary_occ_cap(cfg);
    if (cap <= 0) return true;
    return occurrence > 0 && occurrence <= static_cast<uint32_t>(cap);
}

// Seeds per strand that a larger DNA vote selection keeps (syncmer.h).
inline constexpr int kVoteSeedNestBase = 128;

struct LongReadSeedContext {
    int k = 0;
    const SeedIndex* index = nullptr;

    // Closed-syncmer extraction.
    int syncmer_s = 0;                   // ClosedSyncmerConfig::s
    int syncmer_downsample = 1;          // ClosedSyncmerConfig::downsample

    // Occurrence-aware seed admission. When enabled (DNA), the kernel keeps one
    // lowest-occ seed per read-position strip (occ <= cap). When disabled (RNA),
    // the window votes every extracted seed (capped only by max_query_seeds).
    int max_query_seeds_per_strand = 0;
    // DNA: above kVoteSeedNestBase seeds per strand the selection keeps the
    // kVoteSeedNestBase selection and fills from the larger one (syncmer.h).
    bool nested_vote_seeds = false;
    bool syncmer_occ_aware_enabled = false;  // DNA: true; RNA: false
    LongOccPolicyConfig occ_policy;

    // Positive preset (or explicit override) reference-start diagonal width.
    int vote_diag_bin_width = 64;
    // Adaptive width (DNA presets), in integer arithmetic:
    //   W(L) = clamp(L / vote_diag_slope_den, vote_diag_bin_width, vote_diag_width_max)
    // a length term for indel drift, floored at the base width. A denominator <= 0 drops the
    // length term. False for RNA and for a fixed --dw override.
    bool vote_diag_bin_width_adaptive = false;
    int vote_diag_slope_den = 128;
    int vote_diag_width_max = 2048;
    int min_support = 3;
    int chain_max_candidates_per_window = 0;
    // DNA ratio admission, and the all-chains lane under either admission rule: drain the
    // vote heap in batches, so the exact refine walk is paid per chunk of buckets rather
    // than per bucket. Same output as the sequential drain.
    bool vote_batched_refine = false;
    // DNA: a vote tile that admits no seed votes with its rarest found seed at or under
    // this occurrence (options/dna_profile.h kDnaTileRescueOcc); 0 = none, as on RNA.
    int tile_rescue_occ = 0;
    // DNA, options/dna_profile.h skip_self: the contig that is the query read itself,
    // found by its exact name; -1 when the read is not in the index or the option is off.
    int self_contig = -1;
    const std::vector<std::string>* chr_names = nullptr;

    bool syncmer_occ_aware_active() const {
        return syncmer_occ_aware_enabled;
    }
};

// The occurrence the DNA vote counts for a seed: the index's, less the key's postings on
// the query read's own contig (ctx.self_contig), so 0 for a key that occurs only in the
// read itself.
inline uint32_t vote_seed_occurrence(const LongReadSeedContext& ctx,
                                     const KmerPostingView& view) {
    if (ctx.self_contig < 0 || !view.found()) return view.occurrence;
    const uint32_t self = static_cast<uint32_t>(ctx.self_contig);
    const uint32_t lo =
        view.positions.lower_bound_packed(pack_ref_pos(self, 0));
    const uint32_t hi =
        view.positions.lower_bound_packed(pack_ref_pos(self + 1, 0));
    return view.occurrence - (hi - lo);
}

inline int effective_vote_diag_bin_width(
    const LongReadSeedContext& ctx, int read_span) {
    const int base_width = std::max(1, ctx.vote_diag_bin_width);
    if (!ctx.vote_diag_bin_width_adaptive || read_span <= 0)
        return base_width;
    // read_span > 0 makes the truncating division a floor; integer arithmetic keeps bin
    // boundaries identical on every target.
    const int drift = ctx.vote_diag_slope_den > 0
                          ? read_span / ctx.vote_diag_slope_den
                          : 0;
    return std::max(base_width, std::min(drift, ctx.vote_diag_width_max));
}

}}}  // namespace fa::cpu::lr
