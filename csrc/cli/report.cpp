#include "cli/report.h"

#include "rna/realization/rival_lifecycle.h"
#include "seeding/context.h"
#include "fa_version.h"

#include <sys/resource.h>

#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

namespace fa::cpu::cli {

int64_t peak_rss_bytes() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0)
    return 0;
#ifdef __APPLE__
  return static_cast<int64_t>(ru.ru_maxrss);
#else
  return static_cast<int64_t>(ru.ru_maxrss) * 1024;
#endif
}

namespace {

// User plus system time of every thread of the process.
double cpu_seconds() {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) != 0)
    return 0.0;
  return static_cast<double>(ru.ru_utime.tv_sec + ru.ru_stime.tv_sec) +
         1e-6 * static_cast<double>(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);
}

double gigabytes(int64_t bytes) {
  return static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
}

} // namespace

void write_stats(const RunSummary& summary, std::ostream& out) {
  const auto row = [&](const char* key, const std::string& value) {
    out << "  " << std::left << std::setw(22) << key << value << "\n";
  };
  const auto fmt = [](const char* format, double value) {
    char buf[64];
    std::snprintf(buf, sizeof buf, format, value);
    return std::string(buf);
  };
  const auto i2s = [](int64_t v) { return std::to_string(v); };
  const int64_t unmapped =
      summary.reads > summary.mapped_reads
          ? summary.reads - summary.mapped_reads
          : 0;
  const double mapped_pct =
      summary.reads > 0 ? 100.0 * static_cast<double>(summary.mapped_reads) /
                              static_cast<double>(summary.reads)
                        : 0.0;
  out << "[flashalign] run summary\n";
  row("preset:", summary.preset);
  row("output format:", summary.output_format);
  row("reads:", i2s(summary.reads));
  row("bases:", i2s(summary.bases));
  row("mapped reads:",
      i2s(summary.mapped_reads) + " (" + fmt("%.1f", mapped_pct) + "%)");
  row("unmapped reads:", i2s(unmapped));
  row("records written:", i2s(summary.records));
  row("supplementary:", i2s(summary.supplementary_records));
  row("secondary:", i2s(summary.secondary_records));
  if (summary.rna) {
    row("spliced reads:", i2s(summary.spliced_reads));
    row("junctions:", i2s(summary.junctions));
  }
  row("wall time:", fmt("%.2f s", summary.wall_seconds));
  const double map_s = summary.map_seconds;
  row("throughput:",
      fmt("%.1f reads/s, ",
          map_s > 0.0 ? static_cast<double>(summary.reads) / map_s : 0.0) +
          fmt("%.1f Mbp/s",
              map_s > 0.0 ? static_cast<double>(summary.bases) / map_s / 1e6
                          : 0.0));
  row("peak RSS:", fmt("%.1f GB", gigabytes(summary.peak_rss_bytes)));
}

void write_run_footer(int argc, char** argv, double real_seconds) {
  std::string cmd;
  for (int i = 0; i < argc; ++i) {
    if (i > 0)
      cmd += ' ';
    cmd += argv[i];
  }
  std::fprintf(stderr, "[flashalign] Version: %s\n", FA_VERSION);
  std::fprintf(stderr, "[flashalign] CMD: %s\n", cmd.c_str());
  std::fprintf(stderr,
               "[flashalign] Real time: %.1f sec; CPU: %.1f sec; "
               "Peak RSS: %.1f GB\n",
               real_seconds, cpu_seconds(), gigabytes(peak_rss_bytes()));
}

namespace {

// The resolved configuration as text, each value tagged with its source
// (builtin, preset, index, cli, ...). It needs no aligner, so no reference is
// opened or indexed.
std::string resolved_config_text(const AlignOptions& opt,
                                 const fa::cpu::options::ResolvedOptions& cfg,
                                 fa::cpu::api::BackendMode mode) {
  const auto& mapping = cfg.long_read();
  const auto* rna = cfg.rna();
  const bool is_rna = mode == fa::cpu::api::BackendMode::RnaSplice;
  // (k, s) come from a loaded index, otherwise from the preset.
  const bool seeding_from_index = cfg.index.from_index;

  std::ostringstream os;
  const auto row = [&](const char* key, const std::string& val,
                       const char* src) {
    os << "  " << std::left << std::setw(24) << key;
    // A key or value that fills its column still gets a separating space.
    if (std::string_view(key).size() >= 24)
      os << ' ';
    os << std::setw(16) << val;
    if (val.size() >= 16)
      os << ' ';
    os << "[" << src << "]\n";
    os << std::right; // reset
  };
  const auto i2s = [](long long v) { return std::to_string(v); };
  // Default precision, so 0.8 prints "0.8".
  const auto f2s = [](double value) {
    std::ostringstream text;
    text << value;
    return text.str();
  };
  const auto b2s = [](bool v) { return std::string(v ? "1" : "0"); };
  const auto src = [](bool given, const char* owner) {
    return given ? "explicit" : owner;
  };
  const auto dp_src = [&](const std::optional<int>& o) {
    return src(o.has_value(), "preset");
  };
  // RNA presets leave the tile objective at its struct defaults.
  const char* const tile_owner = is_rna ? "builtin" : "preset";

  row("mode", is_rna ? "rna_long" : "dna_long", "resolved");
  row("preset", opt.preset, opt.preset_source.c_str());
  // asm5, asm10 and asm20 are not qualified on intact chromosomes; ava-ont
  // and ava-hifi are measured on simulated reads only.
  const bool assembly = fa::cpu::api::is_assembly_preset(opt.preset);
  const bool overlap = fa::cpu::api::is_overlap_preset(opt.preset);
  row("preset_status", assembly || overlap ? "experimental" : "supported",
      "preset");
  row("index_source", seeding_from_index ? "faix" : "built", "cli");
  row("k", i2s(cfg.index.k), seeding_from_index ? "index" : "preset");
  row("syncmer_s", i2s(cfg.index.syncmer_s),
      seeding_from_index ? "index" : "preset");
  row("syncmer_downsample", i2s(cfg.index.syncmer_downsample), "builtin");
  row("output_format", opt.format, "cli");
  row("threads", i2s(cfg.common.num_threads),
      src(opt.threads && *opt.threads > 0, "builtin"));

  row("min_support", i2s(cfg.common.min_support),
      src(opt.min_support.has_value(), "preset"));
  row("max_query_seeds", i2s(cfg.common.max_query_seeds_per_strand),
      src(opt.vote_seeds.has_value(), "preset"));
  // Above kVoteSeedNestBase a DNA selection keeps the kVoteSeedNestBase one
  // and fills it (seeding/syncmer.h).
  if (!is_rna)
    row("vote_seed_sampler",
        cfg.common.max_query_seeds_per_strand > fa::cpu::lr::kVoteSeedNestBase
            ? "nested"
            : "original",
        "derived");
  // max_seed_occ is the effective vote occurrence cap; 0 means no filtering.
  fa::cpu::lr::LongOccPolicyConfig occ_cfg;
  occ_cfg.policy = mapping.long_occ_policy;
  occ_cfg.primary_occ_cap = mapping.long_primary_occ_cap;
  occ_cfg.platform_occ_cap = mapping.long_occ_cap;
  const bool occ_explicit = opt.max_vote_occ.has_value();
  // The default DNA rule is max(preset floor, index quantile), capped by an
  // optional ceiling. Until an index has been scanned (e.g. --show-config),
  // the row prints the rule instead of a number.
  const auto& occ_index = mapping.long_occ_index_resolution;
  const bool occ_rule =
      !is_rna && !occ_explicit &&
      mapping.long_occ_policy == fa::cpu::lr::LongOccPolicy::Platform;
  if (occ_rule && !occ_index.resolved) {
    // As resolve_index_occ_cap_into: a ceiling applies when set and not
    // below the floor.
    const bool ceiling = mapping.long_occ_ceiling > 0 &&
                         mapping.long_occ_ceiling >= mapping.long_occ_cap;
    std::ostringstream rule;
    if (ceiling)
      rule << "min(" << mapping.long_occ_ceiling << ", ";
    rule << "max(" << mapping.long_occ_cap << ", index quantile f="
         << fa::cpu::lr::kLongOccQuantileF << ')';
    if (ceiling)
      rule << ')';
    row("max_seed_occ", rule.str(), "rule");
  } else {
    row("max_seed_occ",
        i2s(fa::cpu::lr::effective_long_primary_occ_cap(occ_cfg)),
        occ_index.resolved ? "index" : (occ_explicit ? "explicit" : "preset"));
  }
  row("occ_policy", fa::cpu::lr::long_occ_policy_name(mapping.long_occ_policy),
      occ_explicit ? "explicit" : "preset");
  row("long_occ_cap", i2s(mapping.long_occ_cap),
      occ_index.resolved ? "index" : "preset");
  row("long_primary_occ_cap", i2s(mapping.long_primary_occ_cap),
      occ_explicit ? "explicit" : "builtin");
  // How the index produced the cap.
  if (occ_index.resolved) {
    std::ostringstream detail;
    detail << "f=" << occ_index.quantile_f
           << " n=" << occ_index.distinct_keys
           << " raw=" << occ_index.raw_cap
           << " floor=" << occ_index.floor_cap;
    if (occ_index.ceiling_cap > 0)
      detail << " ceiling=" << occ_index.ceiling_cap;
    detail << " cap=" << occ_index.resolved_cap;
    row("seed_occ_quantile", detail.str(), "index");
  }
  // The vote's empty-tile rescue M; 0 on the RNA presets.
  row("tile_rescue_occ", i2s(mapping.dna_tile_rescue_occ), "preset");
  // The overlap presets' two mechanisms and --dual, printed only under them.
  if (overlap) {
    row("skip_self", b2s(mapping.skip_self), "preset");
    row("all_chains", b2s(mapping.all_chains), "preset");
    row("dual", mapping.dual ? "yes" : "no",
        src(opt.dual.has_value(), "preset"));
  }

  if (!is_rna)
    row("vote_ratio", f2s(mapping.vote_admission_ratio),
        src(opt.dna_vote_admission_ratio.has_value(), "preset"));
  row("vote_diag_bin_width", i2s(mapping.vote_diag_bin_width),
      src(opt.vote_diag_bin_width.has_value(), "preset"));
  // --dw turns the adaptive model off.
  row("vote_diag_adaptive", b2s(mapping.vote_diag_bin_width_adaptive),
      src(opt.vote_diag_bin_width.has_value(), "preset"));
  row("vote_diag_slope_den", i2s(mapping.vote_diag_slope_den),
      src(opt.vote_diag_slope_den.has_value(), "preset"));
  row("vote_diag_width_max", i2s(mapping.vote_diag_width_max),
      src(opt.vote_diag_width_max.has_value(), "preset"));
  row("tile_score.hit", i2s(mapping.query_partition.supported_tile_reward),
      src(opt.tile_supported_reward.has_value(), tile_owner));
  row("tile_score.block", i2s(mapping.query_partition.block_open_cost),
      src(opt.tile_block_open_cost.has_value(), tile_owner));
  row("tile_score.null", i2s(mapping.query_partition.null_tile_cost),
      src(opt.tile_null_cost.has_value(), tile_owner));
  row("tile_score.miss",
      i2s(mapping.query_partition.unsupported_ownership_cost),
      src(opt.tile_unsupported_cost.has_value(), tile_owner));
  row("partition.min_tiles",
      i2s(mapping.query_partition.minimum_supported_tiles_per_non_null_block),
      tile_owner);
  // The splice presets partition into their own 128 tiles.
  row("partition.tiles", i2s(mapping.query_tiles),
      is_rna ? "fixed" : src(opt.tiles.has_value(), "preset"));
  // Ratio admission compares masks on at most kMaxAdmissionQueryTiles tiles.
  if (!is_rna)
    row("partition.admission_tiles",
        mapping.vote_admission_ratio > 0.0
            ? i2s(std::min(mapping.query_tiles,
                           fa::cpu::voting::kMaxAdmissionQueryTiles))
            : std::string("none (count admission)"),
        "derived");
  row("mapq.output_range", "0..60", "fixed");

  // On a splice preset, the fine chain's query gap.
  row("chain_max_gap", i2s(mapping.cigar_local_interval_anchor_chain_max_gap),
      src(opt.dp_max_gap.has_value(), "preset"));
  row("interval_pad", i2s(mapping.cigar_local_interval_anchor_interval_pad),
      "builtin");
  // The bound the vote catalogue runs with.
  if (is_rna)
    row("chain_max_cands", i2s(mapping.chain_max_candidates_per_window),
        "builtin");
  else
    row("chain_max_cands", i2s(fa::cpu::lr::dna_chain_max_candidates(mapping)),
        src(opt.max_cands.has_value(), overlap ? "preset" : "derived"));
  if (!is_rna)
    row("screen_band", i2s(mapping.screen_diag_band),
        src(opt.screen_band.has_value(), "builtin"));

  // lr, lr:hq: -A -B -O -E -z --score-N are the gap-fill row; the end row,
  // which prices every path, is the preset's (dp_end_row). The assembly
  // presets fill under their end row, so the two rows print as one.
  const bool fill = !is_rna;
  row("dp_match", i2s(fill ? mapping.fill_dp_match : mapping.cigar_dp_match),
      dp_src(opt.dp_match));
  row("dp_mismatch",
      i2s(fill ? mapping.fill_dp_mismatch : mapping.cigar_dp_mismatch),
      dp_src(opt.dp_mismatch));
  row("dp_score_n", i2s(fill ? mapping.fill_dp_ambi : mapping.cigar_dp_ambi),
      dp_src(opt.dp_score_n));
  row("dp_gap_open1",
      i2s(fill ? mapping.fill_dp_gap_open1 : mapping.cigar_dp_gap_open1),
      dp_src(opt.dp_gap_open1));
  row("dp_gap_extend1",
      i2s(fill ? mapping.fill_dp_gap_extend1 : mapping.cigar_dp_gap_extend1),
      dp_src(opt.dp_gap_extend1));
  row("dp_gap_open2",
      i2s(fill ? mapping.fill_dp_gap_open2 : mapping.cigar_dp_gap_open2),
      dp_src(opt.dp_gap_open2));
  row("dp_gap_extend2",
      i2s(fill ? mapping.fill_dp_gap_extend2 : mapping.cigar_dp_gap_extend2),
      dp_src(opt.dp_gap_extend2));
  row("dp_zdrop",
      i2s(fill ? mapping.fill_dp_tail_zdrop : mapping.cigar_dp_tail_zdrop),
      dp_src(opt.dp_zdrop));
  row("dp_end_bonus", i2s(mapping.cigar_dp_tail_end_bonus),
      dp_src(opt.dp_end_bonus));
  // On a splice preset dp_bw is the fine chain's band, and -G sets both.
  const bool intron_typed = is_rna && opt.max_intron.has_value();
  row("dp_bw", i2s(mapping.cigar_dp_bw),
      src(opt.dp_bw.has_value() || intron_typed, "preset"));
  row("dp_bw_long", i2s(mapping.cigar_dp_bw_long),
      src(opt.dp_bw_long.has_value() || intron_typed, "preset"));
  row("dp_max_gap", i2s(mapping.cigar_dp_max_gap), dp_src(opt.dp_max_gap));
  row("dp_min_score", i2s(mapping.cigar_dp_min_dp_max),
      dp_src(opt.dp_min_score));
  row("dp_split_min_anchors", i2s(mapping.cigar_dp_split_min_anchors),
      "builtin");
  // A splice preset's inversion Z-drop is splice_inv_zdrop.
  if (fill)
    row("dp_inversion_zdrop", i2s(mapping.fill_dp_inversion_zdrop),
        dp_src(opt.dp_zdrop));
  if (fill && !assembly)
    row("dp_end_row",
        "A" + i2s(mapping.cigar_dp_match) + " B" +
            i2s(mapping.cigar_dp_mismatch) + " N" +
            i2s(mapping.cigar_dp_ambi) + " O" +
            i2s(mapping.cigar_dp_gap_open1) + "," +
            i2s(mapping.cigar_dp_gap_open2) + " E" +
            i2s(mapping.cigar_dp_gap_extend1) + "," +
            i2s(mapping.cigar_dp_gap_extend2) + " z" +
            i2s(mapping.cigar_dp_tail_zdrop) + "," +
            i2s(mapping.cigar_dp_inversion_zdrop),
        "preset");
  if (fill) {
    // A value that is neither typed nor the struct's own is the preset row's.
    const int builtin =
        std::decay_t<decltype(mapping)>{}.min_chain_score;
    row("min_chain_score", i2s(mapping.min_chain_score),
        src(opt.min_chain_score.has_value(),
            mapping.min_chain_score != builtin ? "preset" : "builtin"));
  }

  if (!is_rna) {
    row("dense_diag_min_runs", i2s(mapping.dna_dense_diag_min_runs), "builtin");
    row("tandem_window", i2s(mapping.dna_tandem_window), "builtin");
    row("dp_min_ksw_len", i2s(mapping.cigar_dp_min_ksw_len), "builtin");
    // Unless --max-chain-occ sets it, the pool gate is max_seed_occ, printed
    // by name while that row is still a rule.
    if (mapping.dna_pool_gate_occ == fa::cpu::lr::kDnaPoolGateVoteCap) {
      row("pool_gate_occ",
          occ_rule && !occ_index.resolved
              ? std::string("max_seed_occ")
              : i2s(fa::cpu::lr::dna_chain_occ_thresholds(mapping)
                        .pool_gate_occ),
          "derived");
    } else {
      row("pool_gate_occ", i2s(mapping.dna_pool_gate_occ),
          src(opt.max_chain_occ.has_value(), "builtin"));
    }
  }

  if (is_rna && rna != nullptr) {
    row("rna_min_intron", i2s(rna->min_intron),
        src(opt.min_intron && *opt.min_intron > 0, "preset"));
    row("rna_max_intron", i2s(rna->max_intron),
        src(opt.max_intron && *opt.max_intron > 0, "preset"));
    row("rna_junction_bed",
        rna->junction_bed.empty() ? "off" : rna->junction_bed,
        src(opt.rna_junction_bed.has_value(), "builtin"));
    row("rna_junction_bonus", i2s(rna->splice_junction_bonus),
        src(opt.rna_junction_bonus.has_value(), "preset"));
    const char* strand = "unset";
    switch (rna->strand_mode) {
    case 0:
      strand = "auto";
      break;
    case 1:
      strand = "forward";
      break;
    case 2:
      strand = "reverse";
      break;
    case 3:
      strand = "none";
      break;
    default:
      break;
    }
    row("rna_strand_mode", strand,
        src(opt.splice_strand.has_value(), "builtin"));
    row("splice_transition", i2s(rna->splice_transition), "preset");
    // -z's second value; as in minimap2, a scalar -z sets it too.
    row("splice_inv_zdrop", i2s(rna->splice_inversion_zdrop),
        dp_src(opt.dp_zdrop));
    row("rna_max_loci", i2s(rna->max_locus_chains),
        src(opt.rna_max_loci.has_value(), "builtin"));
  }

  if (!is_rna) {
    // HiFi only: the chain MAPQ margin rule (dna/chain_mapq.h).
    row("dna_chain_mapq.hifi_margin", b2s(mapping.chain_mapq_hifi_margin),
        "builtin");
    // HiFi only: the late inversion probe's local-chain gate
    // (dna/inv_local_chain.h).
    row("dna_inv_probe.local_gate", b2s(mapping.inversion_probe_local_gate),
        "builtin");
    row("dna_rival.pri_ratio", f2s(mapping.pri_ratio),
        src(opt.pri_ratio.has_value(), "builtin"));
    // -N: the alternatives whose whole-query chains enter the ownership
    // selection. Only an installed count prints, so the default text and
    // config_digest stay.
    if (opt.dna_alternative_realize_max)
      row("dna_alternative.realize_max", i2s(mapping.alternative_realize_max),
          "explicit");
  } else if (rna != nullptr) {
    row("rna_rival.pri_ratio", f2s(rna->rival_pri_ratio),
        src(opt.pri_ratio.has_value(), "builtin"));
    row("rna_rival.min_diff", i2s(rna->rival_min_diff),
        src(opt.rna_rival_min_diff.has_value(), "builtin"));
    // Fixed thresholds are printed too, tagged [fixed]. First the
    // co-primary floor, then the chimeric emission floor.
    row("rna_rival.chim_score",
        i2s(fa::cpu::lr::rna::kRnaCoprimaryMinChainScore), "fixed");
    row("rna_rival.chim_anch", i2s(fa::cpu::lr::rna::kRnaCoprimaryMinAnchors),
        "fixed");
    row("rna_rival.chim_qbp", i2s(fa::cpu::lr::rna::kRnaCoprimaryMinQueryBases),
        "fixed");
    row("rna_rival.chimera_max_families",
        i2s(fa::cpu::lr::rna::kRnaChimeraMaxFamilies), "fixed");
    row("rna_rival.chimera_overlap",
        i2s(fa::cpu::lr::rna::kRnaChimeraOverlapNumerator) + "/" +
            i2s(fa::cpu::lr::rna::kRnaChimeraOverlapDenominator),
        "fixed");
    row("rna_rival.chimera_min_qbp",
        i2s(fa::cpu::lr::rna::kRnaChimeraMinQueryBases), "fixed");
    row("rna_rival.realize_max", i2s(rna->rival_realize_max),
        src(opt.rna_max_loci.has_value(), "builtin"));
    row("rna_mapq.qcov_tau", f2s(rna->rna_mapq_qcov_tau), "preset");
    row("rna_mapq.vote_damp", f2s(rna->rna_mapq_disjoint_vote_damp), "preset");
  }

  row("full_read_cigar", b2s(cfg.common.enable_full_read_cigar), "cli");
  const char* cs_mode = "none";
  if (cfg.common.cs == fa::cpu::options::CsMode::Short)
    cs_mode = "short";
  else if (cfg.common.cs == fa::cpu::options::CsMode::Long)
    cs_mode = "long";
  row("cs", cs_mode, src(!opt.cs.empty(), "builtin"));
  row("emit_md", b2s(cfg.common.emit_md), src(opt.emit_md, "builtin"));
  row("emit_eqx", b2s(cfg.common.emit_eqx), src(opt.emit_eqx, "builtin"));
  return os.str();
}

} // namespace

void show_config(const AlignOptions& opt,
                 const fa::cpu::options::ResolvedOptions& cfg,
                 fa::cpu::api::BackendMode mode,
                 std::optional<std::size_t> reference_sequences,
                 std::ostream& out) {
  out << "[flashalign] resolved configuration"
         " (value [source]; source in"
         " {builtin,preset,index,explicit,cli,derived,fixed,resolved})\n";
  out << resolved_config_text(opt, cfg, mode);
  if (reference_sequences) {
    out << "  " << std::left << std::setw(24) << "reference_sequences"
        << *reference_sequences << "\n";
  }
}

} // namespace fa::cpu::cli
