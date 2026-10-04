#include "cli/align_command.h"

#include "cli/parse.h"
#include "cli/output_writer.h"
#include "cli/report.h"

#include "api/aligner.h"
#include "core/types.h"
#include "index/index.h"
#include "dp/ksw2_align.h"
#include "rna/config.h"
#include "seeding/tie_hash.h"
#include "io/fastx.h"
#include "io/reads.h"
#include "io/batch_reader.h"
#include "threading/parallel_for.h"
#include "threading/pipeline.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <sys/stat.h>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace fa::cpu::cli {

namespace {

// Lends idle decode-team threads to the compute window. The read source (and
// its decode team) starts before the compute window exists, so the team calls
// through this hub, which answers "no work" until install() and after
// revoke(). revoke() waits for loaned threads to leave the target, so the
// window can be destroyed while the team still runs; it is called only after
// the pipeline has joined, so no new work can keep it waiting.
class HelperHub {
public:
  bool call(int tid) {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return target_ ? target_(tid) : false;
  }
  void install(std::function<bool(int)> target) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    target_ = std::move(target);
  }
  void revoke() {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    target_ = nullptr;
  }

private:
  std::shared_mutex mutex_;
  std::function<bool(int)> target_;
};

// The smallest -t at which gzip decoding borrows mapping threads.
constexpr int kDecodeLoanMinThreads = 32;

// Threads the decode team borrows from the -t budget: -t/8, clamped to [1, 8].
int decode_loan_threads(int num_threads) {
  const int wanted = num_threads / 8;
  if (wanted < 1)
    return 1;
  return wanted > 8 ? 8 : wanted;
}

// A loan needs every input to be eligible for member-parallel decoding;
// otherwise the loaned threads would sit idle while the serial reader runs.
bool decode_loan_can_be_fed(const std::vector<std::string>& reads_paths) {
  if (reads_paths.empty())
    return false;
  for (const std::string& path : reads_paths) {
    if (!::fa::cpu::io::gz_member_parallel_eligible(path))
      return false;
  }
  return true;
}

// Decompressed bytes the decoder may hold ahead of the reader: 64 MiB per
// thread, clamped to [512 MiB, 4 GiB], unless --io-staging sets it.
uint64_t decode_staging_budget(int num_threads, int io_staging_mib) {
  if (io_staging_mib > 0)
    return static_cast<uint64_t>(io_staging_mib) << 20;
  const uint64_t scaled =
      static_cast<uint64_t>(num_threads < 1 ? 1 : num_threads) * (64ull << 20);
  if (scaled < (512ull << 20))
    return 512ull << 20;
  return scaled > (4ull << 30) ? (4ull << 30) : scaled;
}

// Mapping threads: -t, or all cores but two (for the reader and writer).
int resolve_num_threads(const AlignOptions& opt) {
  return (opt.threads && *opt.threads > 0)
             ? *opt.threads
             : std::max(1, ::fa::cpu::threading::default_thread_count() - 2);
}

// Batches in flight in the compute stage. Forced to 1 at -t 1, so that run
// uses a single OS thread (run_pipeline_serial).
int resolve_batch_window(const AlignOptions& opt, int num_threads) {
  int window = opt.batch_window < 1 ? 1 : opt.batch_window;
  if (num_threads == 1)
    window = 1;
  return window;
}

// The input decode plan, minus the loan callback. The parallel gzip decoder
// borrows threads only with a compute window, inputs that can feed it, and
// -t >= kDecodeLoanMinThreads: below that, serial decoding on the reader
// thread hides behind compute and lending threads is a net loss. The loan is
// at least four threads, so one thread stuck on a long read cannot starve
// the reader.
fa::cpu::io::IoPlan build_io_plan(const AlignOptions& opt, int num_threads,
                                  int batch_window) {
  fa::cpu::io::IoPlan plan;
  plan.parallel_gz = batch_window >= 2 &&
                     num_threads >= kDecodeLoanMinThreads &&
                     decode_loan_can_be_fed(opt.reads_paths);
  if (!plan.parallel_gz)
    return plan;
  plan.gz.workers = static_cast<unsigned>(decode_loan_threads(num_threads));
  plan.gz.staging_budget =
      decode_staging_budget(num_threads, opt.io_staging_mib);
  // The source raises a budget below twice this cap, so keep the cap at half
  // the budget.
  plan.gz.segment_stage_cap =
      std::min<uint64_t>(512ull << 20, plan.gz.staging_budget / 2);
  return plan;
}

// SAM always carries a CIGAR; PAF only with -c (or --cs/--MD).
bool run_needs_cigar(const AlignOptions& opt) {
  return opt.format == "sam" || opt.paf_cigar;
}

// DNA map-only reads no reference bases, so its index load skips unpacking
// the reference. CIGAR realization and RNA need them.
bool run_needs_reference_bases(const AlignOptions& opt) {
  return run_needs_cigar(opt) || fa::cpu::api::is_rna_preset(opt.preset);
}

// Parsed CLI values as typed overrides, plus the CLI's preset-dependent
// runtime defaults.
fa::cpu::api::UserOverrides build_user_overrides(const AlignOptions& opt) {
  fa::cpu::api::UserOverrides u;
  const bool need_cigar = run_needs_cigar(opt);

  // (k, s) is not set here: it comes from the .faix or the preset.
  u.num_threads = resolve_num_threads(opt);
  u.enable_full_read_cigar = need_cigar;
  if (!opt.cs.empty()) {
    u.cs = opt.cs == "long" ? fa::cpu::options::CsMode::Long
                            : fa::cpu::options::CsMode::Short;
  }
  if (opt.emit_md)
    u.emit_md = true;
  if (opt.emit_eqx)
    u.emit_eqx = true;
  if (opt.dna_vote_admission_ratio)
    u.dna_vote_admission_ratio = *opt.dna_vote_admission_ratio;
  if (opt.max_chain_occ)
    u.dna_pool_gate_occ = *opt.max_chain_occ;
  if (opt.tile_supported_reward)
    u.query_tile_supported_reward = *opt.tile_supported_reward;
  if (opt.tile_block_open_cost)
    u.query_tile_block_open_cost = *opt.tile_block_open_cost;
  if (opt.tile_null_cost)
    u.query_tile_null_cost = *opt.tile_null_cost;
  if (opt.tile_unsupported_cost)
    u.query_tile_unsupported_cost = *opt.tile_unsupported_cost;
  if (opt.min_support)
    u.min_support = *opt.min_support;
  if (opt.vote_seeds)
    u.max_query_seeds_per_strand = *opt.vote_seeds;
  if (opt.max_cands)
    u.max_cands = *opt.max_cands;
  if (opt.tiles)
    u.query_tiles = *opt.tiles;
  if (opt.tile_owner)
    u.tile_owner_anchors = *opt.tile_owner == "anchors";
  if (opt.min_chain_score)
    u.min_chain_score = *opt.min_chain_score;
  // --max-vote-occ: 0 turns occurrence filtering off, >0 fixes the cap.
  if (opt.max_vote_occ) {
    if (*opt.max_vote_occ <= 0) {
      u.occ_policy = std::string("off");
    } else {
      u.occ_policy = std::string("fixed");
      u.primary_occ_cap = *opt.max_vote_occ;
    }
  }
  if (opt.vote_diag_bin_width)
    u.vote_diag_bin_width = *opt.vote_diag_bin_width;
  if (opt.vote_diag_slope_den)
    u.vote_diag_slope_den = *opt.vote_diag_slope_den;
  if (opt.vote_diag_width_max)
    u.vote_diag_width_max = *opt.vote_diag_width_max;
  if (opt.dp_match)
    u.dp_match = *opt.dp_match;
  if (opt.dp_mismatch)
    u.dp_mismatch = *opt.dp_mismatch;
  if (opt.dp_score_n)
    u.dp_ambi = *opt.dp_score_n;
  if (opt.dp_gap_open1)
    u.dp_gap_open1 = *opt.dp_gap_open1;
  if (opt.dp_gap_open2)
    u.dp_gap_open2 = *opt.dp_gap_open2;
  if (opt.dp_gap_extend1)
    u.dp_gap_extend1 = *opt.dp_gap_extend1;
  if (opt.dp_gap_extend2)
    u.dp_gap_extend2 = *opt.dp_gap_extend2;
  if (opt.dp_zdrop)
    u.dp_tail_zdrop = *opt.dp_zdrop;
  if (opt.dp_zdrop_inv)
    u.dp_inversion_zdrop = *opt.dp_zdrop_inv;
  if (opt.dp_end_bonus)
    u.dp_tail_end_bonus = *opt.dp_end_bonus;
  if (opt.dp_min_score)
    u.dp_min_dp_max = *opt.dp_min_score;
  if (opt.dp_bw)
    u.dp_bw = *opt.dp_bw;
  if (opt.dp_bw_long)
    u.dp_bw_long = *opt.dp_bw_long;
  if (opt.dp_max_gap)
    u.dp_max_gap = *opt.dp_max_gap;
  // Passed as given; the resolver validates the intron bounds.
  if (opt.min_intron)
    u.rna_min_intron = *opt.min_intron;
  if (opt.max_intron)
    u.rna_max_intron = *opt.max_intron;
  if (opt.rna_junction_bed)
    u.rna_junction_bed = *opt.rna_junction_bed;
  if (opt.rna_junction_bonus)
    u.rna_junction_bonus = *opt.rna_junction_bonus;
  if (opt.pri_ratio)
    u.pri_ratio = *opt.pri_ratio;
  if (opt.rna_rival_min_diff)
    u.rna_rival_min_diff = *opt.rna_rival_min_diff;
  if (opt.rna_max_loci)
    u.rna_max_loci = *opt.rna_max_loci;
  if (opt.dna_alternative_realize_max)
    u.dna_alternative_realize_max = *opt.dna_alternative_realize_max;
  // "auto" (-u b) maps to Unknown, the preset default.
  if (opt.splice_strand) {
    u.rna_strand = static_cast<int>(fa::cpu::lr::rna::parse_strand_mode(
        opt.splice_strand->c_str(), fa::cpu::lr::rna::StrandMode::Unknown));
  }
  return u;
}

// The align target, told apart by its content (api::reference_kind). Of an index only the
// header is read here.
struct ResolvedInputs {
  std::string ref_path;   // FASTA to load; empty when the target is an index
  std::string index_path; // .faix to load
  bool use_index = false; // load index_path instead of building from FASTA
  bool embedded_ref = false; // the reference comes from the index: no FASTA
  // False for an --idx-no-seq index, which supports plain PAF only.
  bool index_has_reference = true;
  // The preset the .faix records, or empty. Used when no -x is given.
  std::string index_preset;
  // The index's seeding and sequence count, from its header.
  int index_k = 0;
  int index_s = 0;
  std::size_t reference_sequences = 0;
  // A multi-part (-I) index, which run_align maps one part at a time.
  bool multipart = false;
};

ResolvedInputs resolve_inputs(const AlignOptions& opt) {
  ResolvedInputs in;
  // Only --show-config runs without a target.
  if (opt.target_path.empty())
    return in;
  if (fa::cpu::api::reference_kind(opt.target_path) ==
      fa::cpu::api::ReferenceKind::Sequences) {
    in.ref_path = opt.target_path;
    return in;
  }
  const fa::cpu::api::IndexProbe probe =
      fa::cpu::api::LongReadAligner::probe_index_header(opt.target_path);
  if (!probe.error.empty())
    throw std::runtime_error(opt.target_path + ": " + probe.error);
  // An --idx-no-seq index still serves plain PAF, which needs only contig
  // names and lengths.
  if (!probe.has_reference && run_needs_cigar(opt)) {
    throw std::runtime_error(
        opt.target_path +
        ": this .faix does not embed its reference (built with --idx-no-seq), "
        "so it cannot realize a CIGAR; map to plain PAF (without "
        "-a/-c/--cs/--MD), or rebuild it with 'flashalign index'");
  }
  in.index_path = opt.target_path;
  in.use_index = true;
  in.embedded_ref = true;
  in.index_has_reference = probe.has_reference;
  in.index_preset = probe.preset;
  in.index_k = probe.k;
  in.index_s = probe.syncmer_s;
  in.reference_sequences = static_cast<std::size_t>(probe.reference_sequences);
  in.multipart = probe.multipart;
  return in;
}

// Refuses SAM text given as a reads file, which kseq would otherwise parse as
// nonsense FASTQ. Only regular uncompressed files are probed: stdin, FIFOs and
// gzip inputs are left alone.
void refuse_sam_text_reads(const std::vector<std::string>& paths) {
  for (const std::string& path : paths) {
    if (path == "-")
      continue;
    struct stat st {};
    if (::stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode))
      continue; // left for the reader to report
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr)
      continue;
    char head[4] = {0, 0, 0, 0};
    const std::size_t got = std::fread(head, 1, sizeof(head), file);
    std::fclose(file);
    if (got != sizeof(head) || head[0] != '@' || head[3] != '\t')
      continue;
    const std::string tag(head + 1, 2);
    if (tag != "HD" && tag != "SQ" && tag != "RG" && tag != "PG" && tag != "CO")
      continue;
    throw std::runtime_error(
        "reads file '" + path +
        "' is SAM text; flashalign reads FASTA/FASTQ (plain, gzip, bgzf) and "
        "unaligned BAM -- convert with 'samtools fastq'");
  }
}

} // namespace

int run_align(const AlignOptions& cli_options) {
  // Only the preset is rewritten, when a .faix supplies it.
  AlignOptions opt = cli_options;
  refuse_sam_text_reads(opt.reads_paths);
  const ResolvedInputs in = resolve_inputs(opt);
  // A multi-part index re-reads the reads once per part, so stdin cannot
  // be used.
  if (in.multipart) {
    for (const std::string& path : opt.reads_paths) {
      if (path == "-") {
        throw std::runtime_error(
            "reads from stdin ('-') cannot be mapped against a multi-part "
            "index, which reads the input once per part; give the reads as "
            "files");
      }
    }
  }
  // The preset is an explicit -x, else the one the .faix records, else lr.
  // An unknown recorded preset is an error rather than a silent lr.
  if (opt.preset_source == "builtin" && !in.index_preset.empty()) {
    if (!fa::cpu::api::preset_is_valid(in.index_preset)) {
      throw std::runtime_error(
          in.index_path + ": this .faix records preset '" + in.index_preset +
          "', which this build does not know; pass an explicit -x, or rebuild "
          "the index with 'flashalign index'");
    }
    opt.preset = in.index_preset;
    opt.preset_source = "index";
  }
  // RNA mapping needs the reference bases even without a CIGAR. Checked
  // after the preset is settled and before --show-config, so a dry run
  // refuses what a real run would.
  if (in.use_index && !in.index_has_reference &&
      fa::cpu::api::is_rna_preset(opt.preset)) {
    throw std::runtime_error(
        in.index_path +
        ": this .faix does not embed its reference (built with --idx-no-seq) "
        "and the spliced presets map from it; rebuild it with 'flashalign "
        "index -x " + opt.preset + "'");
  }
  // --show-config resolves options only: no sequence is read, no aligner is
  // built, and an index is known by its header alone.
  if (opt.show_config) {
    fa::cpu::api::ResolveRequest req;
    req.preset = opt.preset;
    req.user = build_user_overrides(opt);
    std::optional<size_t> reference_sequences;
    if (in.use_index) {
      req.index.has_index = true;
      req.index.k = in.index_k;
      req.index.syncmer_s = in.index_s;
      reference_sequences = in.reference_sequences;
    }
    const fa::cpu::api::ResolvedMapOptions resolved =
        fa::cpu::api::resolve_options(req);
    show_config(opt, resolved.resolved(), resolved.mode(), reference_sequences,
                std::cout);
    return 0;
  }
  // Taken before the reference or index is loaded, for end-to-end wall time.
  const auto program_t0 = std::chrono::steady_clock::now();
  // A multi-part index: its table names every contig and locates each part's image.
  fa::cpu::FaixMultipartTable table;
  if (in.multipart) {
    fa::cpu::FaixLoadStatus status;
    table = fa::cpu::load_faix_multipart_table(in.index_path, &status);
    if (!table.loaded()) {
      throw std::runtime_error(in.index_path + ": " +
                               fa::cpu::faix_load_error_message(status));
    }
  }
  // A .faix target needs no FASTA: the reference comes from the index.
  const bool use_index = in.use_index;
  const bool embedded_ref = in.embedded_ref;
  fa::cpu::io::GenomeLoad genome;
  if (!embedded_ref) {
    genome = fa::cpu::io::load_fasta_genome(in.ref_path);
  }
  // align_batch maps a batch's sequences, given their read-name hashes.
  std::unique_ptr<fa::cpu::api::LongReadAligner> dna_aligner;
  std::function<std::vector<fa::cpu::AlignResult>(
      const std::vector<std::string>&, const std::vector<std::uint32_t>&)>
      align_batch;

  // Fills the output's reference table when no FASTA was loaded.
  auto fill_refs_from = [&](const std::vector<std::string>& names,
                            const std::vector<int64_t>& lens) {
    if (!genome.references.empty())
      return;
    if (names.size() != lens.size()) {
      throw std::runtime_error(
          "flashalign: index reference names/lengths shape mismatch");
    }
    genome.references.reserve(names.size());
    for (size_t i = 0; i < names.size(); ++i) {
      genome.references.emplace_back(names[i], lens[i]);
    }
  };

  fa::cpu::api::BackendMode backend_mode = fa::cpu::api::BackendMode::DnaLong;

  // Index parts to map against: one, or each part of a multi-part index in
  // turn, as minimap2 does without --split-prefix. Reads are re-read per part
  // and output is part-major; each part is its own engine run, with its own
  // MAPQ rivals and occurrence cap. The @SQ table covers the whole reference.
  std::vector<fa::cpu::api::IndexImageSpan> part_spans(1);
  if (in.multipart) {
    part_spans.clear();
    for (const fa::cpu::FaixPartEntry& part : table.parts)
      part_spans.push_back({part.image_offset, part.image_bytes});
  }

  // Resolved once, on the first attach. With an index, (k, s) is read from
  // it; otherwise the resolver takes it from the preset. Every part of a
  // multi-part index has the same seeding.
  std::optional<fa::cpu::api::ResolvedMapOptions> resolved;
  // Loads part `part` (or builds the in-memory index for a FASTA target),
  // installs the resolved options and rebinds align_batch. The caller frees
  // the previous part first.
  auto attach_part = [&](size_t part) {
    fa::cpu::api::ResolveRequest req;
    req.preset = opt.preset;
    req.user = build_user_overrides(opt);
    if (use_index) {
      // The reference unpacks on -t threads, not on every core.
      dna_aligner = std::make_unique<fa::cpu::api::LongReadAligner>(
          fa::cpu::api::LongReadAligner::from_index(
              std::move(genome.sequences), in.index_path,
              opt.min_support.value_or(3),
              /*chain_syncmer_downsample=*/1,
              req.user.num_threads.value_or(0),
              run_needs_reference_bases(opt), part_spans[part]));
      if (!resolved) {
        req.index.has_index = true;
        req.index.k = dna_aligner->config().index.k;
        req.index.syncmer_s = dna_aligner->get_index_syncmer_s();
        resolved = fa::cpu::api::resolve_options(req);
        // As minimap2, note when the index's (k, s) differ from the preset's.
        if (!opt.quiet) {
          const fa::cpu::api::PresetSeeding seeding =
              fa::cpu::api::resolve_preset_seeding(opt.preset);
          const int got_k = dna_aligner->config().index.k;
          const int got_s = dna_aligner->get_index_syncmer_s();
          if (seeding.k != got_k) {
            std::fprintf(
                stderr,
                "[flashalign] note: prebuilt index k=%d overrides preset "
                "%s (k=%d); rebuild the index to change k\n",
                got_k, opt.preset.c_str(), seeding.k);
          }
          if (seeding.syncmer_s != got_s) {
            std::fprintf(stderr,
                         "[flashalign] note: prebuilt index s=%d overrides "
                         "preset %s (s=%d); rebuild the index to change s\n",
                         got_s, opt.preset.c_str(), seeding.syncmer_s);
          }
        }
      }
    } else {
      req.index.has_index = false;
      resolved = fa::cpu::api::resolve_options(req);
      dna_aligner = std::make_unique<fa::cpu::api::LongReadAligner>(
          std::move(genome.sequences), resolved->k(),
          resolved->min_support(), resolved->syncmer_s(),
          resolved->syncmer_downsample(), req.user.num_threads.value_or(0));
    }
    dna_aligner->reconfigure(resolved->resolved());
    backend_mode = resolved->mode();
    auto* p = dna_aligner.get();
    align_batch = [p](const std::vector<std::string>& seqs,
                      const std::vector<std::uint32_t>& name_hashes) {
      return p->align_batch(seqs, &name_hashes);
    };
  };
  attach_part(0);
  // The output's reference table is always the whole reference's.
  if (in.multipart) {
    std::vector<int64_t> lengths;
    lengths.reserve(table.chrom_names.size());
    for (size_t c = 0; c < table.chrom_names.size(); ++c) {
      lengths.push_back(static_cast<int64_t>(table.chr_offsets[c + 1] -
                                             table.chr_offsets[c]));
    }
    fill_refs_from(table.chrom_names, lengths);
  } else {
    fill_refs_from(dna_aligner->chromosome_names(),
                   dna_aligner->chromosome_lengths());
  }

  AlignmentOutputWriter output(opt.output_path, opt.format, genome,
                               opt.command_line, opt.read_group_line);

  // With a window of 1, each batch is aligned and joined in turn, so it costs
  // its slowest read. With two or more, a worker that runs out of reads in
  // batch k starts on k+1; output order is unchanged. Resolved before the
  // reader, since only a window can borrow decode threads.
  const int num_threads = dna_aligner->config().common.num_threads;
  const int batch_window = resolve_batch_window(opt, num_threads);

  // State shared across passes over the reads (one pass per index part).
  // The reader is rebound per pass; the output context is bound once, after
  // the header. batch_bp bounds buffered input, not alignment cost.
  const int64_t batch_bp = opt.batch_bp;
  std::unique_ptr<fa::cpu::io::BatchReader> batch_reader;
  std::optional<AlignmentOutputContext> output_context;
  size_t current_part = 0;
  // Reads seen in the current pass: written by the reader stage, read by the
  // writer stage.
  std::atomic<int64_t> seen{0};
  // The --stats counters, written by the writer stage only.
  RunSummary summary;
  const auto align_t0 = std::chrono::steady_clock::now();

  struct ReadBatch {
    std::vector<fa::cpu::io::FastxRecord> records;
  };
  struct AlignedBatch {
    std::vector<fa::cpu::io::FastxRecord> records;
    std::vector<fa::cpu::AlignResult> results;
  };

  // The records of each in-flight batch, whose sequences were moved into the
  // window. Between collect_windowed() calls its size equals in_flight(),
  // including after a batch throws. Caller thread only.
  std::deque<std::vector<fa::cpu::io::FastxRecord>> windowed_records;
  // Declared last so its workers are joined first on unwind. One per pass.
  std::unique_ptr<fa::cpu::api::WindowedAlignSession> session;

  auto produce = [&]() -> std::optional<ReadBatch> {
    auto records = batch_reader->next(seen);
    if (!records)
      return std::nullopt;
    ReadBatch batch;
    batch.records = std::move(*records);
    return batch;
  };

  // The compute stage runs on the caller thread while the reader and writer
  // stages overlap it. The batch's sequences are moved into the aligner and
  // back rather than copied. Workers see a read's name only as its X31 hash
  // (seeding/tie_hash.h), which breaks exact ties as minimap2's qname does.
  const auto read_name_hashes =
      [](const std::vector<fa::cpu::io::FastxRecord>& records) {
        std::vector<std::uint32_t> hashes;
        hashes.reserve(records.size());
        for (const fa::cpu::io::FastxRecord& record : records)
          hashes.push_back(fa::cpu::lr::tie_name_hash(record.name));
        return hashes;
      };

  // Window of 1: one barrier per batch.
  auto process_barrier = [&](ReadBatch in) -> std::optional<AlignedBatch> {
    std::vector<std::string> seqs;
    seqs.reserve(in.records.size());
    for (auto& rec : in.records)
      seqs.push_back(std::move(rec.seq));
    AlignedBatch done;
    done.results = align_batch(seqs, read_name_hashes(in.records));
    for (size_t i = 0; i < seqs.size(); ++i)
      in.records[i].seq = std::move(seqs[i]);
    done.records = std::move(in.records);
    return std::make_optional(std::move(done));
  };

  // Window of two or more: takes the oldest batch from the window and puts
  // its sequences back into its records. Caller thread only.
  auto collect_windowed = [&]() -> AlignedBatch {
    // Pop before collect(), which drops its batch before rethrowing that
    // batch's exception, so the two queues stay in step.
    AlignedBatch out;
    out.records = std::move(windowed_records.front());
    windowed_records.pop_front();
    fa::cpu::api::AlignedReadBatch done = session->collect();
    for (size_t i = 0; i < done.reads.size(); ++i)
      out.records[i].seq = std::move(done.reads[i]);
    out.results = std::move(done.results);
    return out;
  };

  // Submits a batch and collects the oldest once the window is full, so
  // in_flight() never exceeds the window.
  auto process_windowed = [&](ReadBatch in) -> std::optional<AlignedBatch> {
    std::vector<std::string> seqs;
    seqs.reserve(in.records.size());
    for (auto& rec : in.records)
      seqs.push_back(std::move(rec.seq));
    std::vector<std::uint32_t> name_hashes = read_name_hashes(in.records);
    windowed_records.push_back(std::move(in.records));
    session->submit(std::move(seqs), std::move(name_hashes));
    if (session->in_flight() < static_cast<size_t>(batch_window))
      return std::nullopt;
    return std::make_optional(collect_windowed());
  };

  // At end of input, drains the window oldest-first.
  auto flush_windowed = [&]() -> std::optional<AlignedBatch> {
    if (session->in_flight() == 0)
      return std::nullopt;
    return std::make_optional(collect_windowed());
  };

  fa::cpu::output::SamEmitOptions emit_opts;
  emit_opts.hard_clip_supp = !opt.soft_clip_supp;
  emit_opts.emit_secondary = opt.output_secondary;

  // Writer stage, on its own thread: counters and record output, in input
  // order. Over a multi-part index the counters sum over passes, except input
  // bases, which are counted on the first pass only.
  auto consume = [&](AlignedBatch batch) {
    const auto& records = batch.records;
    const auto& results = batch.results;
    for (size_t i = 0; i < records.size(); ++i) {
      const auto& result = results[i];
      if (current_part == 0)
        summary.bases += static_cast<int64_t>(records[i].seq.size());
      if (result.mapped()) {
        ++summary.mapped_reads;
        const int64_t introns =
            std::count(result.cigar.begin(), result.cigar.end(), 'N');
        if (introns > 0) {
          ++summary.spliced_reads;
          summary.junctions += introns;
        }
      }
      for (const auto& s : result.supplementary)
        if (s.mapped())
          ++summary.supplementary_records;
      // Secondary records and their supplementaries, written only with
      // --secondary yes.
      if (opt.output_secondary && result.mapped()) {
        for (const auto& secondary : result.secondary) {
          if (!secondary.mapped())
            continue;
          ++summary.secondary_records;
          for (const auto& s : secondary.supplementary)
            if (s.mapped())
              ++summary.secondary_records;
        }
      }
      summary.records +=
          write_alignment_records(*output_context, records[i], result);
    }
    if (opt.progress && !opt.quiet) {
      if (in.multipart) {
        std::cerr << "[flashalign] part " << (current_part + 1) << "/"
                  << part_spans.size() << ": aligned " << seen
                  << " reads, wrote " << summary.records << " records\n";
      } else {
        std::cerr << "[flashalign] aligned " << seen << " reads, wrote "
                  << summary.records << " records\n";
      }
    }
  };

  // One pass over the reads per index part. The previous part is freed
  // before the next is loaded, and the header is written on the first pass.
  for (size_t part = 0; part < part_spans.size(); ++part) {
    current_part = part;
    if (part > 0) {
      dna_aligner.reset();
      attach_part(part);
    }
    if (in.multipart && !opt.quiet) {
      const fa::cpu::FaixPartEntry& entry = table.parts[part];
      std::fprintf(stderr,
                   "[flashalign] part %zu/%zu: %llu sequences, %llu bp\n",
                   part + 1, part_spans.size(),
                   static_cast<unsigned long long>(entry.contig_count),
                   static_cast<unsigned long long>(table.part_bp(part)));
    }
    seen.store(0);

    // The gzip decoder borrows some of the -t mapping threads when they are
    // idle rather than adding its own. Its team starts with the reader, before
    // the window exists, so it reaches the window through the hub, which
    // must outlive the team.
    HelperHub helper_hub;
    fa::cpu::io::IoPlan io_plan = build_io_plan(opt, num_threads, batch_window);
    int decode_loan =
        io_plan.parallel_gz ? static_cast<int>(io_plan.gz.workers) : 0;
    if (decode_loan > 0) {
      // Helper tids are the top of the compute width, which the window
      // reserves for them.
      const int helper_tid_base = num_threads - decode_loan;
      HelperHub* hub = &helper_hub;
      io_plan.gz.helper = [hub, helper_tid_base](unsigned ordinal) {
        return hub->call(helper_tid_base + static_cast<int>(ordinal));
      };
    }

    // Opened before the header is written, so a uBAM's @RG/@CO lines can be
    // carried over.
    fa::cpu::io::SequenceReader reader(opt.reads_paths, &io_plan,
                                       /*keep_comment=*/opt.copy_comment);
    // The source may still have fallen back to serial decoding; then drop
    // the loan so the window is not left short of workers. Only the first
    // input is checked.
    if (decode_loan > 0 && !reader.decode_parallel_active())
      decode_loan = 0;
    if (part == 0) {
      output.write_header(reader.source_header_text(), opt.no_header);
      output_context.emplace(output.context(
          /*include_unmapped=*/opt.format == "paf" ? opt.paf_no_hit
                                                   : !opt.sam_hit_only,
          opt.paf_cigar,
          opt.copy_comment, emit_opts));
    }
    batch_reader = std::make_unique<fa::cpu::io::BatchReader>(reader, batch_bp);

    if (batch_window >= 2)
      session = dna_aligner->make_windowed_session(batch_window, decode_loan);
    // From here the decode team's idle threads map reads.
    if (session != nullptr && decode_loan > 0) {
      fa::cpu::api::WindowedAlignSession* target = session.get();
      helper_hub.install([target](int tid) { return target->try_help(tid); });
    }
    // Revokes the loan before the window is torn down, on unwind too.
    struct HelperLoanGuard {
      HelperHub& hub;
      ~HelperLoanGuard() { hub.revoke(); }
    } helper_loan_guard{helper_hub};

    // Release the persistent parallel_for team left over from unpacking the
    // index. The windowed path never uses it; the barrier path respawns it.
    fa::cpu::threading::join_persistent_pool();

    // Reader, compute and writer run concurrently with two batches per
    // hand-off, plus up to the window in the compute stage. Output order is
    // independent of batching, threads and window. -t 1 runs all three stages
    // on this thread.
    if (batch_window >= 2) {
      fa::cpu::threading::run_pipeline<ReadBatch, AlignedBatch>(
          produce, process_windowed, consume, /*depth=*/2, flush_windowed);
    } else if (num_threads == 1) {
      fa::cpu::threading::run_pipeline_serial<ReadBatch, AlignedBatch>(
          produce, process_barrier, consume);
    } else {
      fa::cpu::threading::run_pipeline<ReadBatch, AlignedBatch>(
          produce, process_barrier, consume, /*depth=*/2);
    }
    // Revoke the loan before the window goes.
    helper_hub.revoke();
    session.reset();
    batch_reader.reset();
  }
  output.close();
  // Zero reads usually means an unrecognized input format.
  if (seen == 0) {
    std::cerr << "flashalign: warning: 0 reads read from input"
                 " (empty file, or an unrecognized format?)\n";
  }
  if (opt.progress && !opt.quiet) {
    std::cerr << "[flashalign] done: aligned " << seen << " reads, wrote "
              << summary.records << " records\n";
  }

  if (opt.stats) {
    const auto done = std::chrono::steady_clock::now();
    summary.preset = opt.preset;
    summary.output_format = opt.format;
    summary.rna = backend_mode == fa::cpu::api::BackendMode::RnaSplice;
    summary.reads = seen.load();
    summary.wall_seconds =
        std::chrono::duration<double>(done - program_t0).count();
    summary.map_seconds =
        std::chrono::duration<double>(done - align_t0).count();
    summary.peak_rss_bytes = peak_rss_bytes();
    write_stats(summary, std::cerr);
  }
  return 0;
}

} // namespace fa::cpu::cli
