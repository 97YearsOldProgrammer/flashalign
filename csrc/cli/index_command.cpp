#include "cli/index_command.h"
#include "cli/errors.h"

#include "cli/parse.h"

#include "api/aligner.h"
#include "core/types.h"  // encode_sequence_u8
#include "index/index.h"
#include "index/faix_multipart_writer.h"
#include "io/fastx.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace fa::cpu::cli {

namespace {

struct EncodedReferenceLoad {
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t>> sequences;
};

struct ReferenceMetadata {
    std::vector<std::string> names;
    std::vector<uint64_t> offsets;
    std::unordered_map<std::string, size_t> index_by_name;

    uint64_t total_bp() const {
        return offsets.empty() ? 0 : offsets.back();
    }
};

ReferenceMetadata load_fasta_reference_metadata_sorted(const std::string& path) {
    fa::cpu::io::FastxReader reader(path);
    fa::cpu::io::FastxRecord rec;
    std::vector<std::pair<std::string, uint64_t>> records;
    while (reader.next(rec)) {
        if (rec.seq.empty()) continue;
        records.push_back({
            std::move(rec.name),
            static_cast<uint64_t>(rec.seq.size())});
    }
    if (records.empty()) {
        throw std::runtime_error("reference FASTA has no sequences: " + path);
    }

    std::sort(records.begin(), records.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });

    ReferenceMetadata out;
    out.names.reserve(records.size());
    out.offsets.reserve(records.size() + 1);
    out.index_by_name.reserve(records.size());
    out.offsets.push_back(0);
    for (size_t i = 0; i < records.size(); ++i) {
        if (i > 0 && records[i - 1].first == records[i].first) {
            throw std::runtime_error(
                "duplicate reference sequence name: " + records[i].first);
        }
        out.offsets.push_back(out.offsets.back() + records[i].second);
        out.names.push_back(std::move(records[i].first));
        out.index_by_name.emplace(out.names.back(), i);
    }
    if (out.names.size() > fa::cpu::kFaixMaxContigCount) {
      throw std::runtime_error(
          "the reference has " + std::to_string(out.names.size()) +
          " sequences; the most an index holds is " +
          std::to_string(fa::cpu::kFaixMaxContigCount));
    }
    for (size_t i = 0; i + 1 < out.offsets.size(); ++i) {
      const uint64_t length = out.offsets[i + 1] - out.offsets[i];
      if (length > fa::cpu::kFaixMaxContigBp) {
        throw std::runtime_error("reference sequence " + out.names[i] + " is " +
                                 std::to_string(length) + " bp, longer than the " +
                                 std::to_string(fa::cpu::kFaixMaxContigBp) +
                                 " bp SAM and BAM can address");
      }
    }
    return out;
}

EncodedReferenceLoad load_fasta_genome_encoded_sorted(const std::string& path) {
    fa::cpu::io::FastxReader reader(path);
    fa::cpu::io::FastxRecord rec;
    std::vector<std::pair<std::string, std::vector<uint8_t>>> records;
    while (reader.next(rec)) {
        if (rec.seq.empty()) continue;
        records.push_back({
            std::move(rec.name),
            fa::cpu::encode_sequence_u8(rec.seq)});
    }
    if (records.empty()) {
        throw std::runtime_error("reference FASTA has no sequences: " + path);
    }

    std::sort(records.begin(), records.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });

    EncodedReferenceLoad out;
    out.names.reserve(records.size());
    out.sequences.reserve(records.size());
    for (auto& item : records) {
        out.names.push_back(std::move(item.first));
        out.sequences.push_back(std::move(item.second));
    }
    return out;
}

// -I: greedily cut the name-sorted contigs into parts of at most `budget`
// bases. A contig longer than the budget is a part of its own, never split;
// its index is added to `oversized`.
std::vector<fa::cpu::FaixPartEntry> partition_by_bases(
    const ReferenceMetadata& reference, uint64_t budget,
    std::vector<size_t>& oversized) {
    std::vector<fa::cpu::FaixPartEntry> parts;
    const size_t n = reference.names.size();
    size_t first = 0;
    uint64_t sum = 0;
    for (size_t i = 0; i < n; ++i) {
        const uint64_t len = reference.offsets[i + 1] - reference.offsets[i];
        if (i > first && sum + len > budget) {
            fa::cpu::FaixPartEntry part;
            part.first_contig = first;
            part.contig_count = i - first;
            parts.push_back(part);
            first = i;
            sum = 0;
        }
        if (len > budget) oversized.push_back(i);
        sum += len;
    }
    fa::cpu::FaixPartEntry last;
    last.first_contig = first;
    last.contig_count = n - first;
    parts.push_back(last);
    return parts;
}

// Builds the index of contigs [first, first + count) in one pass over the
// FASTA, with names and, unless --idx-no-seq, the 4-bit packed reference. The
// part is numbered and offset from 0, so it is identical to a single-part index
// of those contigs alone.
fa::cpu::SeedIndex build_reference_run(
    const IndexOptions& opt, const ReferenceMetadata& reference,
    size_t first, size_t count, const fa::cpu::FaixBuildConfig& cfg) {
    std::vector<uint64_t> offsets(count + 1, 0);
    for (size_t j = 0; j < count; ++j) {
        offsets[j + 1] = offsets[j] + (reference.offsets[first + j + 1] -
                                       reference.offsets[first + j]);
    }
    std::vector<uint8_t> packed_reference(
        opt.no_seq ? 0
                   : static_cast<size_t>(
                         fa::cpu::packed_reference_bytes(offsets.back())),
        0);
    fa::cpu::io::FastxReader reader(opt.ref_path);
    std::mutex reader_mutex;
    std::mutex packed_reference_mutex;
    size_t remaining = count;

    fa::cpu::FaixBuildStatus build_status;
    fa::cpu::SeedIndex index = fa::cpu::build_seed_index_streamed(
        offsets,
        cfg,
        [&](size_t& chr_idx, std::vector<uint8_t>& enc) -> bool {
            fa::cpu::io::FastxRecord rec;
            for (;;) {
                {
                    std::lock_guard<std::mutex> lk(reader_mutex);
                    if (remaining == 0) return false;
                    if (!reader.next(rec)) return false;
                }
                if (rec.seq.empty()) continue;
                auto it = reference.index_by_name.find(rec.name);
                if (it == reference.index_by_name.end()) {
                    throw std::runtime_error(
                        "reference sequence disappeared during second pass: "
                        + rec.name);
                }
                const size_t global = it->second;
                if (global < first || global >= first + count) continue;
                chr_idx = global - first;
                const uint64_t expected_len =
                    reference.offsets[global + 1] - reference.offsets[global];
                if (expected_len != static_cast<uint64_t>(rec.seq.size())) {
                    throw std::runtime_error(
                        "reference sequence length changed during second pass: "
                        + rec.name);
                }
                enc = fa::cpu::encode_sequence_u8(rec.seq);
                if (!opt.no_seq) {
                    std::lock_guard<std::mutex> lk(packed_reference_mutex);
                    fa::cpu::pack_reference_4bit_into(
                        enc.data(),
                        static_cast<uint64_t>(enc.size()),
                        offsets[chr_idx],
                        packed_reference);
                }
                {
                    std::lock_guard<std::mutex> lk(reader_mutex);
                    --remaining;
                }
                return true;
            }
        },
        &build_status);
    if (index.empty()) {
        throw std::runtime_error(
            "failed to build index: " +
            fa::cpu::faix_build_error_message(build_status));
    }
    std::vector<std::string> names(reference.names.begin() + first,
                                   reference.names.begin() + first + count);
    if (opt.no_seq) {
        if (!index.set_chromosome_names(std::move(names))) {
            throw std::runtime_error("failed to attach the reference names");
        }
    } else if (!index.set_packed_reference_payload(
                   std::move(names), std::move(packed_reference))) {
        throw std::runtime_error("failed to attach the reference to the index");
    }
    return index;
}

}  // namespace

int run_index(const IndexOptions& opt) {
    // The same seeding resolver as align, so the index matches the preset.
    const fa::cpu::api::PresetSeeding seeding =
        fa::cpu::api::resolve_preset_seeding(opt.preset);
    // 1 <= s <= k; s == k is legal, like minimap2's -w 1.
    if (opt.k && *opt.k < 1) {
        throw UsageError("-k must be at least 1");
    }
    const int build_k = opt.k.value_or(seeding.k);
    if (opt.syncmer_s && (*opt.syncmer_s < 1 || *opt.syncmer_s > build_k)) {
        throw UsageError(
            "-s must be within [1,k] (k=" + std::to_string(build_k) + ")");
    }
    const int chain_syncmer_s =
        opt.syncmer_s ? *opt.syncmer_s : seeding.syncmer_s;

    if (build_k > fa::cpu::kFaixMaxK) {
      throw UsageError(
          "k=" + std::to_string(build_k) + " is too large; the maximum is " +
          std::to_string(fa::cpu::kFaixMaxK));
    }

    if (fa::cpu::api::reference_kind(opt.ref_path) ==
        fa::cpu::api::ReferenceKind::Index) {
        throw std::runtime_error(
            opt.ref_path + " is already a FlashAlign index; `flashalign index` "
            "takes a FASTA or FASTQ file");
    }

    fa::cpu::FaixBuildConfig cfg;
    cfg.k = build_k;
    cfg.syncmer_s = chain_syncmer_s;
    cfg.build_threads = opt.threads.value_or(0);
    // Recorded in the header; align adopts it.
    cfg.preset = opt.preset;

    size_t part_count = 1;
    if (opt.ref_path == "-") {
        // stdin: the parser has already refused -I.
        auto reference = load_fasta_genome_encoded_sorted(opt.ref_path);
        // An empty index means the build failed; build_status says why.
        fa::cpu::FaixBuildStatus build_status;
        fa::cpu::SeedIndex index = fa::cpu::build_seed_index(
            reference.sequences, cfg, &build_status);
        if (index.empty()) {
            throw std::runtime_error(
                "failed to build index: " +
                fa::cpu::faix_build_error_message(build_status));
        }
        const bool attached =
            opt.no_seq
                ? index.set_chromosome_names(std::move(reference.names))
                : index.set_reference_payload(std::move(reference.names),
                                              reference.sequences);
        std::vector<std::vector<uint8_t>>().swap(reference.sequences);
        if (!attached) {
            throw std::runtime_error("failed to attach the reference to the index");
        }
        if (!index.save(opt.out_path))
            throw std::runtime_error("failed to save index: " + opt.out_path);
    } else {
        ReferenceMetadata reference = load_fasta_reference_metadata_sorted(opt.ref_path);
        const size_t chrom_count = reference.names.size();
        // Without -I, or when the reference fits the budget, there is one
        // part and a single-part index is written.
        std::vector<fa::cpu::FaixPartEntry> parts;
        std::vector<size_t> oversized;
        if (opt.batch_bp) {
            parts = partition_by_bases(
                reference, static_cast<uint64_t>(*opt.batch_bp), oversized);
        } else {
            fa::cpu::FaixPartEntry whole;
            whole.first_contig = 0;
            whole.contig_count = chrom_count;
            parts.push_back(whole);
        }
        part_count = parts.size();
        if (!opt.quiet) {
            for (const size_t i : oversized) {
                std::fprintf(stderr,
                    "[flashalign] note: reference sequence %s is %llu bp, "
                    "above the -I budget of %lld; it forms a part of its own\n",
                    reference.names[i].c_str(),
                    static_cast<unsigned long long>(
                        reference.offsets[i + 1] - reference.offsets[i]),
                    static_cast<long long>(*opt.batch_bp));
            }
            if (opt.batch_bp && parts.size() == 1) {
                std::fprintf(stderr,
                    "[flashalign] note: -I %lld: the whole reference (%llu bp) "
                    "fits one part; writing a single-part index\n",
                    static_cast<long long>(*opt.batch_bp),
                    static_cast<unsigned long long>(reference.total_bp()));
            }
        }
        if (parts.size() == 1) {
            fa::cpu::SeedIndex index =
                build_reference_run(opt, reference, 0, chrom_count, cfg);
            if (!index.save(opt.out_path))
                throw std::runtime_error("failed to save index: " + opt.out_path);
        } else {
            // Multi-part: header and global tables first, then each part is
            // built, written and freed in turn, so peak memory is one part's.
#if defined(__GLIBC__)
            // Pin glibc's mmap threshold at its default so large transient
            // buffers are unmapped on free instead of fragmenting the arenas
            // as glibc's dynamic threshold grows.
            mallopt(M_MMAP_THRESHOLD, 128 * 1024);
#endif
            fa::cpu::FaixMultipartWriter::Plan plan;
            plan.k = build_k;
            plan.syncmer_s = fa::cpu::htable_effective_closed_syncmer_s(cfg);
            plan.has_reference = !opt.no_seq;
            plan.preset = cfg.preset;
            plan.chr_offsets = reference.offsets;
            plan.names = reference.names;
            plan.parts = parts;
            fa::cpu::FaixMultipartWriter writer(opt.out_path, std::move(plan));
            if (!writer.opened() || !writer.begin()) {
                throw std::runtime_error("failed to save index: " + opt.out_path);
            }
            for (size_t p = 0; p < parts.size(); ++p) {
                const fa::cpu::FaixPartEntry& part = parts[p];
                const size_t first = static_cast<size_t>(part.first_contig);
                const size_t count = static_cast<size_t>(part.contig_count);
                if (!opt.quiet) {
                    std::fprintf(stderr,
                        "[flashalign] index part %zu/%zu: %zu sequences, "
                        "%llu bp (%s .. %s)\n",
                        p + 1, parts.size(), count,
                        static_cast<unsigned long long>(
                            reference.offsets[first + count] -
                            reference.offsets[first]),
                        reference.names[first].c_str(),
                        reference.names[first + count - 1].c_str());
                }
                {
                    fa::cpu::SeedIndex index = build_reference_run(
                        opt, reference, first, count, cfg);
                    if (!writer.append_part(index)) {
                        throw std::runtime_error(
                            "failed to save index part " +
                            std::to_string(p + 1) + " of " +
                            std::to_string(parts.size()) + ": " +
                            opt.out_path);
                    }
                }
#if defined(__GLIBC__)
                // Return the freed part to the kernel before the next one.
                malloc_trim(0);
#endif
            }
            if (!writer.finish()) {
                throw std::runtime_error("failed to save index: " + opt.out_path);
            }
        }
    }
    if (!opt.quiet) {
        if (part_count > 1) {
            std::cerr << "[flashalign] -I " << *opt.batch_bp
                      << ": multi-part index, " << part_count
                      << " parts; 'align' maps one part at a time\n";
        }
        if (opt.no_seq) {
            std::cerr << "[flashalign] --idx-no-seq: no reference sequence "
                         "stored; this index maps to plain PAF only\n";
        }
    }
    return 0;
}

}  // namespace fa::cpu::cli
