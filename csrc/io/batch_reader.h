// Splits the input into batches of about batch_bp bases, in input order. A
// batch ends with the record that reaches the budget.
#pragma once

#include "reads.h"   // SequenceReader, FastxRecord

#include <atomic>
#include <cstdint>
#include <optional>
#include <vector>

namespace fa { namespace cpu { namespace io {

class BatchReader {
public:
    // batch_bp: bp budget per batch (must be > 0).
    BatchReader(SequenceReader& reader, int64_t batch_bp)
        : reader_(reader), batch_bp_(batch_bp) {}

    // Next batch, or std::nullopt at EOF. `seen` counts the records returned.
    std::optional<std::vector<FastxRecord>> next(std::atomic<int64_t>& seen) {
        if (done_) return std::nullopt;
        std::vector<FastxRecord> batch;
        int64_t batch_bp_seen = 0;
        FastxRecord rec;
        for (;;) {
            if (!reader_.next(rec)) { done_ = true; break; }
            batch_bp_seen += static_cast<int64_t>(rec.seq.size());
            batch.push_back(std::move(rec));
            ++seen;
            if (batch_bp_seen >= batch_bp_) break;
        }
        if (batch.empty()) return std::nullopt;
        return batch;
    }

private:
    SequenceReader& reader_;
    int64_t batch_bp_;
    bool done_ = false;
};

}}}  // namespace fa::cpu::io
