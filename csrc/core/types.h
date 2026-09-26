// Nucleotide encoding and shared alignment types.
#pragma once

#include <flashalign/alignment.hpp>

#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <climits>

namespace fa { namespace cpu {

// A C G T/U -> 0 1 2 3, anything else -> 4.
inline int32_t nuc_encode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return 4;
    }
}

inline std::vector<uint8_t> encode_sequence_u8(const std::string& seq) {
    std::vector<uint8_t> enc(seq.size());
    for (size_t i = 0; i < seq.size(); i++) {
        enc[i] = static_cast<uint8_t>(nuc_encode(seq[i]));
    }
    return enc;
}

using AlignmentOrigin = ::flashalign::AlignmentOrigin;
using AlignResult = ::flashalign::Alignment;

// Supplementary segments to realize beside a primary. Query coordinates are on the
// winning strand, not the forward read.
struct SupplementaryCandidate {
    int diag = 0;    // r - q of the segment
    int q_lo = 0;    // query span
    int q_hi = 0;
};
struct SupplementaryPlan {
    int primary_q_lo = 0;   // primary chain query span
    int primary_q_hi = 0;
    std::vector<SupplementaryCandidate> cands;
};


}}  // namespace fa::cpu
