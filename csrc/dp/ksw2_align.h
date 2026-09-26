// ksw2 helpers for the DP controller (control.h): the score matrix and kernel call counters.
#pragma once

#include "result.h"   // GotohCigarResult

#include "ksw2.h"

#include <cstdint>
#include <cstdlib>    // free
#include <cmath>      // llround
#include <atomic>
#include <string>
#include <vector>

namespace fa { namespace cpu {

// minimap2's m=5 score matrix (A,C,G,T,N), as ksw_gen_simple_mat: +|match| on the diagonal,
// -|mismatch| off it, -|ambi| on the N row and column.
inline void ksw2_simple_mat(int8_t* mat, int match, int mismatch, int ambi) {
    const int8_t a = static_cast<int8_t>(match < 0 ? -match : match);
    const int8_t b = static_cast<int8_t>(mismatch > 0 ? -mismatch : mismatch);
    const int8_t sc_ambi = static_cast<int8_t>(ambi > 0 ? -ambi : ambi);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) mat[i * 5 + j] = (i == j) ? a : b;
        mat[i * 5 + 4] = sc_ambi;
    }
    for (int j = 0; j < 5; ++j) mat[4 * 5 + j] = sc_ambi;
}

inline std::atomic<uint64_t>& ksw2_extz2_call_counter() {
    static std::atomic<uint64_t> calls{0};
    return calls;
}

inline std::atomic<uint64_t>& ksw2_extd2_call_counter() {
    static std::atomic<uint64_t> calls{0};
    return calls;
}

}}  // namespace fa::cpu
