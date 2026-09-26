#pragma once

namespace fa { namespace cpu {

// DP options for the controller (control.h), as minimap2's: scoring (-A/-B/-O/-E/-O2/-E2/
// --score-N) and ksw2 control (-r bands, Z-drop, end bonus, max_sw_mat). Penalties are stored
// positive. -g max_gap is not here: as in minimap2 it sizes the DP window around the kernel.
// The defaults are minimap2's map-ont row; presets override them per call.
struct DpMapOpt {
    // scoring (minimap2 -A/-B/--score-N and the two gap systems -O/-E,-O2/-E2)
    int a = 2;            // -A  match            (>0)
    int b = 4;            // -B  mismatch penalty (stored positive)
    int sc_ambi = 1;      // --score-N            (stored positive)
    int q = 4,  e = 2;    // -O/-E   gap system 1
    int q2 = 24, e2 = 1;  // -O2/-E2 gap system 2

    // control
    int zdrop = 400;      // --z-drop : off-diagonal drop-off that stops/splits DP
    int zdrop_inv = 200;  // minimap2 opt->zdrop_inv: the local-inversion Z-drop.
                          // Gates mm_test_zdrop's reverse-complement probe and
                          // becomes the retry Z-drop on a code-2 verdict.
    int end_bonus = -1;   // flank bonus for reaching the query end (EXTZ_ONLY only)
    int junc_bonus = 9;   // --junc-bonus; inert while junc == nullptr
    int junc_pen   = 5;   // --junc-pen; inert while junc == nullptr
    int bw = 500;         // -r flank band      (narrow; pre ×1.5)
    int bw_long = 20000;  // -r internal band   (wide;   pre ×1.5)
    int64_t max_sw_mat = 100000000; // per-call tlen*qlen cap; over it -> skip (zdropped)

    // Effective bands, as mm_align1 (bw * 1.5 + 1); the internal band is never below the
    // flank band.
    int bw_eff()      const { return static_cast<int>(bw * 1.5 + 1.0); }
    int bw_long_eff() const {
        const int v = static_cast<int>(bw_long * 1.5 + 1.0);
        const int f = bw_eff();
        return v < f ? f : v;
    }
    // Kernel selection, mirroring mm_align_pair(): equal gap systems -> single
    // affine (ksw_extz2_sse); unequal -> dual affine (ksw_extd2_sse).
    bool single_affine() const { return q == q2 && e == e2; }
};

}}  // namespace fa::cpu
