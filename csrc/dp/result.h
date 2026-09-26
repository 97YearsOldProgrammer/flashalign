#pragma once

#include <string>

namespace fa { namespace cpu {

struct GotohCigarResult {
    std::string cigar;
    int ref_offset;      // ref bases skipped at start (semi-global)
    int ref_consumed;    // ref bases consumed (M + D count)
    int matches;         // exact character matches
    int score = 0;       // terminal DP score for accepted tracebacks
    // On a Z-drop the traceback ends at the best-scoring cell and the controller splits
    // there; max_q/max_t are ksw2's zero-based coordinates of that cell.
    bool zdropped = false;
    int max_q = -1;      // raw zero-based query max-cell coordinate
    int max_t = -1;      // raw zero-based target max-cell coordinate
};

}}  // namespace fa::cpu
