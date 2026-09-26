// Umbrella header for the DP layer: the ksw2 adapter and the minimap2-style controller.
#pragma once

#include "../core/cigar.h"  // CIGAR text/stat helpers (output::*) used downstream

#include "result.h"
#include "params.h"
#include "ksw2_align.h"     // ksw2_simple_mat() and the extz2/extd2 call counters
#include "control.h"        // minimap2-style DP controller (dp_extend/dp_fill_gap)
