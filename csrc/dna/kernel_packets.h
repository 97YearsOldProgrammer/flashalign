// The DNA realization geometry between a block's first and last anchors. No
// code here: the planner is plan_verified_geometry() (ordered_anchor_path.h),
// the executors are in realization_controller.cpp and the gap fill is
// dp_fill_gap().
//
// The anchors are the dense chain's path inside the block after
// mm_fix_bad_ends and filter_bad_seeds, in chain order. An interior anchor
// flagged ANCHOR_IGNORE or ANCHOR_TANDEM is skipped as in minimap2; the first
// and last anchors never are. Every corner is a k-mer center, as minimap2's
// mm_adjust_minier.
//
// A region is a maximal run of consecutive retained anchors on one diagonal;
// a long join starts a new one. The gap between two neighbouring anchors of a
// region is the bases between their k-mers (none when they overlap). A region
// of two or more anchors is verified when every gap passes the ungapped
// certificate
//
//   S0 >= (L - 1) * a - 2 * g(1),   g(l) = min(q + e * l, q2 + e2 * l)
//
// with L the gap length and S0 its ungapped score: a path that leaves the
// diagonal keeps at most L - 1 match columns and pays two one-base gaps. A
// verified region is emitted ungapped, with no DP, and cannot Z-drop. It is
// optimal among paths through every anchor center of the region, not among all
// paths between its ends.
//
// A seam joins two neighbouring verified regions with one gap fill from the
// left region's last center to the right region's first. It is minimap2's
// fill: pinned-global DP at bw_long_eff(), or max(q_span, r_span) across a
// long join, with the approximate pass, the Z-drop test and one exact retry.
//
// A stretch is everything else: a region with any failing gap, a lone anchor,
// the material between verified regions that are not neighbours. minimap2's
// loop cuts it into pieces, each ending at the stretch end, at a long join or
// at the first retained anchor where both spans since the previous corner reach
// min_ksw_len. Each piece is one fill, like a seam.
//
// ANCHOR_TANDEM marks an anchor whose key has another posting on the same
// contig within dna_tandem_window. In a tandem array a seed may sit on the
// wrong copy and a corner would pin the path there, so flagged interior anchors
// are skipped and one fill crosses the array. Unlike minimap2's query-side
// MM_SEED_TANDEM this flag is reference-side and dense: a read lying inside an
// array becomes one stretch and one fill.
//
// A Z-drop in a fill splits the block at the fill's maximum cell: the prefix
// is kept and the anchors past the cut form a new segment while at least
// cigar_dp_split_min_anchors remain. A fill over the matrix cap
// (DpMapOpt::max_sw_mat) is returned Z-dropped at its start, so the block
// splits there.
#pragma once
