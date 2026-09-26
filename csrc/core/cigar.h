// CIGAR building, parsing, clipping and replay, shared by the DP and output code.
#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace output {

// Appends `count` of `op`, merged into the last run when the operator matches; a count <= 0
// writes nothing.
void append_cigar_run(std::string& cigar, int count, char op);

// Appends every run of `suffix` through append_cigar_run.
void append_cigar_string(std::string& cigar, const std::string& suffix);

// Puts `prefix` in front of `cigar`, as `append_cigar_string(prefix, cigar)` would, but in
// time linear in `prefix`. Assumes `cigar` has no adjacent same-op runs, so only its leading
// run can merge; the rest is copied as is.
void prepend_cigar_string(std::string& cigar, std::string prefix);

// Packed form: one uint32 per run, `count << 4 | op` with op indexing "MIDNSHP=XB", as in BAM
// and the ksw2 traceback. The packed helpers below produce the same runs as their text
// counterparts, so the DNA realizer can build CIGARs as runs and render text once.
using PackedCigar = std::vector<std::uint32_t>;

constexpr char kPackedCigarOps[] = "MIDNSHP=XB";

constexpr int packed_run_count(std::uint32_t run) {
    return static_cast<int>(run >> 4);
}
constexpr int packed_run_op(std::uint32_t run) {
    return static_cast<int>(run & 0xfu);
}

// Packed append_cigar_run.
void append_packed_run(PackedCigar& cigar, int count, int op);

// Packed append_cigar_string.
void append_packed_cigar(PackedCigar& cigar, const PackedCigar& suffix);

// Packed prepend_cigar_string.
void prepend_packed_cigar(PackedCigar& cigar, PackedCigar prefix);

// Query bases consumed (M/I/S/=/X), or -1 when a run has a zero count.
int packed_cigar_query_consumed(const PackedCigar& cigar);

// Text rendering of a packed CIGAR.
void append_packed_cigar_text(std::string& out, const PackedCigar& cigar);
std::string packed_cigar_text(const PackedCigar& cigar);

struct CigarStats {
    std::string cigar;
    int ref_start_delta = 0;
    int ref_consumed = 0;
    int matches = 0;
    int left_soft = 0;
    int right_soft = 0;
};

// (count, op) pairs; empty when the string is malformed or has a zero count.
std::vector<std::pair<int, char>> parse_cigar_ops(const std::string& cigar);

struct TerminalHardClip {
    std::string cigar;  // a copy with leading/trailing terminal 'S' rewritten to 'H'
    int lead = 0;       // length of the leading terminal clip op (0 if none)
    int trail = 0;      // length of the trailing terminal clip op (0 if none)
};

// Rewrites the terminal soft clips of a supplementary record's CIGAR to hard clips, as
// minimap2 does, and reports their lengths so the caller can trim SEQ/QUAL to match.
TerminalHardClip hard_clip_terminal(const std::string& cigar);

// q_lo/q_hi are on the forward read, while a CIGAR is in reference orientation: on the reverse
// strand the leading clip covers the read's 3' end. The helpers below do this conversion.
struct ForwardQueryInterval {
    int q_lo = 0;
    int q_hi = 0;
    bool valid = false;
};

struct ReferenceOrientedTerminalClips {
    int lead = 0;
    int trail = 0;
};

constexpr int clamp_query_coord(int x, int read_len) {
    return x < 0 ? 0 : (x > read_len ? read_len : x);
}

constexpr ReferenceOrientedTerminalClips
reference_oriented_terminal_clips(int read_len, int forward_q_lo,
                                  int forward_q_hi, bool is_reverse) {
    const int rl = read_len < 0 ? 0 : read_len;
    const int q0 = clamp_query_coord(forward_q_lo, rl);
    const int q1raw = clamp_query_coord(forward_q_hi, rl);
    const int q1 = q1raw < q0 ? q0 : q1raw;
    return is_reverse
        ? ReferenceOrientedTerminalClips{rl - q1, q0}
        : ReferenceOrientedTerminalClips{q0, rl - q1};
}

// Forward-read interval covered by a reference-oriented CIGAR; S and H clips both count.
ForwardQueryInterval forward_query_interval_from_cigar(
    const std::string& cigar, int read_len, bool is_reverse);

// Replaces the terminal clips of a reference-oriented CIGAR with the soft clips for the
// forward interval [forward_q_lo, forward_q_hi). Returns an empty string when the body does
// not consume exactly that interval.
std::string clipped_reference_oriented_cigar(
    const std::string& reference_cigar, int read_len, int forward_q_lo,
    int forward_q_hi, bool is_reverse);

static_assert(reference_oriented_terminal_clips(1000, 100, 600, false).lead == 100,
              "forward segment leading clip");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, false).trail == 400,
              "forward segment trailing clip");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, true).lead == 400,
              "reverse segment leading clip is mirrored");
static_assert(reference_oriented_terminal_clips(1000, 100, 600, true).trail == 100,
              "reverse segment trailing clip is mirrored");

int cigar_query_consumed(const std::string& cigar);

CigarStats cigar_basic_stats(const std::string& cigar);

// Result of replaying a CIGAR against the query (reference orientation) and one reference
// sequence. `valid` is set only when every operation stays in bounds, the clips and body
// cover the whole query, and there is no interior clip.
//
// The counts follow minimap2's mm_update_extra: ambiguous (N) M/I/D residues are left out of
// block_len and matches and counted once in ambiguities, and
// edit_distance = block_len - matches + ambiguities. Consumed lengths include them.
struct CigarReplay {
    bool valid = false;
    int query_consumed = 0;         // M/I/S/H/=/X in full-query coordinates
    int aligned_query_consumed = 0; // M/I/=/X only
    int target_consumed = 0;        // M/D/N/=/X
    int target_start = -1;
    int target_end = -1;
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    int ambiguities = 0;
    int edit_distance = 0;
    int block_len = 0;
    int lead_clip = 0;
    int trail_clip = 0;
    // cs and MD strings without their tag prefixes; empty unless requested.
    std::string cs;
    std::string md;
    // The CIGAR with M/=/X split into =/X runs (--eqx); empty unless requested.
    std::string eqx_cigar;
};

// Optional outputs of replay_cigar; the default asks for none. They follow minimap2's
// write_cs_ds_core, write_MD_core and mm_update_cigar_eqx and never change the counts.
//
//   cs, over the aligned core (clips excluded): ":<len>" (Short) or "=<BASES>" (Long) for an
//   identical run, "*<ref><qry>" for a mismatch, "+<bases>" for I, "-<bases>" for D and
//   "~<t0><t1><len><t[len-2]><t[len-1]>" for N, on the forward reference. An identical run
//   ends at every M/=/X op boundary, so "5M5M" gives ":5:5".
//   MD: running match count; a mismatch writes "<count><REF>", a D writes "<count>^<REFS>",
//   and I and N neither write nor reset the count. Unlike minimap2, the final count is always
//   written, even "0", as the SAM MD grammar requires.
//   =/X: each M/=/X op becomes alternating = and X runs; other ops are copied.
//
// cs, MD and =/X compare base codes, so N against N is a match there although the counts
// treat it as an ambiguity.
struct CigarReplayRequest {
  enum class Cs { None, Short, Long };
  Cs cs = Cs::None;
  bool md = false;
  bool eqx = false;

  bool any() const { return cs != Cs::None || md || eqx; }
};

CigarReplay replay_cigar(const std::string& cigar, const uint8_t* qseq,
                         int qlen, const uint8_t* tseq, int tlen, int ref_start,
                         const CigarReplayRequest& request);

// Counts only.
CigarReplay replay_cigar(const std::string& cigar,
                                const uint8_t* qseq, int qlen,
                                const uint8_t* tseq, int tlen,
                                int ref_start);

}}} // namespace fa::cpu::output
