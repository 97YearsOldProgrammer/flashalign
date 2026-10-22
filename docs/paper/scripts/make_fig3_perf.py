#!/usr/bin/env python3
"""Fig. 3 (v2): speed + peak memory as a 2 x 2 grid of clustered bars.

Grammar
  rows    = metric.  top = Speed (Gbp/hr over 32 threads), bottom = Peak memory
            (GB).  Both columns of a row share one y range, so a row reads
            straight across; only the left column carries the axis label, but
            every panel carries its own numerals (fig-perf's own convention).
  columns = experiment part, left -> right: DNA long reads (HiFi + ONT + ONT
            ultra-long, the third arm of the same panel) and Long-read RNA
            (Kinnex + ONT cDNA).  Column width is len(groups) group-pitches
            (the RNA column gets `pad` extra so its in-frame key and its
            rotated tool names have room), so the group pitch -- and with it
            the bar width -- is one constant across both columns.
  x       = tool, always flash, minibwa, minimap2, rammap, winnowmap2, with the
            names set diagonally (gnuplot's `rotate by -45`) under EVERY panel.
  colour  = readset (figstyle.ARM_FILL), named by a key inside the panel frame.

Nothing is drawn above the panels: no column titles and no panel letters
(the caption's "Left" / "Right" names the columns), no reference and no output-mode sub-line.  REF_LABEL /
MODE_LABEL survive as the fail-loud vocabulary -- a reference or a mode in
perf.tsv with no entry here is still a hard error -- so the caption's claims
("to GRCh38", "ultra-long to T2T-CHM13v2.0", "PAF with CIGAR") stay
contract-checked against the data even though nothing prints them.  The
reference check is per (part, arm): the rows of arm X must carry exactly
PARTS[...]["ref"][X], which is what backs the caption's ultra-long clause.

Colour key: one frameless key INSIDE every panel frame (v1's fig-perf key --
text before a wide sample).  Its corner is not hard-coded; after a first draw
each panel tries "upper right", "upper center", "center right", "upper left"
in turn and keeps the first whose legend box clears (by 2 pt) every bar, every
numeral and every pending slot of that panel.  If none is free the render dies
and asks for more headroom (TOP_FACTOR) rather than shipping a collision.

Placeholder rule: every part draws its full tool roster (PARTS[...]["tools"]).
A roster tool with no rows -- today flash in the DNA part -- is drawn as an
empty dashed "pending" slot with a grey tick label instead of bars, never
silently dropped.  It is derived from the data, so appending the rows to
perf.tsv makes the placeholder disappear on the next render; nothing in this
script has to change.  --no-flash drops MEASURED flash rows: their slot then
says "omitted", and the flag refuses the default --out so the manuscript
figure cannot be overwritten with it.

Values: every bar carries its measured value, rotated 90 deg above the cap,
to one decimal (f"{v:.1f}"; perf.tsv carries two).  The axis top is 1.14 x the
tick ceiling purely to make room for them; the tick range, the pending slot
and the out-of-range contract check all still run to the ceiling (160 / 50).
--no-values suppresses the numerals and changes nothing else, so the two
variants are geometrically identical.  The thread count in the speed label is
read off the rows (one value figure-wide, enforced).

To add a tool: append its rows to perf.tsv (every arm of the part is
required), add its name to TOOL_ORDER and to the roster of every part that
should show it -- an unlisted tool is a hard error by design; the column
grows by one pitch.  To add a part: one PARTS entry (key, panel, tools,
arms, per-arm reference claim) plus its REF_LABEL / MODE_LABEL entries; the
widths and the shared y ranges follow automatically.

Input: data/perf.tsv only.  Single measured pass, 32 threads, no
repetitions -> no error bars are drawn and none exist.
Usage:
  make_fig3_perf.py [--perf TSV] [--out PDF] [--preview PNG]
                       [--no-flash] [--no-values]
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, NullLocator

import figstyle

root = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------- data contract
# One column per part.  `panel`/`arms` select the rows; `ref` is the reference
# set claimed PER ARM and MUST equal what that arm's rows carry (this is what
# backs the caption's "ultra-long to T2T-CHM13v2.0"); `tools` is the roster the
# column promises, in TOOL_ORDER: a roster tool with no rows is drawn as a
# dashed "pending" slot (never silently dropped) and a tool with rows that is
# not on the roster is a hard error.  `pad` is extra column width in group
# pitches (the RNA column needs it for its in-frame key and its tool names).
DNA_TOOLS = ("flash", "minibwa", "minimap2", "rammap", "winnowmap2")
PARTS = [
    dict(key="dna", panel="dna_long", tools=DNA_TOOLS,
         arms=("hifi", "ont", "ul"), pad=0.0,
         ref={"hifi": "grch38_no_alt", "ont": "grch38_no_alt", "ul": "chm13v2"}),
    # the two tools with spliced-alignment lanes; the caption says so
    dict(key="rna", panel="rna", tools=("flash", "minimap2"),
         arms=("pb_kinnex", "ont_cdna"), pad=0.6,
         ref={"pb_kinnex": "grch38_gencode44", "ont_cdna": "grch38_gencode44"}),
]
# The reference / mode vocabulary.  NOTHING is drawn from it any more: it is
# kept as the fail-loud contract behind the caption's prose, so a reference or
# a mode in perf.tsv that this figure has never been told about is a hard
# error instead of a silent caption lie.
REF_LABEL = {"grch38_no_alt": "GRCh38 (no-ALT)",
             "grch38_gencode44": "GRCh38 + GENCODE 44",
             "chm13v2": "T2T-CHM13v2.0"}
MODE_LABEL = {"sam": "SAM output", "cigar": "PAF + CIGAR output"}

TOOL_ORDER = ["flash", "minibwa", "minimap2", "rammap", "winnowmap2"]
HERO = "flash"
DISPLAY = {"flash": "flashalign"}   # tick-label spelling; data key stays "flash"
# Winnowmap2 takes every hardware thread whatever -t says (64 on the
# benchmark host); its rows carry threads=64 as measured and the caption
# says so. Every other row must agree on the one count the speed label claims.
THREAD_EXEMPT = {"winnowmap2"}
ROWS = [  # metric column, y label, (lo, tick ceiling, step) -- same everywhere
    ("gbp_per_hr", "Speed (Gbp/hr over {threads} threads)", (0, 800, 100)),  # threads: from the rows
    ("peak_gb", "Peak memory (GiB)", (0, 70, 10)),  # GiB: the harness reports ru_maxrss / 2^30
]
NUMERIC = ("gbp", "gbp_per_hr", "peak_gb")
MAX_ARMS = max(len(p["arms"]) for p in PARTS)   # sets the one constant bar width

# ------------------------------------------------------------ absolute-inch box
W = 6.9291                      # the paper's \\textwidth, exactly
PH = 1.63                       # panel height, both rows
# top is bare clearance (no column titles, no panel letters); bottom and hspace hold the -45 deg tool names (the
# descender of "winnowmap2" is the binding constraint, checked after the draw).
MARGIN = dict(left=0.42, right=0.10, top=0.10, bottom=0.58,
              wspace=0.26, hspace=0.54)
H = MARGIN["bottom"] + 2 * PH + MARGIN["hspace"] + MARGIN["top"]

FONT = figstyle.BAR_FONT        # 7 pt: y numerals, y labels, key,
                                #       and the diagonal tool names
SMALL = 6.0                     # 6 pt: bar numerals, the pending word
GREY = "#767676"                # pending cue: readable, clearly not a datum
SUBGREY = "#444444"             # bar numerals
TOP_FACTOR = 1.14               # axis top / tick ceiling: headroom for numerals
VAL_GAP = 1.5                   # pt between a bar cap and its numeral

# In-panel key: the corner is searched, not chosen (see place_key).
KEY_LOCS = ["upper right", "upper center", "center right", "upper left"]
KEY_PAD_PT = 2.0                # clearance the key must keep from every datum

ap = argparse.ArgumentParser()
ap.add_argument("--perf", type=Path, default=root / "data" / "perf.tsv")
DEFAULT_OUT = root / "figures/fig3_perf.pdf"
ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
ap.add_argument("--preview", type=Path, default=None)
# parity with make_fig2_mapq.py: --no-flash drops the hero lane. A part whose
# MEASURED flash rows were dropped this way gets an "omitted" slot (not
# "pending"), and the manuscript figure can never be overwritten with it.
ap.add_argument("--no-flash", action="store_true")
# --no-values is a proofing switch: same axes, same geometry, no numerals.
ap.add_argument("--no-values", action="store_true")
a = ap.parse_args()
if a.no_flash and a.out.resolve() == DEFAULT_OUT.resolve():
    raise SystemExit("--no-flash suppresses measured flash lanes: give an explicit "
                     "--out so the manuscript figure is never overwritten")


def die(msg):
    raise SystemExit(f"{a.perf}: {msg}")


# ---------------------------------------------------------------- read + verify
claim = {}
for p in PARTS:
    for arm in p["arms"]:
        if arm not in figstyle.ARM_FILL or arm not in figstyle.ARM_LABEL:
            die(f"arm {arm!r} of part {p['key']!r} has no figstyle colour/label")
        if arm not in p["ref"]:
            die(f"arm {arm!r} of part {p['key']!r} claims no reference")
        if p["ref"][arm] not in REF_LABEL:
            die(f"part {p['key']!r} arm {arm!r} claims reference "
                f"{p['ref'][arm]!r}, which has no REF_LABEL entry")
        if (p["panel"], arm) in claim:
            die(f"arm {arm!r} of panel {p['panel']!r} claimed by two parts")
        claim[(p["panel"], arm)] = p["key"]

for p in PARTS:
    if HERO not in p["tools"] or set(p["tools"]) - set(TOOL_ORDER):
        die(f"part {p['key']!r} roster {list(p['tools'])} must be drawn from "
            f"TOOL_ORDER {TOOL_ORDER} and include {HERO!r}")
    if [t for t in TOOL_ORDER if t in p["tools"]] != list(p["tools"]):
        die(f"part {p['key']!r} roster {list(p['tools'])} is not in TOOL_ORDER order")

lanes = {p["key"]: {} for p in PARTS}
had_flash = {p["key"]: False for p in PARTS}   # flash rows dropped by --no-flash
with open(a.perf, newline="") as fh:
    body = [ln for ln in fh if not ln.startswith("#")]
for r in csv.DictReader(body, delimiter="\t"):
    if r["tool"] not in TOOL_ORDER:
        die(f"unknown tool {r['tool']!r}; known: {TOOL_ORDER}")
    key = claim.get((r["panel"], r["readset"]))
    if key is None:
        die(f"row panel={r['panel']!r} readset={r['readset']!r} belongs to no "
            f"part of this figure; known: {sorted(claim)}")
    if r["tool"] == HERO and a.no_flash:
        had_flash[key] = True
        continue
    if (r["readset"], r["tool"]) in lanes[key]:
        die(f"duplicate row for part {key!r} {r['readset']}/{r['tool']}")
    if r["reference"] not in REF_LABEL:
        die(f"reference {r['reference']!r} has no REF_LABEL entry; known: "
            f"{sorted(REF_LABEL)}")
    if r["mode"] not in MODE_LABEL:
        die(f"mode {r['mode']!r} has no MODE_LABEL entry; known: "
            f"{sorted(MODE_LABEL)}")
    for col in NUMERIC:
        try:
            v = float(r[col])
        except (TypeError, ValueError):
            die(f"{key}/{r['readset']}/{r['tool']}: {col}={r[col]!r} is not a number")
        if v <= 0:
            die(f"{key}/{r['readset']}/{r['tool']}: {col}={r[col]} is not positive")
    lanes[key][(r["readset"], r["tool"])] = r

groups_of, pending_of, omitted_of = {}, {}, {}
for p in PARTS:
    lane = lanes[p["key"]]
    if not lane:
        die(f"no rows for part {p['key']!r} (panel={p['panel']!r}, arms={p['arms']})")
    for col in ("threads", "mode"):
        vals = sorted({r[col] for r in lane.values()
                       if not (col == "threads" and r["tool"] in THREAD_EXEMPT)})
        if len(vals) > 1:
            die(f"part {p['key']!r} mixes {col} values: {vals}")
    # The reference claim is per ARM: hifi/ont go to GRCh38 while the ul arm of
    # the same panel goes to T2T-CHM13v2.0, and the caption says exactly that.
    for arm in p["arms"]:
        refs = sorted({r["reference"] for (rs, _t), r in lane.items() if rs == arm})
        if refs and refs != [p["ref"][arm]]:
            die(f"part {p['key']!r} arm {arm!r}: the caption claims reference "
                f"{p['ref'][arm]!r} but the rows carry {refs}")
    tools = {t for _, t in lane}
    extra = tools - set(p["tools"])
    if extra:
        die(f"part {p['key']!r}: {sorted(extra)} have rows but are not on the "
            f"part's roster {list(p['tools'])}")
    for tool in tools:                       # every tool must cover every arm
        have = {arm for arm, t in lane if t == tool}
        if have != set(p["arms"]):
            die(f"part {p['key']!r}: tool {tool!r} has rows for {sorted(have)} "
                f"but the part's arms are {list(p['arms'])}")
    # The roster is drawn in full: a roster tool with no rows becomes a dashed
    # slot -- "pending" (never measured) or "omitted" (dropped by --no-flash).
    groups_of[p["key"]] = list(p["tools"])
    pending_of[p["key"]] = [t for t in p["tools"] if t not in tools]
    omitted_of[p["key"]] = {t for t in pending_of[p["key"]]
                            if t == HERO and had_flash[p["key"]]}
    for t in pending_of[p["key"]]:
        print(f"    NOTE part {p['key']}: no rows for {t} -> "
              f"{'omitted' if t in omitted_of[p['key']] else 'pending'} slot")
    for metric, _, (lo, hi, _s) in ROWS:     # contract runs to the TICK ceiling
        for (arm, tool), r in sorted(lane.items()):
            v = float(r[metric])
            if not lo <= v <= hi:
                die(f"part {p['key']!r} {tool}/{arm}: {metric}={v} is outside the "
                    f"shared axis [{lo}, {hi}] -- widen the range for ALL parts")

# The thread count is the one claim in the axis label, so it is read off the
# rows too: one value figure-wide, or the shared speed row is meaningless.
threads = sorted({r["threads"] for lane in lanes.values() for r in lane.values()
                  if r["tool"] not in THREAD_EXEMPT})
if len(threads) != 1:
    die(f"rows mix thread counts {threads}; the shared speed row and the caption "
        f"claim one value")
THREADS = threads[0]
for lane in lanes.values():
    for r in lane.values():
        if r["tool"] in THREAD_EXEMPT and r["threads"] != THREADS:
            print(f"    NOTE {r['tool']}/{r['readset']}: threads={r['threads']}, not the "
                  f"label's {THREADS}; the caption carries the caveat")

# ------------------------------------------------------------------- the figure
figstyle.apply_house_style()
fig = plt.figure(figsize=(W, H))
rend = fig.canvas.get_renderer()

avail = W - MARGIN["left"] - MARGIN["right"] - (len(PARTS) - 1) * MARGIN["wspace"]
pitch_total = sum(len(groups_of[p["key"]]) + p["pad"] for p in PARTS)
pitch_in = avail / pitch_total
# One constant slot for the whole figure (gnuplot histogram rule at MAX_ARMS):
# bar width AND in-cluster spacing are the same in every column, so a bar in
# the 2-arm RNA column is exactly as wide as one in the 3-arm DNA column.
SLOT = 1.0 / (MAX_ARMS + 2)
bw = SLOT * figstyle.BAR_WIDTH_FRAC


def cluster_offsets(n_arms):
    """Bar-centre offsets of an n-arm cluster on the figure-wide slot."""
    return [SLOT * (i - (n_arms - 1) / 2.0) for i in range(n_arms)]


def bbox_pad(b, pt):
    """Display-space bbox grown by `pt` printed points on every side."""
    d = pt * fig.dpi / 72.0
    return (b.x0 - d, b.y0 - d, b.x1 + d, b.y1 + d)


def hits(box, b):
    return not (box[2] <= b.x0 or box[0] >= b.x1 or box[3] <= b.y0 or box[1] >= b.y1)


def place_key(ax, arms, obstacles, tag):
    """v1's fig-perf key, inside the frame, at the first corner that is free.

    Placement is measured, not assumed: for each candidate loc the legend is
    realised, drawn, and its window extent tested (2 pt of padding) against
    every bar, every numeral and every pending slot of THIS panel."""
    handles = [Patch(facecolor=figstyle.ARM_FILL[arm], edgecolor="black",
                     linewidth=figstyle.BAR_EDGE_LW, label=figstyle.ARM_LABEL[arm])
               for arm in arms]
    for loc in KEY_LOCS:
        leg = ax.legend(handles=handles, loc=loc, frameon=False,
                        fontsize=figstyle.BAR_FONT, markerfirst=False,
                        handlelength=2.76, handleheight=0.51, handletextpad=0.61,
                        labelspacing=0.05, borderaxespad=0.5)
        fig.canvas.draw()
        box = bbox_pad(leg.get_window_extent(renderer=rend), KEY_PAD_PT)
        if not any(hits(box, o.get_window_extent(renderer=rend)) for o in obstacles):
            print(f"    key {tag}: loc={loc!r}")
            return loc
        leg.remove()
    die(f"panel {tag}: none of {KEY_LOCS} leaves the colour key clear of the "
        f"bars/numerals/pending slots -- raise TOP_FACTOR (now {TOP_FACTOR}) "
        f"to buy headroom, or widen the column's pad")



x0 = MARGIN["left"]
row_axes = [[] for _ in ROWS]     # (part key, ax) per metric row
panels = []                       # (tag, ax, arms, obstacles)
numerals = []                     # every bar's value label, all panels
for ci, p in enumerate(PARTS):
    lane, groups, arms = lanes[p["key"]], groups_of[p["key"]], p["arms"]
    n = len(groups)
    pw = pitch_in * (n + p["pad"])
    offsets = cluster_offsets(len(arms))
    any_row = next(iter(lane.values()))
    print(f"[{p['key']}] groups={groups} arms={list(arms)} "
          f"threads={any_row['threads']} mode={any_row['mode']} "
          f"references={ {arm: p['ref'][arm] for arm in arms} } "
          f"pending={pending_of[p['key']]}")
    for ri, (metric, ylabel, (lo, hi, step)) in enumerate(ROWS):
        y0 = MARGIN["bottom"] + (1 - ri) * (PH + MARGIN["hspace"])
        ax = fig.add_axes([x0 / W, y0 / H, pw / W, PH / H])
        obstacles = []
        for gi, tool in enumerate(groups, start=1):
            if tool in pending_of[p["key"]]:
                # Empty dashed slot spanning the TICK range (not the numeral
                # headroom), as wide as a full cluster. Produced only because
                # the part has no rows for this tool: append them and it is
                # gone.
                half = abs(offsets[0]) + bw / 2
                word = "omitted" if tool in omitted_of[p["key"]] else "pending"
                obstacles.append(ax.add_patch(Rectangle(
                    (gi - half, lo), 2 * half, hi - lo, facecolor="none",
                    edgecolor=GREY, linewidth=0.35, linestyle=(0, (2.2, 1.8)),
                    zorder=1)))
                ax.text(gi, lo + (hi - lo) * 0.5, word, rotation=90,
                        ha="center", va="center", fontsize=SMALL, color=GREY)
                continue
            for arm, off in zip(arms, offsets):
                r = lane[(arm, tool)]
                v = float(r[metric])
                obstacles.append(ax.bar(
                    gi + off, v, width=bw, color=figstyle.ARM_FILL[arm],
                    edgecolor="black", linewidth=figstyle.BAR_EDGE_LW,
                    zorder=3)[0])
                # Rotated, centred on the bar: at a 0.2-pitch cluster spacing
                # the 6 pt glyph column is far narrower than the gap.  One
                # decimal everywhere (perf.tsv carries two): a uniform, stated
                # rounding, not mixed significant figures.  Always created --
                # --no-values HIDES them after the keys are placed, so both
                # variants share one geometry and one key placement.
                numerals.append(ax.annotate(
                    f"{v:.1f}", xy=(gi + off, v), xytext=(0, VAL_GAP),
                    textcoords="offset points", rotation=90, ha="center",
                    va="bottom", fontsize=SMALL, color=SUBGREY, zorder=4,
                    annotation_clip=False))
                obstacles.append(numerals[-1])
                if ri == 0:
                    print(f"    {tool:<11s} {arm:<9s} "
                          f"gbp_per_hr={float(r['gbp_per_hr']):7.2f}  "
                          f"peak_gb={float(r['peak_gb']):6.2f}  {r['source']}")
        # The RNA column's pad keeps its groups centred rather than left-packed.
        ax.set_xlim(0.5 - p["pad"] / 2, n + 0.5 + p["pad"] / 2)
        ax.set_ylim(lo, hi * TOP_FACTOR)
        # Ticks stop at the ceiling: the headroom above it is for the numerals,
        # not a data range, so it must not grow the ladder (a locator would).
        ax.yaxis.set_major_locator(FixedLocator(
            [lo + i * step for i in range(int(round((hi - lo) / step)) + 1)]))
        ax.yaxis.set_minor_locator(NullLocator())
        figstyle.style_bar_axes(ax)
        # fig-perf mirrors the y majors as bare stubs on the spineless right
        # edge; in a grid those stubs would sit beside the next column's
        # numerals and read as stray dashes, so only the last column keeps them.
        ax.tick_params(axis="y", right=(ci == len(PARTS) - 1))
        if ci == 0:                   # axis label left column only, numerals on all
            ax.set_ylabel(ylabel.format(threads=THREADS), fontsize=FONT)
        # gnuplot 'rotate by -45': "\" slant, start-justified at the group
        # centre -- under EVERY panel, both rows (v1's grammar).
        ax.set_xticks(range(1, n + 1))
        ax.set_xticklabels([DISPLAY.get(t, t) for t in groups], rotation=-45, ha="left",
                           rotation_mode="anchor", fontsize=FONT)
        ax.tick_params(axis="x", pad=2.5)
        for lab, tool in zip(ax.get_xticklabels(), groups):
            if tool in pending_of[p["key"]]:
                lab.set_color(GREY)
        row_axes[ri].append((p["key"], ax))
        panels.append((f"{p['key']}/{metric}", ax, arms, obstacles))
    x0 += pw + MARGIN["wspace"]

# Numerals must stay inside their frame: the range contract bounds the BAR at
# the tick ceiling, so the headroom above it is measured after the draw.
fig.canvas.draw()
for tag, ax, _arms, obst in panels:
    top = ax.get_window_extent(renderer=rend).y1
    for o in obst:
        if isinstance(o, Text):
            over = (o.get_window_extent(renderer=rend).y1 - top) / fig.dpi * 72.0
            if over > 0:
                die(f"panel {tag}: numeral {o.get_text()!r} overruns the frame top "
                    f"by {over:.2f} pt -- raise TOP_FACTOR (now {TOP_FACTOR})")

# The key goes in only once every bar and numeral of its panel exists.
locs = [(tag, place_key(ax, arms, obst, tag)) for tag, ax, arms, obst in panels]
if a.no_values:                   # same geometry, same key corners, no numerals
    for t in numerals:
        t.set_visible(False)

# The diagonal names are the one thing that can collide as columns narrow, so
# measure them rather than trust the arithmetic: no name may touch another, run
# off the bottom of the page, or reach into the panel below.
fig.canvas.draw()
lab_boxes = []
for ri, row in enumerate(row_axes):
    boxes = []
    for key, ax in row:
        for lab, tool in zip(ax.get_xticklabels(), groups_of[key]):
            boxes.append((key, tool, lab.get_window_extent(renderer=rend)))
    for (k1, t1, b1), (k2, t2, b2) in [(x, y) for i, x in enumerate(boxes)
                                       for y in boxes[i + 1:]]:
        if not (b1.x1 <= b2.x0 or b2.x1 <= b1.x0 or b1.y1 <= b2.y0 or b2.y1 <= b1.y0):
            die(f"tool names {k1}/{t1} and {k2}/{t2} overlap -- widen "
                f"MARGIN['wspace'] or the column pad")
    for key, tool, b in boxes:
        if b.y0 < 0:
            die(f"tool name {key}/{tool} runs {(-b.y0) / fig.dpi * 72.0:.2f} pt "
                f"off the bottom of the page -- raise MARGIN['bottom'] "
                f"(now {MARGIN['bottom']}) and lower PH")
    if ri + 1 < len(row_axes):
        below = [ax.get_window_extent(renderer=rend) for _k, ax in row_axes[ri + 1]]
        for key, tool, b in boxes:
            for nb in below:
                if not (b.x1 <= nb.x0 or nb.x1 <= b.x0 or b.y1 <= nb.y0 or nb.y1 <= b.y0):
                    die(f"tool name {key}/{tool} of row {ri} reaches into the "
                        f"panel below -- raise MARGIN['hspace'] "
                        f"(now {MARGIN['hspace']}) and lower PH")
    lab_boxes += boxes

drop = min((b.y0 for _k, _t, b in lab_boxes)) / fig.dpi * 72.0
print(f"figure {W:.4f} x {H:.4f} in; group pitch {pitch_in:.4f} in; "
      f"bar {bw * pitch_in:.4f} in; lowest tool name {drop:.2f} pt above the "
      f"page foot; keys " + ", ".join(f"{t.split()[0]}={l}" for t, l in locs))
figstyle.save_figure(fig, a.out, a.preview)
