#!/usr/bin/env python3
"""Fig. 2: MAPQ placement-accuracy curves (minibwa Fig. 2 rendering).

Input (data/points.tsv, from scripts/fig2_points.py): one row per
emitted MAPQ threshold, cumulative from MAPQ 60 downwards.
Usage:
  make_fig2_mapq.py [--points TSV] [--out PDF] [--preview PNG] [--no-flash]
Panels (2 x 3): (a) CHM13 HiFi, (b) CHM13 ONT, (c) HG002 diploid HiFi,
(d) HG002 diploid ONT, (e) RNA human, (f) RNA mouse. The ultra-long cohort has
no panel of its own: it rides in the ONT panel as a second arm, so (b) and (d) draw six curves each.
Marker shape = tool everywhere. Colour: on the four DNA panels colour encodes
the TOOL and is identical across (a)-(d) -- flashalign black star, minimap2
green circle, winnowmap2 sky triangle, minibwa violet square -- and the
ultra-long arm is told apart by a DASHED line, not by colour. The RNA panels
keep minibwa Fig. 2e/2f's colour = tool x read arm, because their two arms are
different platforms rather than two lengths of the same one. One marker per
emitted MAPQ level. The per-arm tool roster is enforced upstream, in
fig2_points.py: minibwa and minimap2 map-ont are HiFi-arm / absent by rule, so
the ONT and ultra-long panels carry flash, minimap2 lr:hq and winnowmap2 only.
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator, NullFormatter, NullLocator

import figstyle

root = Path(__file__).resolve().parent.parent

XLABEL = "Accumulative error rate (#wrong / #aligned)"
YLABEL = "% aligned"
# title, letter, ylim, y major step -- y ranges trimmed per row, minibwa-style.
# The y window and tick step are panel-owned; the x window is NOT hand-set any
# more, it is derived from the points actually drawn in that panel: the decade
# floor of the smallest x and the decade ceiling of the largest (10**floor /
# 10**ceil of log10, i.e. always a 1e-n bound that covers every vertex), so a
# re-measured lane can move the axis but can never fall off it.
# arms: which read-arm suffixes of that points.tsv panel this panel draws. The
# ONT panels take ("(ONT)", "(UL)") -- the ultra-long cohort has no panel of its
# own. None = draw every arm (RNA).
# colour: "tool" resolves every curve through figstyle's bare-label table, so a
# tool has one colour across all four DNA panels and the ultra-long arm is told
# apart by a dashed line instead; "tool_arm" keeps the minibwa Fig. 2e/2f rule
# (colour = tool x read arm) for the RNA panels.
# xlim: hand-set window, or None to derive it (decade floor/ceiling of the drawn
# points -- see the axis-window rule below). The four DNA panels share ONE fixed
# window, DNA_XLIM, so (a)-(d) read straight across; the RNA panels keep the
# derived rule. A fixed window is checked fail-loud: a point outside it aborts
# the render instead of vanishing off the frame.
DNA_XLIM = (1e-8, 1e-1)
# legend: per-panel kwargs merged over LEGEND_BASE. It is a PANELS field, not a
# special case, because where the curves leave a hole is a property of the data
# in that panel. After the draw every legend box is measured against every drawn
# curve of its own panel and the render dies if one touches (LEGEND_PAD, in
# printed points) -- so a re-measured lane cannot quietly slide under the key.
# The six-entry ONT panels (b, d) sit upper-left, in the empty low-error corner;
# (b) tightens its border pad and (d) also shortens its sample line: at the house
# 2.1 em the six-entry box reaches past the flashalign (UL) curve and NONE of the
# nine anchors clears in (d) -- the panel has no empty 0.53 x 0.39 rectangle.
LEGEND_BASE = dict(frameon=False, fontsize=6, handlelength=2.1, handletextpad=0.45,
                   labelspacing=0.05, borderaxespad=0.8, markerfirst=False, numpoints=1)
LEGEND_PAD = 1.5
PANELS = [
    # (title, letter, points.tsv panel, arms or None, xlim or None, ylim, y step, colour, legend)
    ("Simulated HiFi reads to CHM13", "(a)", "Simulated long reads to CHM13", ("(HiFi)",), DNA_XLIM, (90, 100), 2, "tool", {"loc": "lower right"}),
    ("Simulated ONT reads to CHM13", "(b)", "Simulated long reads to CHM13", ("(ONT)", "(UL)"), DNA_XLIM, (88, 100), 2, "tool", {"loc": "upper left", "borderpad": 0.15, "borderaxespad": 0.3}),
    ("Simulated HiFi reads to HG002", "(c)", "Simulated long reads to HG002", ("(HiFi)",), DNA_XLIM, (0, 100), 10, "tool", {"loc": "lower right"}),
    ("Simulated ONT reads to HG002", "(d)", "Simulated long reads to HG002", ("(ONT)", "(UL)"), DNA_XLIM, (0, 100), 10, "tool", {"loc": "upper left", "handlelength": 1.4, "borderpad": 0.15, "borderaxespad": 0.3}),
    ("Simulated long RNA reads to human", "(e)", "Simulated long RNA reads to human", None, None, (75, 100), 5, "tool_arm", {"loc": "lower right"}),
    ("Simulated long RNA reads to mouse", "(f)", "Simulated long RNA reads to mouse", None, None, (70, 100), 5, "tool_arm", {"loc": "lower right"}),
]
COLOUR_KEYS = {"tool", "tool_arm"}
# The tool token in points.tsv is the figstyle registry key; this is the name
# the reader sees. Legends carry NO preset names --
# presets belong in the caption and Methods.
TOOL_LABEL = {"flash": "flashalign"}
# Legend / draw order: tool-major in TOOL_ORDER (hero first), then read arm.
ARM_RANK = {"(HiFi)": 0, "(ONT)": 1, "(UL)": 2}
UL_LINESTYLE = (0, (2.2, 1.3))   # printed points; the ultra-long arm only
# pct_aligned carries eval2roc's 0.5 pseudo-count on the cumulative aligned
# count (bench.mapeval_collect.PSEUDO), so a MAPQ-0 row on a fully placed cohort
# reads 100.00005 %. That is the estimator, not a placement, and it is plotted
# at 100.
Y_CAP = 100.0
TOOL_ORDER = ["flash", "minibwa", "minimap2", "winnowmap2"]

# Absolute-inch layout: 6.9291 in is the paper's true \textwidth, so the figure
# is included at width=\textwidth and prints at exactly these font sizes.
W = 6.9291
# minibwa's own grid: 2 rows x 3 columns, panels ~1.8 in wide at aspect ~1.1.
NROWS, NCOLS = 2, 3
MARGIN = dict(left=0.44, right=0.10, wspace=0.50, bottom=0.40, top=0.26, hspace=0.52)
PW = (W - MARGIN["left"] - MARGIN["right"] - (NCOLS - 1) * MARGIN["wspace"]) / NCOLS
PH = 1.63
H = MARGIN["bottom"] + NROWS * PH + (NROWS - 1) * MARGIN["hspace"] + MARGIN["top"]

ap = argparse.ArgumentParser()
ap.add_argument("--points", type=Path, default=root / "data" / "points.tsv")
ap.add_argument("--out", type=Path, default=None)
ap.add_argument("--preview", type=Path, default=None)
# flash lanes are real measurements and ship by default; --no-flash renders the peer-only twin for comparison.
ap.add_argument("--no-flash", action="store_true")
# paftools-mapeval batching (default): consecutive MAPQ levels with the same
# cumulative wrong count collapse into their last point (labelled by the first
# MAPQ), so a run of zero-error strata is one vertex instead of a vertical tail.
ap.add_argument("--no-batch", action="store_true")
a = ap.parse_args()
if a.out is None:
    a.out = root / "figures" / ("fig2_mapq_noflash.pdf" if a.no_flash else "fig2_mapq.pdf")


def order_key(run):
    """Legend/draw order: tool-major (hero first), then read arm HiFi/ONT/UL."""
    tool = run.split()[0]
    arm = run[run.index(" (") + 1:] if " (" in run else ""
    return (TOOL_ORDER.index(tool), ARM_RANK.get(arm, 0))


curves, criteria = {}, {}
with open(a.points, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["run"].split()[0] == "flash" and a.no_flash:
            continue
        criteria.setdefault(r["panel"], set()).add(r["criterion"])
        # File order is MAPQ 60 -> 0; never sort by x, the doubling-back is real.
        curves.setdefault((r["panel"], r["run"]), []).append(
            (float(r["err_rate"]), min(float(r["pct_aligned"]), Y_CAP),
             int(r["mapq"]), int(r["cum_wrong"])))
if not a.no_batch:
    for key, pts in curves.items():
        out, run_mapq = [], None
        for i, p in enumerate(pts):
            if run_mapq is None:
                run_mapq = p[2]
            if i + 1 == len(pts) or pts[i + 1][3] != p[3]:
                out.append((p[0], p[1], run_mapq, p[3]))
                run_mapq = None
        curves[key] = out
mixed = {p: sorted(c) for p, c in criteria.items() if len(c) > 1}
if mixed:
    raise SystemExit(f"{a.points}: panels mix criterion values: {mixed}")

figstyle.apply_house_style()
fig, axes = plt.subplots(NROWS, NCOLS, figsize=(W, H))
fig.subplots_adjust(left=MARGIN["left"] / W, right=1 - MARGIN["right"] / W,
                    bottom=MARGIN["bottom"] / H, top=1 - MARGIN["top"] / H,
                    wspace=MARGIN["wspace"] / PW, hspace=MARGIN["hspace"] / PH)

for ax in axes.flat[len(PANELS):]:
    ax.remove()

drawn = []   # (letter, title, legend kwargs, legend artist, [Line2D]) for the clash check
for ax, (title, letter, source, arms, xwin, ylim, ystep, colour, legend_kw) in zip(axes.flat, PANELS):
    if colour not in COLOUR_KEYS:
        raise SystemExit(f"{title!r}: unknown colour key {colour!r}; known: {sorted(COLOUR_KEYS)}")
    runs = sorted({r for p, r in curves if p == source
                   and (arms is None or any(r.endswith(a) for a in arms))},
                  key=order_key)
    if not runs:
        raise SystemExit(f"{a.points}: no curves for panel {title!r}")
    # x window: hand-set (the DNA panels, which share one window so they read
    # straight across) or, when xwin is None, the decade floor / ceiling over the
    # points this panel actually draws (post-batching, post---no-flash).
    xs = [p[0] for run in runs for p in curves[(source, run)]]
    if xwin is None:
        xlim = (10.0 ** math.floor(math.log10(min(xs))),
                10.0 ** math.ceil(math.log10(max(xs))))
    else:
        xlim = xwin
        out = [x for x in xs if not (xlim[0] <= x <= xlim[1])]
        if out:
            raise SystemExit(f"{title!r}: {len(out)} points outside the fixed x "
                             f"window {xlim} (first: {out[0]:g}); widen it in PANELS")
    print(f"{letter} {title}  [{next(iter(criteria[source]))}]  "
          f"x=[{xlim[0]:g}, {xlim[1]:g}]{'' if xwin is None else ' (fixed)'}"
          f"  drawn x=[{min(xs):.3g}, {max(xs):.3g}]  y={ylim}")
    lines = []
    for run in runs:
        pts = curves[(source, run)]
        tool = run.split()[0]
        ul = run.endswith("(UL)")
        # DNA panels: the title carries the platform, so the legend is the bare
        # tool name -- except the ultra-long arm, which shares the panel with
        # ONT and must say so. RNA panels keep "tool (arm)" for both arms.
        label = (run.replace(tool, TOOL_LABEL.get(tool, tool), 1) if colour == "tool_arm"
                 else TOOL_LABEL.get(tool, tool) + (" (UL)" if ul else ""))
        # Every vertex must sit inside the panel: with that guaranteed, the line
        # can never leave the frame and clip_on=False only lets marker BODIES
        # straddle the spines (gnuplot behaviour) instead of being bisected.
        # x can no longer be out of range (the window is derived from these
        # very points); the hard fail stays for y, which is panel-owned.
        out = [p for p in pts if not (ylim[0] <= p[1] <= ylim[1])]
        if out:
            raise SystemExit(f"{title!r} / {run!r}: {len(out)} points outside the "
                             f"panel y range (first: {out[0]}); widen PANELS before rendering")
        # "tool": strip the read arm so figstyle resolves the bare-label colour
        # table, giving one colour per tool across every DNA panel. "tool_arm":
        # pass the full label, keeping colour = tool x arm on the RNA panels.
        st = figstyle.style_for_run(run if colour == "tool_arm" else run.split(" (")[0])
        lines += ax.plot([p[0] for p in pts], [p[1] for p in pts], label=label,
                ls=UL_LINESTYLE if (colour == "tool" and ul) else "-",
                color=st["color"], lw=st["lw"], marker=st["marker"], ms=st["ms"],
                mfc="none", mew=st["mew"], zorder=st["zorder"], clip_on=False)
        q60 = next((p for p in pts if p[2] == 60), None)
        tail = f"Q60 err={q60[0]:.4g} kept={q60[1]:.2f}%" if q60 else "no MAPQ 60 row"
        print(f"    {run:<20s} n={len(pts):3d}  {tail}")

    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    # numticks: with a 7-decade window matplotlib's default thins the labels to
    # every other decade and the right-hand bound goes unlabelled; ask for enough
    # ticks that every decade of the fixed DNA window is named. The RNA panels
    # span 3 decades and are unaffected.
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(range(2, 10))))
    ax.xaxis.set_major_formatter(figstyle.log_decade_formatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.tick_params(which="both", top=True, right=True)
    ax.tick_params(labelsize=6)
    ax.set_xlabel(XLABEL, fontsize=6)
    ax.set_ylabel(YLABEL, fontsize=6)
    ax.set_title(title, fontsize=6, pad=3.5)
    figstyle.panel_letter(ax, letter)
    leg = ax.legend(**{**LEGEND_BASE, **legend_kw})
    drawn.append((letter, title, legend_kw, leg, lines))

# The legend must not sit on a curve. Measure it after the draw, in display
# space (a log axis is affine there, so a straight segment stays straight and
# linear sampling is exact), and die naming the panel if anything touches.
fig.canvas.draw()
clashes = []
for letter, title, legend_kw, leg, lines in drawn:
    box = leg.get_window_extent().padded(LEGEND_PAD * fig.dpi / 72.0)
    hits = 0
    for line in lines:
        xy = line.get_transform().transform(line.get_xydata())
        for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:]):
            n_s = max(2, int(math.hypot(x1 - x0, y1 - y0) / 2.0) + 1)
            hits += sum(box.contains(x0 + (x1 - x0) * t / (n_s - 1),
                                     y0 + (y1 - y0) * t / (n_s - 1))
                        for t in range(n_s))
    print(f"    legend {letter} {legend_kw}: "
          f"{'clears' if not hits else f'TOUCHES {hits} sampled points'}")
    if hits:
        clashes.append(f"{letter} {title!r} ({legend_kw})")
if clashes:
    raise SystemExit("legend sits on a curve in: " + "; ".join(clashes)
                     + " -- change that panel's legend field in PANELS "
                       "(another loc, or a shorter handlelength)")

figstyle.save_figure(fig, a.out, a.preview)
