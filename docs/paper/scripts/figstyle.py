#!/usr/bin/env python3
"""House figure style: the minibwa Fig. 2 rendering conventions, in matplotlib.

Every size here comes from vector forensics of minibwa's fig-roc.pdf measured at
its printed scale (0.6001), so the numbers are *printed* points -- valid only if
the figure is included at width=\\textwidth with no further scaling.
Callers must set matplotlib.use("Agg") before importing pyplot.
"""
import math
from pathlib import Path

import matplotlib
from matplotlib import font_manager as fm
from matplotlib.ticker import FuncFormatter

SANS = ["Helvetica", "Helvetica Neue", "Arial", "DejaVu Sans"]


def _has_bold_face(name):
    for f in fm.fontManager.ttflist:
        w = f.weight if isinstance(f.weight, int) else fm.weight_dict.get(f.weight, 400)
        if f.name == name and f.style == "normal" and w >= 600:
            return True
    return False


# macOS ships Helvetica as a .ttc whose bold subface matplotlib never registers,
# so weight='bold' silently renders regular. Demote such families for bold text;
# Arial is metrically identical to Helvetica, so the fallback is invisible.
BOLD_SANS = [f for f in SANS if _has_bold_face(f)] + SANS

# The Okabe-Ito six, in the exact hexes minibwa strokes.
PALETTE = {
    "violet": "#9400D4",
    "green": "#009E73",
    "sky": "#57B5E8",
    "blue": "#0073B2",
    "orange": "#E69E00",
    "red": "#E61F0F",
}

BASE_LW, BASE_MS, BASE_MEW = 0.4, 1.9, 0.32

# Marker shape encodes the tool and is shared between arms; colour encodes the
# curve (tool x arm). This is exactly minibwa Fig. 2e/2f.
_MARKER = {"minibwa": "s", "minimap2": "o", "winnowmap2": "^", "flash": "*"}
_COLOR = {
    "hifi": {"minibwa": PALETTE["violet"], "minimap2": PALETTE["green"],
             "winnowmap2": PALETTE["sky"], "flash": "#000000"},
    "ont": {"minibwa": PALETTE["blue"], "minimap2": PALETTE["orange"],
            "winnowmap2": PALETTE["red"], "flash": "#CC79A7"},
}


def apply_house_style():
    """Set the structural rcParams. Font sizes stay with the individual figure."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": SANS,
        # mathtext for the 10^n tick labels: Helvetica, upright, italic-free.
        "mathtext.fontset": "custom",
        "mathtext.rm": SANS[0],
        "mathtext.it": SANS[0],
        "mathtext.bf": SANS[0],
        "mathtext.sf": SANS[0],
        "mathtext.default": "regular",
        "axes.unicode_minus": False,
        "axes.linewidth": 0.35,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        "xtick.major.width": 0.3, "ytick.major.width": 0.3,
        "xtick.minor.size": 1.0, "ytick.minor.size": 1.0,
        "xtick.minor.width": 0.25, "ytick.minor.width": 0.25,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def style_for_run(run_label):
    """Resolve a points.tsv run label ('minimap2 (ONT)', 'minibwa') to a line style."""
    tool = run_label.split()[0]
    if tool not in _MARKER:
        raise ValueError(f"unknown tool {tool!r} in run label {run_label!r}; "
                         f"known tools: {sorted(_MARKER)}")
    # A bare label is a UL panel, where only one arm exists: use the HiFi colours.
    arm = "ont" if run_label.endswith("(ONT)") else "hifi"
    hero = tool == "flash"
    return {
        "color": _COLOR[arm][tool],
        "marker": _MARKER[tool],
        "lw": BASE_LW * (1.3 if hero else 1.0),
        "ms": BASE_MS * (1.25 if hero else 1.0),
        "mew": BASE_MEW,
        "zorder": 5 if hero else 2,
        "hero": hero,
    }


def log_decade_formatter():
    """10^n major tick labels; \\mathdefault + mathtext.default='regular' keeps
    them in the sans text font rather than the mathtext italic."""
    return FuncFormatter(
        lambda v, _: rf"$\mathdefault{{10^{{{int(round(math.log10(v)))}}}}}$")


def panel_letter(ax, s, dx=-24.0, dy=4.5, fontsize=8.5):
    """Bold letter above-left of the frame, in the y-axis-label column."""
    ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy),
                textcoords="offset points", family=BOLD_SANS, fontweight="bold",
                fontsize=fontsize, ha="left", va="baseline", annotation_clip=False)


# ---- Grouped-bar figures (minibwa fig-perf grammar, printed scale 0.711522) ----
# Arm fills for performance bars. ONT is fig-perf's #D95E0D -- deliberately NOT
# fig-roc/Fig 2's red #E61F0F. One 7 pt size for every text role, no bold.
# Every read set has its own colour: the RNA arms take two Okabe-Ito hues the DNA
# arms do not use, so a Kinnex bar never reads as a HiFi bar.
ARM_FILL = {"hifi": PALETTE["violet"], "ont": "#D95E0D", "ul": PALETTE["sky"],
            "pb_kinnex": PALETTE["green"], "ont_cdna": PALETTE["blue"]}
ARM_LABEL = {"hifi": "HiFi", "ont": "ONT", "ul": "UL",
             "pb_kinnex": "Kinnex", "ont_cdna": "ONT cDNA"}
BAR_EDGE_LW = 0.18   # bar outlines, legend-sample edges, y tick width
BAR_TICK_LEN = 2.2
BAR_FONT = 7
BAR_WIDTH_FRAC = 0.858  # of the unit slot; the rest is the in-group gap


def bar_offsets(n_arms):
    """Clustered-bar centre offsets in group-pitch units, plus the bar width.
    gnuplot histogram rule: unit slot = pitch/(n_arms+2), bars on adjacent
    slots centred in the group, 0.858 slot wide (bars NOT butted)."""
    slot = 1.0 / (n_arms + 2)
    return [slot * (i - (n_arms - 1) / 2.0) for i in range(n_arms)], slot * BAR_WIDTH_FRAC


def style_bar_axes(ax):
    """fig-perf axes grammar: L-spine only (left+bottom, rcParams 0.35 width),
    y majors mirrored as bare inward stubs on the spineless right edge, no
    minors, and no x tick marks at all (tool names are labels only)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", which="major", direction="in", left=True, right=True,
                   length=BAR_TICK_LEN, width=BAR_EDGE_LW, labelsize=BAR_FONT)
    ax.tick_params(axis="y", which="minor", left=False, right=False)
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelsize=BAR_FONT)


def save_figure(fig, pdf_path, preview_png=None, preview_dpi=300):
    """Vector PDF plus an optional raster preview. No bbox_inches: the figure
    width is load-bearing (it must equal \\textwidth exactly)."""
    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path)
    print("wrote", pdf_path)
    if preview_png:
        preview_png = Path(preview_png)
        preview_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview_png, dpi=preview_dpi)
        print("wrote", preview_png)
    return pdf_path
