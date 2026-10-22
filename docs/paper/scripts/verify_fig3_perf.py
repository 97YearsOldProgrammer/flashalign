#!/usr/bin/env python3
"""Mathematical acceptance test for Fig. 3.

Usage: python3 scripts/verify_fig3_perf.py   (after re-rendering fig3_perf.pdf)

Parses figures/fig3_perf.pdf with pdfplumber, recovers every drawn datum
from tick calibration, and matches it against data/perf.tsv.
The layout under test: a 2 x 2 grid, with ultra-long merged into the DNA
column as a third arm/colour, a frameless colour key INSIDE every panel
frame, and the tool names set diagonally under every panel.  Read-only with
respect to the paper folder; all scratch output goes to a fresh temp dir.

State under test: perf.tsv carries 19 rows -> 38 bars, including the five
measured flash lanes (DNA hifi/ont/ul + RNA Kinnex/ONT cDNA), so no tool of
either roster is a placeholder -- there is no "pending" and no "omitted"
slot in the manuscript figure, and the placeholder machinery is exercised
only by the contract renders (7e, 7g).  The speed ladder runs to 800 Gbp/hr
and the memory ladder to 70 GiB (the renderer's ROWS ceilings); winnowmap2
rows carry threads=64 (renderer THREAD_EXEMPT).  The hero tool is DRAWN as
"flashalign" (the renderer's DISPLAY map) while perf.tsv keys its rows
"flash".  28 non-flash bars, 4 placeholder slots in 7e/7g.
"""
import csv, itertools, json, math, re, shutil, subprocess, sys, tempfile
from pathlib import Path

import pdfplumber

PAPER   = Path(__file__).resolve().parent.parent
PDF     = PAPER / "figures/fig3_perf.pdf"
PERF    = PAPER / "data/perf.tsv"
SCRIPT  = PAPER / "scripts/make_fig3_perf.py"
# Every scratch copy of perf.tsv and every test render goes to a fresh temp
# dir, never into the repo; the manuscript figure's mtime is asserted unchanged.
SCRATCH = Path(tempfile.mkdtemp(prefix="verify_fig3_"))
print("scratch dir:", SCRATCH)

# --- from figstyle.ARM_FILL (violet / orange / sky; RNA green / blue) --------
HEX = {"violet": "#9400D4", "orange": "#D95E0D", "sky": "#57B5E8",
       "green": "#009E73", "blue": "#0073B2"}
def rgb(h): return tuple(round(int(h[i:i+2], 16) / 255.0, 4) for i in (1, 3, 5))
COLOR_OF = {rgb(v): k for k, v in HEX.items()}

# column key -> (panel, arms in the part's own order, colour of each arm).
# Ultra-long is the DNA column's third arm, not a column of its own.
COLUMNS = [
    ("DNA", "dna_long", ("hifi", "ont", "ul"),   ("violet", "orange", "sky")),
    ("RNA", "rna",      ("pb_kinnex", "ont_cdna"), ("green", "blue")),
]
# The key words that must appear INSIDE both panels of each column.
KEY_WORDS = {"DNA": {"HiFi", "ONT", "UL"}, "RNA": {"Kinnex", "ONT", "cDNA"}}
N_SWATCH  = {"DNA": 3, "RNA": 2}
TOOL_ORDER = ["flash", "minibwa", "minimap2", "rammap", "winnowmap2"]
HERO = "flash"
# The renderer's DISPLAY map: the hero tool is DRAWN as
# "flashalign" under the panels while perf.tsv still keys its rows "flash".
# Every tick-label comparison below is made in the drawn spelling, every row
# lookup in the tsv key, and TOOL_OF_LABEL is the one place the two meet.
DISPLAY = {"flash": "flashalign"}
LABEL_OF = {t: DISPLAY.get(t, t) for t in TOOL_ORDER}
TOOL_OF_LABEL = {v: k for k, v in LABEL_OF.items()}
# (row key, perf.tsv column, tick ceiling), as in the renderer's ROWS.
ROWS = [("speed", "gbp_per_hr", 800.0), ("memory", "peak_gb", 70.0)]
TOL_FRAC = 0.006                       # 0.6 % of the tick ceiling
PDF_MTIME = PDF.stat().st_mtime        # the manuscript figure must survive every contract test
PERF_MTIME = PERF.stat().st_mtime      # and so must perf.tsv: scratch copies only


def perf_git_lines():
    """`git status --short` lines that mention perf.tsv (may be non-empty).  What check 8
    asserts is that THIS RUN leaves that state exactly as it found it."""
    g = subprocess.run(["git", "status", "--short"], cwd=str(PAPER),
                       capture_output=True, text=True)
    return [l for l in g.stdout.splitlines() if "perf.tsv" in l]


PERF_GIT = perf_git_lines()

results = []
def check(name, ok, detail):
    results.append({"name": name, "pass": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + " :: " + detail)
    return ok

# --------------------------------------------------------------- perf.tsv ---
def read_perf(path):
    with open(path, newline="") as fh:
        body = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(body, delimiter="\t"))

def wanted(rows, drop_flash=False):
    """{(colkey, tool, arm): row} for the rows this figure must draw."""
    out = {}
    for r in rows:
        if drop_flash and r["tool"] == "flash":
            continue
        for ck, panel, arms, _c in COLUMNS:
            if r["panel"] == panel and r["readset"] in arms:
                out[(ck, r["tool"], r["readset"])] = r
    return out

# --------------------------------------------------------------- geometry ---
def norm_color(c):
    if c is None: return None
    if isinstance(c, (int, float)): return (round(float(c), 4),) * 3
    t = tuple(round(float(x), 4) for x in c)
    return t if len(t) == 3 else (round(float(t[0]), 4),) * 3

AXIS_ALIGNED = 1e-6      # |matrix[1]| below this = ordinary horizontal text

def char_angle(c):
    """Rotation of a char's text matrix, in degrees."""
    return math.degrees(math.atan2(c["matrix"][1], c["matrix"][0]))

def is_flat(c):
    """True for ordinary horizontal text.

    NOT the same as pdfplumber's `upright`: pdfminer's test is
    `0 < a*d and b*c <= 0`, which a +/-45 deg rotation also satisfies, so the
    diagonal tool names come back upright=True.  The text matrix is the honest
    witness, so every rotation test here reads the matrix."""
    return c["upright"] and abs(c["matrix"][1]) < AXIS_ALIGNED

def char_pt(c):
    """Type size of a char in points, at any rotation.

    pdfminer's bbox is the font's em box mapped through the text matrix, so
    for a rotation by t the bbox width is adv*|cos t| + size*|sin t|; invert
    that.  (Taking the bbox width raw would score a -45 deg 7 pt 'i' as 4.6.)"""
    ca, sa = abs(c["matrix"][0]), abs(c["matrix"][1])
    if sa < AXIS_ALIGNED:
        return c["size"]
    return ((c["x1"] - c["x0"]) - c["adv"] * ca) / sa

def rotated_groups(page, gap=8.0):
    """Rebuild rotated text.  pdfplumber's reading order for it is unreliable,
    so chars are bucketed by the rotation angle of their text matrix (90 deg
    for the bar numerals, the "pending" word and the y axis labels; -45 deg
    for the diagonal tool names), then chained along their own baseline: same
    perpendicular offset, successive advances closer than `gap`."""
    chars = [c for c in page.chars if c["text"].strip() and not is_flat(c)]
    out = []
    for r, it in itertools.groupby(sorted(chars, key=lambda c: round(char_angle(c))),
                                   lambda c: int(round(char_angle(c)))):
        th = math.radians(r); co, si = math.cos(th), math.sin(th)
        proj = lambda c: c["matrix"][4] * co + c["matrix"][5] * si
        perp = lambda c: -c["matrix"][4] * si + c["matrix"][5] * co
        cs = sorted(it, key=lambda c: (round(perp(c), 1), proj(c)))
        runs, run = [], []
        for c in cs:
            # |.| on the advance: two labels on one -45 deg line (a top-row
            # DNA name and a bottom-row RNA name can share x+y) sort by their
            # 0.1 pt perp difference, so the jump between them may be negative.
            if run and (abs(perp(c) - perp(run[-1])) > 0.5
                        or abs(proj(c) - proj(run[-1])) > gap):
                runs.append(run); run = []
            run.append(c)
        if run: runs.append(run)
        for g in runs:
            out.append(dict(text="".join(c["text"] for c in g), rot=r,
                            x0=min(c["x0"] for c in g), x1=max(c["x1"] for c in g),
                            y0=min(c["y0"] for c in g), y1=max(c["y1"] for c in g),
                            anchor=min(c["matrix"][4] for c in g),
                            size=max(char_pt(c) for c in g)))
    return out

def parse(pdf_path):
    """-> dict with page, panels, bars, ticks, texts."""
    pdf = pdfplumber.open(pdf_path)
    page = pdf.pages[0]
    H, Wd = page.height, page.width

    # panel frames = the white axes background patches (not the page patch)
    frames = [r for r in page.rects
              if r["fill"] and not r["stroke"] and norm_color(r["non_stroking_color"]) == (1.0, 1.0, 1.0)
              and not (r["x0"] < 1 and r["x1"] > Wd - 1)]
    xs = sorted({round(r["x0"], 2) for r in frames})
    ys = sorted({round(r["y0"], 2) for r in frames}, reverse=True)   # speed row on top
    assert len(xs) == 2 and len(ys) == 2, (xs, ys)
    panels = {}
    for ci, (ck, *_r) in enumerate(COLUMNS):
        for ri, (rk, _m, _c) in enumerate(ROWS):
            f = [r for r in frames if abs(r["x0"] - xs[ci]) < .01 and abs(r["y0"] - ys[ri]) < .01][0]
            panels[(ck, rk)] = dict(x0=f["x0"], x1=f["x1"], y0=f["y0"], y1=f["y1"])

    def which(x, y0):
        """A bar is a coloured rect standing ON a panel floor; the colour-key
        swatches share the fills but float inside the panel, well above it."""
        for k, p in panels.items():
            if p["x0"] - .01 <= x <= p["x1"] + .01 and abs(y0 - p["y0"]) < .05:
                return k
        return None

    def inside(x, y):
        for k, p in panels.items():
            if p["x0"] <= x <= p["x1"] and p["y0"] <= y <= p["y1"]:
                return k
        return None

    bars, swatches = [], []
    for r in page.rects:
        c = norm_color(r["non_stroking_color"])
        if not (r["fill"] and r["stroke"] and c in COLOR_OF):
            continue
        xc, yc = (r["x0"] + r["x1"]) / 2, (r["y0"] + r["y1"]) / 2
        k = which(xc, r["y0"])
        if k:
            bars.append(dict(key=k, arm_color=COLOR_OF[c], x0=r["x0"], x1=r["x1"],
                             y0=r["y0"], y1=r["y1"], xc=xc))
        else:
            swatches.append(dict(key=inside(xc, yc), arm_color=COLOR_OF[c], xc=xc, yc=yc))
    pend = [r for r in page.rects if r["stroke"] and not r["fill"]]

    # extract_words on a page stripped of every rotated char: pdfminer would
    # otherwise hand the -45 deg tool names back as "upright" words.
    flat = page.filter(lambda o: o["object_type"] != "char" or is_flat(o))
    up = [dict(text=w["text"], x0=w["x0"], x1=w["x1"], y0=H - w["bottom"], y1=H - w["top"],
               size=w["size"], rot=0)
          for w in flat.extract_words(extra_attrs=["size", "upright"])]
    rot = rotated_groups(page)

    # y tick marks: 2.2 pt inward stubs on each panel's left spine
    ticks = {}
    for k, p in panels.items():
        ticks[k] = sorted(l["y0"] for l in page.lines
                          if abs(l["y1"] - l["y0"]) < .01 and abs(l["x0"] - p["x0"]) < .01
                          and abs((l["x1"] - l["x0"]) - 2.2) < .05
                          and p["y0"] - .05 <= l["y0"] <= p["y1"] + .05)
    return dict(pdf=pdf, page=page, W=Wd, H=H, panels=panels, bars=bars,
                swatches=swatches, pending=pend, up=up, rot=rot, ticks=ticks)

def calibrate(d):
    """value = m*y + c per panel.  The y tick NUMERALS name the ticked values;
    the tick MARKS give their exact y (a numeral's bbox centre carries a
    constant ~0.6 pt alignment offset that would bias every bar)."""
    cal, info = {}, {}
    for k, p in d["panels"].items():
        nums = sorted(((w["y0"] + w["y1"]) / 2, int(w["text"])) for w in d["up"]
                      if re.fullmatch(r"\d+", w["text"]) and w["x1"] <= p["x0"] + .5
                      and w["x1"] > p["x0"] - 30
                      and p["y0"] - 4 <= (w["y0"] + w["y1"]) / 2 <= p["y1"] + 4)
        ys = d["ticks"][k]
        assert len(nums) == len(ys), (k, nums, ys)
        vals = [v for _y, v in nums]
        n = len(ys)
        sy, sv = sum(ys), sum(vals)
        syy = sum(y * y for y in ys); syv = sum(y * v for y, v in zip(ys, vals))
        m = (n * syv - sy * sv) / (n * syy - sy * sy); c = (sv - m * sy) / n
        res = max(abs(m * y + c - v) for y, v in zip(ys, vals))
        # slope from the numeral centres alone, as a cross-check
        ny = [y for y, _v in nums]
        m2 = ((n * sum(y * v for y, v in zip(ny, vals)) - sum(ny) * sv) /
              (n * sum(y * y for y in ny) - sum(ny) ** 2))
        cal[k] = (m, c)
        info[k] = (vals, res, m, m2)
    return cal, info

def tool_labels(d):
    """Diagonal x tick labels -> {(colkey, rowkey): [(tool, x anchor, drawn), ...]}.

    They are set with rotation=-45, ha='left', rotation_mode='anchor', so the
    left edge of the group's bbox IS the tick position.  A label belongs to
    the nearest panel ABOVE it in the same x range -- which is how "under
    every panel, both rows" is tested.  `drawn` is the spelling on the page
    ("flashalign" for the hero); `tool` is the perf.tsv key that label names,
    so every downstream row lookup keeps working off the tsv vocabulary."""
    out = {k: [] for k in d["panels"]}
    for g in d["rot"]:
        if g["text"] not in TOOL_OF_LABEL or g["rot"] == 90:
            continue
        cand = [(p["y0"], k) for k, p in d["panels"].items()
                if p["x0"] - 2 <= g["x0"] <= p["x1"] + 2 and p["y0"] >= g["y1"] - 1]
        if cand:
            out[min(cand)[1]].append((TOOL_OF_LABEL[g["text"]], g["anchor"], g["text"]))
    return {k: sorted(v, key=lambda t: t[1]) for k, v in out.items()}

def recover(d):
    """-> {(colkey, tool, arm, rowkey): [(value, bar), ...]} for every drawn bar."""
    cal, _ = calibrate(d)
    labs = tool_labels(d)
    got, unassigned = {}, []
    for b in d["bars"]:
        ck, rk = b["key"]
        arms, colors = dict((c[0], (c[2], c[3])) for c in COLUMNS)[ck]
        arm = arms[colors.index(b["arm_color"])] if b["arm_color"] in colors else None
        if arm is None: unassigned.append(b); continue
        cand = min(labs[(ck, rk)], key=lambda t: abs(t[1] - b["xc"]))
        m, c = cal[(ck, rk)]
        got.setdefault((ck, cand[0], arm, rk), []).append((m * b["y1"] + c, b))
    return got, unassigned, cal, labs

# ============================== 1. page size ================================
d = parse(PDF)
ok = abs(d["W"] - 498.90) <= 0.1 and d["H"] <= 324
check("1 page size", ok,
      f'width {d["W"]:.4f} pt (target 498.90 +/- 0.1; 6.9291 in = 498.8952), '
      f'height {d["H"]:.2f} pt (<= 324 pt = 4.5 in), 1 page')

# ============================== 2. bars =====================================
rows = read_perf(PERF)
want = wanted(rows)
got, unassigned, cal, labs = recover(d)
_, calinfo = calibrate(d)
prob = []
for k, (vals, res, m, m2) in sorted(calinfo.items(), key=str):
    if res > 0.02 or abs(m - m2) / abs(m) > 2e-3:
        prob.append(f"{k} residual {res:.4f} slope {m:.6f} vs numerals {m2:.6f}")
matched, bad = set(), []
for rk, metric, ceil in ROWS:
    tol = TOL_FRAC * ceil
    for (ck, tool, arm), r in want.items():
        exp = float(r[metric])
        hit = got.get((ck, tool, arm, rk))
        if not hit:
            bad.append(f"no bar for {ck}/{tool}/{arm}/{rk}"); continue
        if len(hit) > 1:
            bad.append(f"{len(hit)} bars for {ck}/{tool}/{arm}/{rk}"); continue
        v = hit[0][0]
        if abs(v - exp) > tol:
            bad.append(f"{ck}/{tool}/{arm}/{rk}: drawn {v:.3f} vs tsv {exp} (tol {tol:.2f})")
        else:
            matched.add((ck, tool, arm, rk))
extra = [k for k in got if k not in matched]
worst, worst_abs = 0.0, ("", 0.0)
for rk, metric, ceil in ROWS:
    for (ck, tool, arm), r in want.items():
        h = got.get((ck, tool, arm, rk))
        if h and len(h) == 1:
            e = abs(h[0][0] - float(r[metric]))
            if e > worst_abs[1]: worst_abs = (f"{ck}/{tool}/{arm}/{rk}", e)
            worst = max(worst, e / ceil)
# order checks: tools left to right, and arms left to right inside a cluster
order_bad = []
for ck, panel, arms, colors in COLUMNS:
    for rk, _m, _c in ROWS:
        seq = [t for t, *_rest in labs[(ck, rk)]]
        exp_seq = [t for t in TOOL_ORDER if t in seq]
        if seq != exp_seq:
            order_bad.append(f"{ck}/{rk} tool order {seq} != {exp_seq}")
        for tool in seq:
            cl = sorted([b for (c2, t2, a2, r2), lst in got.items() if (c2, t2, r2) == (ck, tool, rk)
                         for _v, b in lst], key=lambda b: b["xc"])
            seen = [arms[colors.index(b["arm_color"])] for b in cl]
            if seen and seen != list(arms)[:len(seen)]:
                order_bad.append(f"{ck}/{rk}/{tool} arm order {seen} != {list(arms)}")
ok = (not bad and not extra and not unassigned and not prob and not order_bad
      and len(matched) == 2 * len(want) and len(d["bars"]) == 2 * len(want))
check("2 bars vs perf.tsv", ok,
      f'{len(d["bars"])} filled bars, {len(want)} tsv rows x 2 metrics = {2*len(want)} '
      f'expected; all matched 1:1, worst error {worst*100:.5f}% of the tick ceiling '
      f'(= {worst_abs[1]:.2e} data units at {worst_abs[0]}) '
      f'(tol 0.600%); calibration residual <= '
      f'{max(v[1] for v in calinfo.values()):.4f} in data units; tool order per panel '
      f'{ {f"{k[0]}/{k[1]}": [t for t,*_ in v] for k,v in labs.items()} }; arm order '
      f'hifi<ont<ul, pb_kinnex<ont_cdna confirmed by x'
      + ("" if ok else "; PROBLEMS: " + "; ".join(
          bad + [f"unmatched drawn bar {e}" for e in extra]
          + [f"bar of unknown colour at x={b['xc']:.1f}" for b in unassigned]
          + prob + order_bad)))

# ============================== 3. numerals =================================
num_bad, num_ok = [], 0
for (ck, tool, arm, rk), lst in got.items():
    v, b = lst[0]
    r = want[(ck, tool, arm)]
    exp = f'{float(r[dict((x[0], x[1]) for x in ROWS)[rk]]):.1f}'
    cands = [g for g in d["rot"] if g["rot"] == 90
             and abs((g["x0"] + g["x1"]) / 2 - b["xc"]) <= 1.5
             and g["y0"] >= b["y1"] - 0.5 and g["y0"] <= b["y1"] + 6]
    if not cands:
        num_bad.append(f"{ck}/{tool}/{arm}/{rk}: no rotated numeral within 1.5 pt of x={b['xc']:.2f}")
    elif cands[0]["text"] != exp:
        num_bad.append(f"{ck}/{tool}/{arm}/{rk}: numeral {cands[0]['text']!r} != {exp!r}")
    else:
        num_ok += 1
check("3 printed numerals", not num_bad,
      f"{num_ok}/{len(got)} bars carry a rotated numeral equal to f'{{v:.1f}}' of the "
      f"tsv value, x-centre within 1.5 pt of the bar centre"
      + ("" if not num_bad else "; PROBLEMS: " + "; ".join(num_bad)))

# ============================== 4. placeholders =============================
# Every roster tool of both parts now has measured rows, so the manuscript
# figure must carry NO dashed slot at all -- neither "pending" (never
# measured) nor "omitted" (dropped by --no-flash) -- and the hero must be a
# real tick label under both columns.  The slot machinery itself is still
# exercised, but only by the contract renders 7e (pending) and 7g (omitted).
def pending_slots(d, word="pending"):
    out = []
    for g in d["rot"] + d["up"]:
        if g["text"] == word:
            xc = (g["x0"] + g["x1"]) / 2; yc = (g["y0"] + g["y1"]) / 2
            for k, p in d["panels"].items():
                if p["x0"] <= xc <= p["x1"] and p["y0"] <= yc <= p["y1"]:
                    out.append(k)
    return sorted(out)
slots = pending_slots(d)
omitted = pending_slots(d, "omitted")
hero_cols = sorted({k[0] for k, v in labs.items() if HERO in [t for t, *_ in v]})
hero_drawn = sorted({g for v in labs.values() for t, _x, g in v if t == HERO})
ok = (not slots and not omitted and not d["pending"]
      and hero_cols == ["DNA", "RNA"] and hero_drawn == [LABEL_OF[HERO]])
check("4 no placeholder slots", ok,
      f'"pending" appears {len(slots)} times at {slots} and "omitted" {len(omitted)} '
      f'times at {omitted} (expected 0 + 0: every roster tool of both parts has '
      f'measured rows); {len(d["pending"])} dashed unfilled slot rects (expected 0); '
      f'the hero tool {HERO!r} is a real x tick label, drawn {hero_drawn}, in columns '
      f'{hero_cols}')

# ============================== 5. typography ===============================
sizes = [(char_pt(c), c["text"]) for c in d["page"].chars if c["text"].strip()]
mn = min(sizes)
pf = subprocess.run(["pdffonts", str(PDF)], capture_output=True, text=True).stdout
t3 = [l for l in pf.splitlines()[2:] if l.strip() and "Type 3" in l]
boxes = d["up"] + d["rot"]
ov = []
for a, b in itertools.combinations(boxes, 2):
    dx = min(a["x1"], b["x1"]) - max(a["x0"], b["x0"])
    dy = min(a["y1"], b["y1"]) - max(a["y0"], b["y0"])
    if dx > 0.3 and dy > 0.3:
        ov.append(f'{a["text"]!r}/{b["text"]!r} by {min(dx,dy):.2f} pt')
ok = mn[0] >= 5.9 and not t3 and not ov
check("5 typography", ok,
      f'min char size {mn[0]:.3f} pt (glyph {mn[1]!r}) >= 5.9; fonts: '
      f'{ "; ".join(l.split()[0]+" "+l.split()[1] for l in pf.splitlines()[2:] if l.strip()) }; '
      f'no Type 3: {not t3}; {len(boxes)} text boxes, overlapping pairs > 0.3 pt: '
      f'{len(ov)}{"" if not ov else " -> " + "; ".join(ov[:6])}')

# ============================== 5b. diagonal tool names =====================
name_bad = []
ROSTER = {"DNA": TOOL_ORDER, "RNA": ["flash", "minimap2"]}   # the renderer's PARTS rosters
# compared as DRAWN: the hero's tick label is spelled "flashalign" (DISPLAY)
DRAWN_ROSTER = {ck: [LABEL_OF[t] for t in r] for ck, r in ROSTER.items()}
for ck, panel, arms, colors in COLUMNS:
    for rk, _m, _c in ROWS:
        seen = [g for _t, _x, g in labs[(ck, rk)]]
        if seen != DRAWN_ROSTER[ck]:
            name_bad.append(f"tool names under {ck}/{rk} are {seen}, roster is "
                            f"{DRAWN_ROSTER[ck]}")
rots = sorted({g["rot"] for g in d["rot"] if g["text"] in TOOL_OF_LABEL})
if rots != [-45]:
    name_bad.append(f"tool-name text matrices are at {rots} deg, not -45")
n_names = sum(len(v) for v in labs.values())
check("5b diagonal tool names under every panel", not name_bad,
      f'{n_names} tool names (5 DNA + 2 RNA, per row, x 2 rows = 14 expected), the '
      f'hero drawn as {LABEL_OF[HERO]!r} and matched back to the perf.tsv key '
      f'{HERO!r}; every panel carries its own row of names; every name char rides a rotated text '
      f'matrix at {rots} deg (pdfminer still calls a 45 deg matrix "upright", so the '
      f'matrix, not the flag, is what is tested)'
      + ("" if not name_bad else "; PROBLEMS: " + "; ".join(name_bad)))

# ============================== 6. colour key ===============================
key_bad = []
in_panel = {k: set() for k in d["panels"]}
for w in d["up"]:
    xc, yc = (w["x0"] + w["x1"]) / 2, (w["y0"] + w["y1"]) / 2
    for k, p in d["panels"].items():
        if p["x0"] <= xc <= p["x1"] and p["y0"] <= yc <= p["y1"]:
            in_panel[k].add(w["text"])
sw_per = {k: 0 for k in d["panels"]}
for s in d["swatches"]:
    if s["key"] is None: key_bad.append(f"key swatch at x={s['xc']:.1f} is in no panel")
    else: sw_per[s["key"]] += 1
for ck, panel, arms, colors in COLUMNS:
    for rk, _m, _c in ROWS:
        miss = KEY_WORDS[ck] - in_panel[(ck, rk)]
        if miss: key_bad.append(f"{ck}/{rk} key is missing {sorted(miss)}")
        if sw_per[(ck, rk)] != N_SWATCH[ck]:
            key_bad.append(f"{ck}/{rk} has {sw_per[(ck,rk)]} swatches, want {N_SWATCH[ck]}")
ont_cdna = all(any(abs(a["x1"] - b["x0"]) < 5 and abs(a["y0"] - b["y0"]) < .5
                   for a in d["up"] if a["text"] == "ONT"
                   for b in d["up"] if b["text"] == "cDNA"
                   and d["panels"][(ck, rk)]["x0"] <= b["x0"] <= d["panels"][(ck, rk)]["x1"]
                   and d["panels"][(ck, rk)]["y0"] <= b["y0"] <= d["panels"][(ck, rk)]["y1"])
               for ck in ("RNA",) for rk, _m, _c in ROWS)
if not ont_cdna: key_bad.append('"ONT cDNA" is not an adjacent pair in both RNA panels')
above = [w["text"] for w in d["up"] if w["y0"] > max(p["y1"] for p in d["panels"].values())]
letters = [w["text"] for w in d["up"] if re.fullmatch(r"\([a-h]\)", w["text"])]
if letters:
    key_bad.append(f"panel letters present {letters} (expected none)")
if above:
    key_bad.append(f"text above the top row {sorted(set(above))} (expected no column titles)")
ok = not key_bad and len(d["swatches"]) == 10
check("6 colour key inside every panel", ok,
      f'in-frame key words per panel: '
      f'{json.dumps({f"{k[0]}/{k[1]}": sorted(v & (KEY_WORDS["DNA"] | KEY_WORDS["RNA"])) for k, v in in_panel.items()})}; '
      f'{len(d["swatches"])} key swatches (3+3 DNA, 2+2 RNA = 10), all inside a panel '
      f'frame; "ONT cDNA" adjacent in both RNA panels: {ont_cdna}; no panel letters; nothing '
      f'drawn above the top row (no column titles): {sorted(set(above)) or "confirmed"}'
      + ("" if ok else "; PROBLEMS: " + "; ".join(key_bad)))

# ============================== 7. contract tests ===========================
def run(perf, *flags, out="x.pdf"):
    o = SCRATCH / out
    cmd = [sys.executable, str(SCRIPT), "--perf", str(perf), "--out", str(o),
           "--preview", str(SCRATCH / (Path(out).stem + ".png"))] + list(flags)
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SCRATCH))
    return p.returncode, (p.stdout + p.stderr), o

raw = PERF.read_text()
head = [l for l in raw.splitlines() if l.startswith("#")]
hdr  = [l for l in raw.splitlines() if not l.startswith("#")][0]
body = [l for l in raw.splitlines() if not l.startswith("#")][1:]
def write(name, lines):
    p = SCRATCH / name
    p.write_text("\n".join(head + [hdr] + lines) + "\n")
    return p
sub = []

# (a) unknown tool
ia = next(k for k, l in enumerate(body) if "\tminibwa\t" in l)
f = write("t_a.tsv", body[:ia] + [body[ia].replace("\tminibwa\t", "\tbogus\t")] + body[ia + 1:])
rc, out, _ = run(f, out="a.pdf")
sub.append(("7a unknown tool", rc != 0 and "bogus" in out,
            f"exit {rc}; {out.strip().splitlines()[-1][:150] if out.strip() else '(no output)'}"))
# (b) out of range
r0 = body[0].split("\t"); r0[8] = "999"
f = write("t_b.tsv", ["\t".join(r0)] + body[1:])
rc, out, _ = run(f, out="b.pdf")
last = out.strip().splitlines()[-1] if out.strip() else ""
sub.append(("7b out-of-range value",
            rc != 0 and "999" in out
            and ("axis" in out or f"[0, {ROWS[0][2]:.0f}]" in out),
            f"exit {rc}; {last[:180]}"))
# (c) duplicate row
f = write("t_c.tsv", body + [body[0]])
rc, out, _ = run(f, out="c.pdf")
last = out.strip().splitlines()[-1] if out.strip() else ""
sub.append(("7c duplicate row", rc != 0 and "duplicate" in out.lower(), f"exit {rc}; {last[:180]}"))
# (d) missing winnowmap2/ul row: UL is now an ARM of the DNA part, so a tool
#     that covers hifi+ont but not ul is an incomplete lane set -> hard error.
f = write("t_d.tsv", [l for l in body if not l.startswith("dna_long\tul\twinnowmap2")])
rc, out, o = run(f, out="d.pdf")
last = out.strip().splitlines()[-1] if out.strip() else ""
sub.append(("7d winnowmap2 missing its UL arm -> hard error",
            rc != 0 and "winnowmap2" in out and not o.exists(),
            f"exit {rc}; {last[:190]}"))
# (e) the never-measured placeholder path.  All five flash lanes are measured
#     now, so the only way left to exercise it is to take the rows away: the
#     hero must come back as a dashed "pending" slot (NOT "omitted": nothing
#     was dropped by a flag) in all four panels, the 5 rows x 2 metrics must
#     vanish from the bar count, and every surviving row must still be drawn
#     at its tsv value.  Appending the rows back is what the manuscript figure
#     shows, which is exactly what check 4 asserts.
new = [l for l in body if l.split("\t")[2] != "flash"]
f = write("t_e.tsv", new)
rc, out, o = run(f, out="e.pdf")
det = f"exit {rc}"
okE = rc == 0 and o.exists()
if okE:
    de = parse(o)
    goe, un_e, _c, _l = recover(de)
    ps, om = pending_slots(de), pending_slots(de, "omitted")
    we = wanted(read_perf(f))
    errs = []
    for rk, metric, ceil in ROWS:
        for (ck, tool, arm), r in we.items():
            h = goe.get((ck, tool, arm, rk))
            if not h or len(h) != 1:
                errs.append(f"{ck}/{tool}/{arm}/{rk}: {0 if not h else len(h)} bars")
            elif abs(h[0][0] - float(r[metric])) > TOL_FRAC * ceil:
                errs.append(f"{ck}/{tool}/{arm}/{rk} drawn {h[0][0]:.3f} vs {r[metric]}")
    errs += [f"bar of unknown colour at x={b['xc']:.1f}" for b in un_e]
    errs += [f"drawn bar {k} matches no row" for k in goe if k[:3] not in we]
    exp_slots = sorted(de["panels"])
    okE = (not errs and ps == exp_slots and not om and len(de["pending"]) == 4
           and len(de["bars"]) == 2 * len(we) == 28)
    det = (f'exit 0; flash rows removed -> "pending" {len(ps)} times at {ps} (all four '
           f'panels), "omitted" {len(om)} times, {len(de["pending"])} dashed slot rects '
           f'(4 expected); {len(de["bars"])} bars for {len(we)} rows (28 expected), each '
           f'still recovered at its tsv value'
           + ("" if not errs else "; PROBLEMS: " + "; ".join(errs)))
sub.append(("7e flash rows removed -> pending slots", okE, det))
# (f) --no-values
rc, out, o = run(PERF, "--no-values", out="f.pdf")
okF, detF = rc == 0 and o.exists(), f"exit {rc}"
if okF:
    df = parse(o)
    nums = [g["text"] for g in df["rot"] if re.fullmatch(r"[\d.]+", g["text"])]
    okF = not nums
    detF = (f'exit 0; rotated numeral groups: {len(nums)}{"" if not nums else " -> "+str(nums[:8])}; '
            f'{len(df["bars"])} bars still drawn; "pending" still present '
            f'{len(pending_slots(df))}x')
sub.append(("7f --no-values", okF, detF))
# (g) --no-flash
rc, out, o = run(PERF, "--no-flash", out="g.pdf")
okG, detG = rc == 0 and o.exists(), f"exit {rc}"
if okG:
    dg = parse(o)
    ps = pending_slots(dg); om = pending_slots(dg, "omitted")
    okG = not ps and om == sorted(dg["panels"]) and len(dg["bars"]) == 28
    detG = (f'exit 0; "pending" appears {len(ps)} times at {ps}, "omitted" {len(om)} times at {om} '
            f'(expected 0 + 4: every flash lane is measured now, so --no-flash marks all '
            f'four slots omitted and none pending); {len(dg["bars"])} bars '
            f'(28 = the 14 non-flash rows x 2 metrics)')
sub.append(("7g --no-flash", okG, detG))
# (h) --no-flash with the default --out is refused
ph = subprocess.run([sys.executable, str(SCRIPT), "--no-flash", "--preview", str(SCRATCH / "h.png")],
                    capture_output=True, text=True, cwd=str(SCRATCH))
okH = ph.returncode != 0 and "--out" in (ph.stdout + ph.stderr) and PDF.stat().st_mtime == PDF_MTIME
sub.append(("7h --no-flash refuses the default --out", okH,
            f"exit {ph.returncode}; {(ph.stdout + ph.stderr).strip()[:160]}; figure mtime unchanged: {PDF.stat().st_mtime == PDF_MTIME}"))

for n, o_, dt in sub:
    check(n, o_, dt)

# ============================== 8. git hygiene ==============================
# perf.tsv may carry uncommitted edits, so "perf.tsv is clean" is not the
# right assertion; what this run must prove is that IT changed nothing -- same bytes, same mtime, same
# git status lines as when it started, every scratch copy in the temp dir.
lines = perf_git_lines()
before = PERF.read_text()
same_git = lines == PERF_GIT
same_mtime = PERF.stat().st_mtime == PERF_MTIME
check("8 perf.tsv untouched", raw == before and same_git and same_mtime,
      f'file bytes unchanged during this run: {raw == before}; mtime unchanged: '
      f'{same_mtime}; git status --short in {PAPER} reports {lines or ["(clean)"]} for '
      f'perf.tsv, the same as before the contract renders: {same_git}; every scratch '
      f'copy went to {SCRATCH}')

allpass = all(r["pass"] for r in results)
print("\nALL PASS" if allpass else "\nSOME CHECKS FAILED")
(SCRATCH / "verdict.json").write_text(json.dumps({"checks": results, "all_pass": allpass}, indent=2))
