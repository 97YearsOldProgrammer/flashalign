#!/usr/bin/env python3
"""Emit the mapeval curve-points TSV for the paper's Fig 2 from lane dirs.

Bridges the mapq-curve lane outputs into the repo harness: per-MAPQ Q-rows are
rebuilt from mapeval.json's mapq_hist_correct/mapq_hist_wrong (the json "table"
is the paftools-batched Q-table, too coarse to plot), cross-checked against
mapeval.strata.tsv, then run through bench.mapeval_collect.derive_points (the
eval2roc.pl port: 0.5 pseudo-count, thin-100, EOF flush) so the numbers are
directly comparable to the published minimap2/minibwa figures.

Panels: one per reference, curves labelled "tool (read set)".
Per-arm tool roster: HiFi arms carry flash, minibwa,
minimap2 (map-hifi) and winnowmap2; the ONT and ultra-long arms carry flash,
minimap2 lr:hq and winnowmap2 only -- minimap2's map-ont lanes and minibwa's
ONT/UL lanes are skipped here, so points.tsv holds exactly what Fig. 2 draws.

Usage: python3 fig2_points.py --lanes DIR --bench-repo REPO --out points.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# (ref, readsets, title): one panel group per reference, every curve labelled
# "tool (read set)". The ultra-long arm rides along in the genomic group:
# the renderer draws it inside the ONT panel as a
# second arm, so there is no standalone UL panel and no bare tool label.
PANELS = [
    ("chm13", ("hifi", "ont", "ul"), "Simulated long reads to CHM13"),
    ("hg002dip", ("hifi", "ont", "ul"), "Simulated long reads to HG002"),
    ("human", ("rnahifi", "rnaont"), "Simulated long RNA reads to human"),
    ("mouse", ("rnahifi", "rnaont"), "Simulated long RNA reads to mouse"),
]
READSET_LABEL = {"hifi": "HiFi", "ont": "ONT", "ul": "UL",
                 "rnahifi": "HiFi", "rnaont": "ONT"}
# run order within a panel == legend + style-assignment order
TOOL_ORDER = ["flash", "minibwa", "minimap2", "winnowmap2"]


def parse_lane(name: str):
    """'h-06-mapq-ul-hg002-minibwa' -> (ref, readset, tool) or None."""
    m = re.match(r"[mhce]-\d+-mapq-(.+)$", name)
    if not m:
        return None
    toks = m.group(1).split("-")
    if toks[0] == "ul":
        readset, ref, tool_toks = "ul", toks[1], toks[2:]
    else:
        ref, readset, tool_toks = toks[0], toks[1], toks[2:]
    ref = {"chm13": "chm13", "hg002": "hg002dip", "hg002dip": "hg002dip",
           "human": "human", "mouse": "mouse"}[ref]
    tool = "-".join(tool_toks)
    if tool == "minimap2-ont-hq":
        tool = "minimap2"
    elif tool == "minimap2" and readset in ("ont", "ul"):
        return None  # map-ont arm, not shown -- the ONT/UL arms are lr:hq
    elif tool == "minibwa" and readset in ("ont", "ul"):
        return None  # minibwa is a HiFi-arm peer only
    if tool.startswith("flashalign"):
        tool = "flash"
    return ref, readset, tool


def qrows_from_json(doc, mk_qrow):
    hc, hw = doc["mapq_hist_correct"], doc["mapq_hist_wrong"]
    rows = []
    for q in sorted({int(k) for k in hc} | {int(k) for k in hw}, reverse=True):
        n = hc.get(str(q), 0) + hw.get(str(q), 0)
        if n:
            rows.append(mk_qrow(q, n, hw.get(str(q), 0)))
    return rows


def check_against_strata(rows, strata_path: Path, lane: str):
    """The hist-derived cumulative counts must reproduce every strata row."""
    strata = {int(r["mapq_threshold"]): r
              for r in csv.DictReader(strata_path.open(), delimiter="\t")}
    for q, s in strata.items():
        n = sum(r.n for r in rows if r.mapq >= q)
        w = sum(r.wrong for r in rows if r.mapq >= q)
        if n != int(s["mapped_at_or_above"]) or w != int(s["incorrect_at_or_above"]):
            sys.exit(f"{lane}: hist/strata mismatch at q>={q}: "
                     f"{n}/{w} vs {s['mapped_at_or_above']}/{s['incorrect_at_or_above']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lanes", type=Path, required=True)
    ap.add_argument("--bench-repo", type=Path, required=True,
                    help="checkout providing bench.mapeval_collect")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--thin", type=int, default=None,
                    help="eval2roc -s override (harness default 100)")
    args = ap.parse_args()

    sys.path.insert(0, str(args.bench_repo))
    from bench.mapeval_collect import COLUMNS, DEFAULT_THIN, QRow, derive_points
    thin = args.thin if args.thin is not None else DEFAULT_THIN

    runs = {}  # (ref, tool, readset) -> (rows, total, criterion)
    for lane_dir in sorted(p for p in args.lanes.iterdir() if p.is_dir()):
        parsed = parse_lane(lane_dir.name)
        if parsed is None:
            continue
        ref, readset, tool = parsed
        doc = json.loads((lane_dir / "mapeval.json").read_text())
        rows = qrows_from_json(doc, QRow)
        check_against_strata(rows, lane_dir / "mapeval.strata.tsv", lane_dir.name)
        runs[(ref, tool, readset)] = (rows, doc["counts"]["truth_reads"],
                                      str(doc.get("criterion", "unknown")))

    lines = ["\t".join(COLUMNS)]
    for ref, readsets, panel in PANELS:
        for tool in TOOL_ORDER:
            for readset in readsets:
                if (ref, tool, readset) not in runs:
                    continue
                rows, total, criterion = runs[(ref, tool, readset)]
                label = (tool if readset not in READSET_LABEL
                         else f"{tool} ({READSET_LABEL[readset]})")
                for pt in derive_points(rows, total, thin):
                    lines.append("\t".join((
                        panel, label, criterion, str(pt.mapq),
                        str(pt.cum_aligned), str(pt.cum_wrong),
                        f"{pt.err_rate:.10g}", f"{pt.pct_aligned:.8g}",
                        str(total))))
    args.out.write_text("\n".join(lines) + "\n")
    n_runs = len({(l.split('\t')[0], l.split('\t')[1]) for l in lines[1:]})
    print(f"wrote {args.out}: {len(lines) - 1} points, {n_runs} curves")


if __name__ == "__main__":
    main()
