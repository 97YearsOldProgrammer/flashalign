# FlashAlign: the paper

Source of the FlashAlign preprint. `preprint.pdf` is the compiled paper; `make` rebuilds it
from `main.tex` and `sections/` with the figures in `figures/`.

**Figures.** Fig. 2 (mapping accuracy) is drawn by `scripts/make_fig2_mapq.py` from
`data/points.tsv`, which `scripts/fig2_points.py` builds from the `mapeval` output of each
lane (it imports the harness, `--bench-repo`). Fig. 3 (throughput and memory) is drawn by
`scripts/make_fig3_perf.py` from `data/perf.tsv`; `scripts/verify_fig3_perf.py` checks the
drawn bars against the table.

**Evaluation.** Every number comes from the evaluation harness used for the paper, which is
not part of this release, on the datasets named in the paper, with FlashAlign 0.4.0 and the
peer versions stated in the text. The simulation protocols are described in the paper.

**Citing.** Until the preprint is posted, cite the software: doi:10.5281/zenodo.22985251.
