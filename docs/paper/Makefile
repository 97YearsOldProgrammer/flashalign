# Build the manuscript with tectonic (fetches packages on demand); falls back
# to latexmk if tectonic is absent. `make` also redraws Fig. 2 and Fig. 3 when
# their tables or scripts change.
PDF = main.pdf
TEX = main.tex $(wildcard sections/*.tex) refs.bib

all: $(PDF)

$(PDF): $(TEX) figures/fig2_mapq.pdf figures/fig3_perf.pdf
	@if command -v tectonic >/dev/null; then tectonic --keep-logs main.tex; \
	else latexmk -pdf -interaction=nonstopmode main.tex; fi

figures/fig2_mapq.pdf: data/points.tsv scripts/make_fig2_mapq.py scripts/figstyle.py
	python3 scripts/make_fig2_mapq.py --points data/points.tsv --out $@

figures/fig3_perf.pdf: data/perf.tsv scripts/make_fig3_perf.py scripts/figstyle.py
	python3 scripts/make_fig3_perf.py --perf data/perf.tsv --out $@

clean:
	rm -f main.pdf main.aux main.bbl main.blg main.log main.out main.toc

.PHONY: all clean
