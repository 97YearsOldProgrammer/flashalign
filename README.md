# FlashAlign

Handle high-throughput long-read DNA and RNA data with flashalign.

[![GitHub Downloads](https://img.shields.io/github/downloads/97YearsOldProgrammer/flashalign/total.svg?style=social&logo=github&label=Download)](https://github.com/97YearsOldProgrammer/flashalign/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/flashalign.svg?label=bioconda)](https://anaconda.org/bioconda/flashalign)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/flashalign.svg?label=bioconda%20downloads)](https://anaconda.org/bioconda/flashalign)
[![PyPI](https://img.shields.io/pypi/v/flashalign.svg?style=flat)](https://pypi.org/project/flashalign/)
[![License](https://img.shields.io/github/license/97YearsOldProgrammer/flashalign.svg)](LICENSE)

## Getting started

```sh
git clone https://github.com/97YearsOldProgrammer/flashalign.git
cd flashalign
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target flashalign -j
# long reads against a reference genome
./build/flashalign align -a ref.fa reads.fq > aln.sam
# create an index first and then map
./build/flashalign index -x lr:hq ref.fa              # writes ref.fa.faix
./build/flashalign align -a ref.fa.faix hifi.fq.gz > aln.sam
# use presets
./build/flashalign align -ax lr ref.fa ont.fq.gz > aln.sam                 # Oxford Nanopore genomic reads
./build/flashalign align -ax lr:hq ref.fa hifi.fq.gz > aln.sam             # PacBio HiFi genomic reads
./build/flashalign align -ax splice ref.fa cdna.fq.gz > aln.sam            # spliced long reads (strand unknown)
./build/flashalign align -ax splice:hq -u f ref.fa isoseq.fq.gz > aln.sam  # PacBio Iso-Seq (transcript strand)
./build/flashalign align -ax splice --junc-bed anno.bed ref.fa cdna.fq.gz > aln.sam  # use annotated junctions
./build/flashalign align -cx asm5 ref.fa asm.fa > aln.paf                  # assembly vs reference (experimental)
./build/flashalign align -x ava-ont reads.fq reads.fq > ovl.paf            # all-vs-all read overlap (experimental)
# man page for detailed command line options
man ./flashalign.1
```

## Users' guide

### Installation

For Linux x86-64, a precompiled binary is on the
[release page](https://github.com/97YearsOldProgrammer/flashalign/releases). It runs on any
x86-64 Linux with glibc 2.17 or newer and needs no installation:

```sh
curl -L https://github.com/97YearsOldProgrammer/flashalign/releases/download/v0.1.0/flashalign-0.1.0_x64-linux.tar.bz2 | tar -jxvf -
./flashalign-0.1.0_x64-linux/flashalign version
```

FlashAlign is on [Bioconda](https://anaconda.org/bioconda/flashalign), for Linux and macOS on
x86-64 and ARM:

```sh
conda install -c conda-forge -c bioconda flashalign
```

The Python package (see the [Developers' guide](#developers-guide)) installs with pip. Prebuilt
wheels cover Linux on x86-64 and ARM and macOS on Apple silicon, for Python 3.12 or newer:

```sh
pip install flashalign
```

Otherwise, FlashAlign is built from source with CMake 3.18 or newer, a C and C++17 compiler and
zlib 1.2.9 or newer:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target flashalign -j
./build/flashalign version
```

The result is one binary, `build/flashalign`. The same commands build on x86-64, where the
alignment kernels use SSE4.1, and on ARM, where they go through the bundled sse2neon shim.

### General usage

Without options, `flashalign align` takes a reference and a read file, maps under the `lr`
preset and writes approximate mappings in the PAF format. No base-level alignment is done, so
there is no CIGAR and the coordinates are approximate:

```sh
flashalign align ref.fa reads.fq > approx-mapping.paf
```

You can ask for the CIGAR in the `cg` tag of PAF:

```sh
flashalign align -c ref.fa reads.fq > aln.paf
```

or for base-level alignments in the SAM format:

```sh
flashalign align -a ref.fa reads.fq > aln.sam
```

Reads may be FASTA or FASTQ, plain or gzip'd, and several read files may follow the reference.

A FASTA reference is indexed in memory on every run. To index it once, write the index with
`flashalign index` and give the index in place of the reference:

```sh
flashalign index -x lr:hq ref.fa                     # indexing; writes ref.fa.faix
flashalign align -a ref.fa.faix reads.fq > aln.sam   # alignment
```

The index records the preset it was built with, and `align` maps under that preset unless `-x`
is given. The index owns the seeding: `-k`, `-s` and `-I` are indexing options and can't be
changed during mapping, and an explicit `-x` at mapping time keeps the index's own `-k` and
`-s`.

### Use cases

FlashAlign is tuned to a read type with a preset, `-x`, which sets several options at the
same time. The default is `lr`.

#### Map long genomic reads

```sh
flashalign align -ax lr ref.fa ont.fq.gz > aln.sam       # Oxford Nanopore reads
flashalign align -ax lr:hq ref.fa hifi.fq.gz > aln.sam   # PacBio HiFi reads
```

`lr` is for noisy long reads of ~10% error rate. `lr:hq` is for accurate long reads with an
error rate below 1%, such as PacBio HiFi reads. The two differ in seeding (`-s9` and `-s5`)
and in scoring.

#### Map long mRNA/cDNA reads

```sh
flashalign align -ax splice ref.fa cdna.fq.gz > aln.sam            # Nanopore cDNA
flashalign index -x splice -k 14 ref.fa ref.k14.faix               # Nanopore direct RNA: a k14 index
flashalign align -ax splice -u f ref.k14.faix drna.fq.gz > aln.sam
flashalign align -ax splice:hq -u f ref.fa isoseq.fq.gz > aln.sam  # PacBio Iso-Seq
```

`-u b`, the default, looks for splice sites on both strands; `-u f` looks on the transcript
strand only, for reads already on that strand, such as direct RNA and Iso-Seq. For noisy
Nanopore direct RNA, an index built with `-k 14` finds more junctions, short first and last
exons above all. `splice:hq` differs from `splice` in scoring, in the vote's diagonal bin (48
against 64), in the fine chain's query gap (10,000 against 20,000; `-g` sets it) and in the
MAPQ coverage threshold (0.70 against 0.55).

FlashAlign can take annotated junctions and prefer them during base alignment:

```sh
paftools.js gff2bed anno.gtf > anno.bed
flashalign align -ax splice --junc-bed anno.bed ref.fa cdna.fq.gz > aln.sam
```

`--junc-bed` takes gene annotations in the 12-column BED format, which `paftools.js gff2bed`
converts from GTF or GFF3, or intron positions in 6-column BED with the strand column. A splice
donor or acceptor found in the annotation gets a score bonus, `--junc-bonus` (9 by default).

#### Align an assembly to a reference (experimental)

```sh
flashalign align -cx asm5 ref.fa asm.fa > aln.paf    # an assembly ~0.1% from the reference
```

`asm5`, `asm10` and `asm20` are experimental: they are not qualified on intact chromosomes. They
are for divergences of about 0.1%, 1% and several percent, and take minimap2's scoring for its
presets of the same names, one row for the gap fills and the read ends, with `-r1000,100000`,
`-g10000` and `-S200`. Their placement samples 4,096 vote seeds per strand into a partition of
2,048 query tiles, and, as on every DNA preset, the chains' anchors decide which chain owns each
stretch of the query. They map with an index built with `lr`'s seeding (`-k21 -s9`) and refuse
another.

#### Find overlaps between long reads (experimental)

```sh
flashalign align -x ava-ont reads.fq reads.fq > ovl.paf     # Nanopore read overlap
flashalign align -x ava-hifi reads.fq reads.fq > ovl.paf    # PacBio HiFi read overlap
```

`ava-ont` and `ava-hifi` are experimental: they were measured on simulated reads only. They map a
read set against itself and print every overlap as a map-only PAF record in minimap2's manner:
mapping quality 0 and the tags `tp:A:S`, `cm:i` and `s1:i`, with columns 10 and 11 approximate.
Each read is found in the target by its name and left out of its own mapping, so read names must
be unique. An overlap is printed once, from the shorter of its two reads (minimap2 prints it from the
read whose name sorts first); `--dual=yes` prints it from each read that finds it.
Base-level alignment is not available: `-a`, `-c`, `--cs`, `--MD` and the options that only shape
an alignment are refused, and so are both presets in the Python module.

The seeding is `-k17 -s9` under `ava-ont` and `-k21 -s5` under `ava-hifi`; `flashalign index -x
ava-ont -k INT -s INT reads.fq reads.faix` builds another, and `flashalign align reads.faix
reads.fq` maps on it. The whole read set is one index in memory unless it is built in parts with
`flashalign index -I NUM` ([Multi-part index](#multi-part-index)). On simulated Nanopore reads
`ava-ont` finds 99.8 to 99.9% of the overlaps of 2 kb or more where minimap2 `ava-ont` finds 99.96 to
99.99%, printing a third to two thirds of its false pairs.

### Output

PAF is the default and `-a` writes SAM. `-o FILE` outputs alignments to `FILE` [stdout],
whatever its name. For sorted BAM, pipe SAM to samtools:

```sh
flashalign align -a -x lr:hq ref.fa hifi.fq.gz | samtools sort -o aln.bam
```

Unmapped reads are written to SAM unless `--sam-hit-only` is given, and to PAF only with
`--paf-no-hit`. Secondary alignments are written only with `--secondary=yes`. `--cs`, `--MD` and
`--eqx` add the `cs` tag, the `MD` tag and `=`/`X` CIGAR operators. The PAF columns, the tags and
the `cs` operations are listed under OUTPUT FORMAT in the manual.

### Advanced features

#### The resolved configuration

`--show-config` prints every resolved option with the source of its value, then exits without
loading the reference or the reads; given an index, it reads only the index header:

```sh
flashalign align -x lr:hq --show-config
flashalign align --show-config ref.fa.faix
```

#### Multi-part index

`-I` caps the reference bases loaded into memory for indexing. A reference longer than that is
indexed in parts of whole contigs, and `align` maps the reads against one part at a time,
reading them once per part:

```sh
flashalign index -I 4G ref.fa
flashalign align -a ref.fa.faix reads.fq > aln.sam
```

Mapping quality is incorrect given a multi-part index.

### Getting help

The manual page, [`flashalign.1`](flashalign.1), describes every option, the output format and
the exit status (`man ./flashalign.1`). `flashalign align --help` and `flashalign index --help`
list the common options. Bugs and questions go to the
[issue page](https://github.com/97YearsOldProgrammer/flashalign/issues).

### Citing FlashAlign

The manuscript is not yet public. Until then, please cite the software:
[doi:10.5281/zenodo.22985251](https://doi.org/10.5281/zenodo.22985251).

## Developers' guide

FlashAlign is also a C++ library, with its public headers in
[`include/flashalign/`](include/flashalign/), and a Python package, `flashalign`, that drives
the same library:

```sh
python -m pip install .
python demo.py ref.fa reads.fq
```

The package needs CPython 3.12 or newer and installs no `flashalign` command.
[`demo.py`](demo.py) shows typical use; `help(flashalign.Aligner)` and
[`python/flashalign/_flashalign.pyi`](python/flashalign/_flashalign.pyi) document the API.

## Limitations

- Built for long reads: a read has to carry enough anchors to be collapsed into one or more
  diagonals. Short reads with enough anchors also work, but without a speed advantage over
  minimap2.
- FlashAlign requires SSE4.1 instructions on x86 CPUs or NEON on ARM CPUs. A build without
  them is not provided.
- The assembly presets `asm5`, `asm10` and `asm20` are experimental: they are not qualified on
  intact chromosomes.
- The overlap presets `ava-ont` and `ava-hifi` are experimental: they were measured on simulated
  reads only.

## License

FlashAlign is released under the MIT License; see [`LICENSE`](LICENSE).
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) covers the third-party code.
