# FlashAlign

[![GitHub Downloads](https://img.shields.io/github/downloads/97YearsOldProgrammer/flashalign/total.svg?style=social&logo=github&label=Download)](https://github.com/97YearsOldProgrammer/flashalign/releases)
[![PyPI](https://img.shields.io/pypi/v/flashalign.svg?style=flat)](https://pypi.org/project/flashalign/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/flashalign.svg?label=PyPI%20downloads)](https://pypi.org/project/flashalign/)
[![License](https://img.shields.io/github/license/97YearsOldProgrammer/flashalign.svg)](LICENSE)

## Getting started

```sh
git clone https://github.com/97YearsOldProgrammer/flashalign.git
cd flashalign
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target flashalign -j
# long reads against a reference genome
./build/flashalign align ref.fa reads.fq > aln.sam
# create an index first and then map
./build/flashalign index -x lr:hq ref.fa              # writes ref.fa.faix
./build/flashalign align ref.fa.faix hifi.fq.gz > aln.sam
# use presets
./build/flashalign align -x lr ref.fa ont.fq.gz > aln.sam                 # Oxford Nanopore genomic reads
./build/flashalign align -x lr:hq ref.fa hifi.fq.gz > aln.sam             # PacBio HiFi genomic reads
./build/flashalign align -x splice ref.fa cdna.fq.gz > aln.sam            # spliced long reads (strand unknown)
./build/flashalign align -x splice:hq -u f ref.fa isoseq.fq.gz > aln.sam  # PacBio Iso-Seq (transcript strand)
./build/flashalign align -x splice --junc-bed anno.bed ref.fa cdna.fq.gz > aln.sam  # use annotated junctions
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
preset and writes base-level alignments in the SAM format:

```sh
flashalign align ref.fa reads.fq > aln.sam
```

You can ask for PAF with the CIGAR in the `cg` tag:

```sh
flashalign align -f paf -c ref.fa reads.fq > aln.paf
```

Without `-c`, PAF output is placement only: no base-level alignment is done, so there is no
CIGAR and the coordinates differ from those of the base-level alignment:

```sh
flashalign align -f paf ref.fa reads.fq > approx-mapping.paf
```

Reads may be FASTA or FASTQ, plain or gzip'd, and several read files may follow the reference.

A FASTA reference is indexed in memory on every run. To index it once, write the index with
`flashalign index` and give the index in place of the reference:

```sh
flashalign index -x lr:hq ref.fa                  # indexing; writes ref.fa.faix
flashalign align ref.fa.faix reads.fq > aln.sam   # alignment
```

The index records the preset it was built with, and `align` maps under that preset unless `-x`
is given. The index owns the seeding: `-k`, `-s` and `-I` are indexing options and can't be
changed during mapping, and an explicit `-x` at mapping time keeps the index's own `-k` and
`-s`. Keep one index per seeding: `lr` (`-k21 -s9`), `lr:hq` (`-k21 -s5`), and one index shared
by `splice` and `splice:hq` (`-k15 -s10`).

### Use cases

FlashAlign is tuned to a read type with a preset, `-x`, which sets several options at the
same time. The default is `lr`.

#### Map long genomic reads

```sh
flashalign align -x lr ref.fa ont.fq.gz > aln.sam       # Oxford Nanopore reads
flashalign align -x lr:hq ref.fa hifi.fq.gz > aln.sam   # PacBio HiFi reads
```

`lr` is for noisy long reads of ~10% error rate. `lr:hq` is for accurate long reads with an
error rate below 1%; PacBio HiFi reads belong there, and there is no separate HiFi preset. The
two differ in seeding (`-s9` and `-s5`) and in scoring.

Seeds that occur too often in the reference are ignored. The cutoff follows the reference: it is
the larger of 200 and the occurrence of the top 1.81e-4 fraction of the most frequent seeds in
the index, at most 500 under `lr:hq`, and voting and chaining use the same number. On a human
genome it is 200; on a repeat-rich genome such as wheat it is about 2,200 under `lr` and 500
under `lr:hq`. `--max-vote-occ` sets it for both.

#### Map long mRNA/cDNA reads

```sh
flashalign align -x splice ref.fa cdna.fq.gz > aln.sam            # Nanopore cDNA or direct RNA
flashalign align -x splice:hq -u f ref.fa isoseq.fq.gz > aln.sam  # PacBio Iso-Seq
```

In the splice presets a long deletion is taken as an intron and written as the `N` CIGAR
operator; `-G` sets the maximum intron length (200k by default). By default (`-u b`) the
canonical GT-AG splice sites are looked for on both strands and the transcript strand is
written to the `ts` tag; `-u f` looks on the transcript strand only, for reads that are already
on that strand. `splice:hq` differs from `splice` only in scoring.

FlashAlign can take annotated junctions and prefer them during base alignment:

```sh
paftools.js gff2bed anno.gtf > anno.bed
flashalign align -x splice --junc-bed anno.bed ref.fa cdna.fq.gz > aln.sam
```

`--junc-bed` takes gene annotations in the 12-column BED format, which `paftools.js gff2bed`
converts from GTF or GFF3, or intron positions in 6-column BED with the strand column. A splice
donor or acceptor found in the annotation gets a score bonus, `--junc-bonus` (9 by default).

### Output

SAM is the default; `-f bam` writes BAM and `-f paf` writes PAF. Without `-f`, a `.sam`, `.bam`
or `.paf` extension on `-o` selects the format:

```sh
flashalign align -x lr:hq ref.fa hifi.fq.gz -o aln.bam
```

Unmapped reads are written to SAM and BAM unless `--sam-hit-only` is given; PAF never carries
them. Secondary alignments are written only with `--secondary yes`. `--cs`, `--MD` and `--eqx`
add the `cs` tag, the `MD` tag and `=`/`X` CIGAR operators. The PAF columns, the tags and the
`cs` operations are listed under OUTPUT FORMAT in the manual.

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
flashalign align ref.fa.faix reads.fq > aln.sam
```

Mapping quality is incorrect given a multi-part index.

### Getting help

The manual page, [`flashalign.1`](flashalign.1), describes every option, the output format and
the exit status (`man ./flashalign.1`). `flashalign align --help` and `flashalign index --help`
list the common options. Bugs and questions go to the
[issue page](https://github.com/97YearsOldProgrammer/flashalign/issues).

### Citing FlashAlign

The manuscript is not yet public. Until it is, cite the repository by URL and commit.

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

- Long reads only: a read has to carry enough anchors to be collapsed into one or more
  diagonals.
- FlashAlign requires SSE4.1 instructions on x86 CPUs or NEON on ARM CPUs. A build without
  them is not provided.

## License

FlashAlign is released under the MIT License; see [`LICENSE`](LICENSE).
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) covers the third-party code.
