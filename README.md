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
./build/flashalign index -x splice -k14 ref.fa ref.k14.faix                # noisy Nanopore direct RNA-seq: a k14 index,
./build/flashalign align -ax splice -uf ref.k14.faix drna.fq.gz > aln.sam  #   then map on it
./build/flashalign align -cx asm5 ref.fa asm.fa > aln.paf                  # assembly to assembly/ref alignment (experimental)
./build/flashalign align -x ava-ont reads.fq reads.fq > ovlp.paf           # Nanopore read overlap (experimental)
./build/flashalign align -x ava-hifi reads.fq reads.fq > ovlp.paf          # PacBio HiFi read overlap (experimental)
# man page for detailed command line options
man ./flashalign.1
```

## Users' guide

FlashAlign is a long-read aligner that maps DNA and mRNA sequences against a reference genome.
Typical use cases include: (1) mapping Oxford Nanopore, PacBio HiFi and ultra-long genomic
reads to a reference genome; (2) splice-aware alignment of PacBio Iso-Seq, Nanopore cDNA or
Nanopore direct RNA reads against a reference genome; (3) assembly-to-reference alignment
(experimental); (4) finding overlaps between long reads (experimental).

Presets marked experimental work, but are not yet in their best-tuned form. For long genomic
reads, FlashAlign is nearly twice as fast as the mainstream long-read mappers minimap2 and
minibwa. It is more accurate on simulated long reads and produces biologically meaningful
alignments ready for downstream analyses. For long spliced reads, FlashAlign is more accurate
than minimap2, mapping more reads with fewer errors, at the same or a slightly higher speed.
Detailed evaluations are available from the [preprint](docs/paper/preprint.pdf).

### Installation

For Linux x86-64, a precompiled binary is on the
[release page](https://github.com/97YearsOldProgrammer/flashalign/releases). It runs on any
x86-64 Linux with glibc 2.17 or newer and needs no installation:

```sh
curl -L https://github.com/97YearsOldProgrammer/flashalign/releases/download/v0.4.0/flashalign-0.4.0_x64-linux.tar.bz2 | tar -jxvf -
./flashalign-0.4.0_x64-linux/flashalign version
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

Reads may be FASTA or FASTQ, plain, gzip'd or bgzip'd, or unaligned BAM, and several read files may follow the
reference. Large gzip, bgzip and BAM inputs decode on threads borrowed from `-t`.

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
against 64), in the vote spacing (one seed per 16 query bases against 8, at least 128 per
strand; `--vote-spacing` sets it) and in the MAPQ coverage threshold (0.70 against 0.55); both
take minimap2's query gap, 2000 (`-g`).

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
flashalign align -cx asm5 ref.fa asm.fa > aln.paf    # assembly to assembly/ref alignment
```

`asm5`, `asm10` and `asm20` are for an average divergence of about 0.1%, 1% and several percent,
with minimap2's scoring for its presets of the same names.

#### Find overlaps between long reads (experimental)

```sh
flashalign align -x ava-ont  reads.fq reads.fq > ovlp.paf    # Oxford Nanopore read overlap
flashalign align -x ava-hifi reads.fq reads.fq > ovlp.paf    # PacBio HiFi read overlap
```

The output is PAF without CIGAR, each overlap once, from the shorter read (`--dual=yes` prints it
from both reads); `-a` and `-c` are refused. Read names must be unique.
`ava-ont` seeds with `-k17 -s9` and `ava-hifi` with `-k21 -s5`; a read set mapped more than once is
indexed once with `flashalign index -x ava-ont reads.fq reads.faix`, in parts with `-I` for a large
set ([Multi-part index](#multi-part-index)). The Python module has no overlap presets.

### Output

PAF is the default and `-a` writes SAM. `-o FILE` outputs alignments to `FILE` [stdout],
whatever its name. For sorted BAM, pipe SAM to samtools:

```sh
flashalign align -a -x lr:hq ref.fa hifi.fq.gz | samtools sort -o aln.bam
```

Unmapped reads are written to SAM unless `--sam-hit-only` is given, and to PAF only with
`--paf-no-hit`. Secondary alignments are written only with `--secondary=yes`. `--cs`, `--MD` and
`--eqx` add the `cs` tag, the `MD` tag and `=`/`X` CIGAR operators. The PAF columns, the tags and
the `cs` operations are listed under [OUTPUT FORMAT](https://97yearsoldprogrammer.github.io/flashalign/#OUTPUT_FORMAT) in the manual.

### Advanced features

#### The resolved configuration

`--show-config` prints every resolved option with the source of its value, then exits without
loading the reads or the reference bases; given a single-part index under a DNA or splice preset,
it also loads the index, without its bases, to print the resolved seed occurrence cap:

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

#### The cs optional tag

The `cs` SAM/PAF tag encodes bases at mismatches and INDELs. Like CIGAR, `cs` consists of a
series of operations. Each leading character specifies the operation; the following sequence is
the one involved in the operation.

The `cs` tag is enabled by the command line option `--cs`, which implies `-c`. The following
alignment, for example:

```txt
CGATCGATAAATAGAGTAG---GAATAGCA
||||||   ||||||||||   |||| |||
CGATCG---AATAGAGTAGGTCGAATtGCA
```

is represented as `:6-ata:10+gtc:4*at:3`, where `:[0-9]+` represents an identical block, `-ata`
represents a deletion, `+gtc` an insertion and `*at` indicates reference base `a` is substituted
with a query base `t`. It is similar to the `MD` SAM tag but is standalone and easier to parse.

If `--cs=long` is used, the `cs` string also contains identical sequences in the alignment. The
above example will become `=CGATCG-ata=AATAGAGTAG+gtc=GAAT*at=GCA`. The long form of `cs`
encodes both reference and query sequences in one string. Under a splice preset the `cs` tag
also encodes intron positions and splicing signals with the `~` operation (see
[OUTPUT FORMAT](https://97yearsoldprogrammer.github.io/flashalign/#OUTPUT_FORMAT) in the manual).

#### Working with the PAF format

FlashAlign writes minimap2's PAF: the same 12 columns and the same tag names. The
[`paftools.js`](https://github.com/lh3/minimap2/tree/master/misc) script that comes with
minimap2 therefore reads FlashAlign's PAF and SAM as it reads minimap2's; the commands that
need a CIGAR need `-c` (or `-a`) output.

### Algorithm overview

In the following, FlashAlign command line options have a dash ahead and are highlighted in bold,
with their defaults under the `lr` preset in brackets. The description may help to tune
FlashAlign parameters.

Minimap2 chains a read's anchors colinearly, treating every anchor alike, and other aligners
weight their seeds. FlashAlign instead samples the read's seeds and compresses them into votes for
approximate diagonals; the vote peaks are the candidate placements. Each candidate is then
decompressed in one pass of exact chaining over all of its anchors, organised by a tree, where
minimap2 chains twice, its second pass with a range-maximum query (RMQ).

1. Read **-I** [=*no split*] reference bases, extract (**-k**,**-s**)-closed syncmers
   [=*21*,*9*] and index them. A closed syncmer is a k-mer whose minimal s-mer sits at one of its
   two ends, so about 2/(k-s+1) of the k-mers are seeds.

2. Read **-K** [=*500M*] query bases per batch; the mapping threads draw from two batches at once.
   For each query sequence, do step 3 through 9. Output bytes and order do not depend on **-K**.

3. Look up the syncmers of the query in the index to find the seeds, and ignore, when voting and
   chaining, a seed whose reference occurrence is above **--max-vote-occ** [=*200*; unset, the
   larger of 200 and the occurrence of the top 1.81e-4 fraction of the most frequent seeds]. As
   minimap2 does, drop a seed whose k-mer occurs in the read itself more often than that cap and
   in more than **--q-occ-frac** [=*0.01*] of the read's seeds. A vote tile whose seeds are all
   above the cap still votes with its rarest seed of at most 4095 occurrences.

4. Sample the seeds and vote. Under a DNA preset the read is cut into **-v** [=*128*] equal tiles
   per strand, so the voters are scattered over the read, and each tile votes with one of its
   seeds for a diagonal bin **-W** [=*64*] reference bases wide, widened to L/128 for a read of
   length L, at most **--dw-max** [=*2048*]; a bin with at least **--min-support** [=*3*] votes is
   a vote peak. The splice presets vote with equally spaced seeds instead, one per
   **--vote-spacing** [=*8*] query bases, at least 128 per strand.

5. Admit the vote peaks to the catalogue as candidates: keep every peak whose vote is at least
   **--vote-ratio** [=*0.25*] of the best peak's, at most **--max-cands** [=*16*] per strand. A
   candidate's seed support over the read is counted on the second value of **-v** [=*128*]
   equal query tiles.

6. Under a DNA preset, chain the candidates of each contig and strand once, over the whole read:
   every anchor of the candidates' windows enters one exact chaining pass, the anchors collapsed
   into gap-free diagonal runs whose best predecessors are found exactly, through a k-d tree over
   the runs when the pool is large. The chain links anchors within the band of the second value of **-r**
   [=*500,20000*], across at most that band or **-g** [=*5000*] where that is larger, and a chain
   scoring below **-m** [=*40*] is discarded. Under a splice preset the fine chain takes introns
   up to **-G** [=*200k*] and query gaps up to **-g** [=*2000*].

7. Choose the owner chains as minimap2's `mm_set_parent` does: taking the chains best score
   first, attach each to the first owner whose query span overlaps its own by more than half of
   the shorter, or else make it a new owner; the owners become the primary and supplementary
   records. Keep as secondary alignments the chains attached to an owner that score at least
   **-p** [=*0.8*] times the owner's score, or within 2k of it, at most **-N** [=*5*] per read.

8. Without **-a** or **-c**, write each record from its owner chain's anchors: there is no
   base-level alignment, so the coordinates are approximate. With **-a** or **-c**, align base by
   base: close the gaps between anchors within the second value of **-r**, extend toward the read
   ends within its first value [=*500*] and stop where there are no seeds within **-g**
   [=*5000*]. Truncate where the running score drops too quickly, the Z-drop **-z** [=*800,200*]
   in the gap fills and *400,200* in the read-end extensions, which use `lr`'s own scoring row; a
   drop above the second value is tested for a small inversion, and an alignment needs a peak DP
   score of at least **-S** [=*80*].

9. Estimate the mapping quality against the best rival chain (the `s2` tag), a rival scoring under
   **-m** being weighed as scoring **-m**, and write PAF or SAM. Secondary alignments are written
   only with **--secondary**=*yes* [=*no*], at mapping quality 0; **-y** copies the FASTA/Q
   comments of each read to its records.

### Getting help

The [manual page](https://97yearsoldprogrammer.github.io/flashalign/) describes every option, the
output format and the exit status; it is the HTML form of [`flashalign.1`](flashalign.1)
(`man ./flashalign.1`). `flashalign align --help` and `flashalign index --help`
list the common options. Bugs and questions go to the
[issue page](https://github.com/97YearsOldProgrammer/flashalign/issues).

### Citing FlashAlign

The manuscript is not yet on a preprint server. Until then, please cite the software:
[doi:10.5281/zenodo.22985251](https://doi.org/10.5281/zenodo.22985251).
The preprint is in [`docs/paper/preprint.pdf`](docs/paper/preprint.pdf); the bioRxiv DOI will replace this link.

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
- The assembly and overlap presets are experimental.
- On HiFi reads FlashAlign splits inversion breakpoints less often than minimap2 and calls fewer
  inversions; the Z-drop split is being reworked.
- BAM stores at most 65535 CIGAR operations per record, which an ultra-long read can exceed.
  FlashAlign has no `-L`: it writes the full CIGAR in SAM, and samtools moves a longer CIGAR to
  the `CG:B,I` tag when it converts the SAM to BAM.

## License

FlashAlign is released under the MIT License; see [`LICENSE`](LICENSE).
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) covers the third-party code.
