"""
FlashAlign's compiled extension. Import `flashalign` instead: the package re-exports every name defined here.
"""

from collections.abc import Sequence
import os

# Added by hand: stubgen skips module dunders.
__version__: str


class Config:
    """
    Mapping options, named as `flashalign align --show-config` prints
    them; a dot there (tile_score.hit) is an underscore here.

    -1 (or None, or "") means "take the preset's value"; anything else
    is an explicit setting, which may be rejected with a ValueError.

        fa.Config(preset='lr:hq', threads=8, dp_min_score=50)

    Aligner.config returns the resolved configuration, with every field
    set to the value in use.
    """

    def __init__(self, preset: str = 'lr', **fields) -> None:
        """
        Config(preset='lr', **fields): a Config with those fields set.

        Every keyword must be one of the fields below; an unknown name is a
        ValueError.
        """

    @property
    def preset(self) -> str:
        """Preset name: 'lr', 'lr:hq', 'splice', 'splice:hq' [-x / --preset]."""

    @preset.setter
    def preset(self, arg: str, /) -> None: ...

    @property
    def k(self) -> int:
        """Seed k-mer length [-k]; mapping always uses the index's."""

    @k.setter
    def k(self, arg: int, /) -> None: ...

    @property
    def min_support(self) -> int:
        """Minimum anchor support [--min-support]. -1 = preset-owned."""

    @min_support.setter
    def min_support(self, arg: int, /) -> None: ...

    @property
    def max_query_seeds(self) -> int:
        """Query seeds kept per strand. -1 = preset-owned."""

    @max_query_seeds.setter
    def max_query_seeds(self, arg: int, /) -> None: ...

    @property
    def full_read_cigar(self) -> bool:
        """Realize base-level CIGARs; False maps only (the CLI's plain -f paf)."""

    @full_read_cigar.setter
    def full_read_cigar(self, arg: bool, /) -> None: ...

    @property
    def syncmer_s(self) -> int:
        """Closed-syncmer s [-s]; mapping always uses the index's."""

    @syncmer_s.setter
    def syncmer_s(self, arg: int, /) -> None: ...

    @property
    def syncmer_downsample(self) -> int:
        """Query seed downsampling stride; 1 = none."""

    @syncmer_downsample.setter
    def syncmer_downsample(self, arg: int, /) -> None: ...

    @property
    def vote_diag_bin_width(self) -> int:
        """Vote diagonal bin width [--dw]. -1 = preset-owned."""

    @vote_diag_bin_width.setter
    def vote_diag_bin_width(self, arg: int, /) -> None: ...

    @property
    def long_occ_cap(self) -> int:
        """Seed occurrence cap floor; the index may raise it. -1 = preset-owned."""

    @long_occ_cap.setter
    def long_occ_cap(self, arg: int, /) -> None: ...

    @property
    def long_primary_occ_cap(self) -> int:
        """
        Cap of the fixed occurrence policy that --max-vote-occ selects;
        unused under the presets' policy. -1 = preset-owned.
        """

    @long_primary_occ_cap.setter
    def long_primary_occ_cap(self, arg: int, /) -> None: ...

    @property
    def threads(self) -> int:
        """Worker threads [-t]; 0 = one per available CPU."""

    @threads.setter
    def threads(self, arg: int, /) -> None: ...

    @property
    def cigar_band_frac(self) -> float:
        """CIGAR band as a fraction of the block. -1 = preset-owned."""

    @cigar_band_frac.setter
    def cigar_band_frac(self, arg: float, /) -> None: ...

    @property
    def dp_match(self) -> int:
        """DP matching score [-A]. -1 = preset-owned."""

    @dp_match.setter
    def dp_match(self, arg: int, /) -> None: ...

    @property
    def dp_mismatch(self) -> int:
        """DP mismatch penalty [-B]. -1 = preset-owned."""

    @dp_mismatch.setter
    def dp_mismatch(self, arg: int, /) -> None: ...

    @property
    def dp_score_n(self) -> int:
        """DP ambiguous-base score [--score-N]. -1 = preset-owned."""

    @dp_score_n.setter
    def dp_score_n(self, arg: int, /) -> None: ...

    @property
    def dp_gap_open1(self) -> int:
        """First gap-open penalty [-O]. -1 = preset-owned."""

    @dp_gap_open1.setter
    def dp_gap_open1(self, arg: int, /) -> None: ...

    @property
    def dp_gap_extend1(self) -> int:
        """First gap-extension penalty [-E]. -1 = preset-owned."""

    @dp_gap_extend1.setter
    def dp_gap_extend1(self, arg: int, /) -> None: ...

    @property
    def dp_gap_open2(self) -> int:
        """Second gap-open penalty [-O]. -1 = preset-owned."""

    @dp_gap_open2.setter
    def dp_gap_open2(self, arg: int, /) -> None: ...

    @property
    def dp_gap_extend2(self) -> int:
        """Second gap-extension penalty [-E]. -1 = preset-owned."""

    @dp_gap_extend2.setter
    def dp_gap_extend2(self, arg: int, /) -> None: ...

    @property
    def dp_zdrop(self) -> int:
        """Z-drop score [-z]. -1 = preset-owned."""

    @dp_zdrop.setter
    def dp_zdrop(self, arg: int, /) -> None: ...

    @property
    def dp_end_bonus(self) -> int:
        """Alignment end bonus [--end-bonus]. -2 = preset-owned."""

    @dp_end_bonus.setter
    def dp_end_bonus(self, arg: int, /) -> None: ...

    @property
    def rna_min_intron(self) -> int:
        """Minimum intron length [--min-intron]; splice only."""

    @rna_min_intron.setter
    def rna_min_intron(self, arg: int, /) -> None: ...

    @property
    def rna_max_intron(self) -> int:
        """Maximum intron length [-G]; splice only."""

    @rna_max_intron.setter
    def rna_max_intron(self, arg: int, /) -> None: ...

    @property
    def rna_strand_mode(self) -> int:
        """
        Transcript strand [-u]: -1 unset, 0 auto, 1 forward, 2 reverse, 3 none.
        """

    @rna_strand_mode.setter
    def rna_strand_mode(self, arg: int, /) -> None: ...

    @property
    def tile_score_hit(self) -> int | None:
        """Query-tile supported reward [--tile-score]. None = preset-owned."""

    @tile_score_hit.setter
    def tile_score_hit(self, arg: int | None, /) -> None: ...

    @property
    def tile_score_block(self) -> int | None:
        """Query-tile block-open cost [--tile-score]. None = preset-owned."""

    @tile_score_block.setter
    def tile_score_block(self, arg: int | None, /) -> None: ...

    @property
    def tile_score_null(self) -> int | None:
        """Query-tile null cost [--tile-score]. None = preset-owned."""

    @tile_score_null.setter
    def tile_score_null(self, arg: int | None, /) -> None: ...

    @property
    def tile_score_miss(self) -> int | None:
        """
        Query-tile unsupported-ownership cost [--tile-score]. None = preset-owned.
        """

    @tile_score_miss.setter
    def tile_score_miss(self, arg: int | None, /) -> None: ...

    @property
    def rna_junction_bed(self) -> str:
        """Known-junction BED6/BED12 path [--junc-bed]; splice only."""

    @rna_junction_bed.setter
    def rna_junction_bed(self, arg: str, /) -> None: ...

    @property
    def rna_junction_bonus(self) -> int:
        """Known-junction endpoint bonus [--junc-bonus]. -1 = preset-owned."""

    @rna_junction_bonus.setter
    def rna_junction_bonus(self, arg: int, /) -> None: ...

    @property
    def rna_rival_pri_ratio(self) -> float:
        """Rival-to-primary chain score ratio [-p]. -1 = preset-owned."""

    @rna_rival_pri_ratio.setter
    def rna_rival_pri_ratio(self, arg: float, /) -> None: ...

    @property
    def rna_max_loci(self) -> int:
        """Candidate loci retained per read [-N]. -1 = preset-owned."""

    @rna_max_loci.setter
    def rna_max_loci(self, arg: int, /) -> None: ...

    @property
    def dp_min_score(self) -> int:
        """Minimum DP alignment score [-S]. -1 = preset-owned."""

    @dp_min_score.setter
    def dp_min_score(self, arg: int, /) -> None: ...

    @property
    def dp_bw(self) -> int:
        """
        Chaining and alignment bandwidth [-r]. -1 = preset-owned; a splice
        preset rejects it.
        """

    @dp_bw.setter
    def dp_bw(self, arg: int, /) -> None: ...

    @property
    def dp_bw_long(self) -> int:
        """
        Long-join bandwidth, the second -r value. -1 = preset-owned; a
        splice preset rejects it.
        """

    @dp_bw_long.setter
    def dp_bw_long(self, arg: int, /) -> None: ...

    @property
    def emit_md(self) -> bool:
        """Emit the MD:Z difference string [--MD]."""

    @emit_md.setter
    def emit_md(self, arg: bool, /) -> None: ...

    @property
    def emit_eqx(self) -> bool:
        """Write =/X CIGAR operators instead of M [--eqx]."""

    @emit_eqx.setter
    def emit_eqx(self, arg: bool, /) -> None: ...

    @property
    def cs(self) -> str:
        """
        cs:Z difference string form [--cs]: 'none', 'short' or 'long'.
        Accepts False/True for 'none'/'short'.
        """

    @cs.setter
    def cs(self, arg: object, /) -> None: ...

    def __repr__(self) -> str:
        """The preset, plus every field that differs from the struct default."""

class Hit:
    """
    One alignment record, with mappy's field names.

    Hits are read-only and come only from Aligner.map() and map_batch(),
    in the CLI's record order: the primary, then its supplementary
    segments. Coordinates are 0-based, half-open and on the forward
    strand of both sequences, as in PAF and mappy.
    """

    @property
    def ctg(self) -> str:
        """Target (reference) contig name. PAF column 6."""

    @property
    def ctg_len(self) -> int:
        """Length of `ctg`, in bases. PAF column 7."""

    @property
    def r_st(self) -> int:
        """Target start, 0-based. PAF column 8."""

    @property
    def r_en(self) -> int:
        """Target end, exclusive. PAF column 9."""

    @property
    def q_st(self) -> int:
        """Query start on the forward strand, 0-based. PAF column 3."""

    @property
    def q_en(self) -> int:
        """Query end on the forward strand, exclusive. PAF column 4."""

    @property
    def strand(self) -> int:
        """+1 if the read maps forward, -1 if reverse. PAF column 5."""

    @property
    def mapq(self) -> int:
        """Mapping quality, 0-60. PAF column 12."""

    @property
    def NM(self) -> int:
        """Edit distance over the aligned block (NM:i); -1 without a CIGAR."""

    @property
    def blen(self) -> int:
        """Aligned block length including gaps. PAF column 11."""

    @property
    def mlen(self) -> int:
        """Residue matches in the aligned block. PAF column 10."""

    @property
    def is_primary(self) -> bool:
        """True unless this record is a rival hypothesis (tp:A:S)."""

    @property
    def trans_strand(self) -> int:
        """Transcript strand: +1 / -1 / 0 for unknown (ts:A, RNA only)."""

    @property
    def cigar(self) -> list[tuple[int, int]]:
        """Aligned block as [(length, op), ...]; op indexes "MIDNSHP=XB"."""

    @property
    def cigar_str(self) -> str:
        """The same CIGAR as a string; the PAF cg:Z value (clips dropped)."""

    @property
    def cs(self) -> str:
        """
        cs:Z difference string; "" unless Config.cs or map(cs=...) asks
        for it.
        """

    @property
    def MD(self) -> str:
        """
        MD:Z difference string; "" unless Config.emit_md or map(MD=True)
        asks for it.
        """

    @property
    def score(self) -> int:
        """Alignment score (the AS:i tag), or the chain score without a CIGAR."""

    @property
    def read_len(self) -> int:
        """Length of the query this record was produced from."""

    @property
    def is_supplementary(self) -> bool:
        """True for a supplementary segment of the record before it."""

    @property
    def identity(self) -> float:
        """mlen / blen, or 0.0 without a base-level alignment."""

    @property
    def mismatches(self) -> int:
        """Substituted bases in the aligned block."""

    @property
    def insertions(self) -> int:
        """Inserted query bases in the aligned block."""

    @property
    def deletions(self) -> int:
        """Deleted reference bases in the aligned block."""

    @property
    def ambiguities(self) -> int:
        """Ambiguous (N) bases excluded from mlen and blen (the nn:i tag)."""

    @property
    def is_mapped(self) -> bool:
        """True for every Hit map() returns; unmapped reads yield no Hit."""

    @property
    def secondary(self) -> list[Hit]:
        """
        Rival placement hypotheses, on the primary only.

        These are what MAPQ was computed against, not records of their
        own: the CLI writes secondary rows only under --secondary yes, and
        map() never returns them. Each rival is flattened like the primary
        (the rival, then its supplementary segments). Empty on a
        supplementary Hit and on a rival.
        """

    def __str__(self) -> str:
        """
        mappy's line: the PAF row without its first two columns.

        q_st q_en strand ctg ctg_len r_st r_en mlen blen mapq tp:A ts:A
        cg:Z [cs:Z] [MD:Z], tab-separated. Aligner.paf() gives the full
        PAF row, query name and length included.
        """

    def __repr__(self) -> str:
        """A short one-line summary; str(hit) gives the PAF-style line."""

class Index:
    """
    A closed-syncmer seed index, with its reference embedded.

    An Index is reference-counted, so several Aligners can share one
    copy; Aligner.index returns the one an aligner uses.

        idx = fa.Index.load('ref.faix')       # a .faix written by the CLI
        idx = fa.Index.build_from_fasta('ref.fa', *fa.preset_seeding('lr'))

    k and syncmer_s are fixed when the index is built, and mapping always
    uses the index's pair, never a preset's.
    """

    @staticmethod
    def load(path: str | os.PathLike) -> Index:
        """Load a .faix seed index written by `flashalign index`."""

    @staticmethod
    def build(sequences: Sequence[tuple[str, str]], k: int = 21, syncmer_s: int = 9, threads: int = 0) -> Index:
        """
        Build an index from [(name, seq), ...], embedding the reference.

        Contigs are sorted by name, as the CLI sorts them, so the
        alignments match the CLI's; empty sequences are skipped. The index
        records no preset name.
        """

    @staticmethod
    def build_from_fasta(path: str | os.PathLike, k: int = 21, syncmer_s: int = 9, threads: int = 0) -> Index:
        """
        Build an index from a reference FASTA on disk.

        Plain, gzip or bgzf, or '-' for stdin, read with the CLI's reader.
        Pass *preset_seeding(name) for k and syncmer_s if the aligner will
        run under a preset other than 'lr'.
        """

    @staticmethod
    def is_index_file(path: str | os.PathLike) -> bool:
        """
        Is this path a .faix seed index?

        Decided from the file header, not the name, as `flashalign align`
        decides. False, never an exception, for '-', for a path that
        cannot be opened and for a file that is not a .faix.
        """

    def save(self, path: str | os.PathLike) -> None:
        """Write this index to `path` as a .faix."""

    def seq(self, name: str, start: int = 0, end: int = -1) -> str | None:
        """
        One slice of the embedded reference, or None for an unknown contig.

        Uppercase, with every non-ACGT base as 'N'. end=-1 means the
        contig end, and the range is clamped like mappy's Aligner.seq():
        start<0 becomes 0, end past the contig becomes its length, and
        start>=end gives ''.

        The first call unpacks the whole reference at one byte per base
        (about 3 GB for a human genome) and caches it on the index; later
        calls only copy.
        """

    @property
    def k(self) -> int:
        """Seed k-mer length this index was built with."""

    @property
    def syncmer_s(self) -> int:
        """Closed-syncmer s this index was built with."""

    @property
    def has_reference(self) -> bool:
        """True when the index embeds its reference (seq() needs it)."""

    @property
    def memory_mb(self) -> float:
        """Memory held by the index, in megabytes."""

    @property
    def total_bp(self) -> int:
        """Total reference bases indexed."""

    @property
    def preset(self) -> str:
        """
        The preset name the .faix header records, or ''.

        `flashalign index -x NAME` records it and
        `flashalign align ref.faix reads.fq` maps under it. An
        index built here records none, and Aligner then uses
        'lr' unless a preset is given.
        """

    @property
    def seq_names(self) -> list[str]:
        """Contig names, in index order (sorted by name)."""

    @property
    def seq_lengths(self) -> list[int]:
        """Contig lengths in bases, in seq_names order."""

    @property
    def n_seq(self) -> int:
        """Number of contigs in the index."""

    def __repr__(self) -> str:
        """The seeding, the contig count, the size and any recorded preset."""

class Aligner:
    """
    An aligner over one reference, in the mappy manner.

        a = fa.Aligner('ref.fa', preset='lr:hq')   # or 'ref.faix', or an Index
        for h in a.map(read_seq):
            print(h.ctg, h.r_st, h.r_en, h.strand, h.mapq, h.cigar_str)

    `target` is a path (str or os.PathLike) or an Index. A path is judged
    by its content, not its name: a .faix is loaded, anything else is
    indexed as a reference FASTA with the preset's (k, s).

    `preset` defaults to the preset the .faix records, or 'lr' when it
    records none, as `flashalign align ref.faix reads.fq` does. An
    explicit preset wins over the index's.

    `config` replaces the preset-derived configuration and is used as
    given; `threads` then only sizes an index build from a FASTA, since
    the Config has its own. A `config` and a `preset` that disagree are
    a ValueError.

    Construction raises FileNotFoundError for a missing path, ValueError
    for an unknown preset or a rejected configuration, and RuntimeError
    otherwise. It never returns a broken object, so mappy's
    `if not a: raise` idiom is unnecessary: __bool__ is always True.

    map(), map_batch() and paf() may be called concurrently on one
    Aligner from several Python threads, and release the GIL for the
    whole C++ call. reconfigure() may also run meanwhile; a map() that
    overlaps it uses either the old or the new configuration.
    """

    def __init__(self, target: str | os.PathLike | Index, preset: str | None = None, *, threads: int = 0, config: Config | None = None) -> None:
        """Aligner(target, preset=None, *, threads=0, config=None)"""

    def map(self, seq: str, *, cs: bool | str | None = None, MD: bool | None = None) -> list[Hit]:
        """
        Align one read; the CLI's record set, in the CLI's order.

        Returns [primary, *supplementary segments] as separate Hits, or
        [] for an unmapped read. Rival hypotheses stay on the primary's
        `.secondary`; the CLI writes them as rows only under
        --secondary yes.

        cs and MD request minimap2's difference strings for this call.
        None (the default) keeps Config.cs and Config.emit_md, both off
        unless set; True (or 'short' / 'long') and MD=True turn them on,
        False turns them off. A call whose request differs from the one
        in force reconfigures the aligner first, so a loop that
        alternates cs=True with plain calls reconfigures every time; to
        get cs or MD on every read, set them in the Config.

        The GIL is released for the whole C++ call.
        """

    def map_batch(self, seqs: Sequence[str], *, cs: bool | str | None = None, MD: bool | None = None) -> list[list[Hit]]:
        """
        Align many reads; one inner list per read, in input order.

        Each inner list is what map() would return for that read; cs and
        MD mean the same as for map() and apply to the whole batch. The
        batch runs on `threads` workers and releases the GIL once, so
        this is the fast path for many reads.
        """

    def paf(self, name: str, seq: str, *, cigar: bool = False) -> str:
        """
        The CLI's PAF record(s) for one read, as one string.

        Complete PAF lines, newline-terminated, with the query name and
        length; cigar=True adds the cg:Z tag (the CLI's -c). Under the
        default configuration, cigar=True gives what
        `flashalign align -f paf -c` writes for the read. The CLI's plain
        `-f paf` maps without base-level alignment;
        Config(full_read_cigar=False) with cigar=False gives that output.
        No secondary rows are written, as under the CLI's default. An
        unmapped read gives ''.

        The name seeds the tie-break between equally good loci, as a
        query name does in minimap2; map() and map_batch() have no name
        and seed it from the read length alone, as mappy's map() does.
        """

    def seq(self, name: str, start: int = 0, end: int = -1) -> str | None:
        """
        A slice of the reference, or None for a contig the index lacks.

        Index.seq() by another name; the same clamping, the same uppercase
        ACGT/N alphabet, and the same first-call cost (the whole embedded
        reference is unpacked once and cached).
        """

    def reconfigure(self, config: Config) -> None:
        """
        Install a new configuration on this aligner.

        k and s stay the index's. Safe to call while other threads map;
        a map() that overlaps it uses either the old or the new
        configuration. ValueError for a rejected configuration.
        """

    @property
    def seq_names(self) -> list[str]:
        """Reference contig names, in index order."""

    @property
    def seq_lengths(self) -> list[int]:
        """Reference contig lengths, in seq_names order."""

    @property
    def n_seq(self) -> int:
        """Number of reference contigs."""

    @property
    def k(self) -> int:
        """Seed k-mer length, which is always the index's."""

    @property
    def s(self) -> int:
        """Closed-syncmer s, which is always the index's."""

    @property
    def preset(self) -> str:
        """The preset this aligner resolved its configuration from."""

    @property
    def config(self) -> Config:
        """
        A copy of the resolved configuration from the constructor or
        reconfigure(); a per-call cs/MD request never changes it.
        """

    @property
    def index(self) -> Index:
        """
        The Index this aligner maps against.

        Reference-counted: a second Aligner built on it shares the same
        copy.
        """

    def __bool__(self) -> bool:
        """
        Always True: a failed construction raises instead of returning a
        falsy object.
        """

    def __repr__(self) -> str:
        """
        The preset, the seeding in force, the contig count and the thread count.
        """

class FastxReader:
    """
    An iterator over the records of one FASTA/FASTQ file.

    Made by fastx_read(); yields (name, seq, qual) tuples, or
    (name, seq, qual, comment) with read_comment=True. Reads plain, gzip
    and bgzf files (multi-member gzip included) with the CLI's reader,
    and '-' for stdin; not unaligned BAM. Sequences are uppercased.
    """

    def __iter__(self) -> FastxReader:
        """The reader is its own iterator."""

    def __next__(self) -> tuple:
        """
        The next record, or StopIteration at end of file.

        The GIL is released around each read. `qual` is None for FASTA,
        and `comment` is None when the header carries none.
        """

def fastx_read(path: str | os.PathLike, read_comment: bool = False) -> FastxReader:
    """
    Iterate one FASTA/FASTQ file: (name, seq, qual) per record.

        for name, seq, qual in fa.fastx_read('reads.fq.gz'):
            ...

    `qual` is None for FASTA. read_comment=True yields four-tuples
    (name, seq, qual, comment), with comment None when the header has
    none. Plain, gzip and bgzf are all read (multi-member gzip included),
    and path='-' reads stdin; unaligned BAM is not read here. Sequences
    are uppercased.

    Unlike mappy, a path that cannot be opened raises FileNotFoundError
    rather than returning None.
    """

def revcomp(seq: str) -> str:
    """
    Reverse complement, as mappy's revcomp().

    IUPAC codes are complemented (R<->Y, K<->M, B<->V, D<->H; S, W
    and N map to themselves), case is kept, U becomes A, and any
    other character is kept as is, in reversed position.
    """

def preset_seeding(preset: str) -> tuple[int, int]:
    """
    The (k, syncmer_s) a preset builds its index with.

        fa.Index.build_from_fasta('ref.fa', *fa.preset_seeding('lr:hq'))

    Mapping always uses the index's own seeding, so an index built with
    Index.build_from_fasta for a preset other than 'lr' needs this pair.
    ValueError for an unknown preset.
    """
