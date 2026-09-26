"""FlashAlign: a long-read DNA and RNA aligner with a mappy-style API.

This package drives the same engine as the ``flashalign`` command through its
public C++ API, and ``Aligner.paf()`` writes the same PAF lines as the command
for the same read and configuration.

    import flashalign as fa
    a = fa.Aligner("ref.fa", preset="lr:hq")            # or "ref.faix", or an fa.Index
    for name, seq, qual in fa.fastx_read("reads.fq.gz"):
        for h in a.map(seq):
            print(name, h.ctg, h.r_st, h.r_en, h.strand, h.mapq, h.cigar_str)

``Aligner`` tells a ``.faix`` index from a reference FASTA by its content, not
its name, as the CLI does; a FASTA is indexed with the preset's seeding, and a
``.faix`` is mapped under the preset its header records unless ``preset=`` is
given.

``map()`` returns the primary and its supplementary segments for one read, each
a ``Hit``, or ``[]`` for an unmapped read. ``map_batch()`` maps many reads at
once on ``threads`` workers, with one GIL release for the whole batch.

Errors are raised, never returned as falsy objects: ``FileNotFoundError`` for a
missing path, ``ValueError`` for an unknown preset or a rejected configuration,
and ``RuntimeError`` otherwise.
"""

from ._flashalign import (
    Aligner,
    Config,
    FastxReader,
    Hit,
    Index,
    __version__,
    fastx_read,
    preset_seeding,
    revcomp,
)

__all__ = [
    "Aligner",
    "Config",
    "FastxReader",
    "Hit",
    "Index",
    "__version__",
    "fastx_read",
    "preset_seeding",
    "revcomp",
]
