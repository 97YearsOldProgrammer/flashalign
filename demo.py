import sys
import flashalign as fa

ref, reads = sys.argv[1], sys.argv[2]
a = fa.Aligner(ref, preset="lr:hq")
for name, seq, qual in fa.fastx_read(reads):
    for h in a.map(seq):
        print(name, h.ctg, h.r_st, h.r_en, "+-"[h.strand < 0], h.mapq, h.cigar_str, sep="\t")
