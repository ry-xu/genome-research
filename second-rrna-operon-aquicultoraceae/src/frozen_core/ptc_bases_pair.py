"""Read PTC bases of a divergent 23S by structure-free pairwise alignment (Biopython global, same scoring as
kill_div/ptc_pairwise.py): directly against E. coli 23S, and composed E. coli -> same-genome normal copy -> divergent.

usage: uv run -q --no-project --with biopython python -I src/ptc_bases_pair.py ECOLI.fa NORM.fa DIV.fa
Prints position, E. coli base, normal base, divergent base (direct), divergent base (composed), +-4 nt windows.
"""
import sys
from Bio import Align

POS = [2062, 2251, 2252, 2253, 2449, 2451, 2452, 2453, 2503, 2504, 2583, 2584]


def rd(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper().replace("U", "T")


def mapping(al, a, b):
    m = {}
    aln = al.align(a, b)[0]
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            m[s1 + k] = s2 + k
    return m


def main(ep, npth, dp):
    E, N, D = rd(ep), rd(npth), rd(dp)
    al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
    al.end_gap_score = 0
    en, ed, nd = mapping(al, E, N), mapping(al, E, D), mapping(al, N, D)
    for p in POS:
        q = p - 1
        nb = N[en[q]] if q in en else "-"
        db = D[ed[q]] if q in ed else "-"
        cq = nd.get(en.get(q, -1))
        cb = D[cq] if cq is not None else "-"
        w = D[ed[q] - 4:ed[q] + 5] if q in ed else ""
        print(p, E[q], nb, db, cb, w, sep="\t")


if __name__ == "__main__":
    main(*sys.argv[1:4])
