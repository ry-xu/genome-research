"""Pairwise identity over aligned columns (Biopython local alignment, gaps excluded) for all pairs of given FASTA files."""
import sys, itertools
from Bio import Align
rd = lambda p: "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
S = {p[2:-3]: rd(p) for p in sys.argv[1:]}
for a, b in itertools.combinations(S, 2):
    aln = al.align(S[a], S[b])[0]
    m = n = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            n += 1; m += S[a][s1 + k] == S[b][s2 + k]
    print(f"{a}\t{b}\taligned {n}\tidentity {m/n:.3f}")
