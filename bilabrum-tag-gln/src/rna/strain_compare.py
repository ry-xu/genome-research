"""Pairwise read-codon agreement between 1KP runs at panel codons called in both (are they one strain?).
usage: strain_compare.py A.codons.tsv.gz B.codons.tsv.gz ..."""
import gzip
import sys
from itertools import combinations


def load(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["call"] not in ("mixed", ""):
                d[(r["contig"], r["left0"], r["strand"])] = (r["genome"], r["call"])
    return d


runs = {p.split("/")[-1].split(".")[0]: load(p) for p in sys.argv[1:]}
for a, b in combinations(runs, 2):
    shared = set(runs[a]) & set(runs[b])
    same = sum(runs[a][k][1] == runs[b][k][1] for k in shared)
    da = [k for k in shared if runs[a][k][1] != runs[a][k][0]]
    both = sum(runs[b][k][1] == runs[a][k][1] for k in da)
    print(f"{a} vs {b}: shared codons {len(shared)}, identical read codon {same} ({same / max(len(shared), 1):.4f}); "
          f"of {len(da)} where {a} differs from the genome, {b} has the same codon at {both}")
