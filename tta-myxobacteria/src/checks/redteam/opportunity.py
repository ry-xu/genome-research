"""Mutational opportunity for TTA vs canonical stops, under different single-nucleotide mutation weights.

usage: opportunity.py <panel.tsv> <genes_dir>
Interior codons (residues 22 .. len-21) of complete genes. Weights: uniform; transitions x3; transitions x3 and
G/C->A/T x2. Ratio = opportunity(TTA) / opportunity(stop set), for stops {TAA,TAG,TGA} and {TAA,TAG}.
"""
import sys
from collections import Counter

STOPS = {"TAA", "TAG", "TGA"}
TS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def w(a, b, kappa, sw):
    x = kappa if (a, b) in TS else 1.0
    if a in "GC" and b in "AT":
        x *= sw
    return x


def opp(cu, targets, kappa, sw):
    tot = 0.0
    for c, v in cu.items():
        if c in STOPS or c in targets:
            continue
        for i in range(3):
            for b in "ACGT":
                if b != c[i] and c[:i] + b + c[i + 1:] in targets:
                    tot += v * w(c[i], b, kappa, sw)
    return tot


def main(panel_fn, gdir):
    print("genome\tmodel\tTTA/stop3\tTTA/TAA+TAG")
    for l in open(panel_fn):
        lab, acc, role = l.split("\t")[:3]
        if role != "clade":
            continue
        cu = Counter()
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, c, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                if pb == "0" and pe == "0" and len(nt) % 3 == 0:
                    for i in range(21, len(nt) // 3 - 21):
                        cu[nt[3 * i:3 * i + 3]] += 1
        for name, k, sw in (("uniform", 1, 1), ("ts x3", 3, 1), ("ts x3, GC->AT x2", 3, 2)):
            t = opp(cu, {"TTA"}, k, sw)
            print(f"{lab}\t{name}\t{t / opp(cu, STOPS, k, sw):.3f}\t{t / opp(cu, {'TAA', 'TAG'}, k, sw):.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
