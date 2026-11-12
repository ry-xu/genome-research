"""Contig GC of homologous-continuation events (TTA from stop_background.tsv, STOP from both sides) in clade genomes.

usage: event_gc.py <genes_dir> <stop_background.tsv> <stop_background2.tsv>
Flags events on contigs whose GC is outside the genome's bp-weighted 5-95% contig-GC range.
"""
import sys
from collections import defaultdict


def main(gdir, f1, f2):
    rng, cgc = {}, {}

    def load(lab):
        if lab in rng:
            return
        rows = []
        for line in open(f"{gdir}/{lab}.contigs.tsv"):
            n, L, gc = line.rstrip("\n").split("\t")
            rows.append((float(gc), int(L)))
            cgc[(lab, n)] = (float(gc), int(L))
        rows.sort()
        tot = sum(L for _, L in rows)
        acc, lo, hi = 0, None, None
        for gc, L in rows:
            acc += L
            if lo is None and acc >= 0.05 * tot:
                lo = gc
            if hi is None and acc >= 0.95 * tot:
                hi = gc
        rng[lab] = (lo, hi, sum(g * L for g, L in rows) / tot)

    out = defaultdict(lambda: [0, 0])
    for fn, kcol, side in ((f1, "kind", "5p"), (f2, None, "3p")):
        with open(fn) as f:
            h = next(f).rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(h, line.rstrip("\n").split("\t")))
                if r["role"] != "clade" or r["homologous"] != "1":
                    continue
                kind = r[kcol] if kcol else "STOP"
                if kind == "CTA":
                    continue
                lab = r["genome"]
                load(lab)
                ctg = r["gene"].rsplit("_", 3)[0]
                gc, L = cgc[(lab, ctg)]
                lo, hi, g = rng[lab]
                flag = gc < lo or gc > hi
                out[kind][0] += 1
                out[kind][1] += flag
                if flag:
                    print(f"{lab}\t{kind}\t{side}\t{r['gene']}\t{r.get('codon', '')}\tcontig GC {gc:.3f} (genome {g:.3f}, "
                          f"5-95% {lo:.3f}-{hi:.3f}), len {L}")
    for k, v in out.items():
        print(f"{k}: {v[0]} event rows, {v[1]} on out-of-range-GC contigs")


if __name__ == "__main__":
    main(*sys.argv[1:4])
