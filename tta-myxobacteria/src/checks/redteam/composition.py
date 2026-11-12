"""Per-genome codon position profile of TTA (and comparison Leu codons) in red-team gene calls.

usage: composition.py <genes_dir> <panel.tsv>
For complete genes (no partial end), d = codons upstream of the called stop (d = 1 is the last sense codon).
Bins: 1-10, 11-50, 51-100, 101-150, >150 (first 20 codons excluded). Per genome: GC of genes, counts of TTA, CTA, TTG,
CTG per bin and per 1e5 codons in the bin.
"""
import sys

BINS = [(1, 10), (11, 50), (51, 100), (101, 150), (151, 10**9)]
CODS = ["TTA", "CTA", "TTG", "CTT", "CTG"]


def main(gdir, panel):
    rows = [l.rstrip("\n").split("\t") for l in open(panel)]
    print("label\tgc\tngenes\t" + "\t".join(f"{c}_{a}-{b if b < 1e9 else 'end'}" for c in CODS for a, b in BINS)
          + "\t" + "\t".join(f"ncod_{a}" for a, b in BINS))
    for lab, *_ in rows:
        cnt = {(c, i): 0 for c in CODS for i in range(len(BINS))}
        tot = [0] * len(BINS)
        gcn = ntot = ng = 0
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
                if pb == "1" or pe == "1" or len(nt) % 3:
                    continue
                ng += 1
                gcn += nt.count("G") + nt.count("C")
                ntot += len(nt)
                cods = [nt[k:k + 3] for k in range(0, len(nt) - 3, 3)]  # sense codons, stop dropped
                L = len(cods)
                for k, c in enumerate(cods):
                    if k < 20:
                        continue
                    d = L - k
                    for i, (a, bb) in enumerate(BINS):
                        if a <= d <= bb:
                            tot[i] += 1
                            if c in CODS:
                                cnt[(c, i)] += 1
                            break
        out = [lab, f"{gcn / ntot:.3f}", str(ng)]
        out += [str(cnt[(c, i)]) for c in CODS for i in range(len(BINS))]
        out += [str(t) for t in tot]
        print("\t".join(out))


if __name__ == "__main__":
    main(*sys.argv[1:3])
