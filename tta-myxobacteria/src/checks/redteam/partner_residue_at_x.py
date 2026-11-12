"""Partner residue aligned with the premature stop X for homologous 3'-side events (stop_background2.tsv).

usage: partner_residue_at_x.py <panel.tsv> <genes_dir> <ortho_dir> <stop_background2.tsv>
Prints, per genome and stop codon, the counts of partner residues at the index aligned with X (partner_idx, 0-based).
"""
import sys
from collections import Counter, defaultdict

B = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {}
k = 0
for x in B:
    for y in B:
        for z in B:
            CODE[x + y + z] = AAS[k]
            k += 1


def main(panel_fn, gdir, odir, fn):
    labs = [l.split("\t")[0] for l in open(panel_fn)]
    role = {l.split("\t")[0]: l.split("\t")[2] for l in open(panel_fn)}
    nt = {}
    for lab in labs:
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                r = line.rstrip("\n").split("\t")
                nt[f"{lab}|{r[0]}"] = r[7]
    best = {}
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            if role[r[1].split("|")[0]] == "clade":
                continue
            if r[0] not in best or float(r[9]) > best[r[0]][0]:
                best[r[0]] = (float(r[9]), r[1])
    out = defaultdict(Counter)
    with open(fn) as f:
        next(f)
        for line in f:
            g, rl, gene, cod, j, P0, pid, sc, nl, hom = line.rstrip("\n").split("\t")
            if hom != "1":
                continue
            s = best[f"{g}|{gene}"][1]
            p = "".join(CODE.get(nt[s][i:i + 3], "X") for i in range(0, len(nt[s]) - 3, 3))
            P0 = int(P0)
            win = p[max(0, P0 - 1):P0 + 2]
            out[(g, cod)][p[P0]] += 1
            print(f"{g}\t{cod}\t{gene}\tpartner {s}\tres_at_X {p[P0]}\tneighbourhood {win}")
    print()
    for k2, c in sorted(out.items()):
        print(k2, dict(c))


if __name__ == "__main__":
    main(*sys.argv[1:5])
