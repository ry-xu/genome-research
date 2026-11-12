"""Inspect every interior aligned TTA in clade genes (from ortho_end interior_tta_sites.tsv, R = sister).

usage: interior_sites.py <genes_dir> <ortho_dir> <panel.tsv>
For each site: global protein alignment (BLOSUM62, open -10, extend -0.5) of the clade gene to its best sister RBH
partner and to its RBH partners in every other clade genome. Reports: the partner residue and codon aligned to the
TTA, identity in +-20 aligned columns, the fraction of the clade protein after the TTA that aligns to the partner, and
whether partners carry TTA, a canonical stop or another codon there.
"""
import sys
import warnings
from collections import defaultdict

from Bio import Align
from Bio.Align import substitution_matrices

warnings.filterwarnings("ignore")


def main(gdir, odir, panel_fn):
    panel = [l.rstrip("\n").split("\t") for l in open(panel_fn)]
    labs = [r[0] for r in panel]
    clade = {r[0] for r in panel if r[2] == "clade"}
    nt = {}
    for lab in labs:
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                r = line.rstrip("\n").split("\t")
                nt[f"{lab}|{r[0]}"] = r[7]
    code = {}
    b = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    k = 0
    for x in b:
        for y in b:
            for z in b:
                code[x + y + z] = aas[k]
                k += 1

    def prot(g):
        s = nt[g]
        return "".join(code.get(s[i:i + 3], "X") for i in range(0, len(s) - 3, 3))

    rbh = defaultdict(dict)  # gene -> genome -> partner
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            q, s = line.split("\t")[:2]
            rbh[q][s.split("|")[0]] = s
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    al.mode = "global"
    al.end_gap_score = 0  # free end gaps
    sites = []
    with open(f"{odir}/interior_tta_sites.tsv") as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[0] in clade and r[1] == "sister":
                sites.append((r[0], r[2], int(r[3]), int(r[4]), r[8]))
    sites = sorted(set(sites))
    print("B\tgene\tres\tql\tpartner\tpartner_aa\tpartner_codon\tlocal_id_pm20\tafter_aligned_frac\tsister_or_clade")
    for B, gene, res, ql, best_sis in sites:
        g = f"{B}|{gene}"
        p = prot(g)
        assert nt[g][3 * (res - 1):3 * res] == "TTA"
        partners = [(best_sis, "sister")] + [(rbh[g][o], "clade") for o in sorted(clade) if o != B and o in rbh[g]]
        for part, kind in partners:
            pp = prot(part)
            a = al.align(p, pp)[0]
            # map query index -> target index
            m = {}
            c1, c2 = a.coordinates
            for k in range(len(c1) - 1):
                if c1[k + 1] - c1[k] and c2[k + 1] - c2[k]:
                    for j in range(c1[k + 1] - c1[k]):
                        m[c1[k] + j] = c2[k] + j
            qi = res - 1
            ti = m.get(qi)
            win = [i for i in range(qi - 20, qi + 21) if i in m and i != qi]
            lid = sum(p[i] == pp[m[i]] for i in win) / max(1, len(win))
            after = [i for i in range(qi + 1, len(p)) if i in m]
            af = len(after) / max(1, len(p) - qi - 1)
            if ti is None:
                taa, tcod = "-", "-"
            else:
                taa, tcod = pp[ti], nt[part][3 * ti:3 * ti + 3]
            print(f"{B}\t{gene}\t{res}\t{ql}\t{part}\t{taa}\t{tcod}\t{lid:.2f}\t{af:.2f}\t{kind}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
