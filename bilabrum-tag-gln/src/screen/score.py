"""Per unit and stop codon: in-frame rate, aligned meaning, breadth. Reads stopsense.py scan output.

usage: score.py UNITS.tsv SCANDIR > table.tsv
Columns:
  rate_in    in-frame stop per 1,000 codons over all best-hit spans (codons.npz N_in)
  n          inner aligned columns with this stop and a standard reference residue (Sec U and X excluded)
  sec        inner TGA columns aligned to reference U (selenocysteine), reported only
  rate_al    n per 1,000 inner aligned columns
  prot       distinct reference proteins among the n columns
  top, k     most common reference residue and its count; share = k / n
  prot_top   distinct reference proteins where the stop aligns to top
  bg         top's share among all inner aligned columns (reference composition)
  enr        share / bg
  ident      how often top's own sense codons align to top (the sense-codon baseline)
  sense      share / ident: ~1 when the stop behaves like a codon for top
  second, share2   runner-up residue
  meaning    the amino acid whose own codons' reference-residue profile (in this unit) best explains the stop's
             residues (multinomial likelihood, +0.5 pseudocounts); llr = log-likelihood per column of that profile
             over the unit's reference composition; margin = per-column log-likelihood gap to the next amino acid
  kl, rel    kl = KL(meaning's codon profile || composition), the llr a true codon for the meaning would average;
             rel = llr / kl, roughly the fraction of the stop's occurrences that behave as a codon for the meaning
"""
import csv
import gzip
import os
import sys
from collections import Counter, defaultdict

import numpy as np

B = "TCAG"
CODONS = [a + b + c for a in B for b in B for c in B]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
LET = "ACDEFGHIKLMNPQRSTVWYUXBZJO"
STD = LET[:20]
STOPS = ("TAA", "TAG", "TGA")


def main(units_tsv, scandir):
    cols = ["unit", "stratum", "name", "group", "gc", "codon", "rate_in", "n", "sec", "rate_al", "prot", "top", "k",
            "share", "prot_top", "bg", "enr", "ident", "sense", "second", "share2", "meaning", "llr", "margin",
            "prot_meaning", "kl", "rel", "inner", "bodies", "hits"]
    print("\t".join(cols))
    for u in csv.DictReader(open(units_tsv), delimiter="\t"):
        f = os.path.join(scandir, u["unit"] + ".npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        M, N_in = z["M"], z["N_in"]
        Ms = M[:, :20]
        inner = int(Ms.sum())
        comp = Ms.sum(0) / max(inner, 1)
        aas = sorted(set(AA) - {"*"})
        prof = {}
        for a in aas:
            v = Ms[[j for j in range(64) if AA[j] == a]].sum(0) + 0.5
            prof[a] = np.log(v / v.sum())
        lbg = np.log((Ms.sum(0) + 0.5) / (Ms.sum() + 10))
        prot = defaultdict(set)
        prot_res = defaultdict(set)
        with gzip.open(os.path.join(scandir, u["unit"] + ".stops.tsv.gz"), "rt") as fh:
            next(fh)
            for line in fh:
                contig, s, c, a = line.split("\t")[:4]
                if a in STD:
                    prot[c].add(s)
                    prot_res[(c, a)].add(s)
        for c in STOPS:
            i = CODONS.index(c)
            row = Ms[i]
            n = int(row.sum())
            sec = int(M[i, LET.index("U")])
            rate_in = 1000 * N_in[i] / max(N_in.sum(), 1)
            if n:
                order = np.argsort(-row)
                t, t2 = order[0], order[1]
                share, share2 = row[t] / n, row[t2] / n
                a = STD[t]
                own = [j for j, cc in enumerate(CODONS) if AA[j] == a]
                tot_own = Ms[own].sum()
                ident = Ms[own, t].sum() / tot_own if tot_own else float("nan")
                bg = comp[t]
                ll = sorted(((float(row @ prof[b]), b) for b in aas), reverse=True)
                m = ll[0][1]
                llr = (ll[0][0] - float(row @ lbg)) / n
                pm = np.exp(prof[m])
                kl = float(pm @ (prof[m] - lbg))
                vals = [f"{rate_in:.3f}", n, sec, f"{1000 * n / inner:.3f}", len(prot[c]), a, int(row[t]),
                        f"{share:.3f}", len(prot_res[(c, a)]), f"{bg:.4f}", f"{share / bg:.2f}" if bg else "nan",
                        f"{ident:.3f}", f"{share / ident:.2f}" if ident else "nan", STD[t2], f"{share2:.3f}", m,
                        f"{llr:.3f}", f"{(ll[0][0] - ll[1][0]) / n:.3f}", len(prot_res[(c, m)]), f"{kl:.3f}",
                        f"{llr / kl:.3f}"]
            else:
                vals = [f"{rate_in:.3f}", 0, sec, "0.000", 0, "", 0, "", 0, "", "", "", "", "", "", "", "", "", 0,
                        "", ""]
            print("\t".join(map(str, [u["unit"], u["stratum"], u["name"], u["group"], u["gc"], c] + vals
                                  + [inner, int(z["bodies"]), int(z["hits"])])))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
