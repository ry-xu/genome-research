"""Per-codon tallies of conserved-column consensus residues from records.tsv.

usage: tally_codons.py records.tsv ref_order.json outdir [min_nref=5] [min_frac=0.8]
Writes outdir/codon_table.tsv: one row per (acc, filter, codon) with counts, and outdir/conserved_records.tsv.
filters: mask (panel masked for MAFFT; primary), t4 (panel unmasked), mask_hspagree (MAFFT column == tblastn
HSP column of the anchoring reference residue).
"""
import os
import sys
import csv
import json
from collections import defaultdict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CODONS, T4  # noqa: E402


def consensus(res, nref_total, min_nref, min_frac):
    r = [c for c in res if c not in "-X*"]
    if len(r) < min_nref or len(r) < 0.5 * nref_total:
        return None, len(r), 0.0
    aa, n = Counter(r).most_common(1)[0]
    f = n / len(r)
    return (aa if f >= min_frac else None), len(r), f


def main():
    rec_path, order_path, outdir = sys.argv[1:4]
    min_nref = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    min_frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.8
    os.makedirs(outdir, exist_ok=True)
    order = json.load(open(order_path))
    usage = defaultdict(Counter)
    tally = defaultdict(lambda: defaultdict(Counter))  # (acc, filt) -> codon -> aa counts
    cons_rows = []
    rows = list(csv.DictReader(open(rec_path), delimiter="\t"))
    for r in rows:
        acc, mode, cod = r["acc"], r["mode"], r["codon"]
        nref = len(order[r["gene"]])
        if mode == "mask":
            usage[acc][cod] += 1
        aa, n, f = consensus(r["ref_res"], nref, min_nref, min_frac)
        filt = mode
        if aa is not None:
            tally[(acc, filt)][cod][aa] += 1
            tally[("ALL", filt)][cod][aa] += 1
            if mode == "mask":
                cons_rows.append((acc, r["gene"], r["idx"], r["gpos0"], cod, T4[cod], aa, n, f"{f:.2f}", r["col"], r["hsp_col"], r["ref_res"]))
        if mode == "mask" and r["hsp_col"] == r["col"] and aa is not None:
            tally[(acc, "mask_hspagree")][cod][aa] += 1
            tally[("ALL", "mask_hspagree")][cod][aa] += 1
    for acc in list(usage):
        usage["ALL"].update(usage[acc])
    with open(os.path.join(outdir, "codon_table.tsv"), "w") as fh:
        fh.write("acc\tfilter\tcodon\tt4aa\tusage\tN_conserved\tn_t4aa\tfrac_t4aa\ttop_aa\tn_top\tdist\n")
        for (acc, filt), d in sorted(tally.items()):
            for cod in CODONS:
                c = d.get(cod, Counter())
                N = sum(c.values())
                top, ntop = (c.most_common(1)[0] if N else ("-", 0))
                a4 = T4[cod]
                fh.write(f"{acc}\t{filt}\t{cod}\t{a4}\t{usage[acc][cod]}\t{N}\t{c.get(a4, 0)}\t"
                         f"{(c.get(a4, 0) / N if N else float('nan')):.3f}\t{top}\t{ntop}\t"
                         + ",".join(f"{k}:{v}" for k, v in c.most_common()) + "\n")
    with open(os.path.join(outdir, "conserved_records.tsv"), "w") as fh:
        fh.write("acc\tgene\tidx\tgpos0\tcodon\tt4aa\tcons_aa\tn_ref\tfrac\tcol\thsp_col\tref_res\n")
        for r in cons_rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
