"""Neighbour check: on contigs carrying TTA-at-e*+1 events (runon_events, R = sister), what fraction of the OTHER
genes have an RBH partner in another clade genome, vs the same fraction on contigs carrying canonical "same"-end
events, vs genome-wide. A contaminant contig from a TTA-as-Leu organism would have few clade RBHs.

usage: neighbour_rbh.py <panel.tsv> <genes_dir> <ortho_dir>
"""
import sys
from collections import defaultdict


def main(panel_fn, gdir, odir):
    role = {l.split("\t")[0]: l.split("\t")[2] for l in open(panel_fn)}
    clade = [k for k, v in role.items() if v == "clade"]
    has = defaultdict(set)  # gene -> roles of RBH partners
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            q, s = line.split("\t")[:2]
            has[q].add(role[s.split("|")[0]])
    ev = defaultdict(dict)
    with open(f"{odir}/runon_events.tsv") as f:
        h = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(h, line.rstrip("\n").split("\t")))
            if r["R"] != "sister" or r["B_role"] != "clade":
                continue
            if r["class"] == "runon" and r["w_p1"] == "TTA":
                ev[r["B"]][r["gene"]] = "TTA"
            elif r["class"] == "same":
                ev[r["B"]][r["gene"]] = "STOP"
    print("genome\tclass\tcontigs\tother_genes\tfrac_clade_RBH\tfrac_any_RBH")
    tot = defaultdict(lambda: [0, 0, 0])
    for lab in clade:
        bycontig = defaultdict(list)
        with open(f"{gdir}/{lab}.genes.tsv") as f:
            next(f)
            for line in f:
                gid, ctg = line.split("\t")[:2]
                bycontig[ctg].append(gid)
        for cls in ("TTA", "STOP", "ALL"):
            if cls == "ALL":
                ctgs = set(bycontig)
                ev_genes = set()
            else:
                ev_genes = {g for g, c in ev[lab].items() if c == cls}
                ctgs = {g.rsplit("_", 3)[0] for g in ev_genes}
            others = [g for c in ctgs for g in bycontig[c] if g not in ev_genes]
            n = len(others)
            fc = sum(1 for g in others if "clade" in has[f"{lab}|{g}"]) / max(1, n)
            fa = sum(1 for g in others if has[f"{lab}|{g}"]) / max(1, n)
            tot[cls][0] += n
            tot[cls][1] += fc * n
            tot[cls][2] += fa * n
            print(f"{lab}\t{cls}\t{len(ctgs)}\t{n}\t{fc:.3f}\t{fa:.3f}")
    for cls, (n, c, a) in tot.items():
        print(f"pooled\t{cls}\t-\t{n}\t{c / n:.3f}\t{a / n:.3f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
