"""Per-gene view of one contig: TTA/CTA counts and RBH identity to sister and clade genomes.

usage: contig_inspect.py <panel.tsv> <genes_dir> <ortho_dir> <label> <contig>
"""
import sys
from collections import defaultdict


def main(panel_fn, gdir, odir, lab, ctg):
    role = {l.split("\t")[0]: l.split("\t")[2] for l in open(panel_fn)}
    rb = defaultdict(lambda: {"clade": (0, "-"), "control": (0, "-")})
    with open(f"{odir}/rbh.tsv") as f:
        next(f)
        for line in f:
            r = line.rstrip("\n").split("\t")
            if not r[0].startswith(lab + "|"):
                continue
            grp = role[r[1].split("|")[0]]
            pid = float(r[2])
            if pid > rb[r[0]][grp][0]:
                rb[r[0]][grp] = (pid, r[1].split("|")[0])
    for line in open(f"{gdir}/{lab}.contigs.tsv"):
        n, L, gc = line.rstrip("\n").split("\t")
        if n == ctg:
            print(f"{ctg} length {L} GC {gc}")
    print("gid\tlen_aa\tTTA\tCTA\tGC3\tbest_clade_pid\tbest_sister_pid")
    with open(f"{gdir}/{lab}.genes.tsv") as f:
        next(f)
        for line in f:
            gid, c, b, e, st, pb, pe, nt = line.rstrip("\n").split("\t")
            if c != ctg:
                continue
            cods = [nt[i:i + 3] for i in range(0, len(nt) - 3, 3)]
            gc3 = sum(1 for x in cods if x[2] in "GC") / max(1, len(cods))
            g = f"{lab}|{gid}"
            print(f"{gid}\t{len(cods)}\t{cods.count('TTA')}\t{cods.count('CTA')}\t{gc3:.2f}\t"
                  f"{rb[g]['clade'][0]} {rb[g]['clade'][1]}\t{rb[g]['control'][0]} {rb[g]['control'][1]}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
