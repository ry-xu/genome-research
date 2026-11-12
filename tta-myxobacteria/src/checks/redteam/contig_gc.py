"""Contamination check 1: GC and length of contigs carrying Arm B TTA events vs other Arm B event contigs.

usage: contig_gc.py <panel.tsv> <genes_dir> <armB_events.tsv> [more events files]
Arm B gene ids are <contig>_<begin>_<end>_<strand>; window = e*-1..e*+3, TTA event = TTA anywhere in it.
Per clade genome: genome GC; median contig GC / length for TTA-event contigs and non-TTA-event contigs; Mann-Whitney U
(normal approximation) of contig GC; fraction of TTA-event contigs outside the genome's central 90% GC range (by bp).
"""
import math
import sys
from statistics import median


def mwu(x, y):
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks = {}
    i = 0
    rx = 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            if allv[k][1] == 0:
                rx += r
        i = j
    n1, n2 = len(x), len(y)
    u = rx - n1 * (n1 + 1) / 2
    mu, sd = n1 * n2 / 2, math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (u - mu) / sd if sd else 0
    p = math.erfc(abs(z) / math.sqrt(2))
    return u, z, p


def main(panel, gdir, *evfiles):
    labs = {r.split("\t")[1]: r.split("\t")[0] for r in open(panel).read().splitlines() if r.split("\t")[2] == "clade"}
    ev = {}
    for fn in evfiles:
        with open(fn) as f:
            next(f)
            for line in f:
                r = line.rstrip("\n").split("\t")
                if r[0] in labs:
                    ctg = r[2].rsplit("_", 3)[0]
                    tta = "TTA" in r[8].split(",")
                    ev.setdefault(r[0], []).append((ctg, tta))
    print("acc\tgenome_gc\tn_ctg_TTA\tn_ctg_other\tmedGC_TTA\tmedGC_other\tmedLen_TTA\tmedLen_other\tMWU_z\tMWU_p\t"
          "gc_p05\tgc_p95\tTTA_ctg_outside_90pct")
    for acc, lab in labs.items():
        ctg = {}
        for line in open(f"{gdir}/{lab}.contigs.tsv"):
            n, L, gc = line.split("\t")
            ctg[n] = (int(L), float(gc))
        tot = sum(L for L, _ in ctg.values())
        ggc = sum(L * g for L, g in ctg.values()) / tot
        # bp-weighted GC quantiles
        srt = sorted(ctg.values(), key=lambda x: x[1])
        acc_bp, q05, q95 = 0, None, None
        for L, g in srt:
            acc_bp += L
            if q05 is None and acc_bp >= 0.05 * tot:
                q05 = g
            if q95 is None and acc_bp >= 0.95 * tot:
                q95 = g
        evs = ev.get(acc, [])
        t = {c for c, x in evs if x}
        o = {c for c, x in evs if not x} - t
        gt, go = [ctg[c][1] for c in t], [ctg[c][1] for c in o]
        lt, lo = [ctg[c][0] for c in t], [ctg[c][0] for c in o]
        u, z, p = mwu(gt, go)
        out = sum(1 for g in gt if g < q05 or g > q95)
        print(f"{acc}\t{ggc:.3f}\t{len(t)}\t{len(o)}\t{median(gt):.3f}\t{median(go):.3f}\t{median(lt)}\t{median(lo)}\t"
              f"{z:.2f}\t{p:.3f}\t{q05:.3f}\t{q95:.3f}\t{out}/{len(gt)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
