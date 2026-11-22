"""T1: closest outside hits. For each query rRNA, best hits outside Aquicultoraceae (and outside the order set) with
identity summed over non-overlapping HSPs (>= 1.2 kb aligned for 23S, >= 0.8 kb for 16S), and the query's identity
to its own genome's other copy class.  Usage: python -I src/t1_hits.py
"""
import csv, gzip
from collections import defaultdict

aq = {r["acc"][3:]: r for r in csv.DictReader(open("lab/work/frozen-core/results/lsu/bac/cca/aq/genomes.tsv"), delimiter="\t")}
tax = {}
with open("lab/work/frozen-core/results/lsu/scan23/gtdb_tax.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        tax[p[0]] = p[3]


def agg(fn, ssu):
    hs = defaultdict(list)
    for line in open(fn):
        p = line.rstrip("\n").split("\t")
        q, s, pid, L, qs, qe, ss, se, ev, bs = p[:10]
        t = p[10] if ssu else tax.get(s, "?")
        hs[(q, s)].append((float(bs), int(L), float(pid), int(qs), int(qe), t))
    best = defaultdict(list)
    for (q, s), lst in hs.items():
        lst.sort(reverse=True); took = []
        for h in lst:
            if all(h[4] < t[3] or h[3] > t[4] for t in took):
                took.append(h)
        L = sum(t[1] for t in took); pid = sum(t[1] * t[2] for t in took) / L
        best[q].append((L * pid / 100, L, pid, s, took[0][5]))
    return best


for fn, ssu, minL in [("out/q23_vs_census.tsv", False, 1200), ("out/q16_vs_ssuall.tsv", True, 800)]:
    print(f"== {fn}")
    best = agg(fn, ssu)
    for q in sorted(best):
        rows = [b for b in best[q] if "f__Aquicultoraceae" not in b[4] and b[1] >= minL]
        rows.sort(key=lambda b: -b[2])
        fam = [b for b in best[q] if "f__Aquicultoraceae" in b[4]]
        print(f"{q}  [{len(fam)} in-family subjects; best in-family {max((b[2] for b in fam), default=0):.1f}]")
        for m, L, pid, s, t in rows[:5]:
            tt = ";".join(t.split(" ")[1 if ssu else 0].split(";")[1:6]) if ssu else ";".join(t.split(";")[1:6])
            print(f"   {pid:5.1f}% over {L:5d}  {s.split('~')[0][:22]:22s} {tt}")
