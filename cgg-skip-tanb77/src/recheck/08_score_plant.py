"""Score the planted-frameshift control: was a frameshift of the right direction called within 9 nt of each planted
edit (03_analyze events), and for cgg_gly sites, which residue did 04_fs_detail put opposite the CGG+1 unit?"""
import sys
from collections import Counter, defaultdict
W = "recheck_work"
lab = sys.argv[1]
ev = defaultdict(list)
for line in list(open(f"{W}/out/{lab}.events.tsv"))[1:]:
    f = line.rstrip("\n").split("\t")
    ev[(f[0], f[1])].append((int(f[2]), int(f[3]), f[4], int(f[5])))
cg = defaultdict(list)
for line in list(open(f"{W}/out/{lab}.fs_cgg_sites.tsv"))[1:]:
    f = line.rstrip("\n").split("\t")
    cg[(f[0], f[1])].append((int(f[2]), f[3], f[8], f[7]))
det = Counter(); tot = Counter(); wrongdir = Counter(); dist = Counter(); res = Counter(); res_pairs = Counter()
for line in list(open(f"{W}/out/{lab}.sites.tsv"))[1:]:
    ctg, st, s, codon, raa, t, p0, p = line.rstrip("\n").split("\t")
    p = int(p); tot[t] += 1
    want = "-1" if t == "minus1" else "+1"
    hits = [(g0, g1, d, n) for (g0, g1, d, n) in ev[(ctg, st)] if g0 - 9 <= p <= g1 + 9]
    if any(d == want for (_, _, d, _) in hits):
        det[t] += 1
        dist[(t, min(abs(p - g0) for (g0, g1, d, n) in hits if d == want) // 3)] += 1
    elif hits:
        wrongdir[t] += 1
    if t.startswith("cgg_"):
        near = [x for x in cg[(ctg, st)] if abs(x[0] - p) <= 4]
        if near:
            res[(t, near[0][2] or "no-consensus")] += 1
            for r in near[0][3]: res_pairs[(t, r)] += 1
        else:
            res[(t, "no CGG site")] += 1
for t in tot:
    print(f"{lab} {t}: planted {tot[t]}, detected (right direction, <=9 nt) {det[t]} ({det[t]/tot[t]:.0%}), wrong direction only {wrongdir[t]}")
print("distance (codons) planted edit -> called event:", sorted(dist.items()))
print("cgg_* consensus residue opposite CGG+1 unit:", dict(res), " all clean pairs:", dict(res_pairs.most_common()))
