"""Exploratory (post hoc, not in RULES): how much of the T1/Tf CGG skip signature is recurring conserved sites?
Window = 15 nt + CGG + 15 nt in gene orientation. Single-linkage clusters at Hamming <= 6 over the 33 nt.
Reports clusters by number of UHGV viral genera and vOTUs; also for all T1/Tf in-frame CGG (skip and non-skip)."""
import csv
from collections import defaultdict
CW = "phages_work"
rows = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
cgg = [r for r in rows(f"{CW}/kill1/out/meta.cgg_sites.tsv") if r["set"] in ("T1", "Tf")]
meta = {r["id"]: r for r in rows(f"{CW}/kill1/out/meta_sub.tsv")}
need = {r["id"] for r in cgg}
seqs, name, buf = {}, None, []
for l in open(f"{CW}/data/phages.fna"):
    if l[0] == ">":
        if name in need: seqs[name] = "".join(buf).upper()
        name, buf = l[1:].split()[0], []
    else: buf.append(l.strip())
if name in need: seqs[name] = "".join(buf).upper()
COMP = str.maketrans("ACGTN", "TGCAN")
W = []
for r in cgg:
    g = seqs[r["id"]]; L = len(g)
    b, e, st, k = int(r["begin"]), int(r["end"]), int(r["strand"]), int(r["k"])
    s = g if st == 1 else g.translate(COMP)[::-1]
    p = (b - 1 if st == 1 else L - e) + 3 * k
    w = s[p - 15:p + 18]
    if len(w) == 33: W.append((w, r))
par = list(range(len(W)))
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for i in range(len(W)):
    for j in range(i + 1, len(W)):
        if sum(a != b for a, b in zip(W[i][0], W[j][0])) <= 6:
            par[find(i)] = find(j)
cl = defaultdict(list)
for i in range(len(W)): cl[find(i)].append(i)
def gen(ix): return {meta[W[i][1]["id"]]["vgenus"] for i in ix if meta[W[i][1]["id"]]["vgenus"] != "vGENUS-NA"}
for lab, test in (("skip", lambda r: r["skip"] == "1"), ("non-skip", lambda r: r["skip"] == "0")):
    tot = sum(test(W[i][1]) for i in range(len(W)))
    bins = defaultdict(int); ncl = defaultdict(int)
    for c, ix in cl.items():
        ixs = [i for i in ix if test(W[i][1])]
        if not ixs: continue
        ng = len(gen(ix))
        b = "1 genus or none" if ng <= 1 else ("2-4 genera" if ng <= 4 else "5+ genera")
        bins[b] += len(ixs); ncl[b] += 1
    print(lab, "sites", tot, {b: (bins[b], ncl[b]) for b in sorted(bins)}, "(sites, clusters)")
big = sorted(cl.values(), key=len, reverse=True)[:8]
for ix in big:
    sk = sum(W[i][1]["skip"] == "1" for i in ix)
    print(len(ix), "sites", sk, "skip", len(gen(ix)), "genera", W[ix[0]][0][:15].lower(), W[ix[0]][0][15:18], W[ix[0]][0][18:].lower(), "k", W[ix[0]][1]["k"], "n", W[ix[0]][1]["n"])
print("top clusters: distinct windows, vOTUs, max Hamming to first, skip members' codon at CGG position always CGG by construction")
for ix in big[:6]:
    ws = {W[i][0] for i in ix}
    vo = {meta[W[i][1]["id"]]["votu"] for i in ix}
    mx = max(sum(a != b for a, b in zip(W[ix[0]][0], w)) for w in ws)
    print(len(ix), "windows", len(ws), "vOTUs", len(vo), "max_ham", mx)
