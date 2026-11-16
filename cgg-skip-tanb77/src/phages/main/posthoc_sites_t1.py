"""Post hoc (descriptive), after the Pfam-check agent 19:09: effective-n view of the CGG skip signature.
[T1-only variant of posthoc_sites.py: Tf excluded from T.] In-frame CGGs (kill1 meta.cgg_sites.tsv) of T1+Tf and R phages at GC < 32 (gc_bin <= 30), each as a 33-nt window
(15 + CGG + 15, gene orientation, as the Pfam-check agent's recur.py). Single-linkage clusters at Hamming <= 6 over T and R together.
Reports (a) clusters called skip (>= half of the group's members skip) per group, (b) skip/O for CGGs in clusters
seen in one viral genus only, (c) clusters holding both T and R members (same conserved site, two host sides)."""
import csv, sys, numpy as np
H = int(sys.argv[1]) if len(sys.argv) > 1 else 6
from collections import defaultdict
CW = "phages_work"
rows = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
grp = {"T1": "T", "R": "R"}
cgg = [r for r in rows(f"{CW}/kill1/out/meta.cgg_sites.tsv") if r["set"] in grp and int(r["gc_bin"]) <= 30]
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
W, R = [], []
for r in cgg:
    g = seqs[r["id"]]; L = len(g)
    b, e, st, k = int(r["begin"]), int(r["end"]), int(r["strand"]), int(r["k"])
    s = g if st == 1 else g.translate(COMP)[::-1]
    p = (b - 1 if st == 1 else L - e) + 3 * k
    w = s[p - 15:p + 18]
    if len(w) == 33 and w[15:18] == "CGG": W.append(w); R.append(r)
print("windows:", len(W), "dropped:", len(cgg) - len(W))
A = np.frombuffer("".join(W).encode(), dtype=np.uint8).reshape(len(W), 33)
n = len(W); par = np.arange(n)
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for i0 in range(0, n, 256):
    d = (A[i0:i0 + 256, None, :] != A[None, :, :]).sum(2)
    for a, b in zip(*np.nonzero(d <= H)):
        i = i0 + a
        if b > i:
            ra, rb = find(i), find(b)
            if ra != rb: par[ra] = rb
cl = defaultdict(list)
for i in range(n): cl[find(i)].append(i)
G = lambda i: grp[R[i]["set"]]; S = lambda i: R[i]["skip"] == "1"
gen = lambda i: meta[R[i]["id"]]["vgenus"] if meta[R[i]["id"]]["vgenus"] != "vGENUS-NA" else R[i]["id"]
print("clusters:", len(cl))
for g in ("T", "R"):
    sites = [i for i in range(n) if G(i) == g]
    print(f"\n== {g}: CGGs {len(sites)}, skip {sum(S(i) for i in sites)}")
    cs = [[i for i in m if G(i) == g] for m in cl.values()]; cs = [c for c in cs if c]
    sk = sum(1 for c in cs if sum(S(i) for i in c) * 2 >= len(c))
    print(f"(a) clusters with {g} members: {len(cs)}; called skip: {sk} = {sk/len(cs):.3f}")
    one = [c for c in cs if len({gen(i) for i in c}) == 1]
    O1 = sum(len(c) for c in one); S1 = sum(S(i) for c in one for i in c)
    sk1 = sum(1 for c in one if sum(S(i) for i in c) * 2 >= len(c))
    print(f"(b) one-genus clusters: {len(one)} clusters, skip clusters {sk1} = {sk1/len(one):.3f}; CGGs skip/O {S1}/{O1} = {S1/O1:.3f}")
    many = [c for c in cs if len({gen(i) for i in c}) >= 5]
    print(f"    clusters in >= 5 genera: {len(many)}, CGGs {sum(len(c) for c in many)}, skip {sum(S(i) for c in many for i in c)}")
mixed = [m for m in cl.values() if {G(i) for i in m} == {"T", "R"}]
print(f"\n(c) clusters with both T and R members: {len(mixed)}")
for m in sorted(mixed, key=len, reverse=True):
    t = [i for i in m if G(i) == "T"]; r = [i for i in m if G(i) == "R"]
    print(f"  T skip {sum(S(i) for i in t)}/{len(t)} ({len({gen(i) for i in t})} genera) | R skip {sum(S(i) for i in r)}/{len(r)} ({len({gen(i) for i in r})} genera) | {W[m[0]]}")
