"""T3 summary: matched-divergence comparison, region profiles, and site intactness against the census consensus.
Usage: python -I src/t3_summary.py 23S|16S
"""
import sys, csv, statistics as st, math
sys.path.insert(0, "src")
from alnlib import read_sto, read_cols, short

kind = sys.argv[1]
rows = list(csv.DictReader(open(f"out/t3_{kind}_pairs.tsv"), delimiter="\t"))
for r in rows:
    for k in list(r):
        if k not in ("cat", "pair"):
            r[k] = float(r[k])
div = [r for r in rows if r["cat"] == "div-div"]
lo, hi = min(r["pid"] for r in div) - 0.02, max(r["pid"] for r in div) + 0.02
fun = [r for r in rows if r["cat"] == "functional" and lo <= r["pid"] <= hi]
nn = [r for r in rows if r["cat"] == "norm-norm fam"]
regs = [k[4:] for k in rows[0] if k.startswith("rel_")]


def med(xs):
    xs = [x for x in xs if x == x]
    return st.median(xs) if xs else float("nan")


print(f"== {kind}: div-div pid {lo+0.02:.3f}-{hi-0.02:.3f}; matched functional pairs (pid {lo:.3f}-{hi:.3f}): {len(fun)}")
for lab, S in [("div-div", div), ("matched functional", fun), ("norm-norm family", nn)] + \
        [(r["cat"], [r]) for r in rows if r["cat"] not in ("div-div", "functional", "norm-norm fam")]:
    print(f"{lab:30s} n={len(S):3d} pid~{med([r['pid'] for r in S]):.3f} rate_inv {med([r['rate_inv>=0.99'] for r in S]):.4f} "
          f"(max {max(r['rate_inv>=0.99'] for r in S):.4f}) rate_var {med([r['rate_var<0.80'] for r in S]):.3f} "
          f"inv_depl {med([r['inv_depl'] for r in S]):.3f} (range {min(r['inv_depl'] for r in S):.3f}-{max(r['inv_depl'] for r in S):.3f})")
print("region relative rates (median):")
print("  " + " ".join(f"{g:>8s}" for g in regs))
for lab, S in [("div-div", div), ("functional", fun), ("norm-norm", nn)]:
    print(f"  {lab:12s}" + " ".join(f"{med([r['rel_'+g] for r in S]):8.2f}" for g in regs))


# Spearman rho between div-div median profile and functional median profile
def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(o):
        r[i] = k
    return r


a = [med([r['rel_' + g] for r in div]) for g in regs]
b = [med([r['rel_' + g] for r in fun]) for g in regs]
ra, rb = rank(a), rank(b)
n = len(a)
rho = 1 - 6 * sum((x - y) ** 2 for x, y in zip(ra, rb)) / (n * (n * n - 1))
print(f"Spearman rho(div-div profile, functional profile) over {n} domains = {rho:.2f}")
