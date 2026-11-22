"""Deviation from census consensus at conserved columns, split by what the nucleotide touches in E. coli 70S (7K00).
Classes (priority order): ligand (tRNA/mRNA), bridge (other subunit), protein (r-protein only), free (none).
Expected per class from the copy's own rate in census-conservation bins (0.95-0.99, 0.99-0.999, >=0.999), so classes
are compared at equal conservation.  Also: changes shared by all divergent copies (ancestral) vs not.
Usage: python -I src/t3_classes.py 23S|16S
"""
import sys, csv, math
sys.path.insert(0, "src")
from alnlib import read_sto, read_cols, short

kind = sys.argv[1]
aln = read_sto("out/aln23.sto" if kind == "23S" else "out/aln16.sto")
cols = read_cols(kind)
A = {(short(n) if "|" in n else n): s for n, s in aln.items()}
src = open("src/t3_profile.py").read()
exec(src[src.index('if kind == "23S":'):src.index("epos = ")])
ctl = {"23S": {"tt_div": "tt_23S_NZ_CP041345.1_1989407-", "Ars_decay": "ar_23S_JBQMRP010000030.1_14720+",
               "Ars2_decay": "as_23S_JBHEQW010000092.1_13747+"},
       "16S": {"tt_div": "tt_16S_NZ_CP041345.1_1992611-", "tb_typeB": "tb_16S_NC_014165.1_2893291-",
               "Ars_decay": "ar_16S_JBQMRP010000030.1_13534+"}}[kind]
seqs = {**{f"div_{g}": A[k] for g, k in DIV.items() if k in A}, **{f"norm_{g}": A[k] for g, k in OWN.items() if k in A},
        **{k: A[v] for k, v in ctl.items() if v in A}}

con = {}
for r in csv.DictReader(open("out/contacts_7k00.tsv"), delimiter="\t"):
    if r["kind"] == kind:
        p = int(r["ecoli_pos"])
        con[p] = "ligand" if r["ligand"] == "1" else "bridge" if r["bridge"] == "1" else "protein" if r["protein"] == "1" else "free"
epos = [r["ecoli_pos"] for r in cols]; f = [r["f_cons"] for r in cols]; cons = [r["cons"].replace("T", "U") for r in cols]
BINS = [(0.95, 0.99), (0.99, 0.999), (0.999, 1.01)]
part = [r["partner"] for r in cols]
STRAT = len(sys.argv) > 2 and sys.argv[2] == "paired"
def b_of(x, i=None):
    for k, (lo, hi) in enumerate(BINS):
        if lo <= x < hi:
            return k + (3 if (STRAT and i is not None and part[i] >= 0) else 0)
    return None
idx = [i for i in range(len(cols)) if f[i] >= 0.95 and epos[i] in con]
CL = ["ligand", "bridge", "protein", "free"]
print(f"== {kind}: conserved columns (f>=0.95) with a structure position: {len(idx)}; by class: " +
      ", ".join(f"{c} {sum(con[epos[i]] == c for i in idx)}" for c in CL))
print(f"{'copy':12s} " + " ".join(f"{c:>18s}" for c in CL) + "   (obs/exp, ratio)")
agg = {c: [0, 0.0] for c in CL}
for nm, s in seqs.items():
    # per-bin rate over all classes
    bn = [0] * 6; bd = [0] * 6
    for i in idx:
        if s[i] in "-.N":
            continue
        k = b_of(f[i], i); bn[k] += 1; bd[k] += s[i] != cons[i]
    rate = [bd[k] / bn[k] if bn[k] else 0 for k in range(6)]
    cells = []
    for c in CL:
        o = e = 0
        for i in idx:
            if con[epos[i]] != c or s[i] in "-.N":
                continue
            o += s[i] != cons[i]; e += rate[b_of(f[i], i)]
        cells.append(f"{o:4d}/{e:6.1f} {o / e if e else float('nan'):5.2f}")
        if nm.startswith("div_"):
            agg[c][0] += o; agg[c][1] += e
    print(f"{nm:12s} " + " ".join(f"{x:>18s}" for x in cells))
print("all div copies: " + ", ".join(f"{c} {agg[c][0]}/{agg[c][1]:.1f} = {agg[c][0]/agg[c][1]:.2f}" for c in CL))

# ancestral (all div copies deviate, same base) vs lineage-specific, per class
dv = [seqs[k] for k in seqs if k.startswith("div_")]
print("\nper class, columns where divergent copies deviate from consensus:")
for c in CL:
    cc = [i for i in idx if con[epos[i]] == c and all(s[i] not in "-.N" for s in dv)]
    allsame = sum(1 for i in cc if len({s[i] for s in dv}) == 1 and dv[0][i] != cons[i])
    anydev = sum(1 for i in cc if any(s[i] != cons[i] for s in dv))
    print(f"  {c:8s} cols {len(cc):4d}  any copy deviates {anydev:4d} ({anydev/len(cc):.2f})  all copies share the same change {allsame:4d} ({allsame/len(cc):.2f})")

# per r-protein contact set (div copies pooled), same expectation
pn = {}
for r in csv.DictReader(open("out/contacts_7k00.tsv"), delimiter="\t"):
    if r["kind"] == kind and r["proteins"]:
        pn[int(r["ecoli_pos"])] = r["proteins"].split(",")
print("\nper r-protein: pooled divergent copies obs/exp at its contact columns (f>=0.95)")
res = {}
for nm, s in seqs.items():
    if not nm.startswith("div_"):
        continue
    bn = [0] * 6; bd = [0] * 6
    for i in idx:
        if s[i] in "-.N":
            continue
        k = b_of(f[i], i); bn[k] += 1; bd[k] += s[i] != cons[i]
    rate = [bd[k] / bn[k] if bn[k] else 0 for k in range(6)]
    for i in idx:
        if s[i] in "-.N" or epos[i] not in pn:
            continue
        for p in pn[epos[i]]:
            a = res.setdefault(p, [0, 0.0, set()]); a[0] += s[i] != cons[i]; a[1] += rate[b_of(f[i], i)]; a[2].add(i)
for p, (o, e, cc) in sorted(res.items(), key=lambda x: x[1][0] / max(x[1][1], 1e-9)):
    if e >= 3:
        print(f"  {p:5s} cols {len(cc):3d} obs {o:3d} exp {e:5.1f} ratio {o/e:.2f}")
