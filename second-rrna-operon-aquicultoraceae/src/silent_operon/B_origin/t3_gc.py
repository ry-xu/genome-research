"""GC content of divergent vs normal copies at census-paired and unpaired columns, and share of G-C among canonical pairs.
Usage: python -I src/t3_gc.py 23S|16S"""
import sys
sys.path.insert(0, "src")
from alnlib import read_sto, read_cols, short
kind = sys.argv[1]
aln = read_sto("out/aln23.sto" if kind == "23S" else "out/aln16.sto"); cols = read_cols(kind)
A = {(short(n) if "|" in n else n): s for n, s in aln.items()}
src = open("src/t3_profile.py").read(); exec(src[src.index('if kind == "23S":'):src.index("epos = ")])
ctl = {"23S": {"tt_div": "tt_23S_NZ_CP041345.1_1989407-", "tt_norm": "tt_23S_NZ_CP041345.1_533890+", "Ars_decay": "ar_23S_JBQMRP010000030.1_14720+", "Ars_norm": "ar_23S_JBQMRP010000004.1_315056+"},
       "16S": {"tt_div": "tt_16S_NZ_CP041345.1_1992611-", "tt_norm": "tt_16S_NZ_CP041345.1_531966+", "Ars_decay": "ar_16S_JBQMRP010000030.1_13534+", "Ars_norm": "ar_16S_JBQMRP010000004.1_313238+"}}[kind]
seqs = {**{f"div_{g}": A[k] for g, k in DIV.items() if k in A}, **{f"norm_{g}": A[k] for g, k in OWN.items() if k in A}, **{k: A[v] for k, v in ctl.items() if v in A}}
part = [r["partner"] for r in cols]
print(f"== {kind}: GC paired | GC unpaired | GC-pair share of canonical pairs | GU share")
for nm, s in seqs.items():
    P = [s[i] for i in range(len(s)) if part[i] >= 0 and s[i] not in "-.N"]
    U = [s[i] for i in range(len(s)) if part[i] < 0 and s[i] not in "-.N"]
    gcp = sum(c in "GC" for c in P) / len(P); gcu = sum(c in "GC" for c in U) / len(U)
    pr = [s[i] + s[part[i]] for i in range(len(s)) if part[i] > i and s[i] not in "-.N" and s[part[i]] not in "-.N"]
    can = [p for p in pr if p in ("AU", "UA", "GC", "CG", "GU", "UG")]
    print(f"  {nm:10s} {gcp:.3f} | {gcu:.3f} | {sum(p in ('GC','CG') for p in can)/len(can):.3f} | {sum(p in ('GU','UG') for p in can)/len(can):.3f}")
