"""For each CAT tRNA (tRNAscan -O calls): (1) matches over the full reference length (global alignment, end gaps free),
(2) neighbouring tRNAs (tRNAscan -O order) to place it in the conserved slots: Met-e before trnL(TAA), Ile-CAU between
trnW and trnQ in Saccharina/Desmarestia, Met-f elsewhere."""
import re
from Bio import SeqIO
from Bio.Align import PairwiseAligner
B = "lab/work/organelle-reading/checks/chor_wgs/dbcheck/"
exec(open(B + "d6_helpers.py").read())
for g in G:
    rows = parse_tscan(B + f"trna/{g}.tscanO.txt"); ss = parse_ss(B + f"trna/{g}.tscanO.ss")
    order = [(int(p[2]), p[4] + "-" + p[5]) for p in rows]
    print(f"-- {NAME[g]} tRNA order (tRNAscan -O, begin:type-anticodon):")
    print("   " + " ".join(f"{b}:{t}" for b, t in order))
    for i, p in enumerate(rows):
        if p[5] != "CAT":
            continue
        b, e, s = ss[int(p[1])]
        prev = order[i - 1][1] if i > 0 else "(start)"; nxt = order[i + 1][1] if i + 1 < len(order) else "(end)"
        res = []
        for rk, rs in REF.items():
            aln = al.align(s, rs)[0]
            m = sum(s[i2] == rs[j2] for (s1, e1), (s2, e2) in zip(*aln.aligned) for i2, j2 in zip(range(s1, e1), range(s2, e2)))
            res.append((m / len(rs), m, len(rs), rk))
        res.sort(reverse=True)
        print(f"   CAT {b}-{e} (len {len(s)}) neighbours: {prev} | {nxt}")
        print("      " + "; ".join(f"{rk} {m}/{L}={f:.3f}" for f, m, L, rk in res))
