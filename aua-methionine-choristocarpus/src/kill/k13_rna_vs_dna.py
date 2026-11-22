"""Kill check 7b: are the k12 mito reads transcripts or leftover DNA? Compare read depth over annotated CDS, rRNA, tRNA and
intergenic NC_085326 positions (DNA gives roughly flat depth; RNA is gene-biased), and depth spread across CDS (per-gene mean).
usage: uv run python -I checks/chor_wgs/kill/k13_rna_vs_dna.py SAM GENBANK"""
import re, sys
from collections import defaultdict
from Bio import SeqIO
SAM, GB = sys.argv[1:3]
r = SeqIO.read(GB, "genbank"); L = len(r.seq)
dep = [0] * L
for l in open(SAM):
    if l.startswith("@"): continue
    f = l.split("\t"); flag, pos, mq, cig = int(f[1]), int(f[3]) - 1, int(f[4]), f[5]
    if flag & (4 | 256 | 2048) or mq < 20: continue
    rp = pos
    for n, op in re.findall(r"(\d+)([MIDNSHX=])", cig):
        n = int(n)
        if op in "M=XD":
            for i in range(rp, min(rp + n, L)): dep[i] += 1
        if op in "M=XDN": rp += n
cls = ["intergenic"] * L; gene = {}
for ft in r.features:
    if ft.type not in ("CDS", "rRNA", "tRNA"): continue
    name = (ft.qualifiers.get("gene") or ft.qualifiers.get("product") or ["?"])[0]
    for p in ft.location:
        cls[p] = ft.type
        if ft.type == "CDS": gene[p] = name
tot = defaultdict(lambda: [0, 0])
for i in range(L): tot[cls[i]][0] += dep[i]; tot[cls[i]][1] += 1
for k, (s, n) in sorted(tot.items()): print(f"{k}: {n} positions, mean depth {s / n:.2f}")
g = defaultdict(lambda: [0, 0])
for p, nm in gene.items(): g[nm][0] += dep[p]; g[nm][1] += 1
m = sorted((s / n, nm) for nm, (s, n) in g.items())
print(f"CDS genes {len(m)}: per-gene mean depth min {m[0][0]:.2f} ({m[0][1]}), median {m[len(m) // 2][0]:.2f}, max {m[-1][0]:.2f} ({m[-1][1]}); max/min {m[-1][0] / max(m[0][0], 0.01):.1f}")
print("lowest 3:", [(nm, round(v, 2)) for v, nm in m[:3]], "highest 3:", [(nm, round(v, 2)) for v, nm in m[-3:]])
