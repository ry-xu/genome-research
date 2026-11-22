"""Structure-free pairing test for the divergent 23S, with a site-rate-matched neutral null.
1. Census columns -> positions in DASXVE normal 23S (from its cmalign row; the normal copy aligns confidently).
2. Divergent copy placed on the normal copy by Biopython global pairwise alignment (no structure used).
3. For census base pairs (columns.tsv partner) with both residues present in the divergent copy: fraction canonical
   (WC/GU), compared with the normal copy and with neutral nulls that mutate the normal copy at the same overall identity:
   (a) uniform rate; (b) per-column rate proportional to census variability (1 - f_cons), floor 0.001.
4. Changed-pair breakdown: pairs where div differs from normal at one side vs both sides, fraction canonical.
usage: uv run --with biopython python -I pairing_free.py NORM_ROWNAME DIV_FASTA [N_SIM]"""
import random, sys
from Bio import Align
CAN = {"AU", "UA", "GC", "CG", "GU", "UG"}
normrow, divfa = sys.argv[1], sys.argv[2]
nsim = int(sys.argv[3]) if len(sys.argv) > 3 else 200
cols = {}
for l in open("../../columns.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "col": continue
    cols[int(f[0])] = (f[4], float(f[5]), int(f[7]))
row = ""
for l in open("../aq/lsu_loci.sto"):
    if l.startswith(normrow + " ") : row += l.split()[1]
row = row.upper().replace("T", "U")
nseq = "".join(c for c in row if c.isalpha())
col2n, k = {}, 0
for i, c in enumerate(row):
    if c.isalpha(): col2n[i] = k; k += 1
rd = lambda p: "".join(l.strip() for l in open(p) if not l.startswith(">")).upper().replace("T", "U")
dseq = rd(divfa)
al = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
al.end_gap_score = 0
aln = al.align(nseq, dseq)[0]
n2d = {}
for (a, b), (c, d) in zip(*aln.aligned):
    for t in range(b - a): n2d[a + t] = c + t
pairs = [(i, j) for i, (cons, fc, j) in cols.items() if j > i and i in col2n and j in col2n]
def frac(seqget):
    n = can = 0
    for i, j in pairs:
        a, b = seqget(i), seqget(j)
        if a is None or b is None: continue
        n += 1; can += (a + b) in CAN
    return n, can / n
norm_get = lambda c: nseq[col2n[c]]
div_get = lambda c: dseq[n2d[col2n[c]]] if col2n[c] in n2d else None
# identity over mapped census columns
mapped = [c for c in col2n if col2n[c] in n2d]
ident = sum(norm_get(c) == div_get(c) for c in mapped) / len(mapped)
print(f"div {divfa}: census columns in normal {len(col2n)}, mapped to div {len(mapped)}, identity over mapped columns {ident:.3f}")
print(f"normal copy: pairs {frac(norm_get)[0]}, f_canonical {frac(norm_get)[1]:.3f}")
nd, fd = frac(div_get)
print(f"divergent copy (pairwise placement): pairs {nd}, f_canonical {fd:.3f}")
s1 = s1c = s2 = s2c = 0
for i, j in pairs:
    a, b = div_get(i), div_get(j)
    if a is None or b is None: continue
    da, db = a != norm_get(i), b != norm_get(j)
    if da and db: s2 += 1; s2c += (a + b) in CAN
    elif da or db: s1 += 1; s1c += (a + b) in CAN
print(f"changed pairs vs normal: one side {s1} (canonical {s1c/s1:.3f}), both sides {s2} (canonical {s2c/s2:.3f})")
# nulls
rng = random.Random(7)
var = {c: max(0.001, 1 - cols[c][1]) for c in mapped}
target = 1 - ident
def sim(rate):
    out = []
    for _ in range(nsim):
        m = {}
        for c in mapped:
            b = norm_get(c)
            if rng.random() < rate[c]: b = rng.choice([x for x in "ACGU" if x != b])
            m[c] = b
        n = can = 0
        for i, j in pairs:
            if i in m and j in m:
                n += 1; can += (m[i] + m[j]) in CAN
        out.append(can / n)
    out.sort(); return out
uni = {c: target for c in mapped}
s = target / (sum(var.values()) / len(var))
prop = {c: min(1.0, var[c] * s) for c in mapped}
for nm, r in [("uniform", uni), ("variability-proportional", prop)]:
    o = sim(r)
    print(f"null {nm} (mean sub rate {sum(r.values())/len(r):.3f}): f_canonical median {o[len(o)//2]:.3f}, max {o[-1]:.3f} of {nsim}; observed {fd:.3f}")
# paired-column sub rate in the observed div vs null expectation
pc = set(x for p in pairs for x in p)
obs_pc = sum(norm_get(c) != div_get(c) for c in mapped if c in pc) / sum(1 for c in mapped if c in pc)
obs_uc = sum(norm_get(c) != div_get(c) for c in mapped if c not in pc) / sum(1 for c in mapped if c not in pc)
exp_pc = sum(prop[c] for c in mapped if c in pc) / sum(1 for c in mapped if c in pc)
print(f"substitution rate at paired columns {obs_pc:.3f} vs unpaired {obs_uc:.3f}; variability-null expects paired {exp_pc:.3f}")
