"""C2 (RULES.md): Altair's scan.py logic (W=31, MINCOV 30, gap>=15 inside window+-10, loci = best dG per 2W stretch),
with ViennaRNA DNA parameters (Mathews 2004) instead of RNA. Same two SAMs.
usage: agents/altair/rnaenv/bin/python -I c2.py out_prefix sam.zst ..."""
import sys, re, subprocess, numpy as np, RNA
RNA.params_load_DNA_Mathews2004()
GENOME = "/Users/ryan/lab/work/archaeal-introns/culture/data/ref/genome.fna"
G = "".join(l.strip() for l in open(GENOME) if l[0] != ">").upper(); N = len(G)
MINCOV, MAPQ, GAPMIN, PAD, FLANK, W = 30, 20, 15, 5, 10, 31
CIG = re.compile(r"(\d+)([MIDNSHP=X])"); RC = str.maketrans("ACGT", "TGCA")
def dg(s): return min(RNA.fold(s)[1], RNA.fold(s.translate(RC)[::-1])[1])
out, sams = sys.argv[1], sys.argv[2:]
blocks, gaps = [], []
for sam in sams:
    p = subprocess.Popen(["zstd", "-dc", sam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        if line[0] == "@": continue
        f = line.split("\t", 6); flag = int(f[1])
        if flag & 0x904 or int(f[4]) < MAPQ: continue
        pos = int(f[3]) - 1; bs = pos
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "M=XD" and not (op == "D" and n >= GAPMIN): pos += n
            elif op == "N" or op == "D":
                blocks.append((bs, pos)); gaps.append((pos, n)); pos += n; bs = pos
        blocks.append((bs, pos))
    p.wait()
blocks = np.array(blocks, np.int64); gaps = np.array(gaps, np.int64).reshape(-1, 2)
print("blocks", len(blocks), "gaps>=15", len(gaps), flush=True)
order = np.argsort(gaps[:, 0]); gs = gaps[order, 0]; ge_s = (gaps[:, 0] + gaps[:, 1])[order]
S = W + 2 * PAD
d = np.zeros(N + 1, np.int64); ln = blocks[:, 1] - blocks[:, 0]; b = blocks[ln >= S]
np.add.at(d, b[:, 0], 1); np.add.at(d, b[:, 1] - S + 1, -1); cov = np.cumsum(d)[:N]
def joined(a):
    lo = np.searchsorted(gs, a - FLANK); hi = np.searchsorted(gs, a + W + FLANK)
    return int(np.sum(ge_s[lo:hi] <= a + W + FLANK))
# detector check windows: L02 (3 placements) and 859933, scored by the same code even if through < MINCOV
for name, a in [("L02", 4487605), ("L02+1", 4487606), ("L02+2", 4487607), ("S859", 859933), ("S859_W31mid", 859940)]:
    print(f"CHECK {name} a={a} dG_DNA={dg(G[a:a+W]):.2f} through={int(cov[a-PAD])} joined={joined(a)}", flush=True)
xs = np.nonzero(cov >= MINCOV)[0]; a_all = xs + PAD; a_all = a_all[a_all + W + FLANK < N]
rows = []
for a in a_all:
    s = G[a:a + W]
    if "N" in s: continue
    rows.append((a, dg(s), int(cov[a - PAD]), joined(a)))
r = np.array(rows, dtype=float); np.save(f"{out}_W31dna.npy", r)
taken = np.zeros(N, bool); loci = []
for k in np.argsort(r[:, 1]):
    a = int(r[k, 0])
    if taken[a]: continue
    taken[max(0, a - W):a + W] = True; loci.append(r[k])
L = np.array(loci); np.save(f"{out}_W31dna_loci.npy", L)
print("windows", len(r), "loci", len(L), flush=True)
