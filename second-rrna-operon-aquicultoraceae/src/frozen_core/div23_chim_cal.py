import collections, random, subprocess, sys, tempfile, os
sys.path.insert(0, "src/frozen_core")
from div23_chim import win, run
random.seed(7)
F = "results/lsu/scan23/census_lsu_all.fna"
seqs, name = {}, None
for l in open(F):
    if l[0] == ">": name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l.strip())
seqs = {k: "".join(v) for k, v in seqs.items() if len("".join(v)) >= 2700}
tax = {}
for l in open("results/lsu/scan23/gtdb_tax.tsv"):
    f = l.rstrip("\n").split("\t")
    tax[f[0]] = f[-1].split(";")
fam = collections.defaultdict(lambda: collections.defaultdict(list))
for k in seqs:
    a = k.split("|")[0]
    a = a[3:] if a[:3] in ("GB_", "RS_") else a
    t = tax.get(a) or tax.get(k.split("|")[0])
    if t and len(t) >= 6: fam[t[4]][t[5]].append(k)
fams = [f for f in fam if len(fam[f]) >= 2]
res = []
with tempfile.TemporaryDirectory() as td:
    while len(res) < 300:
        f = random.choice(fams); g1, g2 = random.sample(list(fam[f]), 2)
        a, b = random.choice(fam[f][g1]), random.choice(fam[f][g2])
        open(td + "/q", "w").write(f">q\n{seqs[a]}\n"); open(td + "/s", "w").write(f">s\n{seqs[b]}\n")
        out = subprocess.run(["blastn", "-task", "blastn", "-query", td + "/q", "-subject", td + "/s", "-outfmt",
                              "6 pident length nident qseq sseq", "-evalue", "1e-20"], capture_output=True, text=True).stdout
        hs = [l.split("\t") for l in out.splitlines() if int(l.split("\t")[1]) >= 300]
        if not hs: continue
        ident = sum(int(h[2]) for h in hs) / sum(int(h[1]) for h in hs)
        ws = [win(h[3], h[4].strip()) for h in hs]
        res.append((ident, max(run(w) for w in ws), min(x for w in ws for x in w) if any(ws) else 100))
for lo, hi in ((0.80, 0.85), (0.85, 0.88), (0.88, 0.90), (0.90, 0.92), (0.92, 0.95), (0.95, 1.0)):
    b = [r for r in res if lo <= r[0] < hi]
    if b:
        fp = sum(1 for r in b if r[1] >= 600 and r[2] <= 90)
        print(f"{lo:.2f}-{hi:.2f}\tn={len(b)}\tchimera-call={fp}\trun>=600={sum(1 for r in b if r[1]>=600)}\tmax_run={max(r[1] for r in b)}")
