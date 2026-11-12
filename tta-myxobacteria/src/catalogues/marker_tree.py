"""Quick placement check: concatenated core-marker ML tree (mafft + iqtree3 LG+G4 -fast, SH-aLRT) of the 10 clade
genomes, candidates and outgroup panel genomes; reports whether each candidate falls inside the clade split.
usage: uv run --no-project --with treeswift python -I src/catalogues/marker_tree.py WORKDIR OUTGROUPS_CSV cand.faa...
env: MAX_MARKERS (default 150), THREADS for iqtree3 (default 1). The tree is rooted on the last outgroup.
"""
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen

ROOT = Path(__file__).resolve().parents[2]
MDB = ROOT / "data/catalogues/markerdb"
work = Path(sys.argv[1]); work.mkdir(parents=True, exist_ok=True)
outg = sys.argv[2].split(",")
cands = sys.argv[3:]
seqs = {}
for h, s in screen.read_fa(MDB / "markers.faa"):
    g, _, m = h.split("|")
    seqs[(g, m)] = s
taxa = list(screen.CLADE) + outg
cn = [screen.gname(c) for c in cands]
q = work / "cands.faa"
allp = {}
with open(q, "w") as o:
    for c in cands:
        for h, s in screen.read_fa(c):
            o.write(f">{h}\n{s}\n"); allp[h] = s
best = {}
for qq, s, pid, bs in screen.diamond(q, MDB / "markers", 25, threads=1):
    g, m = qq.split("|")[0], s.split("|")[2]
    if (g, m) not in best or bs > best[(g, m)][0]:
        best[(g, m)] = (bs, qq)
for (g, m), (bs, qq) in best.items():
    seqs[(g, m)] = allp[qq]
taxa += cn
ms = sorted({m for (_, m) in seqs})
pres = {m: sum((t, m) in seqs for t in taxa) for m in ms}
use = sorted([m for m in ms if pres[m] >= 0.9 * len(taxa)], key=lambda m: -pres[m])[:int(os.environ.get("MAX_MARKERS", 150))]
print("taxa", len(taxa), "markers used", len(use), flush=True)
aln = defaultdict(str)
for m in use:
    f = work / f"m{m}.faa"
    with open(f, "w") as o:
        for t in taxa:
            if (t, m) in seqs:
                o.write(f">{t}\n{seqs[(t, m)]}\n")
    r = subprocess.run(["/opt/homebrew/bin/mafft", "--auto", "--thread", "1", "--quiet", str(f)], capture_output=True, text=True)
    a = dict(screen.read_fa_text(r.stdout)) if hasattr(screen, "read_fa_text") else None
    if a is None:
        a, name = {}, None
        for line in r.stdout.splitlines():
            if line.startswith(">"):
                name = line[1:].strip(); a[name] = ""
            else:
                a[name] += line.strip()
    L = len(next(iter(a.values())))
    for t in taxa:
        aln[t] += a.get(t, "-" * L)
cf = work / "concat.faa"
with open(cf, "w") as o:
    for t in taxa:
        o.write(f">{t}\n{aln[t]}\n")
print("alignment length", len(aln[taxa[0]]), flush=True)
subprocess.run(["/opt/homebrew/bin/iqtree3", "-s", str(cf), "-m", "LG+G4", "-T", os.environ.get("THREADS", "1"), "-fast", "--alrt", "1000",
                "--prefix", str(work / "tree"), "-redo", "-quiet"], check=True)
import treeswift
t = treeswift.read_tree_newick(str(work / "tree.treefile"))
lab = {l.label: l for l in t.traverse_leaves()}
t.reroot(lab[outg[-1]])  # root on a far outgroup (all panel outgroups are outside the clade in R232)
t = treeswift.read_tree_newick(t.newick())
m = t.mrca(set(screen.CLADE))
inside = {l.label for l in m.traverse_leaves()}
print("clade MRCA support (SH-aLRT):", m.label, "| leaves under MRCA:", len(inside))
print("outgroups inside clade MRCA (should be none):", sorted(inside & set(outg)))
for c in cn:
    print(c, "INSIDE" if c in inside else "OUTSIDE")
# sister/basal check: is there a split containing exactly the 10 refs (+ inside candidates) without c?
print(t.newick())
