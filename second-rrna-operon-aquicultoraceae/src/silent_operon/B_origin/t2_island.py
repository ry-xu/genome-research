"""T2 island signals: GC of the +-10 kb flanks (rRNA genes excluded) against the genome's 20-kb windows (z), GC of the
rRNA genes themselves, and tRNA genes (aragorn) within the operon and +-10 kb.
Usage: python -I src/t2_island.py
"""
import csv, subprocess, statistics as st, tempfile, os
from Bio import SeqIO

ARAG = "~/dev/tools/bin/aragorn1.2.41"
units = [r for r in csv.DictReader(open("out/units.tsv"), delimiter="\t") if r["family"] == "Aquicultoraceae"]
loci = {r["id"]: r for r in csv.DictReader(open("out/loci.tsv"), delimiter="\t")}
gen = SeqIO.index("data/aq65.fna", "fasta")
gc = lambda s: (s.count("G") + s.count("C")) / max(1, len(s) - s.count("N"))
win = {}
for acc in {u["acc"] for u in units}:
    ws = []
    for k in gen:
        if k.startswith(acc + "|"):
            s = str(gen[k].seq).upper()
            for i in range(0, len(s) - 20000 + 1, 20000):
                ws.append(gc(s[i:i + 20000]))
    win[acc] = ws
print("unit       genus         label  flankGC  genome(mean,sd,n)   z     rRNA genes GC                  tRNAs (operon / flanks)")
for u in sorted(units, key=lambda x: (x["genus"], x["label"])):
    s = str(gen[u["contig"]].seq).upper()
    a, b = int(u["start"]), int(u["end"])
    fl = s[max(0, a - 1 - 10000):a - 1] + s[b:b + 10000]
    ws = win[u["acc"]]
    mu, sd = (st.mean(ws), st.stdev(ws)) if len(ws) > 2 else (float("nan"), float("nan"))
    z = (gc(fl) - mu) / sd if sd == sd and len(fl) > 2000 else float("nan")
    rg = []
    for m in u["members"].split(","):
        L = loci[m]; rg.append(f"{L['kind']}:{gc(s[int(L['start'])-1:int(L['end'])]):.3f}")
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">r\n{s[max(0, a - 1 - 10000):b + 10000]}\n"); tf = fh.name
    out = subprocess.run([ARAG, "-t", "-gcbact", "-w", tf], capture_output=True, text=True).stdout
    os.remove(tf)
    off = max(0, a - 1 - 10000)
    tin, tfl = [], []
    for line in out.splitlines():
        p = line.split()
        if len(p) >= 5 and p[1].startswith("tRNA-"):
            c = p[2].strip("c[]").split(",")
            lo = off + int(c[0]); hi = off + int(c[1])
            (tin if (lo >= a and hi <= b) else tfl).append(p[1])
    print(f"{u['unit']:10s} {u['genus']:13s} {u['label']:5s}  {gc(fl):.3f}   {mu:.3f},{sd:.3f},{len(ws):3d}   {z:5.1f}  {' '.join(rg):30s} {','.join(tin) or '-'} / {','.join(tfl) or '-'}")
