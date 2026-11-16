"""Context of CGG +1 skips pooled over a set: the base after CGG (N of CGGN) against Gly codon usage, and the codon
before the CGG (the P-site codon if the peptidyl-tRNA slips) against background codon usage.
usage: python -I src/analysis/skip_context.py <set.tsv> <outdir> [exclude_acc,...]
Background: codons in frame-continuous DIAMOND hits (n_fs == 0) of the same genomes."""
import csv
import gzip
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from candidate import read_fasta  # noqa: E402
from summarize import flag  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")
from paths import GENOMES as ATLAS  # noqa: E402


def rc(s):
    return s.translate(COMP)[::-1]


setf, outdir = sys.argv[1], Path(sys.argv[2])
excl = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else set()
bg, pg, pc, nb, nn = Counter(), Counter(), Counter(), Counter(), 0
for r in csv.DictReader(open(setf), delimiter="\t"):
    g = r["gca"]
    if g in excl or not (outdir / f"{g}.fs.tsv.gz").exists():
        continue
    fs = list(csv.DictReader(gzip.open(outdir / f"{g}.fs.tsv.gz", "rt"), delimiter="\t"))
    flag(fs)
    sites = [x for x in fs if x["clean"] and x["type"] == "+1" and x["old_codon"] == "CGG"]
    if not sites:
        continue
    seqs = read_fasta(str(ATLAS / g[:8] / f"{g}.fna.gz"))
    for h in csv.DictReader(gzip.open(outdir / f"{g}.hits.tsv.gz", "rt"), delimiter="\t"):
        if h["n_fs"] != "0":
            continue
        s = seqs[h["contig"]][int(h["gstart"]) - 1:int(h["gend"])].upper()
        s = rc(s) if h["strand"] == "-" else s
        cods = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
        bg.update(cods)
        pg.update(a for a, b in zip(cods, cods[1:]) if b[:2] == "GG")
    for x in sites:
        c = x["ctx"]
        pc[c[17:20].upper()] += 1
        nb[c[23].upper()] += 1
        nn += 1
tot = sum(v for k, v in bg.items() if set(k) <= set("ACGT"))
gly = {b: bg["GG" + b] for b in "ACGT"}
gt = sum(gly.values())
print(f"sites {nn}; background codons {tot}; CGG in background {bg['CGG']} ({1e3 * bg['CGG'] / tot:.3f}/1000)")
print("N of CGGN:   " + "  ".join(f"{b} {nb[b]} ({100 * nb[b] / nn:.0f}%)" for b in "ACGT"))
print("GGN usage:   " + "  ".join(f"GG{b} {100 * gly[b] / gt:.0f}%" for b in "ACGT"))
print("codon before CGG, obs/exp (exp from background usage), top and bottom by log ratio with obs+exp >= 10:")
rows = []
for k in sorted(bg):
    if not set(k) <= set("ACGT"):
        continue
    e = nn * bg[k] / tot
    if pc[k] + e >= 10:
        rows.append((math.log2((pc[k] + 0.5) / (e + 0.5)), k, pc[k], e))
rows.sort()
for t in rows[::-1][:10]:
    print(f"   {t[1]} obs {t[2]} exp {t[3]:.1f} log2 {t[0]:+.2f}")
print("   ...")
for t in rows[:6]:
    print(f"   {t[1]} obs {t[2]} exp {t[3]:.1f} log2 {t[0]:+.2f}")
pt = sum(v for k, v in pg.items() if set(k) <= set("ACGT"))
print("codon before CGG against codons before in-frame GGN (Gly context), obs+exp >= 10:")
rows = []
for k in sorted(pg):
    if not set(k) <= set("ACGT"):
        continue
    e = nn * pg[k] / pt
    if pc[k] + e >= 10:
        rows.append((math.log2((pc[k] + 0.5) / (e + 0.5)), k, pc[k], e))
rows.sort()
for t in rows[::-1][:12]:
    print(f"   {t[1]} obs {t[2]} exp {t[3]:.1f} log2 {t[0]:+.2f}")
print("   ...")
for t in rows[:8]:
    print(f"   {t[1]} obs {t[2]} exp {t[3]:.1f} log2 {t[0]:+.2f}")
chi = sum((pc[k] - nn * bg[k] / tot) ** 2 / (nn * bg[k] / tot) for k in bg if set(k) <= set("ACGT") and bg[k])
print(f"chi2 (61-63 df) {chi:.0f}")
