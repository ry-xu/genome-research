"""Enolase gene barcode: Cafeteria (CAG) vs Bilabrum (TAG) at human ENOA Q166 (32 of 32 1KP RNA reads TAG).
usage (from the folder root): uv run python -I src/fig_bilabrum_enoa.py OUT.png ;
Needs the two genomes in data/dl/ (not included; see src/README.md); codons from the reviewer's tblastn (standard code) HSPs of human ENOA (P06733); one locus per genome (contig, start, end)."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "review"))
from stops_in_hsps import read_fasta, walk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

R = "data/review"
SRC = [("caf", "data/dl/GCA_008330645.1_CrBVI_genomic.fna", f"{R}/caf_human76.tsv", ("VLTN01000067.1", 93900, 95700)),
       ("bil", "data/dl/GCA_965178125.1_puBilCfla1.1_genomic.fna",
        f"{R}/bil_human76.tsv", ("OZ237754.1", 81600, 83000))]
SITE = 166  # human ENOA Gln; Bilabrum TAG, Cafeteria CAG


def codons(genome, hits, locus):
    seqs = read_fasta(genome)
    out = []  # (x, codon, qpos)
    for line in open(hits):
        f = line.rstrip("\n").split("\t")
        s = f[1].split("|")[1]
        lo_, hi_ = sorted((int(f[8]), int(f[9])))
        if "ENOA_HUMAN" not in f[0] or s != locus[0] or lo_ < locus[1] or hi_ > locus[2]:
            continue
        h = dict(qstart=int(f[6]), sstart=int(f[8]), sframe=int(f[13]), qseq=f[15], sseq=f[16])
        last, k = int(f[6]) - 1, 0
        for col, qp, a, b, cod, p in walk(h, seqs[s]):
            if qp is None:
                k += 1
                out.append((last + 0.5 + 0.1 * k, cod, None))
            else:
                last, k = qp, 0
                out.append((qp, cod, qp))
    seen, res = set(), []
    for x, c, q in sorted(out):
        if q is not None and q in seen:
            continue
        seen.add(q)
        res.append((x, c, q))
    return res


TINT = ["#FDF1E8", "#FAE0CE", "#F7CFB2", "#F4C19C"]
ORANGE = "#E5650E"
INK, MUTED, RULE = "#1F1F1F", "#8A8A85", "#9A9A95"


def tint(c):
    return TINT[sum("ACGT".index(b) for b in c) % 4] if set(c) <= set("ACGT") else "#FFFFFF"


data = {k: codons(g, h, c) for k, g, h, c in SRC}
lo = min(x for d in data.values() for x, _, _ in d)
hi = max(x for d in data.values() for x, _, _ in d) + 1
W, H = 12, 6.75
fig = plt.figure(figsize=(W, H), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.axis("off")
plt.rcParams["font.family"] = "Helvetica"
BX0, BX1 = 0.9, 11.1
ax.text(0.6, 0.62, "A stop codon in the middle of an essential gene (enolase)", fontsize=27, fontweight="bold",
        color=INK, va="center")
rows = [("caf", ("Cafeteria", ", a relative"), "standard code: TAG = stop", 1.35),
        ("bil", ("Bilabrum", ""), "reads TAG as glutamine", 3.95)]
for key, name, sub, y in rows:
    d = data[key]
    t = ax.text(BX0, y, name[0], fontsize=21, fontweight="bold", color=INK, va="center", style="italic")
    if name[1]:
        ax.annotate(name[1], xycoords=t, xy=(1, 0.5), va="center", fontsize=21, fontweight="bold", color=INK)
    ax.text(BX0, y + 0.47, sub, fontsize=16, color=MUTED, va="center")
    # zoom: 7 codons centred on the site
    i = next(j for j, (x, c, q) in enumerate(d) if q == SITE)
    win = d[i - 3:i + 4]
    bw, gap, zx0, zy = 0.88, 0.1, 3.95, y + 0.42
    for j, (x, c, q) in enumerate(win):
        xx = zx0 + j * (bw + gap)
        site = q == SITE
        hot = c in ("TAG", "TAA", "TGA")
        ax.add_patch(Rectangle((xx, zy), bw, 0.6, facecolor=ORANGE if hot else tint(c),
                               edgecolor=INK if site else "none", lw=1.6 if site else 0))
        ax.text(xx + bw / 2, zy + 0.31, c, fontsize=23, family="Menlo", ha="center", va="center",
                color="white" if hot else INK, fontweight="bold" if site else "normal")
    by = y + 1.33
    bh = 0.55
    for x, c, q in d:
        xx = BX0 + (x - lo) / (hi - lo) * (BX1 - BX0)
        w = (BX1 - BX0) / (hi - lo)
        hot = c in ("TAG", "TAA", "TGA")
        ax.add_patch(Rectangle((xx, by), w * (1.0 if hot else 1.0), bh, facecolor=ORANGE if hot else tint(c),
                               edgecolor="none", lw=0))
        if hot:
            ax.add_patch(Rectangle((xx - w, by), 3 * w, bh, facecolor=ORANGE, edgecolor="none"))
    ax.add_patch(Rectangle((BX0, by), BX1 - BX0, bh, facecolor="none", edgecolor="#DDDDDA", lw=0.8))
    sx = BX0 + (SITE - lo + 0.5) / (hi - lo) * (BX1 - BX0)
    ax.add_patch(Rectangle((sx - 0.08, by - 0.06), 0.16, bh + 0.12, facecolor="none", edgecolor="#6A6A66", lw=1.0))
    zr = zx0 + 7 * (bw + gap) - gap
    ax.plot([zx0 - 0.05, sx - 0.08], [zy + 0.66, by - 0.06], color=RULE, lw=0.9)
    ax.plot([zr + 0.05, sx + 0.08], [zy + 0.66, by - 0.06], color=RULE, lw=0.9)
ax.text(BX1, 6.27, "Enolase is a core enzyme of glycolysis. Boxed: glutamine (Q) in human enolase; all 32 RNA reads from a related Bilabrum strain carry TAG there.",
        fontsize=11, color=MUTED, ha="right", va="center")
ax.text(BX1, 6.52, "The other orange bar in Bilabrum is a second TAG, opposite lysine (K) in human enolase; all 36 RNA reads there carry TAG too.",
        fontsize=11, color=MUTED, ha="right", va="center")
fig.savefig(sys.argv[1], dpi=300, facecolor="white")
for k, d in data.items():
    print(k, len(d), "codons; stops:", [(round(x, 1), c, q) for x, c, q in d if c in ("TAG", "TAA", "TGA")])
