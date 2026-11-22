"""DNA-vs-RNA bar figure (figures/div_rrna_bars.png): second vs normal 23S read pairs of Ca. Aquicultor secundus.
usage (from the folder root): uv run python -I src/fig_div_rrna_bars.py OUT.png
Read pairs are recounted from the PAFs in data/expression_23s/ by counts() of src/frozen_core/fig_div23.py
(identity >= 0.98, >= 60 nt, Crystal Geyser loci, per read pair): DNA SRR6051694 (28 May 2015); RNA SRR11558543,
SRR11558542, SRR11558425 (loci.paf). The two later runs are added as constants; their counts are in
data/silent_operon/rna/ (the PAFs are too large to include). In the lab this script loaded
work/frozen-core/src/fig_div23.py; only that path was changed."""
import importlib.util, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

spec = importlib.util.spec_from_file_location("fd", "src/frozen_core/fig_div23.py")
fd = importlib.util.module_from_spec(spec); spec.loader.exec_module(fd)
dna = fd.counts("SRR6051694.loci.paf")
# Two more Crystal Geyser RNA runs counted by the C_rna counter (data/silent_operon/rna/*.counts.txt) (same rule; reproduces the
# SRR11558543 count): SRR11558439 (CG05, 25 May 2015) 0 : 74 and SRR6241672 (CG13, not rRNA-depleted) 0 : 43.
NEW_DIV, NEW_NORM = 0, 74 + 43
rna = [fd.counts(p) for p in ("SRR11558543.loci.paf", "SRR11558542.loci.paf", "loci.paf")]
ROWS = [("DNA", "1 metagenome", dna["div"], dna["norm"]),
        ("RNA", "5 metatranscriptomes", sum(c["div"] for c in rna) + NEW_DIV, sum(c["norm"] for c in rna) + NEW_NORM)]
print(ROWS)

ORANGE, INK, MUTED, GREY = "#E8650A", "#1a1a1a", "#8a8580", "#c9c5c0"
plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
fig = plt.figure(figsize=(9.2, 2.35), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.20, 0.20, 0.62, 0.78])
ax.set_xlim(0, 1); ax.set_ylim(1.45, -0.85); ax.axis("off")
for i, (lab, sub, d, n) in enumerate(ROWS):
    f = d / (d + n)
    if d:
        ax.barh(i, f, color=ORANGE, height=0.62, lw=0)
    ax.barh(i, 1 - f, left=f + (0.004 if d else 0), color=GREY, height=0.62, lw=0)
    ax.text(-0.03, i - 0.1, lab, ha="right", va="center", fontsize=16, color=INK)
    ax.text(-0.03, i + 0.24, sub, ha="right", va="center", fontsize=11, color=MUTED)
    ax.text(1.03, i, f"{d} : {n}", ha="left", va="center", fontsize=15, color=INK)
# direct colour labels over the DNA bar, no legend
d0 = ROWS[0][2] / (ROWS[0][2] + ROWS[0][3])
ax.text(d0 / 2, -0.52, "second 23S", ha="center", va="bottom", fontsize=13, color=ORANGE, fontweight="bold")
ax.text(d0 + (1 - d0) / 2, -0.52, "normal 23S", ha="center", va="bottom", fontsize=13, color=MUTED)
fig.text(0.20, 0.03, "23S read pairs, Ca. Aquicultor secundus, Crystal Geyser, Utah.\n"
         "DNA 28 May 2015; RNA from 5 libraries, one not rRNA-depleted.", fontsize=10.5, color=MUTED, ha="left")
fig.savefig(sys.argv[1], facecolor="white")
