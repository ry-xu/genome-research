"""Codon-row figure for the mito-proline post (outbox/2026-11-11-mito-proline-codons).
usage: uv run --no-project --with matplotlib python -I fig_proline_rows.py OUT.png

Mitochondrial Cox2, profile mito_Fun_cox2 columns 100-108 (Clover's census profiles,
work/new-stops-euk/data/mito_profiles.hmm), aligned the census way (best hmmsearch domain of each gbmitoorf ORF,
match states, posterior '*' at all 27 codons). Each 27-nt stretch occurs once in its genome, forward strand, in
work/proline-codon-yeasts/data/genomes/all.fa; first base given below. Column 104 is in the census set of conserved
proline positions with CUU/CUA in all four (work/new-stops-euk/results/cun_pro/audit/*.census_ctn_p.txt,
'CTA mito_Fun_cox2:104'). Amino acids: mold mitochondrial code (table 4); CUA shown as P in the four.
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORANGE, INK, MUTED, PALE = "#E8650A", "#1a1a1a", "#8a8580", "#fbe3d1"
# (name, accession, first base of the 27 nt, DNA as in the genome)
ROWS = [("Wickerhamiella brachini", "BK074528.1", 6574, "TTAATAGCATTCCCATCATTTATATTA"),
        ("Wickerhamiella siamensis", "BK074538.1", 31423, "TTAATCGCATTCCTATCATTTATCTTA"),
        ("Wickerhamiella infanticola", "BK074532.1", 10766, "TTAATTGCATTCCTATCATTTATATTA"),
        ("Wickerhamiella verensis", "BK074542.1", 22393, "TTAATTGCATTCCTATCATTTATATTA")]
# The Sugiyamaella pair (BK074251.1 S. americana relative; BK074262.1 S. pinicola, CTA at 1086+12) was cut for the tweet.
SUB = {"BK074528.1": "relative", "BK074251.1": "relative"}
B = "TCAG"
AA4 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # NCBI table 4, TCAG order


def aa(c, reassigned):
    if reassigned and c in ("CTT", "CTA"):
        return "P"
    return AA4[16 * B.index(c[0]) + 4 * B.index(c[1]) + B.index(c[2])]


plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
ys, y = [], 0.0
for r in ROWS:
    if r is None:
        y -= 0.45
        continue
    ys.append((y, r)); y -= 1.0
top, bot = 0.55, y + 0.35
fig = plt.figure(figsize=(12, 0.62 * (top - bot) + 1.1), dpi=200); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.0, 0.17, 1.0, 0.80])
ax.set_xlim(-15.5, 36.5); ax.set_ylim(bot, top); ax.axis("off")
X = [4 * k for k in range(9)]  # codon k centred at x = 4k + 1
ax.add_patch(plt.Rectangle((X[4] - 0.6, bot + 0.05), 3.2, top - bot - 0.1, color=PALE, lw=0))
for y, (name, acc, pos, dna) in ys:
    changed = acc not in SUB
    for k in range(9):
        c = dna[3 * k:3 * k + 3]
        hot = changed and c in ("CTT", "CTA")
        col, w = (ORANGE, "bold") if hot else (INK, "normal")
        ax.text(X[k] + 1, y + 0.1, c.replace("T", "U"), ha="center", va="center", family="Menlo", fontsize=21,
                color=col, fontweight=w)
        ax.text(X[k] + 1, y - 0.3, aa(c, changed), ha="center", va="center", fontsize=12,
                color=ORANGE if hot else MUTED, fontweight=w)
    if acc in SUB:
        ax.text(-1.6, y + 0.12, name, ha="right", va="center", fontsize=15, style="italic", color=INK)
        ax.text(-1.6, y - 0.28, SUB[acc], ha="right", va="center", fontsize=11, color=MUTED)
    else:
        ax.text(-1.6, y - 0.05, name, ha="right", va="center", fontsize=15, style="italic", color=INK)
fig.text(0.075, 0.035, "Part of mitochondrial Cox2, aligned. Orange: CUA, normally leucine, where the relative has CCA (proline).\n"
         "In four yeast genomes, 286 of 290 CUU/CUA codons at strongly conserved positions sit at proline positions.",
         fontsize=10.5, color=MUTED, ha="left", va="bottom", linespacing=1.5)
fig.savefig(sys.argv[1], facecolor="white")
