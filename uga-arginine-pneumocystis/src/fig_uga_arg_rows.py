"""Codon-row figure for the Pneumocystis UGA = arginine post (figures/uga_arg_rows.png).
usage: uv run python -I src/fig_uga_arg_rows.py OUT.png   (from the folder root)

Mitochondrial cox1 CDS codons 430-438 (0-based) of P. jirovecii NC_020331.1 and Pneumocystis sp. 'macacae'
MT726214.1 (host Macaca mulatta; the closest known relative of P. jirovecii, Cisse et al. 2021), translated with the
mold mitochondrial code (NCBI table 4): both read G L Q G M P R [R / UGA] I, the MPRR motif. Each 27-nt stretch
occurs once in its genome, forward strand, in data/genomes/fig_genomes.fa (checked below); first base given.
The UGA is P. jirovecii cox1 codon index 437, one of the 23 strict-arginine sites (data/sites/isolates.tsv,
JX499143.1 cox1 437; data/rna/pjir_*.tsv). Pairwise alignment to cow COX1 puts it at R439 of the MPRR pair
(R438/R439 bind the heme propionates); arginine in 99.5% of 6,188 genomes at this column (pre-post check).
The UGA's amino acid is shown as R, as the claim is "most likely".

Re-written for this folder from the description in data/checks/tweet_check.md, after the drafting agent's own copy
(adapted from the proline figure script) was not kept; it redraws the posted figure.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORANGE, INK, MUTED, PALE = "#E8650A", "#1a1a1a", "#8a8580", "#fbe3d1"
# (name, subtitle, accession, first base of the 27 nt, DNA as in the genome)
ROWS = [("Pneumocystis of macaques", "closest known relative", "MT726214.1", 5374, "GGACTTCAAGGTATGCCTAGAAGAATT"),
        ("Pneumocystis jirovecii", "human", "NC_020331.1", 17546, "GGATTACAGGGTATGCCCAGATGAATT")]
RELATIVE = {"MT726214.1"}
HOT = 7  # the codon shown highlighted (UGA in P. jirovecii)
B = "TCAG"
AA4 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # NCBI table 4, TCAG order


def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def aa(c, reassigned):
    if reassigned and c == "TGA":
        return "R"
    return AA4[16 * B.index(c[0]) + 4 * B.index(c[1]) + B.index(c[2])]


genomes = read_fasta(Path(__file__).resolve().parent.parent / "data" / "genomes" / "fig_genomes.fa")
for _, _, acc, pos, dna in ROWS:
    g = genomes[acc]
    assert g.count(dna) == 1 and g.find(dna) + 1 == pos, (acc, pos)

plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
ys, y = [], 0.0
for r in ROWS:
    ys.append((y, r)); y -= 1.0
top, bot = 0.424, y + 0.475  # the bracket and its label sit above the axes (not clipped)
fig = plt.figure(figsize=(12, 2.515), dpi=200); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.0, 0.274, 1.0, 0.566])
ax.set_xlim(-15.5, 36.5); ax.set_ylim(bot, top); ax.axis("off")
X = [4 * k for k in range(9)]  # codon k centred at x = 4k + 1
ax.add_patch(plt.Rectangle((X[HOT] - 0.6, bot + 0.05), 3.2, top - bot - 0.1, color=PALE, lw=0))
# bracket over the arginine pair (codons 6 and 7)
bx0, bx1, by = X[6] - 0.4, X[7] + 2.4, 0.51
ax.plot([bx0, bx0, bx1, bx1], [by - 0.1, by, by, by - 0.1], color=INK, lw=1.2, clip_on=False)
ax.text((bx0 + bx1) / 2, by + 0.28, "arginine pair that binds heme", ha="center", va="center", fontsize=15, color=INK)
for y, (name, sub, acc, pos, dna) in ys:
    changed = acc not in RELATIVE
    for k in range(9):
        c = dna[3 * k:3 * k + 3]
        hot = changed and c == "TGA"
        col, w = (ORANGE, "bold") if hot else (INK, "normal")
        ax.text(X[k] + 1, y + 0.1, c.replace("T", "U"), ha="center", va="center", family="Menlo", fontsize=21,
                color=col, fontweight=w)
        ax.text(X[k] + 1, y - 0.3, aa(c, changed), ha="center", va="center", fontsize=12,
                color=ORANGE if hot else MUTED, fontweight=w)
    genus, _, rest = name.partition(" ")
    label = f"$\\it{{{genus}}}$ {rest}" if rest == "of macaques" else f"$\\it{{{genus}}}$ $\\it{{{rest}}}$"
    ax.text(-1.6, y + 0.12, label, ha="right", va="center", fontsize=15, color=INK)
    ax.text(-1.6, y - 0.28, sub, ha="right", va="center", fontsize=11, color=MUTED)
fig.text(0.075, 0.035, "Part of mitochondrial Cox1, aligned. Orange: UGA, tryptophan in the code these fungi are filed "
         "under, where the relative has AGA (arginine).\n"
         "Arginine is at this position in 99.5% of 6,188 mitochondrial genomes.",
         fontsize=10.5, color=MUTED, ha="left", va="bottom", linespacing=1.5)
fig.savefig(sys.argv[1], facecolor="white")
