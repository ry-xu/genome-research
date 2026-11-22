"""Figure for the write-up: Shikimatogenerans tree (rooted on the Bostrichidae symbiont, cladogram, ultrafast bootstrap >= 70 shown)
beside each genome's TGA share of Trp codons (TGA / (TGA + TGG), data/symbionts/symbiont_profile.tsv) and its RF2 state.
usage: uv run python -I src/throscid_figure.py [treefile] [out.png]   (from the folder root; drew figures/trixagus_tree_rf2.png)"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "data" / "throscidae" / "concat.faa.treefile"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "figures" / "trixagus_tree_rf2.png"
# label -> accession, host, country, RF2 state (README "RF2 is lost in each reader"; hosts from BioSample)
# "front gone": no RF2 homology before Tcar residue 124; "TGA inside": in-frame TGA at codon 129 (README, 18:11)
TIPS = {"Tser_sericeus": ("GCA_040210545.2", "Trixagus sericeus", "USA", "deleted"),
        "Tder_dermestoides": ("GCA_040211365.2", "Trixagus dermestoides", "Germany", "deleted"),
        "Ttur_turgidus": ("GCA_040409785.2", "Trixagus turgidus", "Japan", "deleted"),
        "Tcar_carinifrons": ("GCA_040211755.2", "Trixagus carinifrons", "Belgium", "SPF lost"),
        "Tduv_duvalii": ("GCA_040210955.2", "Trixagus duvalii", "Germany", "SPF lost, front gone"),
        "Tmey": ("GCA_037094695.1", "Trixagus meybohmi", "Belgium", "SPF lost, TGA inside"),
        "Aulonothroscus_AUS03": ("GCA_040212085.2", "Aulonothroscus sp. AUS03", "Australia", "SPF lost"),
        "S_bostrichidophilus": ("GCA_029851105.1", "Dinoderus porcellus", "Germany", "SPF intact")}
READER, OTHER = "#eb6834", "#2a78d6"  # reference palette slots 2 and 1 (light)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def main():
    share = {}
    for r in csv.DictReader(open(HERE / "data" / "symbionts" / "symbiont_profile.tsv"), delimiter="\t"):
        tga, tgg = float(r["TGA_per_k"]), float(r["TGG_per_k"])
        share[r["acc"]] = 100 * tga / (tga + tgg)
    tree = Phylo.read(TREE, "newick")
    tree.root_with_outgroup("S_bostrichidophilus")
    tree.ladderize(reverse=True)
    tips = tree.get_terminals()
    y = {t: i for i, t in enumerate(tips)}
    depth = {}

    def place(clade, d):  # cladogram: x = depth in nodes; tips aligned at max depth afterwards
        depth[clade] = d
        for c in clade.clades:
            place(c, d + 1)
    place(tree.root, 0)
    maxd = max(depth[t] for t in tips)

    def ypos(c):
        if c.is_terminal():
            return y[c]
        return sum(ypos(k) for k in c.clades) / len(c.clades)

    fig = plt.figure(figsize=(10, 4.4), facecolor=SURF)
    ax = fig.add_axes([0.01, 0.1, 0.20, 0.78])
    bx = fig.add_axes([0.53, 0.1, 0.27, 0.78])
    for a in (ax, bx):
        a.set_facecolor(SURF)
    for c in tree.find_clades(order="preorder"):
        x0, y0 = (maxd if c.is_terminal() else depth[c]), ypos(c)
        for k in c.clades:
            xk = maxd if k.is_terminal() else depth[k]
            ax.plot([x0, x0], [y0, ypos(k)], color=INK2, lw=1.4, solid_capstyle="round")
            ax.plot([x0, xk], [ypos(k), ypos(k)], color=INK2, lw=1.4, solid_capstyle="round")
        if not c.is_terminal() and c is not tree.root and c.confidence is not None and c.confidence >= 70:
            ax.text(depth[c] - 0.08, y0 - 0.12, f"{c.confidence:.0f}", ha="right", va="bottom", fontsize=8, color=INK2)
    ax.set_xlim(-0.2, maxd + 0.1)
    ax.set_ylim(len(tips) - 0.5, -0.5)
    ax.axis("off")
    for t in tips:
        acc, host, country, rf2 = TIPS[t.name]
        reader = rf2 == "deleted"
        fig.text(0.215, 0.1 + 0.78 * (len(tips) - 0.5 - y[t]) / len(tips), f"{host}  ({country})", fontsize=9.5,
                 color=INK, va="center", style="italic",
                 weight="bold" if reader else "normal")
        v = share[acc]
        bx.barh(y[t], v, height=0.55, color=READER if reader else OTHER, edgecolor=SURF, linewidth=0)
        bx.text(v + 1.5, y[t], f"{v:.0f}%", va="center", fontsize=9, color=INK)
        fig.text(0.815, 0.1 + 0.78 * (len(tips) - 0.5 - y[t]) / len(tips), rf2, fontsize=9,
                 color=INK if reader else INK2, va="center", weight="bold" if reader else "normal")
    bx.set_ylim(len(tips) - 0.5, -0.5)
    bx.set_xlim(0, 100)
    bx.set_yticks([])
    bx.set_xticks([0, 25, 50, 75, 100])
    bx.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=8.5, color=INK2)
    bx.grid(axis="x", color=GRID, lw=0.8)
    bx.set_axisbelow(True)
    for s in ("top", "right", "left"):
        bx.spines[s].set_visible(False)
    bx.spines["bottom"].set_color(GRID)
    bx.tick_params(axis="x", length=0)
    fig.text(0.53, 0.93, "TGA share of Trp codons", fontsize=10, color=INK, weight="bold")
    fig.text(0.815, 0.93, "RF2 gene", fontsize=10, color=INK, weight="bold")
    fig.text(0.01, 0.93, "Ca. Shikimatogenerans, by host beetle", fontsize=10, color=INK, weight="bold")
    fig.text(0.01, 0.02, "Tree: 99 proteins, IQ-TREE LG+F+G4, rooted on the Dinoderus (Bostrichidae) symbiont; not to scale; "
             "ultrafast bootstrap >= 70 shown. TGA share = TGA / (TGA + TGG).", fontsize=7.5, color=INK2)
    fig.savefig(OUT, dpi=200, facecolor=SURF)
    print(OUT)


if __name__ == "__main__":
    main()
