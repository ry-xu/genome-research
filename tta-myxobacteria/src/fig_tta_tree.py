"""Outbox figure: the 16-genome Polyangiaceae panel on the GTDB R232 tree, with TTA gene ends and the TTA tRNA.

Counts are the adversarial review's close-relative test (src/checks/redteam/ortho/end_codon_summary.tsv, rows with
sister genomes as the reference): genes that run on past the point where a sister ortholog stops at a canonical
stop, and how many of them have TTA as the next codon. tRNA-Leu(UAA) from src/checks/machinery/trna_panel.tsv.

usage: uv run --no-project --with biopython,matplotlib,pyrodigal,numpy,requests python -I src/fig_tta_tree.py <out.png>
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from blind_panel import family_subtree  # noqa: E402

INK, ORANGE, TREE, MUTED = "#1c1c1a", "#e8650a", "#b9b8b2", "#8a8a84"
plt.rcParams["font.family"] = ["Helvetica", "Arial", "DejaVu Sans"]


def main():
    key = {r["acc"]: r for r in csv.DictReader(open(ROOT / "results/blind_panel_key.tsv"), delimiter="\t")}
    ends = {}
    for r in csv.DictReader(open(ROOT / "src/checks/redteam/ortho/end_codon_summary.tsv"), delimiter="\t"):
        if r["B"] != "B" and r["R"] == "sister":
            ends[r["B"].split("_", 1)[1]] = (int(r["TTA_p1"]), int(r["runon"]))
    trna = {r["acc"]: int(r["Leu_TAA"]) for r in csv.DictReader(open(ROOT / "src/checks/machinery/trna_panel.tsv"),
                                                                 delimiter="\t")}
    t = family_subtree()
    for tip in t.get_terminals():
        if tip.name not in key:
            t.prune(tip)
    t.ladderize(reverse=True)
    tips = t.get_terminals()
    assert len(tips) == 16
    clade_root = t.common_ancestor([tp for tp in tips if key[tp.name]["role"] == "clade"])
    in_clade = {id(c) for c in clade_root.find_clades()}

    y = {id(tp): i for i, tp in enumerate(tips)}
    height = {}

    def place(c):
        if c.is_terminal():
            height[id(c)] = 0
            return
        for ch in c.clades:
            place(ch)
        height[id(c)] = 1 + max(height[id(ch)] for ch in c.clades)
        y[id(c)] = (y[id(c.clades[0])] + y[id(c.clades[-1])]) / 2

    place(t.root)
    H = height[id(t.root)]
    x = {k: (H - h) for k, h in height.items()}

    fig = plt.figure(figsize=(12, 8.4), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 12)
    ax.set_ylim(8.4, 0)
    ax.axis("off")
    x0, x1 = 0.5, 3.4  # tree span
    row0, dy = 1.75, 0.36

    def X(c):
        return x0 + (x1 - x0) * x[id(c)] / H

    def Y(c):
        return row0 + dy * y[id(c)]

    def style(c):
        return (ORANGE, 2.4) if id(c) in in_clade else (TREE, 1.6)

    for c in t.find_clades():
        for ch in c.clades:
            col, lw = style(c)
            ax.plot([X(c), X(c)], [Y(c), Y(ch)], color=col, lw=lw, solid_capstyle="butt")
            col, lw = style(ch)
            ax.plot([X(c), X(ch)], [Y(ch), Y(ch)], color=col, lw=lw, solid_capstyle="butt")
    ax.plot([x0 - 0.15, x0], [Y(t.root)] * 2, color=TREE, lw=1.6)

    cx = [3.75, 6.6, 9.2]
    ax.text(0.45, 0.72, "Genes end at TTA in one clade, not its relatives", fontsize=27, fontweight="bold", color=INK, va="center")
    for xx, h in zip(cx, ["Genus", "Ends at TTA", "tRNA for TTA"]):
        ax.text(xx, 1.33, h, fontsize=15, color=MUTED, va="center")
    for tp in tips:
        k = key[tp.name]
        clade = k["role"] == "clade"
        n_tta, n = ends[tp.name.split("_")[2].split(".")[0]]
        yy = Y(tp)
        w = "bold" if clade else "normal"
        ax.text(cx[0], yy, k["genus"][3:], fontsize=13.5, color=INK if clade else MUTED, fontweight=w, va="center")
        ax.text(cx[1], yy, f"{n_tta} of {n}", fontsize=13.5, color=INK if clade else MUTED, fontweight=w, va="center")
        lost = trna[tp.name] == 0
        ax.text(cx[2], yy, "not found" if lost else "found", fontsize=13.5, color=ORANGE if lost else MUTED,
                fontweight="bold" if lost else "normal", va="center")
    ax.text(0.5, 8.0, "Ends at TTA: of the genes that run on past the point where a close relative's gene stops, "
            "how many have TTA as the next codon.   Tree: GTDB R232, topology only.", fontsize=10.5, color=MUTED, va="center")
    fig.savefig(sys.argv[1], dpi=300, facecolor="white")
    print("->", sys.argv[1])


if __name__ == "__main__":
    main()
