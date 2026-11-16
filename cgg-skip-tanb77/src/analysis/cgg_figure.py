"""Figure and its table: skipped CGG (share of clean +1 breaks) and in-frame CGG (per 1,000 codons), by lineage and
CGG-reader status (codon-hunt's ARAGORN screen: no tRNA with anticodon CCG or TCG = no reader; genomes missing from the
screen are left out).
usage: python -I src/analysis/cgg_figure.py <out.tsv> <out.png>
inputs: data/frozen_test/cgg_test_{T,cag698r,moll,C50}.tsv (src/analysis/cgg_test.py) and data/frozen_test/inframe_{tanb77,cag698,moll,clos}.tsv
(src/analysis/inframe_cgg.py)"""
import csv
import statistics as st
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[2] / "data" / "frozen_test"  # cgg_test_*.tsv and inframe_*.tsv
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import TRNA_SCREEN as SCREEN  # noqa: E402
GROUPS = [("TANB77", "T", "tanb77"), ("CAG-698", "cag698r", "cag698"), ("Mollicutes", "moll", "moll"),
          ("Other low-GC\nClostridia", "C50", "clos")]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
COL = {"no reader": "#eb6834", "reader": "#2a78d6"}


def main():
    out_tsv, out_png = sys.argv[1], sys.argv[2]
    scr = {r["acc"]: {x.split(":")[0] for x in r["anticodons"].split(",") if x}
           for r in csv.DictReader(open(SCREEN), delimiter="\t")}
    rows = []
    for name, ct, inf in GROUPS:
        brk = {r["gca"]: r for r in csv.DictReader(open(HERE / f"cgg_test_{ct}.tsv"), delimiter="\t")}
        agg = {}
        for r in csv.DictReader(open(HERE / f"inframe_{inf}.tsv"), delimiter="\t"):
            ac = scr.get(r["gca"])
            if ac is None or r["gca"] not in brk:
                continue
            k = "reader" if ac & {"CCG", "TCG"} else "no reader"
            a = agg.setdefault(k, {"n": 0, "cgg": 0, "p1": 0, "per_k": []})
            a["n"] += 1
            a["cgg"] += int(brk[r["gca"]]["cgg"])
            a["p1"] += int(brk[r["gca"]]["p1"])
            a["per_k"].append(float(r["CGG_per_k"]))
        for k in ("no reader", "reader"):
            a = agg[k]
            rows.append([name.replace("\n", " "), k, a["n"], a["cgg"], a["p1"], a["cgg"] / a["p1"] if a["p1"] else 0,
                         st.median(a["per_k"])])
    w = csv.writer(open(out_tsv, "w"), delimiter="\t", lineterminator="\n")
    w.writerow(["lineage", "reader", "genomes", "cgg_breaks", "clean_plus1_breaks", "cgg_share", "inframe_cgg_per_k_median"])
    w.writerows([r[:5] + [f"{r[5]:.3f}", f"{r[6]:.3f}"] for r in rows])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=True, gridspec_kw={"wspace": 0.08})
    fig.patch.set_facecolor(SURF)
    ys, labels = [], []
    for i, (name, _, _) in enumerate(GROUPS):
        for j, k in enumerate(("no reader", "reader")):
            ys.append(i * 2.6 + j * 1.0)
            labels.append((name, k))
    for ax in (a1, a2):
        ax.set_facecolor(SURF)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
    for y, r in zip(ys, rows):
        c = COL[r[1]]
        a1.barh(y, 100 * r[5], height=0.7, color=c, edgecolor=SURF, linewidth=2)
        a1.text(100 * r[5] + 1.5, y, f"{r[3]:,} of {r[4]:,}", va="center", fontsize=8.5, color=INK2)
        a2.plot(r[6], y, "o", ms=8, color=c, markeredgecolor=SURF, markeredgewidth=2)
        a2.text(r[6] * 1.25, y, f"{r[6]:.3f}  (n {r[2]})", va="center", fontsize=8.5, color=INK2)
    a1.set_yticks(ys)
    a1.set_yticklabels([f"{n}  ·  {k}" if k == "no reader" else k for n, k in labels], fontsize=9)
    for t, (n, k) in zip(a1.get_yticklabels(), labels):
        t.set_color(INK if k == "no reader" else INK2)
    a1.invert_yaxis()
    a1.set_xlim(0, 100)
    a1.set_xlabel("CGG share of clean +1 frame breaks (%)")
    a1.set_title("Skipped CGG", loc="left", fontsize=11, fontweight="bold", color=INK)
    a2.set_xscale("log")
    a2.set_xlim(0.006, 12)
    a2.set_xticks([0.01, 0.1, 1, 10])
    a2.set_xticklabels(["0.01", "0.1", "1", "10"])
    a2.minorticks_off()
    a2.set_xlabel("In-frame CGG per 1,000 codons (median genome, log scale)")
    a2.set_title("CGG read as a codon", loc="left", fontsize=11, fontweight="bold", color=INK)
    fig.text(0.01, -0.06, "Reader: a tRNA with anticodon CCG or UCG (ARAGORN). Orange: no reader; blue: reader. "
             "Breaks pooled over genomes; n = genomes.", fontsize=8, color=INK2)
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor=SURF)


if __name__ == "__main__":
    main()
