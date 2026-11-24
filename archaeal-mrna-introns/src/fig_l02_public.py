"""Build the README figure (figures/asgard_b35_l02.png) from data/asgard_b35/.

Usage (from the folder root):
  uv run python -I src/fig_l02_public.py
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "asgard_b35"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "figures" / "asgard_b35_l02.png"

P, D = 4487605, 31  # 0-based start and length on CP104013.1
GENES = [("Gene 1 (UYP47557)", 4487083, 4487611), ("Gene 2 (UYP47558)", 4487702, 4488440)]
LO, HI = P - 260, P + D + 260

# Exact sequences from experimental-tests/asgard-b35-l02/sequences.txt (plus strand = mRNA sense)
UP_DNA = "AATAACATGATC"                         # 12 nt
INS_DNA = "AATTAAGAATAGGATTTCCCTATCTTAGATT"     # 31 nt; positions 3..5 are TAA stop codon
DN_DNA = "AAACCGTGTACT"                         # 12 nt
UP_RNA = UP_DNA.replace("T", "U")
DN_RNA = DN_DNA.replace("T", "U")


def load_rna():
    joined, through = [], []
    with open(DATA / "rna_mol_templates.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["site"] != "L02" or r["place"] != "0":
                continue
            t = (int(r["start"]), int(r["end"]))
            if r["class"] == "joined":
                joined.append(t + (int(r["gap_start"]), int(r["gap_len"])))
            else:
                through.append(t)
    return joined, through


def load_dna():
    joined, through = [], []
    with open(DATA / "fig_l02_dna_templates.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            t = (int(r["start"]), int(r["end"]))
            if r["class"] == "joined":
                joined.append(t)
            else:
                through.append(t)
    return joined, through


def main():
    rj, rt = load_rna()
    dj, dt = load_dna()
    assert (len(rj), len(rt), len(dj), len(dt)) == (385, 26, 0, 46), (len(rj), len(rt), len(dj), len(dt))
    anyread = next(
        int(r.split("\t")[4])
        for r in open(DATA / "rna_mol_gaps.tsv")
        if r.split("\t")[:3] == ["CP104013.1", str(P), str(D)]
    )
    assert anyread == 383, anyread

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon

    ORANGE = "#E8650A"
    ORANGE_DARK = "#B94800"
    STOP_RED = "#C8281E"
    GREY = "#ABA9A2"
    INK = "#141414"
    MUTED = "#686660"
    BAND = "#FDE7D7"
    BAND_LIGHT = "#FEF2E8"
    RULE = "#E2DFD8"

    plt.rcParams.update({"font.family": "Helvetica", "font.size": 13})

    rows_r = sorted(rj, key=lambda t: (t[0], t[1])) + sorted(rt, key=lambda t: (t[0], t[1]))
    rows_d = sorted(dt, key=lambda t: (t[0], t[1]))

    fig = plt.figure(figsize=(15.2, 9.8), dpi=220)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045, 0.962,
        "Asgard archaeon B35: 31 bases present in DNA, missing from 94% of RNA",
        fontsize=22, fontweight="bold", color=INK, va="top",
    )

    # ---- Panel A: Compact 2-row DNA vs RNA sequence zoom ----
    ax_seq = fig.add_axes([0.045, 0.768, 0.91, 0.145])
    ax_seq.axis("off")

    L_up, L_ins, L_dn = len(UP_DNA), len(INS_DNA), len(DN_DNA)
    L_tot = L_up + L_ins + L_dn  # 55 bases

    ax_seq.set_xlim(-9.8, L_tot + 0.8)
    ax_seq.set_ylim(-0.92, 1.85)

    ax_seq.add_patch(Rectangle((L_up - 0.5, -0.78), L_ins, 2.52, color=BAND, lw=0, zorder=0))

    y_dna = 0.98
    y_rna = 0.08
    ax_seq.add_patch(
        Rectangle((L_up + 3 - 0.45, y_dna - 0.30), 2.9, 0.60,
                  facecolor="#FAD7D4", edgecolor=STOP_RED, lw=1.1, zorder=2)
    )
    ax_seq.text(L_up + 4, y_dna + 0.44, "stop", ha="center", va="bottom",
                fontsize=11, fontweight="bold", color=STOP_RED)

    ax_seq.text(-1.2, y_dna, "DNA (genome)", ha="right", va="center", fontsize=13.5, fontweight="bold", color=INK)
    ax_seq.text(-1.2, y_rna, "RNA (94% of reads)", ha="right", va="center", fontsize=13.5, fontweight="bold", color=ORANGE_DARK)

    full_dna = UP_DNA + INS_DNA + DN_DNA
    for i, ch in enumerate(full_dna):
        in_ins = L_up <= i < L_up + L_ins
        in_stop = L_up + 3 <= i < L_up + 6
        col = STOP_RED if in_stop else (ORANGE_DARK if in_ins else INK)
        wt = "bold" if (in_stop or in_ins) else "normal"
        ax_seq.text(i, y_dna, ch, ha="center", va="center", family="Menlo",
                    fontsize=11.5, color=col, fontweight=wt, zorder=3)

    for i, ch in enumerate(UP_RNA):
        ax_seq.text(i, y_rna, ch, ha="center", va="center", family="Menlo", fontsize=11.5, color=INK, zorder=3)
    for j, ch in enumerate(DN_RNA):
        i = L_up + L_ins + j
        ax_seq.text(i, y_rna, ch, ha="center", va="center", family="Menlo", fontsize=11.5, color=INK, zorder=3)

    ax_seq.plot([L_up - 0.1, L_up + L_ins - 0.9], [y_rna, y_rna],
                color=ORANGE, lw=1.6, ls=(0, (4, 3)), zorder=2)
    ax_seq.text(
        L_up + L_ins / 2 - 0.5, y_rna, "  31 bases missing  ",
        ha="center", va="center", fontsize=12, fontweight="bold", color=ORANGE_DARK,
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=ORANGE, lw=1.0), zorder=4,
    )

    ax_seq.text(
        L_up + L_ins / 2 - 0.5, -0.50,
        "skips a stop codon and joins the genes in frame (matches ~66 of 107 relatives)",
        ha="center", va="center", fontsize=11.8, color=ORANGE_DARK,
    )

    # ---- Panel B: Gene labels + RNA pileup + DNA pileup ----
    left_b, right_b = 0.045, 0.745
    gs = fig.add_gridspec(
        3, 1,
        height_ratios=[32, len(rows_r), len(rows_d) * 2.1],
        hspace=0.075,
        left=left_b, right=right_b, top=0.722, bottom=0.135,
    )
    axg, axr, axd = (fig.add_subplot(gs[i]) for i in range(3))

    for ax in (axg, axr, axd):
        ax.set_xlim(LO, HI)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_yticks([])
        ax.tick_params(axis="x", length=0, labelbottom=False)

    seq_pos = ax_seq.get_position()
    x_a0 = seq_pos.x0 + seq_pos.width * ((L_up - 0.5) - (-9.8)) / ((L_tot + 0.8) - (-9.8))
    x_a1 = seq_pos.x0 + seq_pos.width * ((L_up + L_ins - 0.5) - (-9.8)) / ((L_tot + 0.8) - (-9.8))
    y_a_bot = seq_pos.y0 + seq_pos.height * (-0.78 - (-0.92)) / (1.85 - (-0.92))

    g_pos = axg.get_position()
    x_b0 = g_pos.x0 + g_pos.width * (P - LO) / (HI - LO)
    x_b1 = g_pos.x0 + g_pos.width * (P + D - LO) / (HI - LO)
    y_b_top = g_pos.y0 + g_pos.height * (0.30 - (-1.05)) / (1.55 - (-1.05))

    fig.patches.append(
        Polygon(
            [(x_a0, y_a_bot), (x_a1, y_a_bot), (x_b1, y_b_top), (x_b0, y_b_top)],
            closed=True, facecolor=BAND, edgecolor="none", alpha=0.45,
            transform=fig.transFigure, zorder=0,
        )
    )

    axg.set_ylim(-1.05, 1.55)
    for name, s, e in GENES:
        a, b = max(s, LO), min(e, HI)
        axg.text((a + b) / 2, -0.30, name, ha="center", va="center", fontsize=12.5, color=MUTED, zorder=5)
    axg.add_patch(Rectangle((P, -0.95), D, 1.25, color=ORANGE, lw=0, zorder=3))
    axg.text(P + D / 2, 0.55, "31 bases", ha="center", va="bottom", fontsize=12, fontweight="bold", color=ORANGE_DARK)

    def draw_pileup(ax, rows, gapped):
        nrow = len(rows)
        ax.set_ylim(nrow + 0.5, -0.5)
        ax.add_patch(Rectangle((P, -0.5), D, nrow + 1, color=BAND_LIGHT, lw=0, zorder=0))
        ax.axvline(P, color="#F3CBB0", lw=0.8, zorder=1)
        ax.axvline(P + D, color="#F3CBB0", lw=0.8, zorder=1)
        for i, t in enumerate(rows):
            s, e = max(t[0], LO), min(t[1], HI)
            if gapped and len(t) == 4:
                gs_, gl = t[2], t[3]
                ax.plot([s, gs_], [i, i], color=ORANGE, lw=0.9, solid_capstyle="butt", zorder=2)
                ax.plot([gs_ + gl, e], [i, i], color=ORANGE, lw=0.9, solid_capstyle="butt", zorder=2)
            else:
                ax.plot([s, e], [i, i], color=GREY, lw=1.0, solid_capstyle="butt", zorder=2)

    draw_pileup(axr, rows_r, True)
    draw_pileup(axd, rows_d, False)

    axr.annotate(
        "Each line = 1 sequenced molecule\n(sorted by start position)",
        xy=(P - 62, 300),
        xytext=(LO + 10, 295),
        fontsize=11.5,
        color=MUTED,
        va="center",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=RULE, lw=0.9),
        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0, connectionstyle="arc3,rad=0.08"),
        zorder=6,
    )

    xr = HI + 12

    def draw_side_label(ax, y_mid, line1, line2, col1, col2):
        y0, y1 = ax.get_ylim()  # inverted: y0 > y1
        ax_h = ax.get_position().height
        dy = (y0 - y1) * (0.0135 / ax_h)
        ax.text(xr, y_mid - dy, line1, color=col1, fontweight="bold", va="center", fontsize=13.5, clip_on=False)
        ax.text(xr, y_mid + dy, line2, color=col2, fontweight="normal", va="center", fontsize=13.0, clip_on=False)

    draw_side_label(axr, len(rj) * 0.48,
                    f"RNA: {len(rj)} molecules", "missing the 31 bases (94%)",
                    ORANGE_DARK, ORANGE_DARK)
    draw_side_label(axr, len(rj) + len(rt) * 0.55,
                    f"RNA: {len(rt)} molecules", "with the 31 bases (6%)",
                    MUTED, MUTED)
    draw_side_label(axd, len(rows_d) * 0.50,
                    f"DNA: {len(dt)} fragments", "all with the 31 bases (100%)",
                    INK, MUTED)

    axr.axhline(len(rj) - 0.5, color="white", lw=1.8, zorder=5)

    axd.tick_params(axis="x", labelbottom=True, labelsize=12, colors=MUTED, length=4, width=0.8)
    rel_ticks = [-200, -100, 0, 100, 200]
    axd.set_xticks([P + r for r in rel_ticks])
    axd.set_xticklabels([f"{r:+d}" if r != 0 else "0" for r in rel_ticks])
    axd.set_xlabel("Position relative to the 31-base stretch (bases)", fontsize=12, color=MUTED, labelpad=6)

    fig.text(
        0.045, 0.062,
        f"Ca. Lokiarchaeum ossiferum B35 (CP104013.1, position 0 = 4,487,605; RuBisCO-family gene; Schleper lab culture). "
        f"{anyread} of {len(rj)} gapped RNA molecules lack the exact 31 bases predicted from DNA.\n"
        "RNA: SRR27313660 & SRR32148620 (first 10M read pairs each, unstranded; identical-end pairs merged). "
        "DNA: SRR22273278 (first 8M pairs, 46 fragments).",
        fontsize=10.5, color=MUTED, va="top",
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
