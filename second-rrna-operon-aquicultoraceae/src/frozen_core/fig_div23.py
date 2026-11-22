"""Figure for the divergent 23S post: PTC bases (left) and DNA vs RNA reads on the Crystal Geyser loci (right).

usage: uv run -q --no-project --with matplotlib python -I src/frozen_core/fig_div23.py OUT.png
Read-pair (fragment) counts are recomputed here from the PAFs in data/expression_23s (identity >= 0.98, alignment >= 60 nt,
Crystal Geyser population loci only). Conservation comes from data/ptc/columns.tsv. Divergent-copy bases are
those of the Crystal Geyser (Aquicultor_A) copy, read by structure-free pairwise alignment with src/frozen_core/ptc_bases_pair.py;
each listed base agrees whether mapped directly from E. coli or via either DASXVE01 copy. The 2449-2453 loop is
aligner-dependent (a 1-nt gap moves), so it is shown as one row: the UAACA motif is absent under every alignment.
"""
import collections
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# repo layout (run from the folder root); in the lab: ROOT = work/frozen-core/results/lsu/bac/, EXPR = ROOT + "cca/expr/"
ROOT = "data/ptc/"
EXPR = "data/expression_23s/"
CG = ["PFKS", "PFIW", "MNXI", "JAACTV", "PFTP", "PFNG", "PEXG", "PFFU"]
DIV = {2062: "C", 2251: "U", 2252: "U", 2253: "U", 2503: "U", 2504: "A", 2583: "A", 2584: "A"}
LOOP = (2449, 2453, "UAACA", "absent")
INK, MUTED, NORM_C, DIV_C, GRID = "#1f2328", "#59636e", "#8c959f", "#c4432b", "#d0d7de"


def counts(paf):
    """Fragments, not alignment lines: fastq-dump without --split wrote each mate pair as one 302-nt read, which
    minimap2 aligns as about two local hits. A fragment counts for the class of its passing CG-locus alignments."""
    cls = {l.split("\t")[0]: l.split("\t")[1] for l in open(EXPR + "loci_class_all.tsv")}
    frag = collections.defaultdict(set)
    for line in open(EXPR + paf):
        f = line.split("\t")
        t, m, a = f[5], int(f[9]), int(f[10])
        if t in cls and a >= 60 and m / a >= 0.98 and any(x in t.split("|")[1] for x in CG):
            frag[f[0]].add(cls[t])
    c = collections.Counter("both" if len(v) > 1 else next(iter(v)) for v in frag.values())
    assert not c["both"], c
    return c


def main(out):
    cons = {}
    for line in open(ROOT + "columns.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[1].isdigit() and (int(f[1]) in DIV or LOOP[0] <= int(f[1]) <= LOOP[1]):
            cons[int(f[1])] = (f[2], float(f[5]))
    runs = [("DNA  2015-05-28", "SRR6051694.loci.paf"), ("RNA  2015-05-25", "SRR11558543.loci.paf"),
            ("RNA  2015-05-26", "SRR11558542.loci.paf"), ("RNA  2015-05-27", "loci.paf")]
    rc = [(lab, counts(p)) for lab, p in runs]
    for lab, c in rc:
        print(lab, dict(c))
    rows = [(str(p), *cons[p], DIV[p]) for p in sorted(DIV) if p < LOOP[0]]
    lo = min(cons[p][1] for p in range(LOOP[0], LOOP[1] + 1))
    rows.append((f"{LOOP[0]}-{LOOP[1] % 100}", LOOP[2], lo, LOOP[3]))
    rows += [(str(p), *cons[p], DIV[p]) for p in sorted(DIV) if p > LOOP[1]]

    fig = plt.figure(figsize=(12.5, 5.6), dpi=150)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.02, 0.07, 0.40, 0.72])
    ax.axis("off")
    ax.set_xlim(0, 4.2)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    fig.text(0.02, 0.93, "Peptidyl-transferase-centre bases", fontsize=12, color=INK, weight="bold")
    for x, h in ((0.45, "E. coli\nposition"), (1.55, "bacteria (census)\nbase, share"), (3.15, "divergent copy,\nCrystal Geyser")):
        fig.text(0.02 + 0.40 * x / 4.2, 0.835, h, ha="center", va="center", fontsize=9.5, color=MUTED)
    for i, (pos, nt, fc, d) in enumerate(rows):
        loop = "-" in pos
        ax.text(0.45, i, pos, ha="center", va="center", fontsize=11, color=INK)
        ax.text(1.55, i, f"{nt}  {'≥' if loop else ''}{fc * 100:.2f}%", ha="center", va="center", fontsize=10.5,
                color=MUTED, family="monospace" if loop else None)
        ax.text(3.15, i, d, ha="center", va="center", fontsize=10.5 if loop else 13, family=None if loop else "monospace",
                weight="bold", color=DIV_C if d != nt else INK)
        ax.plot([0, 4.0], [i + 0.5, i + 0.5], color=GRID, lw=0.6)

    bx = fig.add_axes([0.58, 0.22, 0.30, 0.58])
    for i, (lab, c) in enumerate(rc):
        n = c["div"] + c["norm"]
        d = c["div"] / n
        if d > 0:
            bx.barh(i, d, color=DIV_C, height=0.55)
        bx.barh(i, 1 - d, left=d + (0.004 if d > 0 else 0), color=NORM_C, height=0.55)
        bx.text(1.03, i, f"{c['div']} : {c['norm']}", va="center", fontsize=10.5, color=INK)
    bx.set_yticks(range(len(rc)), [lab for lab, _ in rc], fontsize=10.5, color=INK)
    bx.invert_yaxis()
    bx.set_xlim(0, 1)
    bx.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", "25%", "50%", "75%", "100%"], fontsize=9.5, color=MUTED)
    for sp in ("top", "right", "left"):
        bx.spines[sp].set_visible(False)
    bx.spines["bottom"].set_color(GRID)
    bx.tick_params(length=0)
    bx.text(1.03, -0.7, "div : normal", fontsize=9.5, color=MUTED)
    fig.text(0.47, 0.93, "Read pairs on the same Crystal Geyser 23S loci, identity ≥98%", fontsize=12, color=INK, weight="bold")
    fig.text(0.47, 0.875, "share of read pairs on the divergent copy (red) vs the normal copy (grey)", fontsize=9.5, color=MUTED)
    fig.text(0.47, 0.06, "One species, 8 genomes from Crystal Geyser, Utah.\nDNA = metagenome; RNA = viral-fraction metatranscriptomes.",
             fontsize=8.5, color=MUTED, linespacing=1.4)
    fig.savefig(out, facecolor="white")

if __name__ == "__main__":
    main(sys.argv[1])
