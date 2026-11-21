"""Two- or three-row SRL figure for outbox/2026-11-14-srl-gaaa.
usage: uv run --no-project --with matplotlib python -I fig_srl_rows.py OUT.png [--pfreud]
Sequences (E. coli 2652-2667) were read from the genomes: E. coli K-12 GCF_000005845.2 (7/7 copies),
Z. mobilis 18 complete assemblies (61/61 copies, agents/altair/kills/srl_zymo/dl), P. freudenreichii 53 complete
assemblies (work/curation/pfreud/srl_counts.tsv, every copy GCAA)."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORANGE, INK, MUTED, PALE = "#E8650A", "#1a1a1a", "#8a8580", "#fbe3d1"
rows = [("E. coli", "almost all bacteria", "CUAGUACGAGAGGACC"),
        ("Zymomonas mobilis", "all 18 complete genomes, every copy", "UUAGUACGAAAGGACC")]
if "--pfreud" in sys.argv:
    rows.append(("Propionibacterium freudenreichii", "all 53 complete genomes, every copy", "UUAGUACGCAAGGACC"))
ref = rows[0][2]
plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
H = 0.9 + 0.95 * len(rows)
fig = plt.figure(figsize=(9.2, H), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.01, 0.16, 0.98, 0.83])
n = len(rows)
ax.set_xlim(-8.6, 15.6); ax.set_ylim(-0.85, n + 0.4); ax.axis("off")
ax.add_patch(plt.Rectangle((6.5, -0.5), 4, n, color=PALE, lw=0))
for i, p in ((0, 2652), (7, 2659), (10, 2662), (15, 2667)):
    ax.text(i, -0.78, str(p), ha="center", fontsize=11, color=MUTED)
# one bracket over the shaded loop, with a plain label (no per-site callouts, no title)
yb = n - 0.35
ax.plot([6.6, 6.6, 10.4, 10.4], [yb, yb + 0.18, yb + 0.18, yb], color=INK, lw=1.2, solid_capstyle="butt")
ax.text(8.5, yb + 0.32, "the loop ricin attacks", ha="center", va="bottom", fontsize=13, color=INK)
for k, (name, sub, s) in enumerate(rows):
    y = n - 1 - k
    for i, (b, b0) in enumerate(zip(s, ref)):
        diff = b != b0 and 7 <= i <= 10
        ax.text(i, y, b, ha="center", va="center", family="Menlo", fontsize=22,
                color=ORANGE if diff else INK, fontweight="bold" if diff else "normal")
    ax.text(-1.0, y + 0.14, name, ha="right", va="center", fontsize=16, style="italic", color=INK)
    ax.text(-1.0, y - 0.26, sub, ha="right", va="center", fontsize=11, color=MUTED)
fig.text(0.075, 0.045, "23S ribosomal RNA, E. coli positions 2652-2667. Orange: the letter that differs from E. coli.",
         fontsize=10.5, color=MUTED, ha="left")
fig.savefig(sys.argv[1], facecolor="white")
