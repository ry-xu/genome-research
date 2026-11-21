"""Figure for outbox/2026-11-14-srl-gaaa. usage: python src/fig_srl.py data/fig_srl OUT.png
All numbers from data/fig_srl/{loops,phyla}.tsv (src/fig_srl_numbers.py). Needs matplotlib (run: uv run --no-project --with matplotlib python -I;
frozen-core has no matplotlib)."""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = sys.argv[1]
loops = list(csv.DictReader(open(f"{d}/loops.tsv"), delimiter="\t"))
phyla = list(csv.DictReader(open(f"{d}/phyla.tsv"), delimiter="\t"))
ORANGE, INK, MUTED, LIGHT, PALE = "#E8650A", "#1a1a1a", "#8a8580", "#d4d0cb", "#fbe3d1"
plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
fig = plt.figure(figsize=(15, 7.4), dpi=200); fig.patch.set_facecolor("white")
axA = fig.add_axes([0.03, 0.12, 0.50, 0.66]); axB = fig.add_axes([0.70, 0.12, 0.27, 0.66])
# A: E. coli 2652-2667, loop 2659-2662 (string index 7-10) shaded; bases differing from E. coli in orange
ref = loops[0]["modal_2652_2667"]
axA.axis("off"); axA.set_xlim(-11.5, 16.5); axA.set_ylim(-0.7, len(loops) + 0.6)
axA.add_patch(plt.Rectangle((6.55, -0.5), 4, len(loops) + 0.15, color=PALE, lw=0))
for i, p in ((0, 2652), (7, 2659), (10, 2662), (15, 2667)):
    axA.text(i, len(loops) + 0.05, str(p), ha="center", fontsize=11, color=MUTED)
for k, r in enumerate(loops):
    y = len(loops) - 1 - k; s = r["modal_2652_2667"]
    for i, (b, b0) in enumerate(zip(s, ref)):
        diff = b != b0
        axA.text(i, y, b, ha="center", va="center", family="Menlo", fontsize=19,
                 color=ORANGE if diff else INK, fontweight="bold" if diff else "normal")
    g, mg, n = int(r["genomes"]), int(r["metagenome_derived"]), int(r["genomes_with_modal"])
    sub = "" if k == 0 else (f"{g} genomes, {mg} from metagenomes" if mg else f"{g} isolate genomes (2 strains)")
    axA.text(-1.2, y + (0.13 if sub else 0), r["group"], ha="right", va="center", fontsize=15, color=INK,
             style="italic" if "Zymomonas" in r["group"] or "coli" in r["group"] else "normal")
    if sub:
        axA.text(-1.2, y - 0.27, sub, ha="right", va="center", fontsize=11, color=MUTED)
axA.set_title("The loop reads GAAA instead of GAGA", loc="left", x=0.0, fontsize=18, fontweight="bold", color=INK, pad=10)
# B: carriers by GTDB phylum, metagenome-derived vs not
y = list(range(len(phyla)))[::-1]
for yi, r in zip(y, phyla):
    n, i = int(r["carriers"]), int(r["cultured_isolates"])
    axB.barh(yi, n - i, height=0.62, color=LIGHT)
    if i: axB.barh(yi, i, left=n - i, height=0.62, color=ORANGE)
    axB.text(n + 0.5, yi, str(n), va="center", fontsize=12, color=INK)
    axB.text(-0.6, yi, r["phylum"], va="center", ha="right", fontsize=12, color=INK)
axB.set_xlim(0, 28); axB.set_ylim(-0.6, len(phyla) - 0.4); axB.axis("off")
axB.set_title("Genomes per GTDB phylum", loc="left", x=-0.5, fontsize=18, fontweight="bold", color=INK, pad=10)
for yy, c, lab in ((-0.10, ORANGE, "cultured isolate"), (-0.16, LIGHT, "metagenome-assembled or symbiont")):
    axB.add_patch(plt.Rectangle((-0.5, yy - 0.012), 0.04, 0.03, color=c, transform=axB.transAxes, clip_on=False))
    axB.text(-0.43, yy, lab, color=INK, fontsize=11, va="center", transform=axB.transAxes)
fig.suptitle(f"Sarcin-ricin loop of 23S rRNA: GAAA in every full-length copy, {sum(int(r['carriers']) for r in phyla)} bacterial genomes",
             x=0.03, ha="left", fontsize=22, fontweight="bold", color=INK, y=0.97)
fig.text(0.03, 0.015, "E. coli numbering, 2652-2667. Shaded: the tetraloop 2659-2662.\nRows show the most common sequence in each group; orange letters differ\nfrom E. coli. Groups are GTDB orders, except Zymomonas (genus).",
         fontsize=11, color=MUTED, ha="left")
fig.savefig(sys.argv[2], facecolor="white")
