"""Outbox figure: where each genome's ATA codons sit (conserved Met vs Ile positions).
Counts are parsed from b_read.out and b_read_context.out, nothing typed by hand.
usage (from the folder root): uv run python -I src/fig_aua_bars.py OUT.png
"""
import os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "codons")
ORANGE, GREY, INK, MUTED, BG = "#eb6834", "#c3c2b7", "#0b0b0b", "#52514e", "#fcfcfb"

counts = {}
for f in ("b_read.out", "b_read_context.out"):
    for line in open(os.path.join(HERE, f)):
        m = re.match(r"(\S+) drop \S+ ATA:M=(\d+)/I=(\d+)", line)
        if m:
            counts[m.group(1)] = (int(m.group(2)), int(m.group(3)))

# (label, sublabel, key); None = gap
ROWS = [("Choristocarpus", "brown alga", "cten_WGS"),
        ("Discosporangium", "its sister genus", "dmes_WGS"),
        ("Phaeostrophion", "brown alga", "NEG_Phaeostrophion"),
        ("Protohalopteris", "brown alga", "NEG_Protohalopteris"),
        ("Sargassum", "brown alga", "NEG_Sargassum"),
        ("Saccharina", "brown alga (kelp)", "NEG_Saccharina"),
        ("Petalonia", "brown alga", "NEG_Petalonia"),
        None,
        ("Sphaerosorus", "yellow-green alga\n(class reads AUA as Met)", "CTX_Sphaerosorus")]

plt.rcParams.update({"font.family": "Helvetica", "font.size": 13})
fig, ax = plt.subplots(figsize=(9, 5.6), dpi=200)
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
y = 0
for r in ROWS:
    if r is None:
        y += 0.6; continue
    name, sub, key = r
    m, i = counts[key]
    s = m / (m + i)
    hl = key == "cten_WGS"
    ax.barh(-y, s, color=ORANGE, height=0.62, edgecolor=BG, linewidth=1.5)
    ax.barh(-y, 1 - s, left=s, color=GREY, height=0.62, edgecolor=BG, linewidth=1.5)
    ax.text(-0.03, -y + 0.07, name, ha="right", va="center", color=INK, fontsize=13,
            fontstyle="italic", fontweight="bold" if hl else "normal")
    ax.text(-0.03, -y - (0.42 if "\n" in sub else 0.27), sub, ha="right", va="center", color=MUTED, fontsize=9.5, linespacing=1.0)
    ax.text(1.02, -y, f"{m} of {m + i}", ha="left", va="center", color=INK, fontsize=12,
            fontweight="bold" if hl else "normal")
    y += 1
top = 0.75
ax.text(0.0, top, "at methionine positions", color=ORANGE, fontsize=12, fontweight="bold", va="bottom")
ax.text(1.0, top, "at isoleucine positions", color=MUTED, fontsize=12, fontweight="bold", va="bottom", ha="right")
ax.text(1.02, top, "AUA at Met", color=MUTED, fontsize=10, va="bottom")
ax.set_xlim(0, 1); ax.set_ylim(-y + 0.4, 1.0)
ax.axis("off")
fig.text(0.5, 0.955, "Where each mitochondrial genome uses the codon AUA", ha="center", fontsize=15, color=INK)
fig.subplots_adjust(left=0.30, right=0.86, top=0.90, bottom=0.03)
fig.savefig(sys.argv[1], facecolor=BG)
print({k: counts[k] for k in [r[2] for r in ROWS if r]})
