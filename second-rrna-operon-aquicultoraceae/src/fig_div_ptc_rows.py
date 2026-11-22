"""Two-row letter figure: normal vs second (divergent) 23S of Ca. Aquicultor secundus (GCA_009993245.1),
an aligned stretch around the P loop (E. coli G2251-G2253 = GGG; UUU in the second copy).
Sequences come from data/silent_operon/origin/aln23.sto (cmalign; work/silent-operon/B_origin/out/aln23.sto in the lab), rows 23S_013 (norm) and 23S_012 (div).
usage (from the folder root): uv run python -I src/fig_div_ptc_rows.py OUT.png"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALN = "data/silent_operon/origin/aln23.sto"
seqs = {}
for l in open(ALN):
    if l.startswith(("#", "//")) or not l.strip():
        continue
    p = l.split()
    if len(p) == 2:
        seqs[p[0]] = seqs.get(p[0], "") + p[1]
up = lambda s: s.upper().replace("T", "U").replace(".", "-")
N = up(seqs["23S_013|GCA_009993245.1|Aquicultor_A|norm"])
D = up(seqs["23S_012|GCA_009993245.1|Aquicultor_A|div"])
cols = [i for i in range(len(N)) if not (N[i] == "-" and D[i] == "-")]
N = "".join(N[i] for i in cols); D = "".join(D[i] for i in cols)
k = next(i for i in range(len(N) - 3) if N[i:i + 3] == "GGG" and D[i:i + 3] == "UUU")
a, b = k - 30, k + 33
n, d = N[a:b], D[a:b]
ORANGE, INK, MUTED, PALE = "#E8650A", "#1a1a1a", "#8a8580", "#fbe3d1"
plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
L = len(n)
fig = plt.figure(figsize=(13, 1.85), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0.0, 0.17, 1.0, 0.82]); ax.axis("off")
ax.set_xlim(-11.5, L + 0.5); ax.set_ylim(-0.55, 1.95)
ax.add_patch(plt.Rectangle((30 - 0.5, -0.45), 3, 1.9, color=PALE, lw=0))
yb = 1.55
ax.plot([29.6, 29.6, 32.4, 32.4], [yb - 0.1, yb, yb, yb - 0.1], color=INK, lw=1.1)
ax.text(31, yb + 0.07, "holds the tRNA", ha="center", va="bottom", fontsize=11.5, color=INK)
for y, name, s in ((1.0, "normal 23S", n), (0.0, "second 23S", d)):
    ax.text(-0.8, y, name, ha="right", va="center", fontsize=14, color=INK)
    for i, ch in enumerate(s):
        diff = y == 0.0 and ch != n[i]
        ax.text(i, y, ch, ha="center", va="center", family="Menlo", fontsize=12.5,
                color=ORANGE if diff else INK, fontweight="bold" if diff else "normal")
same = sum(1 for x, y2 in zip(n, d) if x == y2 and x != "-")
fig.text(0.07, 0.04, f"Ca. Aquicultor secundus 23S, aligned stretch at the catalytic centre (it differs more here than across the whole gene, which is ~79% identical). "
         f"Orange: differences ({L - same} of {L}).", fontsize=10, color=MUTED, ha="left")
fig.savefig(sys.argv[1], facecolor="white")
print(n); print(d)
