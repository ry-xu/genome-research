"""image.png: one real CGG skip site, rpoC of GCA_963931445.1 (TANB77, CAG-508, CALXVA01; chicken caecum).
All letters come from site.json, written by site.py from the genome FASTA and the B. subtilis genome.
usage: uv run --no-project --with matplotlib python -I src/image.py"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(f"{HERE}/site.json"))
OR, INK, MUTE, FAINT, LIGHT, BG = "#E8650A", "#1c1c1a", "#8a8a84", "#c9c8c2", "#F8DCC6", "#ffffff"
B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {x + y + z: AA[16 * i + 4 * j + k] for i, x in enumerate(B) for j, y in enumerate(B) for k, z in enumerate(B)}
THREE = dict(F="Phe", L="Leu", S="Ser", Y="Tyr", C="Cys", W="Trp", P="Pro", H="His", Q="Gln", R="Arg", I="Ile", M="Met",
             T="Thr", N="Asn", K="Lys", V="Val", A="Ala", D="Asp", E="Glu", G="Gly", X="?")

NL, NR = 6, 6                                    # codons shown left of the C and from the glycine on
up = [S["dna_up"][i:i + 3] for i in range(0, len(S["dna_up"]), 3)][-NL:]
dn = [S["dna_dn"][i:i + 3] for i in range(0, len(S["dna_dn"]), 3)][:NR]
prot = [CODE[c] for c in up] + [CODE[c] for c in dn]
assert "".join(prot) == S["prot_window"][10 - NL:10 + NR], (prot, S["prot_window"])
assert S["c"] == "C" and dn[0] == "GGA" and CODE[dn[0]] == "G"
REF = "Bacillus subtilis"
ref = list(S["refs"][REF]["window"][10 - NL:10 + NR])
old = [S["c"] + S["dna_dn"][0:2]] + [S["dna_dn"][i:i + 3] for i in range(2, 12, 3)]   # CGG AAA TGT TAA
old_aa = []
for c in old:
    old_aa.append("stop" if CODE[c] == "*" else THREE[CODE[c]])
    if CODE[c] == "*": break
old = old[:len(old_aa)]
assert old[0] == "CGG"

W, H, DPI = 16.0, 9.0, 100                       # 1600 x 900 px
fig = plt.figure(figsize=(W, H), dpi=DPI); fig.patch.set_facecolor(BG)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
plt.rcParams["font.family"] = "Helvetica"
MONO = "Menlo"

# horizontal layout: label column, then codon cells; the skipped C gets a narrow cell
X0, CW, GAP, CWC = 3.30, 0.86, 0.08, 0.46
xs, x = [], X0
for k in range(NL):
    xs.append(x); x += CW + GAP
xc = x; x += CWC + GAP
for k in range(NR):
    xs.append(x); x += CW + GAP
cx = [v + CW / 2 for v in xs]
LX = 0.55
y_dna, y_prot, y_ref = 5.55, 4.35, 2.85
ch = 0.70

ax.text(LX, 8.25, "CGG: skip the C, read GGA as glycine", fontsize=38, fontweight="bold", color=INK, va="center")
ax.text(LX, 7.55, "RNA polymerase gene (rpoC) of a TANB77 gut bacterium from a chicken", fontsize=21, color=MUTE, va="center")

# row labels
ax.text(LX, y_dna, "DNA", fontsize=24, fontweight="bold", color=INK, va="center")
ax.text(LX, y_prot, "Protein", fontsize=24, fontweight="bold", color=INK, va="center")
ax.text(LX, y_ref, "Same protein in\nBacillus subtilis", fontsize=22, fontweight="bold", color=INK, va="center",
        linespacing=1.25)

# DNA codons
for k, c in enumerate(up + dn):
    gly = k == NL
    ax.add_patch(plt.Rectangle((xs[k], y_dna - ch / 2), CW, ch, facecolor=OR if gly else LIGHT, edgecolor="none"))
    ax.text(cx[k], y_dna, c, ha="center", va="center", fontsize=25, family=MONO, color="white" if gly else INK,
            fontweight="bold" if gly else "normal")
# the skipped C: grey, no fill, struck through
ax.add_patch(plt.Rectangle((xc, y_dna - ch / 2), CWC, ch, facecolor="none", edgecolor=FAINT, lw=1.5, ls=(0, (3, 3))))
ax.text(xc + CWC / 2, y_dna, S["c"], ha="center", va="center", fontsize=25, family=MONO, color=FAINT)
# marker: "C skipped" above it, straight arrow into the glycine codon
ya = y_dna + ch / 2 + 0.38
ax.text(xc + CWC / 2, ya + 0.32, "skipped", ha="center", va="center", fontsize=19, color=INK)
ax.annotate("", xy=(xs[NL] + 0.12, ya), xytext=(xc + 0.04, ya),
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=2, mutation_scale=22, shrinkA=0, shrinkB=0))
ax.plot([xc + 0.04, xc + 0.04], [ya - 0.12, ya + 0.0], color=INK, lw=2, solid_capstyle="butt")

# protein of this bacterium, 3-letter, directly under each codon
for k, a in enumerate(prot):
    gly = k == NL
    ax.text(cx[k], y_prot, THREE[a], ha="center", va="center", fontsize=26, color=OR if gly else INK,
            fontweight="bold" if gly else "normal")
# connector from codon to amino acid
for k in range(len(prot)):
    ax.plot([cx[k], cx[k]], [y_dna - ch / 2 - 0.1, y_prot + 0.3], color=OR if k == NL else FAINT, lw=1.6 if k == NL else 1.0)

# relative: identical residues in ink with a bar, others muted
nid = 0
for k, (a, b) in enumerate(zip(prot, ref)):
    same = a == b; nid += same
    if same:
        ax.plot([cx[k], cx[k]], [y_prot - 0.3, y_ref + 0.32], color=OR if k == NL else INK, lw=2.4 if k == NL else 1.6)
    ax.text(cx[k], y_ref, THREE[b], ha="center", va="center", fontsize=26,
            color=(OR if k == NL else INK) if same else MUTE, fontweight="bold" if k == NL else "normal")

# what the old frame would do
y_old = 1.55
ax.text(LX, y_old, "Without the skip:", fontsize=20, color=INK, fontweight="bold", va="center")
t = ax.text(LX + 2.75, y_old, "  ".join(old), fontsize=21, family=MONO, color=INK, va="center")
bb = t.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
ax.text(bb.x1 + 0.25, y_old, "=  " + " ".join(old_aa) + ": the protein is cut short", fontsize=20, color=INK, va="center")

ax.text(LX, 0.62, f"Genome GCA_963931445.1 (genus CALXVA01, order TANB77), from GTDB. Bars: same amino acid in B. subtilis "
        f"RpoC ({nid} of {len(prot)}).\nIllumina DNA reads from the same sample: 7 of 7 keep the C.", fontsize=16, color=MUTE,
        va="center", linespacing=1.5)

out = os.path.join(os.path.dirname(HERE), "image.png")
fig.savefig(out, facecolor=BG, dpi=DPI)
print(out, "identity", nid, "of", len(prot), "old frame", old, old_aa, "ref", "".join(ref), "prot", "".join(prot))
