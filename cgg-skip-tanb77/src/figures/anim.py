"""ribosome.gif: the rpoC skip site of GCA_963931445.1 as an animation, in the style of the TTA post's
01_ribosome_final_aa_v2 (same palette, two-subunit ribosome, smootherstep motion, amino-acid beads).
The ribosome steps 3-3-3-3, stops on CGG (no tRNA reads it), hops one letter to GGA, reads glycine, and continues in
the new frame. Letters come from site.json (site.py).
usage: uv run --no-project --with matplotlib python -I src/anim.py <frames_dir> <dpi> <fps> [still t ...]"""
import json, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(f"{HERE}/site.json"))
OR, INK, LIGHT, RIB, BG, MUTE, FAINT = "#E8650A", "#1c1c1a", "#F8DCC6", "#DAD9D3", "#ffffff", "#8a8a84", "#c9c8c2"
B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {x + y + z: AA[16 * i + 4 * j + k] for i, x in enumerate(B) for j, y in enumerate(B) for k, z in enumerate(B)}
THREE = dict(F="Phe", L="Leu", S="Ser", Y="Tyr", C="Cys", W="Trp", P="Pro", H="His", Q="Gln", R="Arg", I="Ile", M="Met",
             T="Thr", N="Asn", K="Lys", V="Val", A="Ala", D="Asp", E="Glu", G="Gly")
CAPTION = "RNA polymerase gene rpoC · CALXVA01, order TANB77, chicken gut (GCA_963931445.1) · real sequence"

NL, NR = 4, 4
up = S["dna_up"][-3 * NL:]; dn = S["dna_dn"][:3 * NR]
LET = up + S["c"] + dn                                   # 25 letters
CI = 3 * NL                                              # index of the skipped C
assert LET[CI:CI + 3] == "CGG" and LET[CI + 1:CI + 4] == "GGA"
POS = [0, 3, 6, 9, CI, CI + 1, CI + 4, CI + 7, CI + 10]  # ribosome start letter at each stop
READ = {p: LET[p:p + 3] for p in POS if p != CI}         # codons read (CGG is not)
AAS = {p: THREE[CODE[c]] for p, c in READ.items()}
assert AAS[CI + 1] == "Gly"
GROUP = [(0, 3), (3, 6), (6, 9), (9, 12), (12, 13), (13, 16), (16, 19), (19, 22), (22, 25)]  # boxes drawn

FW, FH = 12.8, 7.2
LW, GAPW = 0.37, 0.20
xl, x = [], 0.0
gi = {}
for g, (a, b) in enumerate(GROUP):
    for i in range(a, b): xl.append(x + (i - a) * LW + LW / 2); gi[i] = g
    x += (b - a) * LW + GAPW
X0 = (FW - (x - GAPW)) / 2
xl = [X0 + v for v in xl]
YL, CH = 3.35, 0.62

# timeline: dwell at each stop, then move to the next
T0, DW, MV, DW_CGG, MV_HOP, HOLD = 0.4, 0.26, 0.26, 1.25, 0.45, 2.2
sched, t = [], T0
for k, p in enumerate(POS):
    dw = DW_CGG if p == CI else DW
    mv = 0 if k == len(POS) - 1 else (MV_HOP if p == CI else MV)
    sched.append((p, t, t + dw, t + dw + mv)); t += dw + mv
TOTAL = t + HOLD


def ease(u):
    u = min(max(u, 0.0), 1.0); return u ** 3 * (u * (6 * u - 15) + 10)


def xc(p):                    # centre of the 3 letters from p
    return (xl[p] + xl[p + 2]) / 2


def state(t):
    if t < T0: return xc(POS[0]), float(POS[0]), set()
    done = set()
    for k, (p, a, b, c) in enumerate(sched):
        if t < a: break
        if p != CI and t >= a + 0.6 * (b - a): done.add(p)
        if t < b: return xc(p), float(p), done
        if t < c:
            u = ease((t - b) / (c - b)); return xc(p) + u * (xc(POS[k + 1]) - xc(p)), p + u * (POS[k + 1] - p), done
    return xc(POS[-1]), float(POS[-1]), done


def draw(ax, t):
    ax.clear(); ax.set_xlim(0, FW); ax.set_ylim(0, FH); ax.axis("off")
    xr, pos, done = state(t)
    hopped = pos > CI + 0.5
    # DNA boxes
    for g, (a, b) in enumerate(GROUP):
        l, r = xl[a] - LW / 2 + 0.02, xl[b - 1] + LW / 2 - 0.02
        if a == CI:
            if hopped:
                ax.add_patch(Rectangle((l, YL - CH / 2), r - l, CH, facecolor="none", edgecolor=FAINT, lw=1.5, ls=(0, (3, 3)), zorder=2))
            else:
                ax.add_patch(Rectangle((l, YL - CH / 2), r - l, CH, facecolor=LIGHT, edgecolor="none", zorder=2))
            continue
        gly = a == CI + 1 and (CI + 1) in done
        ax.add_patch(Rectangle((l, YL - CH / 2), r - l, CH, facecolor=OR if gly else LIGHT, edgecolor="none", zorder=2))
    for i, ch in enumerate(LET):
        gly = gi[i] == gi[CI + 1] and (CI + 1) in done
        col = FAINT if (i == CI and hopped) else ("white" if gly else INK)
        ax.text(xl[i], YL, ch, ha="center", va="center", fontsize=25, family="Menlo", color=col, zorder=6,
                fontweight="bold" if gly else "normal")
    # ribosome
    ls_bot, ls_h, ls_w = YL + CH / 2 - 0.1, 0.9, 1.75
    ss_top, ss_h, ss_w = YL - CH / 2 + 0.1, 0.44, 1.35
    ax.add_patch(FancyBboxPatch((xr - ls_w / 2, ls_bot), ls_w, ls_h, boxstyle="round,pad=0,rounding_size=0.36",
                                facecolor=RIB, edgecolor=INK, lw=1.4, zorder=5))
    ax.add_patch(FancyBboxPatch((xr - ss_w / 2, ss_top - ss_h), ss_w, ss_h, boxstyle="round,pad=0,rounding_size=0.2",
                                facecolor=RIB, edgecolor=INK, lw=1.4, zorder=5))
    # labels under the ribosome
    ylab = ss_top - ss_h - 0.36
    p_cgg = sched[POS.index(CI)]
    if p_cgg[1] <= t < p_cgg[2]:
        ax.text(xc(CI), ylab, "CGG: no tRNA reads it", ha="center", va="center", fontsize=21, color=INK)
    elif p_cgg[2] <= t < p_cgg[3]:
        ax.text(xc(CI) + 0.2, ylab, "skip 1 letter", ha="center", va="center", fontsize=21, color=INK)
    elif hopped and t < sched[POS.index(CI + 1)][3] + 0.15:
        ax.text(xc(CI + 1), ylab, "GGA = glycine", ha="center", va="center", fontsize=21, color=OR, fontweight="bold")
    if t >= sched[-1][2] + 0.2:
        al = ease((t - sched[-1][2] - 0.2) / 0.5)
        l, r = xl[CI + 1] - LW / 2, xl[len(LET) - 1] + LW / 2
        yb = ylab - 0.12
        ax.plot([l, r], [yb + 0.3] * 2, color=INK, lw=1.6, alpha=al, solid_capstyle="butt")
        ax.plot([l, l], [yb + 0.3, yb + 0.45], color=INK, lw=1.6, alpha=al); ax.plot([r, r], [yb + 0.3, yb + 0.45], color=INK, lw=1.6, alpha=al)
        ax.text((l + r) / 2, yb - 0.12, "the rest of the gene: one letter over", ha="center", va="center", fontsize=21, color=INK, alpha=al)
    # protein: beads above each read codon, joined by a line into the ribosome
    top = ls_bot + ls_h; by = top + 0.4
    bw, bh = 0.86, 0.46
    read = sorted(done)
    xs = [xc(p) for p in read]
    tail = xr
    x = min(xs[0] if xs else xr, xr) - 0.6
    ax.plot([x, xs[0] - bw / 2 if xs else xr], [by] * 2, color=INK, lw=1.8, solid_capstyle="butt", zorder=3)
    for k in range(len(xs) - 1):
        ax.plot([xs[k] + bw / 2, xs[k + 1] - bw / 2], [by] * 2, color=INK, lw=1.8, solid_capstyle="butt", zorder=3)
    last = xs[-1] + bw / 2 if xs else x
    if tail > last:
        ax.plot([last, tail], [by] * 2, color=INK, lw=1.8, solid_capstyle="butt", zorder=3)
    under = any(abs(v - xr) < bw / 2 for v in xs)
    ax.plot([xr, xr], [top, by - (bh / 2 if under else 0)], color=INK, lw=1.8, solid_capstyle="butt", zorder=3)
    for p, v in zip(read, xs):
        gly = p == CI + 1
        ax.add_patch(FancyBboxPatch((v - bw / 2, by - bh / 2), bw, bh, boxstyle="round,pad=0,rounding_size=0.12",
                                    facecolor=OR if gly else LIGHT, edgecolor=INK, lw=1.8, zorder=4))
        ax.text(v, by, AAS[p], ha="center", va="center", fontsize=18, color="white" if gly else INK,
                fontweight="bold" if gly else "normal", zorder=4)
    # title removed 16 Nov at Ryan's request
    ax.text(FW / 2, 0.45, CAPTION, ha="center", va="center", fontsize=13.5, color=MUTE)


def main(out, dpi, fps, stills=()):
    os.makedirs(out, exist_ok=True)
    plt.rcParams["font.family"] = "Helvetica"
    fig = plt.figure(figsize=(FW, FH), dpi=dpi); fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    if stills:
        for t in stills:
            draw(ax, t if t >= 0 else TOTAL - 1 / fps); fig.savefig(f"{out}/still_{t:.2f}.png", facecolor=BG)
        print(LET, AAS, f"TOTAL {TOTAL:.2f}", [(p, round(a, 2), round(b, 2), round(c, 2)) for p, a, b, c in sched]); return
    n = int(round(TOTAL * fps))
    for i in range(n):
        draw(ax, i / fps); fig.savefig(f"{out}/f{i:04d}.png", facecolor=BG)
    print(f"{n} frames, {TOTAL:.2f} s at {fps} fps, {int(FW * dpi)}x{int(FH * dpi)}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), [float(x) for x in sys.argv[4:]])
