"""Figure for outbox/2026-11-14-asgard-protein-intron: molecules over the L02 site, RNA above, DNA below.

  uv run --no-project --with matplotlib python -I src/fig_l02.py   (run from culture/)
      -> out/fig_l02_dna_templates.tsv, ../../../outbox/2026-11-14-asgard-protein-intron/image.png

RNA: the L02 molecules of out/rna_mol_templates.tsv (place 0; written by src/molecules.py for SRR27313660 and
SRR32148620), the table behind the post's 385 / 26. DNA: fragments of out/sam/SRR22273278.sam.zst classified at L02
with the same rules as src/molecules.py (junctions.parse, MIN_ID, FLANK, FLANK_MM), drawn as their templates (outer
ends over mapped mates, clipped bases added back). Counts are asserted against the post before drawing.
"""
import collections, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from junctions import parse, MIN_ID, FLANK, FLANK_MM, CLIP

P, D = 4487605, 31                      # L02 removal, 0-based start and length (data/sites.tsv)
GENES = [("UYP47557", 4487083, 4487611), ("UYP47558", 4487702, 4488440)]   # GFF CDS, 0-based half-open, + strand
LO, HI = P - 260, P + D + 260
OUT = os.path.join(HERE, "..", "..", "..", "outbox", "2026-11-14-asgard-protein-intron", "image.png")


def rna():
    joined, through = [], []
    with open(os.path.join(HERE, "out", "rna_mol_templates.tsv")) as f:
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


def dna():
    cls, tmpl = {}, {}
    path = os.path.join(HERE, "out", "sam", "SRR22273278.sam.zst")
    for pass_ in (1, 2):
        proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        for line in proc.stdout:
            if line[0] == "@":
                continue
            f = line.split("\t", 10)
            if f[2] != "CP104013.1":
                continue
            pos0, cigar, fr = int(f[3]) - 1, f[5], f[0].rsplit(".", 1)[0]
            end, ident, gaps, ops = parse(cigar, pos0)
            if pass_ == 1:
                if ident < MIN_ID or end <= P - FLANK or pos0 >= P + D + FLANK:
                    continue
                if len(gaps) == 1:
                    gs, gl, before = gaps[0]
                    if gs >= P - FLANK and gs + gl <= P + D + FLANK and before >= FLANK and len(ops) - before >= FLANK \
                            and ops[before - FLANK:before].count("x") <= FLANK_MM and ops[before:before + FLANK].count("x") <= FLANK_MM:
                        cls[fr] = "joined"
                elif not gaps and ((pos0 <= P - FLANK and end >= P + FLANK) or (pos0 <= P + D - FLANK and end >= P + D + FLANK)):
                    cls.setdefault(fr, "through")
            elif fr in cls:
                cl = CLIP.match(cigar)
                us, ue = pos0 - int(cl.group(1) or 0), end + int(cl.group(2) or 0)
                t = tmpl.get(fr)
                tmpl[fr] = (us, ue) if t is None else (min(t[0], us), max(t[1], ue))
        if proc.wait() != 0:
            sys.exit("zstd failed")
    with open(os.path.join(HERE, "out", "fig_l02_dna_templates.tsv"), "w") as fo:
        fo.write("fragment\tclass\tstart\tend\n")
        for fr in sorted(cls):
            fo.write(f"{fr}\t{cls[fr]}\t{tmpl[fr][0]}\t{tmpl[fr][1]}\n")
    return [tmpl[f] for f in cls if cls[f] == "joined"], [tmpl[f] for f in cls if cls[f] == "through"]


def main():
    rj, rt = rna()
    dj, dt = dna()
    print(f"RNA molecules: without the stretch {len(rj)}, through {len(rt)}; DNA fragments: without {len(dj)}, through {len(dt)}")
    assert (len(rj), len(rt), len(dj), len(dt)) == (385, 26, 0, 46), "counts differ from the post"
    exact = sum(1 for t in rj if (t[2], t[3]) == (P, D))
    anyread = next(int(r.split("\t")[4]) for r in open(os.path.join(HERE, "out", "rna_mol_gaps.tsv"))
                   if r.split("\t")[:3] == ["CP104013.1", str(P), str(D)])        # templates with a read at P+D
    print(f"RNA molecules with exactly {P}+{D}: {anyread} in at least one read, {exact} in the read molecules.py kept")
    assert anyread == 383, "differs from the post"

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyArrow

    ORANGE, GREY, INK, MUTED, BAND = "#eb6834", "#b4b2aa", "#0b0b0b", "#6f6e69", "#fbe3d6"
    plt.rcParams.update({"font.family": "Helvetica", "font.size": 15})
    rows_r = sorted(rj, key=lambda t: (t[0], t[1])) + sorted(rt, key=lambda t: (t[0], t[1]))
    rows_d = sorted(dt, key=lambda t: (t[0], t[1]))
    fig = plt.figure(figsize=(15, 9.4), dpi=200)
    gs = fig.add_gridspec(3, 1, height_ratios=[34, len(rows_r), len(rows_d) * 1.6], hspace=0.06,
                          left=0.06, right=0.80, top=0.86, bottom=0.10)
    axg, axr, axd = (fig.add_subplot(gs[i]) for i in range(3))
    for ax in (axg, axr, axd):
        ax.set_xlim(LO, HI)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_yticks([])
        ax.tick_params(axis="x", length=0, labelbottom=False)
    # genes
    axg.set_ylim(-1, 1)
    for name, s, e in GENES:
        a, b = max(s, LO), min(e, HI)
        axg.add_patch(FancyArrow(a, -0.35, b - a, 0, width=0.7, head_width=0.7, head_length=10 if e <= HI else 0,
                                 length_includes_head=True, color="#d8d6cf", lw=0, zorder=4))
        axg.text((a + b) / 2, -0.35, name, ha="center", va="center", fontsize=14, color=INK, zorder=5)
    axg.add_patch(Rectangle((P, -1), D, 1.3, color=ORANGE, lw=0, zorder=3))
    axg.text(P + D / 2, 0.62, "31 bases", ha="center", va="center", fontsize=14, color=INK)

    def draw(ax, rows, nrow, gapped):
        ax.set_ylim(nrow + 0.5, -0.5)
        ax.add_patch(Rectangle((P, -0.5), D, nrow + 1, color=BAND, lw=0, zorder=0))
        for i, t in enumerate(rows):
            s, e = max(t[0], LO), min(t[1], HI)
            if gapped and len(t) == 4:
                gs_, gl = t[2], t[3]
                ax.plot([s, gs_], [i, i], color=ORANGE, lw=1.1, solid_capstyle="butt")
                ax.plot([gs_ + gl, e], [i, i], color=ORANGE, lw=1.1, solid_capstyle="butt")
            else:
                ax.plot([s, e], [i, i], color=GREY, lw=1.1, solid_capstyle="butt")

    draw(axr, rows_r, len(rows_r), True)
    draw(axd, rows_d, len(rows_d), False)
    axd.tick_params(axis="x", labelbottom=True, labelsize=13, colors=MUTED)
    ticks = [P - 200, P - 100, P, P + D + 100, P + D + 200]
    axd.set_xticks(ticks)
    axd.set_xticklabels([f"{x:,}" for x in ticks])
    # labels to the right
    xr = HI + 8
    axr.text(xr, len(rj) / 2, f"RNA: {len(rj)} molecules\nlack the 31 bases", color=INK, va="center", fontsize=16, clip_on=False)
    axr.text(xr, len(rj) + len(rt) / 2, f"{len(rt)} read through", color=MUTED, va="center", fontsize=16, clip_on=False)
    axd.text(xr, len(rows_d) / 2, f"DNA: {len(dt)} fragments\nread through,\n0 lack them", color=MUTED, va="center", fontsize=16,
             clip_on=False)
    axg.text(xr, 0, "annotated genes", color=MUTED, va="center", fontsize=14, clip_on=False)
    fig.text(0.06, 0.95, "Asgard archaeon B35: 31 bases missing from a protein gene's RNA, not its DNA",
             fontsize=24, fontweight="bold", color=INK, va="top")
    fig.text(0.06, 0.905, "Ca. Lokiarchaeum ossiferum B35, CP104013.1; the shaded 31 bases hold the first gene's stop codon. "
             "Removing them joins the two genes in frame.", fontsize=15, color=MUTED, va="top")
    fig.text(0.06, 0.065, "RNA lines are molecules (read pairs with the same outer ends merged); DNA lines are fragments; "
             f"clipped to the window.\n{anyread} of the {len(rj)} RNA molecules lack exactly the predicted 31 bases. "
             "RNA: SRR27313660 and SRR32148620, first 10 M read pairs each; neither library is stranded.\n"
             "DNA: SRR22273278, first 8 M read pairs. Positions 0-based on CP104013.1.",
             fontsize=12.5, color=MUTED, va="top")
    fig.savefig(OUT, facecolor="white")
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
