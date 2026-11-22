"""Figure: in Porphyridium aerugineum mitochondria, what codon the second isolate (CCMP1948) has where the first
(SAG 1380-2) has AGG, ATG or another arginine codon (counts from codon_pair.py, results/agg_met/isolate_pair.txt).
usage: python fig_agg.py out.png"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, GREY, INK, INK2, ORANGE, LIGHT = '#2a78d6', '#a3a29c', '#0b0b0b', '#52514e', '#eb6834', '#dddcd6'
# rows: SAG 1380-2 codon class -> (same codon, ATG<->AGG swap, another Arg codon, anything else)
ROWS = [('AGG\n(an arginine codon\nin the standard code)', 37, 17, 0, 4),
        ('ATG\n(methionine)', 42, 16, 0, 13),
        ('AGA, CGN\n(the other arginine\ncodons)', 41, 0, 34, 3)]
SEGS = [('same codon', GREY), ('ATG ↔ AGG', ORANGE), ('another arginine codon', BLUE), ('anything else', LIGHT)]


def main(out):
    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=200)
    fig.patch.set_facecolor('#fcfcfb'); ax.set_facecolor('#fcfcfb')
    for i, (lab, *v) in enumerate(ROWS):
        y, n, left = len(ROWS) - 1 - i, sum(v), 0
        for k, (c, (name, col)) in enumerate(zip(v, SEGS)):
            if c:
                ax.barh(y, c / n, left=left, height=0.56, color=col, edgecolor='#fcfcfb', lw=2)
                if c / n > 0.06:
                    ax.text(left + c / n / 2, y, str(c), ha='center', va='center', fontsize=10,
                            color='#fcfcfb' if col in (BLUE, ORANGE) else INK)
            left += c / n
        ax.text(1.015, y, f'{n} sites', va='center', fontsize=9, color=INK2)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([r[0] for r in ROWS][::-1], fontsize=9.5, color=INK)
    # direct segment labels above the top bar
    top = len(ROWS) - 1
    v = ROWS[0][1:]; n = sum(v)
    ax.text(v[0] / n / 2, top + 0.42, 'same codon', ha='center', fontsize=9, color=INK2)
    ax.text((v[0] + v[1] / 2) / n, top + 0.42, 'ATG ↔ AGG swap', ha='center', fontsize=9, color=ORANGE)
    v = ROWS[2][1:]; n = sum(v)
    ax.text((v[0] + v[2] / 2) / n, 0 - 0.47, 'another arginine codon', ha='center', va='top', fontsize=9, color=BLUE)
    ax.text(1 - v[3] / n / 2, 0 - 0.47, 'other', ha='right', va='top', fontsize=9, color=INK2)
    ax.text(0.995, top - 0.36, '0 to another arginine codon', ha='right', va='top', fontsize=9, color=BLUE)
    ax.text(-0.02, top + 0.62, 'codon in isolate SAG 1380-2', ha='right', va='bottom', fontsize=9, color=INK2)
    ax.set_xlim(0, 1); ax.set_ylim(-0.75, top + 0.75)
    ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
    ax.set_xlabel('share of aligned sites, by the codon in the other isolate (CCMP1948)', color=INK2)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color('#c3c2b7')
    ax.tick_params(colors=INK2, left=False)
    fig.suptitle('Two isolates of the red alga Porphyridium aerugineum: in mitochondrial genes, AGG trades\n'
                 'places with ATG (methionine), never with an arginine codon', x=0.01, ha='left', color=INK, fontsize=11.5)
    fig.text(0.01, 0.01, 'Codon at each aligned position of 14 shared mitochondrial genes, isolate SAG 1380-2 '
             '(NC_062302.1) vs CCMP1948 (MT483995.1).', fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 0.97, 1))
    fig.savefig(out, facecolor=fig.get_facecolor())


if __name__ == '__main__':
    main(*sys.argv[1:])
