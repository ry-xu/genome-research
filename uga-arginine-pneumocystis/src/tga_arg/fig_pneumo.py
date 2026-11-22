"""Figure: per Pneumocystis mitogenome, UGA codons at conserved Arg positions (bars) next to the regular Arg codons
there; UGA at conserved Trp is 0 in all. usage: python fig_pneumo.py pneumo_cols.tsv out.png"""
import sys, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, GREY, INK, INK2, BG = '#2a78d6', '#c9c8c1', '#0b0b0b', '#52514e', '#fcfcfb'
NAME = {'NC_020331.1': 'P. jirovecii (human)', 'NC_013660.2': 'P. carinii (rat)', 'NC_020332.1': 'P. murina (mouse)',
        'NC_060318.1': 'P. wakefieldiae (rat)', 'NC_060319.1': 'P. oryctolagi (rabbit)', 'NC_060320.1': 'P. canis (dog)',
        'MT726214.1': "P. sp. 'macacae' (macaque)", 'MT726212.1': "P. sp. 'ludovicianus' (prairie dog)", 'PP484879.1': 'P. sp. SQ178 (squirrel)'}


def main(tsv, out):
    rows = {r['acc']: r for r in csv.DictReader(open(tsv), delimiter='\t')}
    order = [a for a in NAME if a in rows][::-1]
    fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=200)
    fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
    for i, a in enumerate(order):
        r = rows[a]; u, s = int(r['TGA_strictR']), int(r['Arg_strictR'])
        ax.barh(i, u, height=.62, color=BLUE)
        ax.barh(i, s, left=u, height=.62, color=GREY)
        ax.text(u + s + .6, i, f'{u}', va='center', color=BLUE if u else INK2, fontsize=9, fontweight='bold' if u else None)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([NAME[a] for a in order], fontsize=9, color=INK, style='italic')
    ax.set_xlim(0, 56)
    ax.set_xlabel('codons at conserved arginine positions', color=INK2)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color('#c3c2b7'); ax.tick_params(colors=INK2, left=False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=BLUE, label='UGA (tryptophan in table 4, the code these fungi are filed under)'),
                       Patch(color=GREY, label='regular arginine codons (CGN, AGA, AGG)')],
              loc='lower left', bbox_to_anchor=(0, 1.0), ncol=1, frameon=False, fontsize=8.5, labelcolor=INK2)
    ax.set_title('UGA replaces arginine codons in human Pneumocystis', loc='left', pad=34,
                 color=INK, fontsize=12)
    fig.text(0.01, 0.01, 'Conserved = one amino acid at >= 90% of aligned mitochondrial genomes (weighted by family). '
             '\nUGA sits at conserved tryptophan positions 0 times in all 9 genomes, where tryptophan is UGG (35-37 per genome).',
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, facecolor=fig.get_facecolor(), bbox_inches='tight', pad_inches=0.15)


if __name__ == '__main__':
    main(*sys.argv[1:])
