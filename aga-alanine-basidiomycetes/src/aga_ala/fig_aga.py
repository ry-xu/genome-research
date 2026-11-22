"""Figure: per fungal mitochondrial genome, AGA codons at strict Ala vs strict Arg columns, coloured by whether the
genome has a tRNA(UCU) with the alanine identity pair G3:U70 (ARAGORN). usage: python fig_aga.py trna_vs_usage.tsv out.png"""
import sys, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BLUE, GREY, INK, INK2 = '#2a78d6', '#a3a29c', '#0b0b0b', '#52514e'


def main(tsv, out):
    rows = list(csv.DictReader(open(tsv), delimiter='\t'))
    rng = np.random.default_rng(1)
    x = np.array([int(r['AGA_at_Arg']) for r in rows], float)
    y = np.array([int(r['AGA_at_Ala']) for r in rows], float)
    g = np.array([int(r['G3U70']) > 0 for r in rows])
    jx, jy = x + rng.uniform(-.3, .3, len(x)), y + rng.uniform(-.22, .22, len(y))
    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=200)
    fig.patch.set_facecolor('#fcfcfb'); ax.set_facecolor('#fcfcfb')
    ax.scatter(jx[~g], jy[~g], s=9, c=GREY, alpha=.35, lw=0, label='no such tRNA found')
    ax.scatter(jx[g], jy[g], s=46, c=BLUE, edgecolors='#fcfcfb', lw=1.2, zorder=3,
               label='has a tRNA for AGA with the alanine mark G3:U70')
    ax.legend(loc='upper right', frameon=False, labelcolor=INK2, fontsize=9, markerscale=1.2)
    ax.set_xlim(-1.2, 50); ax.set_ylim(-0.8, 10)
    ax.set_xlabel('AGA codons at conserved arginine positions (per genome)', color=INK2)
    ax.set_ylabel('AGA codons at conserved alanine positions', color=INK2)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color('#c3c2b7')
    ax.tick_params(colors=INK2)
    na = sum(r['class'] == 'ala' for r in rows); nag = sum(r['class'] == 'ala' and int(r['G3U70']) > 0 for r in rows)
    nr = sum(r['class'] == 'arg' for r in rows); nrg = sum(r['class'] == 'arg' and int(r['G3U70']) > 0 for r in rows)
    ax.annotate(f'{na} genomes put AGA only at alanine positions\n(Trichosporon and relatives, Tremella, Tilletia,\nCystobasidium and a few others);\n'
                f'{nag} of them carry a tRNA for AGA with the alanine\nidentity mark (G3:U70)', xy=(0.6, 6.5), xytext=(6, 7.2),
                color=BLUE, fontsize=10, va='center', arrowprops=dict(arrowstyle='-', color=BLUE, lw=1))
    ax.annotate(f'{nr:,} genomes put AGA at arginine positions (3 or more)\nand never at alanine; {nrg} carry that tRNA',
                xy=(20, 0.25), xytext=(17, 2.6), color=INK2, fontsize=10, va='center',
                arrowprops=dict(arrowstyle='-', color=INK2, lw=1))
    ax.set_title('Fungal mitochondria: where AGA, an arginine codon in fungi, sits in proteins', loc='left',
                 color=INK, fontsize=12)
    fig.text(0.01, 0.01, f'{len(rows):,} fungal mitochondrial genomes; conserved = one amino acid at >= 90% of aligned '
             'mitochondrial genomes (weighted by family). Points jittered.', fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, facecolor=fig.get_facecolor())


if __name__ == '__main__':
    main(*sys.argv[1:])
