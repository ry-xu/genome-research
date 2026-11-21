"""Figure: which codon spells proline at each strict Pro column, per genome, next to the placement tree.

usage: python pro_grid.py accs.txt names.fa tree.nwk|- outgroup out.png census.npz [census2.npz ...]
accs.txt: one accession per line (the genomes to show); names.fa: FASTA whose headers are '>acc Genus species';
tree '-' orders rows by name. Columns: strict Pro columns where >= 60% of the shown genomes have a site
and at least half of those sites are CCN or CTN (drops columns where the shown clade has another residue). Cells:
CCN, CTN, other codon, or no site. Also writes out.png's .tsv with the per-genome counts shown.
env HIGHLIGHT: comma-separated accessions whose labels are drawn bold in the CTN colour; TITLE: figure title.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
sys.path.insert(0, 'src/codes_census')
from sitescan_org import load, AA, CODONS

accf, namef, treef, outgroup, out = sys.argv[1:6]
D = load(sys.argv[6:])
want = [l.strip() for l in open(accf) if l.strip()]
names = {}
for l in open(namef):
    if l.startswith('>'):
        a, _, n = l[1:].strip().partition(' ')
        names[a] = n
gi = {a: i for i, a in enumerate(D['acc'])}
want = [a for a in want if a in gi]
if treef != '-':
    from Bio import Phylo
    tree = Phylo.read(treef, 'newick')
    tree.root_with_outgroup(outgroup)
    tree.ladderize()
    order = [c.name for c in tree.get_terminals()]
    want = [a for a in order if a in gi]
else:
    tree = None
    want.sort(key=lambda a: (names.get(a, ''), a))

P = AA.index('P')
cls = np.full(64, 2)
for i, c in enumerate(CODONS):
    if c.startswith('CC'):
        cls[i] = 0
    elif c.startswith('CT'):
        cls[i] = 1
col, sg, sc = D['col'], D['sg'], D['sc']
pro = D['strict'] & (D['cons'] == P)
G = np.array([gi[a] for a in want])
sel = np.isin(sg, G) & pro[col]
cols = np.unique(col[sel])
have = np.zeros((len(G), len(cols)), bool)
grid = np.full((len(G), len(cols)), 3)
ci = {c: j for j, c in enumerate(cols)}
ri = {g: r for r, g in enumerate(G)}
for g, c, k in zip(sg[sel], col[sel], sc[sel]):
    grid[ri[g], ci[c]] = cls[k]
    have[ri[g], ci[c]] = True
spelled = (grid <= 1).sum(0)
keep = (have.mean(0) >= 0.6) & (spelled >= 0.5 * have.sum(0))
grid = grid[:, keep]
order = np.lexsort((cols[keep],))
grid = grid[:, order]

with open(out.rsplit('.', 1)[0] + '.tsv', 'w') as fh:
    fh.write('acc\tname\tCCN\tCTN\tother\tno_site\n')
    for r, a in enumerate(want):
        n = np.bincount(grid[r], minlength=4)
        fh.write(f'{a}\t{names.get(a, "")}\t{n[0]}\t{n[1]}\t{n[2]}\t{n[3]}\n')

H = 0.22 * len(want) + 1.4
fig = plt.figure(figsize=(13, H))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.1], wspace=0.02)
cmap = ListedColormap(['#3b6fb6', '#e07b1a', '#9a9a9a', '#f2f2f2'])
ax2 = fig.add_subplot(gs[1])
ax2.imshow(grid, aspect='auto', cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
ax2.set_yticks([])
ax2.set_xticks([])
prots = {D['models'][m].split('_')[-1] for m in (cols[keep] // D['K']).tolist()}  # cox1 has two profile models
ax2.set_xlabel(f'{grid.shape[1]} of the {len(cols)} positions in {len(prots)} mitochondrial proteins where '
               f'≥ 90% of surveyed mitogenomes (weighted by family) have proline,\n'
               f'shown where ≥ 60% of these {len(want)} genomes have an aligned codon and at least half of those codons are CCN or CUN',
               fontsize=9)
for s in ax2.spines.values():
    s.set_visible(False)
ax1 = fig.add_subplot(gs[0])
labels = [names.get(a, a) for a in want]
HL = {names.get(a, a) for a in os.environ.get('HIGHLIGHT', '').split(',') if a}
if tree is not None:
    lab = dict(zip(want, labels))
    for c in tree.get_terminals():
        c.name = lab.get(c.name, c.name)
    for c in tree.get_nonterminals():
        c.confidence = None
        c.name = None
    Phylo.draw(tree, axes=ax1, do_show=False, label_func=lambda c: c.name if c.is_terminal() else None)
    ax1.set_ylim(len(want) + 0.5, 0.5)
    ax1.axis('off')
    for t in ax1.texts:
        t.set_fontsize(8)
        t.set_fontstyle('italic')
        if t.get_text().strip() in HL:
            t.set_fontweight('bold'); t.set_color('#c0610a')
    # widen the tree panel's x range so the longest tip label ends inside it, not on the grid
    fig.canvas.draw()
    rend, W = fig.canvas.get_renderer(), ax1.get_window_extent().width
    x0, R = ax1.get_xlim()
    for t in ax1.texts:
        R = max(R, x0 + (t.get_position()[0] - x0) * W / (W - t.get_window_extent(rend).width - 6))
    ax1.set_xlim(x0, R)
else:
    ax1.set_ylim(len(want) - 0.5, -0.5)
    for r, l in enumerate(labels):
        ax1.text(0.98, r, l, ha='right', va='center', fontsize=8, fontstyle='italic')
    ax1.axis('off')
from matplotlib.patches import Patch
ax2.legend(handles=[Patch(color='#3b6fb6', label='CCN (proline in every code)'),
                    Patch(color='#e07b1a', label='CUU or CUA (leucine in the standard code)'),
                    Patch(color='#9a9a9a', label='other codon'), Patch(color='#f2f2f2', label='no aligned site')],
           loc='upper center', bbox_to_anchor=(0.5, 1.0 + 1.0 / H), ncol=2, fontsize=8, frameon=False)
if os.environ.get('TITLE'):
    fig.suptitle(os.environ['TITLE'], x=0.12, y=0.975, ha='left', fontsize=13)
fig.savefig(out, dpi=150, bbox_inches='tight')
print(grid.shape, file=sys.stderr)
