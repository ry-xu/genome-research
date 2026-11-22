"""Re-join ARAGORN tRNA(UCU) calls with the AGA classes; list all G3:U70 UCU carriers among fungi."""
import numpy as np, csv
from collections import defaultdict
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
acc, tax = C['acc'], C['tax']
ala = set(acc[np.load(A + 'audit/ala_class_idx.npy')]); arg = set(acc[np.load(A + 'audit/arg_class_idx.npy')])
fung = set(acc[C['fungi']])
tr = defaultdict(list)
allacc = set()
for r in csv.DictReader(open(A + 'fungi.aragorn.tsv'), delimiter='\t'):
    allacc.add(r['acc'])
    if r['anticodon'].upper() == 'TCT':
        tr[r['acc']].append(r)
print('genomes in aragorn table', len(allacc), 'fungal census genomes', len(fung), 'census genomes absent from aragorn table', len(fung - allacc))
gu = {a for a, rs in tr.items() if any(r['pair3_70'] == 'G:T' for r in rs)}
print('ala class with G3U70 UCU:', len(ala & gu), 'of', len(ala), '; without:', sorted(ala - gu))
print('arg class with G3U70 UCU:', len(arg & gu), 'of', len(arg), '; arg-class genomes with no UCU tRNA at all:', len(arg - set(tr)))
print('arg class absent from aragorn table:', len(arg - allacc))
ti = {a: t for a, t in zip(acc, tax)}
print('\nall fungal G3:U70 UCU carriers (acc, class, taxon, pair, N73, stem5, stem3):')
for a in sorted(gu, key=lambda a: ti.get(a, '')):
    cl = 'ala' if a in ala else 'arg' if a in arg else 'other'
    for r in tr[a]:
        if r['pair3_70'] == 'G:T':
            print(a, cl, ';'.join(ti.get(a, '?').split(';')[5:8]), ti.get(a, '?').split(';')[-1], r['pair3_70'], r['N73'], r['stem5'], r['stem3'], r['stem_pairs'])
