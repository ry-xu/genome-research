"""Per-genome AGA usage (all aligned / strict by residue), CGN, AGG and tRNA(UCU) calls for a taxon."""
import sys, csv, numpy as np
from collections import Counter, defaultdict
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
acc, tax = C['acc'], C['tax']; sg, sc, col = C['sg'], C['sc'], C['col']
st = C['strict'][col]; cc = C['cons'][col]
tr = defaultdict(list)
for r in csv.DictReader(open(A + 'fungi.aragorn.tsv'), delimiter='\t'):
    if r['anticodon'] == 'TCT': tr[r['acc']].append(r['pair3_70'] + '/' + r['N73'])
G = len(acc)
tot = np.bincount(sg, minlength=G)
aga = sc == COD.index('AGA'); cgn = np.isin(sc, [COD.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG')])
n_aga = np.bincount(sg[aga], minlength=G); n_cgn = np.bincount(sg[cgn], minlength=G)
n_agg = np.bincount(sg[sc == COD.index('AGG')], minlength=G)
idx = np.flatnonzero(aga & st)
by = defaultdict(Counter)
for i in idx: by[sg[i]][AA[cc[i]]] += 1
for taxon in sys.argv[1].split(','):
    gs = [g for g in range(G) if taxon in tax[g].split(';')]
    print(f'== {taxon} ({len(gs)} genomes) acc genus aligned AGA AGG CGN AGA/kcod strictAGA tRNA(UCU)')
    for g in sorted(gs, key=lambda g: tax[g]):
        print(acc[g], tax[g].split(';')[-1][:18], tot[g], n_aga[g], n_agg[g], n_cgn[g], f'{1000 * n_aga[g] / max(tot[g], 1):.1f}', dict(by[g]), tr.get(acc[g], []))
