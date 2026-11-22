"""Pseudo-replication check: distinct strict columns with AGA (by consensus residue) per lineage, and the class
contrast recomputed counting each (column, codon class) once per lineage."""
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
acc, tax = C['acc'], C['tax']; sg, sc, col = C['sg'], C['sc'], C['col']; K = int(C['K']); models = C['models']
st = C['strict'][col]; cc = C['cons'][col]
ala = np.load(A + 'audit/ala_class_idx.npy')
LIN = ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales')
lin = {g: next(k for k in LIN if k in tax[g].split(';')) for g in ala}
AGA = COD.index('AGA'); SYN = [COD.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGG')]
m = np.isin(sg, ala) & st
out = defaultdict(lambda: defaultdict(set))
for i in np.flatnonzero(m & (sc == AGA)):
    out[lin[sg[i]]][AA[cc[i]]].add(int(col[i]))
tot = Counter()
for l in LIN:
    d = out[l]
    ng = sum(1 for g in ala if lin[g] == l)
    print(f'{l} ({ng} genomes): distinct strict columns with AGA by residue:', {k: len(v) for k, v in sorted(d.items(), key=lambda x: -len(x[1]))})
    for k, v in d.items(): tot[k] += len(v)
allc = defaultdict(set)
for l in LIN:
    for k, v in out[l].items(): allc[k] |= v
print('sum over lineages:', dict(tot), '; distinct over all 32:', {k: len(v) for k, v in allc.items()})
# per-lineage distinct-column 2x2 (AGA@A, AGA@R, syn@A, syn@R)
for l in LIN:
    gs = [g for g in ala if lin[g] == l]
    mm = np.isin(sg, gs) & st
    def dc(codmask, res): return len(set(col[mm & codmask & (cc == AA.index(res))].tolist()))
    t = [dc(sc == AGA, 'A'), dc(sc == AGA, 'R'), dc(np.isin(sc, SYN), 'A'), dc(np.isin(sc, SYN), 'R')]
    print(l, 'distinct columns AGA@A, AGA@R, syn@A, syn@R:', t, 'p=%.2g' % fisher_exact([[t[0], t[1]], [t[2], t[3]]], 'greater')[1])
# list the Ala columns per lineage with model names
for l in LIN:
    print(l, sorted(f'{models[k // K]}:{k % K}' for k in out[l]['A']))
