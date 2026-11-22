"""Recompute the remaining README lead numbers: all-columns consensus (32 and same-class Arg genomes), commonest
non-GCN codon at strict Ala columns per class, T. asahii CP116791.1 AGA counts."""
import numpy as np
from collections import Counter
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
acc, tax = C['acc'], C['tax']; sg, sc, col = C['sg'], C['sc'], C['col']; strict, cons = C['strict'], C['cons']
G = len(acc); st = strict[col]; cc = cons[col]
iA, iR = AA.index('A'), AA.index('R'); AGA = COD.index('AGA')
SYN = [COD.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGG')]
GCN = [COD.index(c) for c in ('GCT', 'GCC', 'GCA', 'GCG')]; CGN = [COD.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG')]
ala = np.load(A + 'audit/ala_class_idx.npy'); arg = np.load(A + 'audit/arg_class_idx.npy')
def allcols(gs, name):
    m = np.isin(sg, gs); print(f'== {name}: {len(gs)} genomes')
    for k, cods in (('AGA', [AGA]), ('CGN', CGN), ('GCN', GCN)):
        s = m & np.isin(sc, cods); n = int(s.sum()); c = Counter(AA[i] for i in cc[s])
        print(f'  {k}: {n} sites; ' + ', '.join(f'{a} {v / n:.3f}' for a, v in c.most_common(6)) + f' | R {c["R"]} sites')
allcols(ala, 'ala-class 32')
CL = ('Tremellomycetes', 'Exobasidiomycetes', 'Cystobasidiomycetes')
same = [g for g in arg if any(k in tax[g].split(';') for k in CL)]
print('same-class arg genomes:', [(acc[g], tax[g].split(';')[-1]) for g in same])
allcols(same, 'same-class arg-class')
a = np.bincount(sg[(sc == AGA) & st & (cc == iA)], minlength=G); r = np.bincount(sg[(sc == AGA) & st & (cc == iR)], minlength=G)
both = [g for g in np.flatnonzero(C['fungi']) if r[g] >= 3 and any(k in tax[g].split(';') for k in CL)]
print('same-class genomes with >=3 AGA@strict-R (any Ala):', [(acc[g], tax[g].split(';')[-1], a[g], r[g]) for g in both])
if set(both) != set(same): allcols(both, 'same-class >=3 AGA@R incl. Ala-users')
print('== non-GCN codons at strict Ala columns, all genomes per class (top 4) and AGA count')
fung = np.flatnonzero(C['fungi'])
for cl in CL + ('Ustilaginomycetes', 'Agaricomycetes', 'Microbotryomycetes'):
    gs = [g for g in fung if cl in tax[g].split(';')]
    s = np.isin(sg, gs) & st & (cc == iA) & ~np.isin(sc, GCN)
    c = Counter(COD[i] for i in sc[s]); print(f'  {cl} ({len(gs)} g): {c.most_common(4)}  AGA={c["AGA"]}')
for lin in ('Trichosporonales', 'Tremellales', 'Agaricales', 'Ustilaginales'):
    gs = [g for g in fung if lin in tax[g].split(';')]
    s = np.isin(sg, gs) & st & (cc == iA) & ~np.isin(sc, GCN)
    c = Counter(COD[i] for i in sc[s]); print(f'  {lin} ({len(gs)} g): {c.most_common(4)}  AGA={c["AGA"]}')
for q in ('CP116791.1',):
    g = [i for i in range(G) if acc[i] == q]
    for gi in g:
        s = (sg == gi) & (sc == AGA)
        print(q, tax[gi].split(';')[-1], 'AGA aligned', int(s.sum()), 'strict', Counter(AA[i] for i in cc[s & st]))
