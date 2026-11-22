import numpy as np
from collections import Counter
S = 'scratch/'
C = np.load(S + 'mycache.npz')
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
sg, sc, col, strict, cons = C['sg'], C['sc'], C['col'], C['strict'], C['cons']
tax = C['tax'].tolist(); ala = np.load(S + 'ala32.npy'); st = strict[col]; cc = cons[col]
A = AA.index('A'); AGA = COD.index('AGA'); GCN = [COD.index(c) for c in ('GCT','GCC','GCA','GCG')]
LIN = ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales')
ta = tn = ga = gn = 0
for l in LIN:
    gs = [g for g in ala if l in tax[g].split(';')]
    m = np.isin(sg, gs) & st
    x = m & (sc == AGA)
    a = len(set(col[x & (cc == A)].tolist())); n = len(set(col[x & (cc != A)].tolist()))
    y = m & np.isin(sc, GCN)
    a2 = len(set(col[y & (cc == A)].tolist())); n2 = len(set(col[y & (cc != A)].tolist()))
    ta += a; tn += n; ga += a2; gn += n2
    print(f'{l}: AGA distinct Ala {a}, non-Ala {n} | GCN distinct Ala {a2}, non-Ala {n2}')
print(f'sum per lineage: AGA Ala share {ta/(ta+tn):.2f} ({ta}/{ta+tn}); GCN {ga/(ga+gn):.2f} ({ga}/{ga+gn})')
