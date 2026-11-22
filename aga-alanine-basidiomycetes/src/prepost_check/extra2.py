import numpy as np, csv
from collections import Counter
from scipy.stats import fisher_exact
S = 'scratch/'
C = np.load(S + 'mycache.npz')
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
sg, sc, col, strict, cons = C['sg'], C['sc'], C['col'], C['strict'], C['cons']
acc, tax = C['acc'].tolist(), C['tax'].tolist(); K = int(C['K']); models = C['models'].tolist(); G = len(acc)
ai = {a: i for i, a in enumerate(acc)}
ala = np.load(S + 'ala32.npy'); st = strict[col]; cc = cons[col]
A, R, Gly = AA.index('A'), AA.index('R'), AA.index('G'); AGA = COD.index('AGA'); GCN = [COD.index(c) for c in ('GCT','GCC','GCA','GCG')]
print('strict columns total', strict.sum(), '; strict Ala columns', int((strict & (cons == A)).sum()), '; strict Arg', int((strict & (cons == R)).sum()))
g = ai['QB252587.1']
m = (sg == g) & (sc == AGA)
print('Sporobolomyces QB252587.1: AGA aligned', m.sum(), 'at strict', (m & st).sum(), Counter(AA[i] for i in cc[m & st]), 'all-col consensus', Counter(AA[i] for i in cc[m]).most_common(4))
spor = [i for i in range(G) if 'Sporidiobolales' in tax[i].split(';')]
mm = np.isin(sg, spor) & (sc == AGA)
print('Sporidiobolales (21 genomes): AGA aligned', mm.sum(), 'at strict', (mm & st).sum(), Counter(AA[i] for i in cc[mm & st]))
# non-Ala strict sites of AGA in 32: by column and lineage
m32 = np.isin(sg, ala) & (sc == AGA) & st & (cc != A)
print('AGA at non-Ala strict columns: distinct cols', len(set(col[m32].tolist())), Counter((AA[c], models[k // K].split('_')[-1] + ':' + str(k % K)) for c, k in zip(cc[m32], col[m32])).most_common(12))
a = np.isin(sg, ala) & st
x = [[int((a & (sc == AGA) & (cc == Gly)).sum()), int((a & (sc == AGA) & (cc == A)).sum())], [int((a & np.isin(sc, GCN) & (cc == Gly)).sum()), int((a & np.isin(sc, GCN) & (cc == A)).sum())]]
print('Gly vs Ala, AGA vs GCN at strict cols:', x, 'p', fisher_exact(x)[1])
nonA = [[int((a & (sc == AGA) & (cc != A)).sum()), int((a & (sc == AGA) & (cc == A)).sum())], [int((a & np.isin(sc, GCN) & (cc != A)).sum()), int((a & np.isin(sc, GCN) & (cc == A)).sum())]]
print('non-Ala vs Ala, AGA vs GCN:', nonA, 'GCN Ala share', nonA[1][1] / sum(nonA[1]))
# figure range check: max AGA@Arg in fungi
fung = np.array(['Fungi' in t.split(';') for t in tax])
agaR = np.bincount(sg[(sc == AGA) & st & (cc == R)], minlength=G)
print('max AGA@Arg (strict, sites) among fungi', agaR[fung].max(), '; genomes > 50:', int((agaR[fung] > 50).sum()))
