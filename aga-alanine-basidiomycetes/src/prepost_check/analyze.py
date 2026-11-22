import sys, numpy as np, csv
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from math import comb
C = np.load('scratch/mycache.npz')
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
sg, sc, col, strict, cons = C['sg'], C['sc'], C['col'], C['strict'], C['cons']
acc, tax, table = C['acc'].tolist(), C['tax'].tolist(), C['table']; K = int(C['K']); models = C['models'].tolist()
G = len(acc)
fungi = np.array(['Fungi' in t.split(';') for t in tax])
print('fungal genomes in census:', fungi.sum())
st = strict[col]; cc = cons[col]
A, R = AA.index('A'), AA.index('R')
AGA = COD.index('AGA'); AGG = COD.index('AGG'); CGN = [COD.index(c) for c in ('CGT','CGC','CGA','CGG')]; GCN=[COD.index(c) for c in ('GCT','GCC','GCA','GCG')]
SYN = CGN + [AGG]
isAGA = sc == AGA
a_ala = np.bincount(sg[isAGA & st & (cc == A)], minlength=G)
a_arg = np.bincount(sg[isAGA & st & (cc == R)], minlength=G)
ala = [g for g in range(G) if fungi[g] and a_ala[g] >= 1 and a_arg[g] == 0]
arg = [g for g in range(G) if fungi[g] and a_arg[g] >= 3 and a_ala[g] == 0]
both = [g for g in range(G) if fungi[g] and a_ala[g] >= 1 and a_arg[g] >= 1]
print('ala class', len(ala), '; arg class', len(arg), '; both', [(acc[g], a_ala[g], a_arg[g], tax[g].split(";")[-1]) for g in both if a_ala[g]>=2])
print('genomes with any AGA at strict Ala (fungi):', int((a_ala[fungi] >= 1).sum()))
# species names from records
org = {}
for f in ('mito', 'gbmito', 'gbmitoorf'):
    for r in csv.reader(open(f'lab/work/new-stops-euk/data/{f}.records.tsv'), delimiter='\t'):
        org[r[0]] = r[1]
LIN = ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales')
lin = {g: next((k for k in LIN if k in tax[g].split(';')), 'other') for g in ala}
print('per lineage genomes:', Counter(lin.values()))
names = sorted(org.get(acc[g], '?') for g in ala)
print('organism names (%d unique):' % len(set(names)))
for g in sorted(ala, key=lambda g: (lin[g], org.get(acc[g]))):
    print('  ', acc[g], lin[g], org.get(acc[g]), 'AGA@A', a_ala[g], 'table', table[g])
print('per-genome AGA@Ala range', min(a_ala[ala]), max(a_ala[ala]), 'total', a_ala[ala].sum())
m32 = np.isin(sg, ala)
sel = m32 & isAGA & st
print('AGA at strict columns in 32:', sel.sum(), Counter(AA[i] for i in cc[sel]).most_common())
print('  Ala share', (cc[sel] == A).mean())
dcols = set(col[sel & (cc == A)].tolist())
print('distinct strict Ala columns with AGA over all 32:', len(dcols))
perlin = {l: set(col[np.isin(sg, [g for g in ala if lin[g]==l]) & isAGA & st & (cc == A)].tolist()) for l in LIN}
print('distinct per lineage:', {l: len(v) for l, v in perlin.items()}, 'sum', sum(len(v) for v in perlin.values()))
# per lineage distinct-column Fisher
for l in LIN:
    gs = [g for g in ala if lin[g] == l]
    mm = np.isin(sg, gs) & st
    def dc(cm, r): return len(set(col[mm & cm & (cc == r)].tolist()))
    t = [dc(isAGA, A), dc(isAGA, R), dc(np.isin(sc, SYN), A), dc(np.isin(sc, SYN), R)]
    # also synonyms at AGA's Ala columns specifically
    syn_at_aga_cols = len(set(col[mm & np.isin(sc, SYN) & np.isin(col, list(perlin[l]))].tolist()))
    p = fisher_exact([[t[0], t[1]], [t[2], t[3]]], alternative='greater')[1]
    # site-level pooled (not distinct)
    s = [int((mm & isAGA & (cc == A)).sum()), int((mm & isAGA & (cc == R)).sum()), int((mm & np.isin(sc, SYN) & (cc == A)).sum()), int((mm & np.isin(sc, SYN) & (cc == R)).sum())]
    ps = fisher_exact([[s[0], s[1]], [s[2], s[3]]], alternative='greater')[1]
    print(f'{l} n={len(gs)} distinct [AGA@A, AGA@R, syn@A, syn@R] = {t} p={p:.2g}; syn at AGA\'s own Ala cols={syn_at_aga_cols}; site-level {s} p={ps:.2g}')
# AGG in the 32
selg = m32 & (sc == AGG) & st
print('AGG at strict cols in 32:', selg.sum(), 'distinct', len(set(col[selg].tolist())), Counter(AA[i] for i in cc[selg]))
# all columns
for name, cods in (('AGA', [AGA]), ('CGN', CGN), ('GCN', GCN), ('AGG', [AGG])):
    s = m32 & np.isin(sc, cods)
    c = Counter(AA[i] for i in cc[s]); n = s.sum()
    s2 = s & ~st
    c2 = Counter(AA[i] for i in cc[s2]); n2 = s2.sum()
    print(f'all cols {name}: n={n}', ', '.join(f'{a} {v/n:.3f}' for a, v in c.most_common(6)), f'R={c["R"]/max(n,1):.3f} ({c["R"]})',
          f'| non-strict only n={n2}:', ', '.join(f'{a} {v/n2:.3f}' for a, v in c2.most_common(5)), f'R={c2["R"]/max(n2,1):.3f}')
# same-class arg genomes: classes Tremellomycetes/Exobasidiomycetes/Cystobasidiomycetes with AGA >=3 strict Arg
CL = ('Tremellomycetes', 'Exobasidiomycetes', 'Cystobasidiomycetes')
sc_arg = [g for g in arg if any(k in tax[g].split(';') for k in CL)]
print('same-class arg genomes:', [(acc[g], tax[g].split(';')[-1]) for g in sc_arg])
s = np.isin(sg, sc_arg) & isAGA
c = Counter(AA[i] for i in cc[s]); print('  AGA all cols', s.sum(), {k: round(v / s.sum(), 3) for k, v in c.most_common(3)})
s = np.isin(sg, sc_arg) & np.isin(sc, CGN)
c = Counter(AA[i] for i in cc[s]); print('  CGN all cols', s.sum(), {k: round(v / s.sum(), 3) for k, v in c.most_common(3)})
# background: across all fungal genomes, arginine codons at strict Ala columns (selection-free rate)
fg = np.flatnonzero(fungi)
mf = np.isin(sg, fg) & st & (cc == A)
non32 = mf & ~m32
for name, cods in (('AGA', [AGA]), ('CGN', CGN), ('AGG', [AGG])):
    x = non32 & np.isin(sc, cods)
    print(f'outside the 32: {name} at strict Ala: {x.sum()} sites in {len(set(sg[x].tolist()))} genomes')
# all non-GCN codons at strict Ala in the 32
x = m32 & st & (cc == A) & ~np.isin(sc, GCN)
print('non-GCN codons at strict Ala in 32:', Counter(COD[i] for i in sc[x]).most_common(8))
x = m32 & st & np.isin(sc, GCN)
print('GCN at strict cols in 32:', x.sum(), Counter(AA[i] for i in cc[x]).most_common(5))
np.save('scratch/ala32.npy', np.array(ala))
np.save('scratch/arg.npy', np.array(arg))
