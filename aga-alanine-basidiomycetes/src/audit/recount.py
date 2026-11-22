"""Recompute the README lead numbers from the cached census (independent of codon_contrast.py / trna_vs_usage.py)."""
import numpy as np, csv
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
AA = 'ACDEFGHIKLMNPQRSTVWY'
B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
acc, tax, src = C['acc'], C['tax'], C['src']
sg, sc, col = C['sg'], C['sc'], C['col']
strict, cons = C['strict'], C['cons']
G = len(acc)
st = strict[col]; cc = cons[col]
AGA = COD.index('AGA'); SYN = [COD.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGG')]
iA, iR = AA.index('A'), AA.index('R')
def cnt(mask): return np.bincount(sg[mask], minlength=G)
a = cnt((sc == AGA) & st & (cc == iA)); r = cnt((sc == AGA) & st & (cc == iR))
sa = cnt(np.isin(sc, SYN) & st & (cc == iA)); sr = cnt(np.isin(sc, SYN) & st & (cc == iR))
fung = np.flatnonzero(C['fungi'])
# compare with author's per-genome table
auth = {}
with open(A + 'aga_ala_Fungi.tsv') as f:
    next(f)
    for row in csv.DictReader(f, delimiter='\t'):
        if row['acc'] != 'pooled':
            auth.setdefault(row['acc'], []).append(tuple(int(row[k]) for k in ('codon@Y', 'codon@STD', 'syn@Y', 'syn@STD')))
dups = Counter(acc[fung])
print('fungal genomes', len(fung), 'distinct acc', len(dups), 'acc in >1 census file:', [(k, v) for k, v in dups.items() if v > 1][:10])
mism = 0
for g in fung:
    mine = (a[g], r[g], sa[g], sr[g])
    if tuple(mine) not in auth.get(acc[g], []):
        mism += 1
print('genomes whose 4 counts differ from aga_ala_Fungi.tsv:', mism)
ala = [g for g in fung if a[g] >= 1 and r[g] == 0]
arg = [g for g in fung if r[g] >= 3 and a[g] == 0]
print('ala-class', len(ala), 'arg-class', len(arg), 'ala-class distinct acc', len(set(acc[ala])))
# lineage table
def lin(t):
    for k in ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales'):
        if k in t.split(';'): return k
    return 'other:' + ';'.join(t.split(';')[4:7])
L = defaultdict(lambda: [0, 0, 0, []])
for g in ala:
    x = L[lin(tax[g])]; x[0] += 1; x[1] += a[g]; x[2] += r[g]; x[3].append(acc[g])
for k, v in L.items(): print(k, v[:3])
for g in fung:
    if a[g] >= 1 and r[g] >= 1 and r[g] < 3: print('both', acc[g], tax[g].split(';')[-1], a[g], r[g])
# class contrasts
for cl in ('Tremellomycetes', 'Exobasidiomycetes', 'Cystobasidiomycetes', 'Ustilaginomycetes', 'Agaricomycetes', 'Microbotryomycetes'):
    gs = [g for g in fung if cl in tax[g].split(';')]
    t = np.array([a[gs].sum(), r[gs].sum(), sa[gs].sum(), sr[gs].sum()])
    p = fisher_exact([[t[0], t[1]], [t[2], t[3]]], alternative='greater')[1]
    print(f'{cl}: {len(gs)} genomes, AGA@A {t[0]} AGA@R {t[1]} syn@A {t[2]} syn@R {t[3]} p={p:.2g}; '
          f'AGA@R from ala-class genomes: {sum(r[g] for g in gs if a[g] >= 1)}')
np.save(A + 'audit/ala_class_idx.npy', np.array(ala)); np.save(A + 'audit/arg_class_idx.npy', np.array(arg))
with open(A + 'audit/ala_class.tsv', 'w') as o:
    o.write('acc\tsrc\tlineage\tspecies_tax\tAGA@A\tAGA@R\tsyn@A\tsyn@R\n')
    for g in sorted(ala, key=lambda g: (lin(tax[g]), acc[g])):
        o.write(f'{acc[g]}\t{src[g]}\t{lin(tax[g])}\t{tax[g].split(";")[-1]}\t{a[g]}\t{r[g]}\t{sa[g]}\t{sr[g]}\n')
