"""Are AGA-at-strict-Ala columns shared across the three putative independent origins more than chance?
Groups: Tremellomycetes (Trichosporonales+Tremellales), Tilletiales, Pucciniomycotina (Cystobasidiales+Sporidiobolales).
Null: each group's observed number of distinct AGA@Ala columns drawn uniformly (or weighted by the group's Ala-codon
occupancy) from the strict Ala columns where that group has any Ala codon (GCN or AGA) aligned."""
import numpy as np
from collections import Counter
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
tax = C['tax']; sg, sc, col = C['sg'], C['sc'], C['col']; K = int(C['K']); models = C['models']
st = C['strict'][col]; cc = C['cons'][col]; iA = AA.index('A')
ala = np.load(A + 'audit/ala_class_idx.npy')
GRP = {'Tremellomycetes': ('Trichosporonales', 'Tremellales'), 'Tilletiales': ('Tilletiales',),
       'Pucciniomycotina': ('Cystobasidiales', 'Sporidiobolales')}
AGA = COD.index('AGA'); ALA = [COD.index(c) for c in ('GCT', 'GCC', 'GCA', 'GCG', 'AGA')]
obs, pool, wts = {}, {}, {}
for gname, lins in GRP.items():
    gs = [g for g in ala if any(l in tax[g].split(';') for l in lins)]
    m = np.isin(sg, gs) & st & (cc == iA)
    obs[gname] = set(col[m & (sc == AGA)].tolist())
    occ = Counter(col[m & np.isin(sc, ALA)].tolist())
    pool[gname] = np.array(sorted(occ)); wts[gname] = np.array([occ[c] for c in pool[gname]], float)
    print(f'{gname}: {len(gs)} genomes, {len(obs[gname])} AGA@Ala columns, pool {len(pool[gname])} strict Ala columns')
names = list(GRP)
def stat(sets):  # columns present in >= 2 groups
    c = Counter(x for s in sets for x in s); return sum(1 for v in c.values() if v >= 2)
o = stat([obs[n] for n in names])
shared = Counter(x for n in names for x in obs[n])
print('observed columns in >= 2 groups:', o, sorted(f'{models[k // K]}:{k % K}' for k, v in shared.items() if v >= 2))
rng = np.random.default_rng(1)
for mode in ('uniform', 'occupancy-weighted'):
    sims = []
    for _ in range(20000):
        sets = []
        for n in names:
            p = None if mode == 'uniform' else wts[n] / wts[n].sum()
            sets.append(set(rng.choice(pool[n], size=len(obs[n]), replace=False, p=p).tolist()))
        sims.append(stat(sets))
    sims = np.array(sims)
    print(f'{mode} null: mean {sims.mean():.2f}, 99.9th pct {np.percentile(sims, 99.9):.0f}, P(>= {o}) = {(sims >= o).mean():.2g} (20,000 draws)')
# hotspot check: are AGA@Ala columns less purely Ala (more S/T/G minority) across all 22,698 genomes than other pool columns?
freq = C['freq']; print('freq shape', freq.shape)
allpool = set().union(*[set(pool[n].tolist()) for n in names]); agac = set().union(*obs.values())
def prof(cs):
    f = freq[sorted(cs)]; f = f / f.sum(1, keepdims=True)
    return {a: round(float(f[:, AA.index(a)].mean()), 3) for a in 'ASGTCV'}
print('mean residue share, AGA@Ala columns   (n=%d):' % len(agac), prof(agac))
print('mean residue share, other pool columns (n=%d):' % len(allpool - agac), prof(allpool - agac))
from scipy.stats import mannwhitneyu
fa = freq[:, iA] / freq.sum(1)
x, y = fa[sorted(agac)], fa[sorted(allpool - agac)]
print('Ala share median %.3f vs %.3f, Mann-Whitney two-sided p = %.2g' % (np.median(x), np.median(y), mannwhitneyu(x, y).pvalue))
