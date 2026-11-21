import sys, numpy as np
from collections import Counter
sys.path.insert(0, 'src/codes_census')
sys.path.insert(0, 'src')
from sitescan_org import load, AA, CODONS
accf, tsv = sys.argv[1], sys.argv[2]
D = load(sys.argv[3:])
want = [l.split()[0] for l in open(accf) if l.strip()]
gi = {a: i for i, a in enumerate(D['acc'])}
print('missing', [a for a in want if a not in gi])
want = [a for a in want if a in gi]
P = AA.index('P'); L = AA.index('L')
cls = np.full(64, 2)
for i, c in enumerate(CODONS):
    if c.startswith('CC'): cls[i] = 0
    elif c.startswith('CT'): cls[i] = 1
col, sg, sc = D['col'], D['sg'], D['sc']
pro = D['strict'] & (D['cons'] == P)
G = np.array([gi[a] for a in want])
sel = np.isin(sg, G) & pro[col]
cols = np.unique(col[sel])
have = np.zeros((len(G), len(cols)), bool); grid = np.full((len(G), len(cols)), 3)
ci = {c: j for j, c in enumerate(cols)}; ri = {g: r for r, g in enumerate(G)}
for g, c, k in zip(sg[sel], col[sel], sc[sel]):
    grid[ri[g], ci[c]] = cls[k]; have[ri[g], ci[c]] = True
spelled = (grid <= 1).sum(0)
keep = (have.mean(0) >= 0.6) & (spelled >= 0.5 * have.sum(0))
print('strict-Pro cols with any of the 49:', len(cols), 'kept', keep.sum())
K = D['K']; models = D['models']
km = [models[c // K] for c in cols[keep]]
print('models of kept cols', Counter(km))
prot = Counter(m.split('_')[-1] for m in km); print('distinct proteins', len(prot), prot)
g2 = grid[:, keep]
ref = {l.split('\t')[0]: l.rstrip('\n').split('\t')[2:] for l in open(tsv) if not l.startswith('acc')}
bad = 0
for r, a in enumerate(want):
    n = np.bincount(g2[r], minlength=4).tolist()
    if [str(x) for x in n] != ref[a]: bad += 1; print('MISMATCH', a, n, ref[a])
print('rows mismatching tsv:', bad, 'of', len(want))
# all strict Pro columns: how many models/proteins in whole census set for the four
for a in ['BK074532.1','BK074538.1','BK074542.1','BK074262.1']:
    m = (sg == gi[a]) & pro[col] & np.isin(sc, [i for i,c in enumerate(CODONS) if c.startswith('CT')])
    print(a, 'CTN@strictP', m.sum(), 'models', Counter(models[c//K] for c in col[m]).most_common())
    mL = (sg == gi[a]) & D['strict'][col] & (D['cons'][col] == L)
    print('   at strict L codons:', Counter(CODONS[k] for k in sc[mL]).most_common())
# genome-level usage of CTC/CTG in the four (all aligned sites)
for a in ['BK074532.1','BK074538.1','BK074542.1','BK074262.1']:
    m = (sg == gi[a])
    cc = Counter(CODONS[k] for k in sc[m])
    print(a, 'aligned-site usage CTT/CTC/CTA/CTG', [cc[c] for c in ('CTT','CTC','CTA','CTG')], 'CCN', [cc[c] for c in ('CCT','CCC','CCA','CCG')])
