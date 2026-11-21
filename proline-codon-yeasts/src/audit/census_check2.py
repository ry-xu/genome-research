import sys, numpy as np, csv
from collections import Counter
AA = 'ACDEFGHIKLMNPQRSTVWY'
CODONS = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
d = np.load(sys.argv[1])
sg, sc, col, res = d['sg'], d['sc'], d['col'], d['res']
freq, tot, strict, cons = d['freq'], d['tot'], d['strict'], d['cons']
acc = d['acc'].tolist(); tax = d['tax'].tolist(); models = d['models'].tolist(); K = int(d['K'])
gi = {a: i for i, a in enumerate(acc)}
FOUR = ['BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1']
P, L = AA.index('P'), AA.index('L')
CTN = [CODONS.index(c) for c in CODONS if c.startswith('CT')]
CCN = [CODONS.index(c) for c in CODONS if c.startswith('CC')]
TTR = [CODONS.index(c) for c in ('TTA', 'TTG')]
def counts(st, cn, g):
    m = sg == g
    ok = st[col]
    out = []
    for y, grp in ((P, CCN), (P, CTN), (L, CTN), (L, TTR)):
        out.append(int((m & ok & (cn[col] == y) & np.isin(sc, grp)).sum()))
    return out
print('acc  P:CCN P:CTN L:CTN L:TTR (as in pro_leu)')
for a in FOUR + ['BK074539.1', 'BK074529.1', 'BK074524.1']:
    print(a, counts(strict, cons, gi[a]))
# per-codon at strict P for the four
for a in FOUR:
    m = (sg == gi[a]) & strict[col] & (cons[col] == P)
    print(a, 'codons at strict P:', {CODONS[k]: v for k, v in Counter(sc[m].tolist()).most_common()})
# leave-out: remove the four genomes' own votes from freq
w = 1 / 105.
f2 = freq.copy()
for a in FOUR:
    m = (sg == gi[a]) & (res >= 0)
    np.add.at(f2, (col[m], res[m]), -w)
f2 = np.maximum(f2, 0)
t2 = f2.sum(1)
s2 = (t2 > 0) & (f2.max(1) >= 0.9 * np.maximum(t2, 1e-12))
c2 = f2.argmax(1)
print('leave-four-out consensus:')
for a in FOUR + ['BK074539.1']:
    print(a, counts(s2, c2, gi[a]))
# how much do the four move columns? columns whose strict/cons changes
ch = (s2 != strict) | (c2 != cons)
print('columns changed by removing the four:', int(ch.sum()))
# Strict P columns: dominated by whom? weight share of Dipodascomycetes at the four's CTN@P columns
for a in FOUR[:2]:
    m = (sg == gi[a]) & strict[col] & (cons[col] == P) & np.isin(sc, CTN)
    cols = np.unique(col[m])
    print(a, 'n distinct CTN@P cols', len(cols), 'median tot weight', float(np.median(tot[cols])), 'min share', float((freq[cols].max(1) / tot[cols]).min()))
# dump CTN@strictP keys for comparing with reads tsv
for a in FOUR:
    m = (sg == gi[a]) & strict[col] & (cons[col] == P) & np.isin(sc, CTN)
    keys = sorted({(CODONS[c], f'{models[k // K]}:{k % K}') for c, k in zip(sc[m], col[m])})
    with open(f'{sys.argv[2]}/{a}.census_ctn_p.txt', 'w') as fh:
        for c, mk in keys:
            fh.write(f'{c}\t{mk}\n')
