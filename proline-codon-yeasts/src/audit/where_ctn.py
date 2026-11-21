import sys, numpy as np
from collections import Counter
sys.path.insert(0, 'src/codes_census')
sys.path.insert(0, 'src')
from sitescan_org import load, AA, CODONS
D = load(sys.argv[1:])
gi = {a: i for i, a in enumerate(D['acc'])}
P = AA.index('P'); L = AA.index('L')
share = D['freq'].max(1) / np.maximum(D['tot'], 1e-12)
ctn = np.array([c.startswith('CT') for c in CODONS])
for a in ['BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1', 'BK074539.1']:
    m = (D['sg'] == gi[a]) & ctn[D['sc']]
    cols = D['col'][m]
    cons = D['cons'][cols]; sh = share[cols]
    leu_share = D['freq'][cols, L] / np.maximum(D['tot'][cols], 1e-12)
    pro_share = D['freq'][cols, P] / np.maximum(D['tot'][cols], 1e-12)
    print(a, 'CTN aligned', m.sum(), 'consensus residue:', Counter(AA[c] for c in cons).most_common(6))
    for thr in (0.5, 0.7, 0.8, 0.9):
        print(f'   cons=L with share>={thr}:', int(((cons == L) & (sh >= thr)).sum()), f' cons=P share>={thr}:', int(((cons == P) & (sh >= thr)).sum()))
    print('   max Leu share at a CTN site:', round(float(leu_share.max()), 3), ' sites with Leu share>=0.5:', int((leu_share >= 0.5).sum()))
