import sys, os, numpy as np
sys.path.insert(0, 'src/codes_census')
from sitescan_org import load, AA, CODONS
taxa = sys.argv[1].split(','); D = load(sys.argv[2:])
col = D['col']; ok = D['strict'][col]; cons = D['cons'][col]
P, L = AA.index('P'), AA.index('L')
grp = {'CCN': [CODONS.index(c) for c in CODONS if c.startswith('CC')], 'CTN': [CODONS.index(c) for c in CODONS if c.startswith('CT')],
       'TTR': [CODONS.index(c) for c in ('TTA', 'TTG')]}
print('acc\tgenus_species\tP:CCN\tP:CTN\tL:CTN\tL:TTR\tall_aligned:CCN\tall_aligned:CTN')
for g, a in enumerate(D['acc']):
    t = D['tax'][g].split(';')
    if not any(x in t for x in taxa): continue
    m = D['sg'] == g
    r = [((m & ok & (cons == y)) & np.isin(D['sc'], grp[k])).sum() for y, k in ((P, 'CCN'), (P, 'CTN'), (L, 'CTN'), (L, 'TTR'))]
    r += [(m & np.isin(D['sc'], grp[k])).sum() for k in ('CCN', 'CTN')]
    print(a, t[-1], *r, sep='\t')
