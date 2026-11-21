"""CUU/CUA codons of the four CUN=Pro genomes at strict census columns, by the column's consensus residue, using the
census loader itself (sitescan_org.load). Recount of the pre-post checker's 286-of-291 (prepost/cun-pro_work/).
usage: python cun_by_res.py census.npz [census2.npz ...]"""
import sys, os
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

FOUR = {'BK074532.1': 'W. infanticola', 'BK074538.1': 'W. siamensis', 'BK074542.1': 'W. verensis',
        'BK074262.1': 'S. pinicola'}
D = load(sys.argv[1:])
gi = {a: i for i, a in enumerate(D['acc'])}
cun = np.isin(D['sc'], [i for i, c in enumerate(CODONS) if c in ('CTT', 'CTA')])
tot = Counter()
for a, name in FOUR.items():
    m = (D['sg'] == gi[a]) & cun
    st = D['strict'][D['col'][m]]
    c = Counter(AA[x] for x in D['cons'][D['col'][m]][st])
    tot.update(c)
    print(a, name, 'CUU/CUA at strict columns', int(st.sum()), dict(c.most_common()))
print('all four', sum(tot.values()), dict(tot.most_common()))
