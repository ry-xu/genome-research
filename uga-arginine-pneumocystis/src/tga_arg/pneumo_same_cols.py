"""At each census column where P. jirovecii (NC_020331.1) has TGA at a strict Arg column, the codon each other
Pneumocystis genome has there. usage: python pneumo_same_cols.py out.tsv census.npz [...]"""
import sys, os
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

out, *paths = sys.argv[1:]
D = load(paths)
col, sc, sg = D['col'], D['sc'], D['sg']
acc = list(D['acc'])
pj = acc.index('NC_020331.1')
pn = [g for g in range(len(acc)) if 'Pneumocystis' in D['tax'][g].split(';') and g != pj]
R = AA.index('R')
m = (sg == pj) & (sc == CODONS.index('TGA')) & D['strict'][col] & (D['cons'][col] == R)
cols = sorted(set(col[m].tolist()))
tot = Counter()
with open(out, 'w') as fh:
    fh.write('col\t' + '\t'.join(acc[g] for g in pn) + '\n')
    for c in cols:
        row = []
        for g in pn:
            k = np.where((sg == g) & (col == c))[0]
            x = CODONS[sc[k[0]]] if len(k) else '-'
            row.append(x); tot[x] += 1
        fh.write(f'{c}\t' + '\t'.join(row) + '\n')
print('jirovecii TGA strict-R columns', len(cols), 'other genomes', len(pn))
print(sorted(tot.items(), key=lambda t: -t[1]))
