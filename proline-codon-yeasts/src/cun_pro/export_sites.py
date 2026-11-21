"""Write every aligned CUU/CUA codon of the four CUN=Pro genomes with its census column, the column's consensus residue,
its weighted share and whether the column is strict (share >= 0.9). Same loader and consensus as cun_by_res.py.
usage: python -I src/cun_pro/export_sites.py OUT.tsv mito.census.npz gbmito.census.npz gbmitoorf.census.npz
Made for this repository so that the 286-of-290 count can be redone from data/census/cun_sites.tsv
(src/check_numbers.py) without the census files."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

FOUR = {'BK074532.1': 'W. infanticola', 'BK074538.1': 'W. siamensis', 'BK074542.1': 'W. verensis',
        'BK074262.1': 'S. pinicola'}
D = load(sys.argv[2:])
gi = {a: i for i, a in enumerate(D['acc'])}
K = D['K']
share = D['freq'].max(1) / np.maximum(D['tot'], 1e-12)
want = [CODONS.index(c) for c in ('CTT', 'CTA')]
with open(sys.argv[1], 'w') as f:
    f.write('acc\tspecies\tcodon\tmodel\tcolumn\tconsensus\tshare\tstrict\n')
    for a, name in FOUR.items():
        m = (D['sg'] == gi[a]) & np.isin(D['sc'], want)
        rows = sorted(zip(D['col'][m].tolist(), D['sc'][m].tolist()))
        for col, sc in rows:
            f.write(f"{a}\t{name}\t{CODONS[sc]}\t{D['models'][col // K]}\t{col % K}\t{AA[D['cons'][col]]}\t"
                    f"{share[col]:.3f}\t{int(D['strict'][col])}\n")

# Per genome: all codons at strict proline columns (the "91-93 positions" in the write-up), into OUT's folder
P = AA.index('P')
ccn = [i for i, c in enumerate(CODONS) if c.startswith('CC')]
with open(os.path.join(os.path.dirname(sys.argv[1]) or '.', 'strict_pro_codons.tsv'), 'w') as f:
    f.write('acc\tspecies\tstrict_pro_sites\tCUU_CUA\tCCN\tother\n')
    for a, name in FOUR.items():
        m = (D['sg'] == gi[a]) & D['strict'][D['col']] & (D['cons'][D['col']] == P)
        n, cu, cc = int(m.sum()), int((m & np.isin(D['sc'], want)).sum()), int((m & np.isin(D['sc'], ccn)).sum())
        f.write(f'{a}\t{name}\t{n}\t{cu}\t{cc}\t{n - cu - cc}\n')
