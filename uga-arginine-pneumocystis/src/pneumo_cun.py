"""What do CUN codons face in Pneumocystis mitochondria? (Ma et al. 2013 report an extra tRNA-Thr(UAG) in P. jirovecii;
a UAG anticodon reads CUN, and CUN = Thr is the yeast mitochondrial change, table 3.)
For every Pneumocystis genome in the RefSeq and GenBank mito censuses: codons CTN, TTR (Leu) and ACN (Thr) at profile
columns whose most likely emission has probability >= 0.6, tallied by that consensus residue.
usage: python -I src/pneumo_cun.py data/codes/mito.census.npz data/codes/gbmito.census.npz"""
import sys
from collections import Counter
import numpy as np

AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
GROUPS = {'CTN': ['CTT', 'CTC', 'CTA', 'CTG'], 'TTR': ['TTA', 'TTG'], 'ACN': ['ACT', 'ACC', 'ACA', 'ACG']}
tot = Counter()
print('acc\tgroup\tn\tconsensus residues')
for f in sys.argv[1:]:
    d = np.load(f)
    sg, sc, sm, sk, emax, earg = (d[k] for k in ('site_genome', 'site_codon', 'site_model', 'site_col', 'emax', 'earg'))
    for g in [i for i, t in enumerate(d['taxonomy']) if 'Pneumocystis' in str(t)]:
        m = sg == g
        e, a = emax[sm[m], sk[m]], earg[sm[m], sk[m]]
        for name, cods in GROUPS.items():
            sel = np.isin(sc[m], [COD.index(x) for x in cods]) & (e >= 0.6)
            cnt = Counter(AA[x] for x in a[sel])
            tot.update({(name, k): v for k, v in cnt.items()})
            print(f"{d['acc'][g]}\t{name}\t{sel.sum()}\t" + ' '.join(f'{k}{v}' for k, v in cnt.most_common()))
print('total\t' + ' '.join(f'{n}:{k}{v}' for (n, k), v in sorted(tot.items())))
