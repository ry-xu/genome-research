"""Audit: nucleotide and amino-acid identity of the 14 shared CDSs between SAG 1380-2 and CCMP1948 (MAFFT globalpair).
usage: python -I nt_identity.py A.gb B.gb"""
import sys, importlib.util, os
spec = importlib.util.spec_from_file_location('pair_nt', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pair_nt.py'))
p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)
from Bio.Seq import Seq
a, _ = p.cds(sys.argv[1]); b, _ = p.cds(sys.argv[2])
tot = {'nt': [0, 0], 'aa': [0, 0]}
for g in sorted(set(a) & set(b)):
    x, y = p.mafft_pair(a[g][0], b[g][0], 'nt')
    m = sum(1 for i, j in zip(x, y) if i == j != '-'); n = sum(1 for i, j in zip(x, y) if i != '-' and j != '-')
    tot['nt'][0] += m; tot['nt'][1] += n
    pa = str(Seq(a[g][0]).translate(table=4)); pb = str(Seq(b[g][0]).translate(table=4))
    x, y = p.mafft_pair(pa.replace('*', 'X'), pb.replace('*', 'X'), 'aa')
    m2 = sum(1 for i, j in zip(x, y) if i == j != '-'); n2 = sum(1 for i, j in zip(x, y) if i != '-' and j != '-')
    tot['aa'][0] += m2; tot['aa'][1] += n2
    print(f'{g}\tnt {m}/{n} {m/n:.3f}\taa {m2}/{n2} {m2/n2:.3f}')
for k, (m, n) in tot.items():
    print(f'TOTAL {k}: {m}/{n} = {m/n:.3f}')
