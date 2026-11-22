"""Decoy: for every codon, where does it sit at strict columns in the 32 ala-class genomes (pooled per lineage)?
A codon whose strict sites mostly vote for one non-standard residue would be a decoy 'hit'. Also: codon usage."""
import numpy as np
from collections import Counter
from Bio.Data import CodonTable
A = 'lab/work/new-stops-euk/results/aga_ala/'
C = np.load(A + 'audit/cache.npz', allow_pickle=True)
T4 = CodonTable.unambiguous_dna_by_id[4].forward_table
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
acc, tax = C['acc'], C['tax']
sg, sc, col = C['sg'], C['sc'], C['col']
st = C['strict'][col]; cc = C['cons'][col]
ala = np.load(A + 'audit/ala_class_idx.npy')
lin = {}
for g in ala:
    p = tax[g].split(';')
    lin[g] = next(k for k in ('Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales') if k in p)
groups = {'all32': set(ala)}
for g, l in lin.items(): groups.setdefault(l, set()).add(g)
out = open(A + 'audit/decoy.tsv', 'w')
out.write('group\tcodon\tT4\taligned\tstrict\tstd_share\ttop_nonstd\ttop_nonstd_n\ttop_nonstd_share\tvotes\n')
for name, gs in groups.items():
    m = np.isin(sg, list(gs))
    for ci, c in enumerate(COD):
        mm = m & (sc == ci)
        n_al = int(mm.sum()); ms = mm & st
        cnt = Counter(AA[x] for x in cc[ms]); n = sum(cnt.values())
        aa = T4.get(c, '*')
        ns = [(k, v) for k, v in cnt.most_common() if k != aa]
        top, tn = ns[0] if ns else ('-', 0)
        out.write(f'{name}\t{c}\t{aa}\t{n_al}\t{n}\t{cnt.get(aa, 0) / max(n, 1):.2f}\t{top}\t{tn}\t{tn / max(n, 1):.2f}\t{dict(cnt.most_common(5))}\n')
out.close()
# print the interesting part
import csv
rows = list(csv.DictReader(open(A + 'audit/decoy.tsv'), delimiter='\t'))
print('all32: codons by aligned count; strict sites; share at standard residue; top non-standard')
for r in sorted([r for r in rows if r['group'] == 'all32'], key=lambda r: int(r['aligned'])):
    if int(r['aligned']) < 2000 or r['codon'] == 'AGA':
        print(r['codon'], r['T4'], r['aligned'], r['strict'], r['std_share'], r['top_nonstd'], r['top_nonstd_n'], r['top_nonstd_share'], r['votes'])
print('\nany codon/group with >=4 strict sites where std share < 0.5:')
for r in rows:
    if int(r['strict']) >= 4 and float(r['std_share']) < 0.5:
        print(r['group'], r['codon'], r['T4'], r['strict'], r['std_share'], r['votes'])
