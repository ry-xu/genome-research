"""Align each shared CDS of two GenBank records (protein level, table 4) and tabulate the codon in B at every AGG of A.
usage: python codon_pair.py A.gb B.gb CODON"""
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

def cds(path):
    out = {}
    rec = SeqIO.read(path, 'genbank')
    for f in rec.features:
        if f.type != 'CDS': continue
        g = (f.qualifiers.get('gene') or f.qualifiers.get('product') or ['?'])[0].lower()
        nt = str(f.extract(rec.seq)).upper()
        out[g] = nt
    return out

def tr(nt):
    return str(Seq(nt[:len(nt)//3*3]).translate(table=4)).replace('*', 'X')

A, B, codon = sys.argv[1], sys.argv[2], sys.argv[3]
a, b = cds(A), cds(B)
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.mode = -10, -0.5, 'global'
tot = Counter(); same_nt = diff_nt = 0
for g in sorted(set(a) & set(b)):
    pa, pb = tr(a[g]), tr(b[g])
    aln = al.align(pa, pb)[0]
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1):
            i, j = s1 + k, s2 + k
            ca, cb = a[g][3*i:3*i+3], b[g][3*j:3*j+3]
            if ca == cb: same_nt += 1
            else: diff_nt += 1
            if ca == codon:
                tot[cb] += 1
                print(f'{g}\t{i+1}\t{ca}\t{cb}\t{pa[max(0,i-3):i+4]}\t{pb[max(0,j-3):j+4]}')
print('# shared genes', sorted(set(a) & set(b)), file=sys.stderr)
print(f'# codons identical at aligned positions: {same_nt}/{same_nt+diff_nt}', file=sys.stderr)
print('# codon in B at A', codon, dict(tot.most_common()), file=sys.stderr)
