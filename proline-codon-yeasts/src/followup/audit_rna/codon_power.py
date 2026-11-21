"""(c)/(f) what the W. versatilis test can see at the codons that matter: CUN codons (by codon) with cp2 at RNA
depth>=20 (caller's pileup), Sc-aligned exon U-cp2 count, and ONT-variant sites with RNA depth>=20 not called.
usage: codon_power.py ref.fa genes.tsv rna.bam rdd.tsv ont_check.txt"""
import sys, csv
sys.path.insert(0, 'src/followup')  # run from the folder root
import rdd
from collections import OrderedDict, Counter
fa, genes, bam, rddtsv = sys.argv[1:5]
g = rdd.fasta1(fa)
rows = OrderedDict()
for r in csv.DictReader(open(genes), delimiter='\t'):
    rows.setdefault(r['gene'], []).append((int(r['pos1']), int(r['strand'])))
pr = rdd.pileup(fa, bam); dep = {p: sum(c.values()) for p, c in pr.items()}
allc, tested = Counter(), Counter()
for name, rs in rows.items():
    for i in range(0, len(rs) - 2, 3):
        cod = ''.join(g[p - 1] for p, _ in rs[i:i + 3])  # all genes + strand here
        allc[cod] += 1
        if dep.get(rs[i + 1][0], 0) >= 20: tested[cod] += 1
for c in ('CTT', 'CTA', 'CTC', 'CTG', 'TTA', 'TTG', 'CCT', 'CCA', 'CCC', 'CCG'):
    print(f'{c}: in ORFs {allc[c]}, cp2 at RNA depth>=20 {tested[c]}')
print('U at cp2, codons with cp2 at depth>=20:', sum(v for k, v in tested.items() if k[1] == 'T'))
