"""Two-residue test: does codon X sit at columns of residue B (vs its table-4 residue A) more than its synonyms do?
Fisher one-sided on [[X@B, X@A], [syn@B, syn@A]] using usage_at_aa_columns.tsv from codon_stats.py.
usage: pair_test.py <usage_at_aa_columns.tsv> <codon> <A> <B> <syn1,syn2,...>"""
import sys, csv, collections
from scipy.stats import fisher_exact
f, X, A, B, syn = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5].split(',')
cnt = collections.defaultdict(int)
for r in csv.DictReader(open(f), delimiter='\t'):
    cnt[(r['acc'], r['def'], r['aa'], r['codon'])] += int(r['count_at_cons_cols_of_aa'])
keys = sorted({(a, d) for a, d, _, _ in cnt})
print(f'acc\tdef\t{X}@{B}\t{X}@{A}\tsyn@{B}\tsyn@{A}\tp_one_sided')
for acc, d in keys:
    xb, xa = cnt[(acc, d, B, X)], cnt[(acc, d, A, X)]
    sb = sum(cnt[(acc, d, B, s)] for s in syn); sa = sum(cnt[(acc, d, A, s)] for s in syn)
    p = fisher_exact([[xb, xa], [sb, sa]], alternative='greater')[1] if xb + xa else float('nan')
    print(f'{acc}\t{d}\t{xb}\t{xa}\t{sb}\t{sa}\t{p:.2e}')
