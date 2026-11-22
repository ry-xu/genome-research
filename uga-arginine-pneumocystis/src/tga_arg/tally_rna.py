"""Sum reads_codons.py tables over runs: per site, reads showing the genomic codon vs each other triplet.
usage: python tally_rna.py CODON table.tsv [table2.tsv ...]
Prints per-site rows for CODON (consensus, strict, reads, same, top others) and totals by consensus residue, plus the
same totals for every codon in the tables (genome-wide mismatch background)."""
import sys, ast
from collections import Counter, defaultdict


def main(codon, *tables):
    site = {}
    for t in tables:
        for line in open(t):
            f = line.rstrip('\n').split('\t')
            if f[0] == 'gene':
                continue
            k = (f[0], int(f[1]), f[2])
            r = site.setdefault(k, {'cons': f[4], 'strict': f[6], 'col': f[3], 'reads': 0, 'same': 0, 'other': Counter()})
            r['reads'] += int(f[7]); r['same'] += int(f[8]); r['other'].update(ast.literal_eval(f[9]))
    tot = defaultdict(lambda: [0, 0, Counter()])
    for (g, i, c), r in sorted(site.items()):
        if c == codon:
            print(g, i, c, r['col'], r['cons'], 'strict' if r['strict'] == '1' else '', r['reads'], r['same'],
                  dict(r['other'].most_common(3)), sep='\t')
        key = (c, r['cons'] if c == codon else '*')
        tot[key][0] += r['reads']; tot[key][1] += r['same']; tot[key][2].update(r['other'])
    print('-- totals (codon, consensus): reads, same, top other triplets')
    for k in sorted(tot):
        n, s, o = tot[k]
        print(k, n, s, f'{s / n:.3f}' if n else '-', dict(o.most_common(4)), sep='\t')


if __name__ == '__main__':
    main(*sys.argv[1:])
