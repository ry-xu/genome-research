"""List conserved-column occurrences of a codon with distance (codons) to the enclosing segment ends.
usage: inspect_codon.py <ortho_dir> <acc> <codon> [min_frac=0.8]"""
import sys, csv, os
d, acc, codon = sys.argv[1:4]
mf = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
segs = [l.rstrip('\n').split('\t') for l in open(os.path.join(d, f'{acc}.segments.tsv'))]
rows = list(csv.DictReader(open(os.path.join(d, f'{acc}.columns.tsv')), delimiter='\t'))
for r in rows:
    if r['codon'] != codon or int(r['nres']) < 6 or int(r['ntop']) < mf * int(r['nres']):
        continue
    x = int(r['coord'])
    edge = None
    for g, st, lo, hi, *_ in segs:
        lo, hi = int(lo), int(hi)
        if g == r['gene'] and lo - 3 <= x <= hi + 3:
            edge = min(x - lo, hi - x) // 3
    print(acc, r['gene'], r['strand'], x, codon, '->', r['top'], f"{r['ntop']}/{r['nres']}", r['residues'],
          'col', r['col'], 'codons_from_segment_edge', edge)
