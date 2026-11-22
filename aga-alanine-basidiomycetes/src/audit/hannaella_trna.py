"""Check README: Hannaella NC_039747 has a tRNA(UCU) only at ARAGORN threshold_scale 0.8, with G3:U70 and A73.
Uses the project's acceptor parse (org_trnas.acceptor). argv: fasta"""
import sys, pyaragorn
sys.path.insert(0, 'lab/work/new-stops-euk/src/codes_census')
from org_trnas import acceptor
seq = ''.join(l.strip() for l in open(sys.argv[1]) if l[0] != '>')
for ts in (1.0, 0.9, 0.8):
    try: f = pyaragorn.RNAFinder(4, tmrna=False, linear=False, threshold_scale=ts)
    except TypeError as e: print('no threshold_scale', e); break
    n = 0
    for g in f.find_rna(seq):
        if g.type != 'tRNA': continue
        n += 1
        if g.anticodon.upper() == 'TCT' or g.amino_acid == 'Arg':
            s = g.sequence() if callable(g.sequence) else g.sequence
            print(ts, g.amino_acid, g.anticodon.upper(), g.begin, g.end, g.strand, f'{g.energy:.1f}', acceptor(s), s)
    print(ts, 'tRNAs', n)
