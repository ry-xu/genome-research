"""Consensus residue under each codon at ALL aligned columns (no strict filter), pooled over a genome list.
usage: python all_columns.py accs.txt census.npz [census2.npz ...]
For AGA, CGN (Arg synonyms) and GCN (Ala) prints the share of sites whose column consensus is each residue, and the
mean consensus share of R and A at the codon's columns (how Arg-like vs Ala-like its positions are)."""
import sys, os
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS


def main(accs, *paths):
    want = set(open(accs).read().split())
    D = load(paths)
    gi = [g for g, a in enumerate(D['acc']) if a in want]
    m = np.isin(D['sg'], gi)
    sc, col = D['sc'][m], D['col'][m]
    cons = D['cons'][col]
    groups = {'AGA': ['AGA'], 'CGN': ['CGT', 'CGC', 'CGA', 'CGG'], 'AGG': ['AGG'], 'GCN': ['GCT', 'GCC', 'GCA', 'GCG']}
    print(f'{len(gi)} genomes')
    for k, cods in groups.items():
        sel = np.isin(sc, [CODONS.index(c) for c in cods])
        n = int(sel.sum())
        c = Counter(AA[i] for i in cons[sel])
        top = ', '.join(f'{a} {v / n:.2f}' for a, v in c.most_common(8)) if n else ''
        rk = ' '.join(f'{a}={c[a] / n:.3f}' for a in 'RKA') if n else ''
        print(f'{k}: {n} sites; column consensus: {top} | {rk}')


if __name__ == '__main__':
    main(*sys.argv[1:])
