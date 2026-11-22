"""Codons of other genomes at the columns where one genome uses a codon.
usage: python same_column.py ACC CODON OTHER_ACC[,OTHER_ACC...] census.npz [census2.npz ...]
For each site where ACC uses CODON: model:column, consensus, strict, then each other genome's codon(s) there."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS


def main(acc, codon, others, *paths):
    D = load(paths)
    others = others.split(',')
    gi = {a: i for i, a in enumerate(D['acc'])}
    sg, sc, col = D['sg'], D['sc'], D['col']
    K = D['K']
    by = {}
    for a in others:
        k = np.flatnonzero(sg == gi[a])
        d = {}
        for c, x in zip(col[k], sc[k]):
            d.setdefault(int(c), []).append(CODONS[x])
        by[a] = d
    k = np.flatnonzero((sg == gi[acc]) & (sc == CODONS.index(codon)))
    print('model:col', 'cons', 'strict', *others, sep='\t')
    for c in sorted(set(int(x) for x in col[k])):
        print(f"{D['models'][c // K]}:{c % K}", AA[D['cons'][c]], int(D['strict'][c]),
              *[','.join(by[a].get(c, ['-'])) for a in others], sep='\t')


if __name__ == '__main__':
    main(*sys.argv[1:])
