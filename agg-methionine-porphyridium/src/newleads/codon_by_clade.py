"""Per-genome counts of chosen codons at strict columns by consensus residue, for genomes whose taxonomy matches a text.
usage: python codon_by_clade.py TAXON_TEXT CODON[,CODON...] RES[,RES...] census.npz [census2.npz ...]
Prints one row per matching genome: acc, table, last 3 taxonomy ranks, then codon@res counts (strict columns) and the
codon's total aligned sites."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS


def main(text, codons, resids, *paths):
    D = load(paths)
    codons, resids = codons.split(','), resids.split(',')
    sel = [g for g, t in enumerate(D['tax']) if text in t]
    sg, sc, col = D['sg'], D['sc'], D['col']
    st, cc = D['strict'][col], D['cons'][col]
    m = np.isin(sg, sel)
    sgm, scm, stm, ccm = sg[m], sc[m], st[m], cc[m]
    hdr = ['acc', 'tab', 'taxon'] + [f'{c}@{r}' for c in codons for r in resids] + [f'{c}_n' for c in codons]
    print('\t'.join(hdr))
    for g in sel:
        k = sgm == g
        row = [D['acc'][g], str(D['table'][g]), ';'.join(D['tax'][g].split(';')[-3:])]
        for c in codons:
            ci = [x for x in CODONS if all(a in (b, "N") for b, a in zip(x, c))]
            ci = [CODONS.index(x) for x in ci]
            for r in resids:
                row.append(str(int(np.sum(k & np.isin(scm, ci) & stm & (ccm == AA.index(r))))))
        for c in codons:
            row.append(str(int(np.sum(k & np.isin(scm, [CODONS.index(x) for x in CODONS if all(a in (b, "N") for b, a in zip(x, c))])))))
        print('\t'.join(row))


if __name__ == '__main__':
    main(*sys.argv[1:])
