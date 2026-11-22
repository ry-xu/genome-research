"""Codon-specific contrast: does a codon sit at strict columns of residue Y more than its standard synonyms do?
usage: python codon_contrast.py TAXON|ACC CODON Y census.npz [census2.npz ...]
For every genome whose accession equals TAXON or whose taxonomy contains it, counts at strict columns:
  a = CODON at Y columns, b = CODON at STD(CODON) columns, c = standard synonyms at Y, d = synonyms at STD columns,
and a one-sided Fisher p for a/b > c/d. Real Y<->STD substitutions move all synonyms alike; a meaning change moves
only CODON. Also prints the pooled 2x2 over the genomes listed."""
import sys, os
import numpy as np
from scipy.stats import fisher_exact
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sitescan_org import load, AA, CODONS, STD


def main(taxon, codon, y, *paths):
    D = load(paths)
    ci, std = CODONS.index(codon), STD[codon]
    syn = [CODONS.index(c) for c in CODONS if STD.get(c) == std and c != codon]
    yi, si = AA.index(y), AA.index(std)
    strict, cons = D['strict'], D['cons']
    ok = strict[D['col']]
    cy = ok & (cons[D['col']] == yi)
    cs = ok & (cons[D['col']] == si)
    isc = D['sc'] == ci
    iss = np.isin(D['sc'], syn)
    rows = [g for g, a in enumerate(D['acc']) if a == taxon or taxon in D['tax'][g].split(';')]
    print(f'{codon} (std {std}) vs synonyms {[CODONS[s] for s in syn]} at strict {y} / {std} columns; {len(rows)} genomes')
    print('acc\ttable\tgenus\tcodon@Y\tcodon@STD\tsyn@Y\tsyn@STD\tp')
    tot = np.zeros(4, int)
    G = len(D['acc'])
    cnt = np.stack([np.bincount(D['sg'][a & b], minlength=G) for a, b in ((isc, cy), (isc, cs), (iss, cy), (iss, cs))], 1)
    for g in rows:
        n = cnt[g]
        tot += n
        p = fisher_exact([[n[0], n[1]], [n[2], n[3]]], alternative='greater')[1] if n.sum() else 1.0
        print(f"{D['acc'][g]}\t{D['table'][g]}\t{D['tax'][g].split(';')[-1]}\t" + '\t'.join(map(str, n)) + f'\t{p:.1e}')
    p = fisher_exact([[tot[0], tot[1]], [tot[2], tot[3]]], alternative='greater')[1] if tot.sum() else 1.0
    print('pooled\t\t\t' + '\t'.join(map(str, tot)) + f'\t{p:.1e}')


if __name__ == '__main__':
    main(*sys.argv[1:])
