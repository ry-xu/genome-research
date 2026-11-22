"""Selection-free AGA lineage test (pre-post check fix 6): every census genome in each order, whatever its AGA use.
Distinct strict columns per order: AGA at Ala / Arg columns vs its synonyms (CGN, AGG) at Ala / Arg columns;
one-sided Fisher. Census loader: sitescan_org.load.
usage: python lineage_nosel.py census.npz [census2.npz ...]"""
import sys, os
import numpy as np
from scipy.stats import fisher_exact
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

ORDERS = ['Trichosporonales', 'Tremellales', 'Tilletiales', 'Cystobasidiales', 'Sporidiobolales']
D = load(sys.argv[1:])
A, R = AA.index('A'), AA.index('R')
aga = CODONS.index('AGA'); syn = [CODONS.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGG')]
st = D['strict'][D['col']]; cc = D['cons'][D['col']]
for o in ORDERS:
    gs = [i for i, t in enumerate(D['tax']) if o in t.split(';')]
    m = np.isin(D['sg'], gs) & st
    def ncol(codons, res):
        k = m & np.isin(D['sc'], codons) & (cc == res)
        return len(np.unique(D['col'][k]))
    t = [[ncol([aga], A), ncol([aga], R)], [ncol(syn, A), ncol(syn, R)]]
    p = fisher_exact(t, alternative='greater')[1]
    print(f'{o}\tgenomes {len(gs)}\tAGA@Ala {t[0][0]} AGA@Arg {t[0][1]} syn@Ala {t[1][0]} syn@Arg {t[1][1]}\tp {p:.2g}')
