"""Per genome in a taxon: placed TGA and TGG codons (all columns), and TGA/TGG at strict Trp columns, from the census.
usage: python clade_tga.py TAXON census.npz [...]"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

taxon, *paths = sys.argv[1:]
D = load(paths)
col, sc, sg = D['col'], D['sc'], D['sg']
keep = np.array([taxon in t.split(';') for t in D['tax']])
m = keep[sg]
col, sc, sg = col[m], sc[m], sg[m]
ok, cons = D['strict'][col], D['cons'][col]
W = AA.index('W')
n = len(D['acc'])
tga, tgg = sc == CODONS.index('TGA'), sc == CODONS.index('TGG')
cnt = lambda x: np.bincount(sg[x], minlength=n)
a, b, c, d = cnt(tga), cnt(tgg), cnt(tga & ok & (cons == W)), cnt(tgg & ok & (cons == W))
print('acc\ttable\tTGA\tTGG\tTGA_strictW\tTGG_strictW\ttaxon')
for g in np.where(keep)[0]:
    t = D['tax'][g].split(';')
    print(f"{D['acc'][g]}\t{D['table'][g]}\t{a[g]}\t{b[g]}\t{c[g]}\t{d[g]}\t{';'.join(t[-4:])}")
