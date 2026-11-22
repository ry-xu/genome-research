"""Per genome of a taxon: TGA, TGG and Arg codons (CGN, AGR) at strict Arg / strict Trp columns of the census, plus
TGA at all placed sites by column consensus. usage: python pneumo_cols.py TAXON out.tsv census.npz [census2.npz ...]"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

taxon, out, *paths = sys.argv[1:]
D = load(paths)
col, sc, sg = D['col'], D['sc'], D['sg']
ok, cons = D['strict'][col], D['cons'][col]
R, W = AA.index('R'), AA.index('W')
arg = np.isin(sc, [CODONS.index(c) for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG')])
tga, tgg = sc == CODONS.index('TGA'), sc == CODONS.index('TGG')
with open(out, 'w') as fh:
    fh.write('acc\torganism_tax\tTGA_strictR\tTGA_strictW\tTGA_all\tTGA_allR\tTGG_strictW\tTGG_strictR\tArg_strictR\tArg_strictW\n')
    for g, a in enumerate(D['acc']):
        if taxon not in D['tax'][g].split(';'):
            continue
        m = sg == g
        v = [(tga & ok & (cons == R)), (tga & ok & (cons == W)), tga, (tga & (cons == R)), (tgg & ok & (cons == W)),
             (tgg & ok & (cons == R)), (arg & ok & (cons == R)), (arg & ok & (cons == W))]
        fh.write(f"{a}\t{D['tax'][g].split(';')[-1]}\t" + '\t'.join(str(int((x & m).sum())) for x in v) + '\n')
print(open(out).read())
