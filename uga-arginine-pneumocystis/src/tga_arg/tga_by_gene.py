"""TGA per gene model for given accessions, with column consensus (strict or not), from raw census npz files.
usage: python tga_by_gene.py acc[,acc...] census.npz [...]"""
import sys, os
from collections import Counter, defaultdict
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'codes_census'))
from sitescan_org import load, AA, CODONS

accs, *paths = sys.argv[1:]
accs = accs.split(',')
D = load(paths)
TGA, TGG = CODONS.index('TGA'), CODONS.index('TGG')
for p in paths:
    z = np.load(p, allow_pickle=True)
    A = list(z['acc']); models = z['models']
    for a in accs:
        if a not in A: continue
        g = A.index(a); m = z['site_genome'] == g
        out = defaultdict(Counter)
        for mod, cod in zip(z['site_model'][m], z['site_codon'][m]):
            if cod in (TGA, TGG): out[str(models[mod])][CODONS[cod]] += 1
        print(a, str(z['taxonomy'][g]).split(';')[-1], {k: dict(v) for k, v in sorted(out.items())})
