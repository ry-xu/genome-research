"""ARAGORN (pyaragorn, codon-hunt's settings) on the pig TANB77 representatives: tRNAs with 6- or 8-nt anticodon loops.
usage: python -I aragorn_loops.py <fna...>"""
import sys
from collections import Counter

import pyaragorn

f = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
cnt = Counter()
for path in sys.argv[1:]:
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    for c, s in seqs.items():
        for g in f.find_rna("".join(s)):
            cnt[g.anticodon_length] += 1
            if g.anticodon_length != 3:
                print(path.split("/")[-1], c, g.begin, g.end, g.strand, g.amino_acids, g.anticodon, g.energy,
                      g.sequence(), sep="\t")
print("anticodon lengths:", dict(cnt))
