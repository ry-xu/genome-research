"""ARAGORN (pyaragorn, code 11, linear) tally over genomes: anticodon-loop lengths, Arg/Gly anticodons, any CCG/TCG
anticodon, per-genome presence. usage: python -I aragorn_tally.py <fna...>"""
import sys
from collections import Counter

import pyaragorn

f = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
loops, ac, n_g, odd = Counter(), Counter(), 0, []
for path in sys.argv[1:]:
    n_g += 1
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    seen = set()
    for c, s in seqs.items():
        for g in f.find_rna("".join(s)):
            loops[g.anticodon_length] += 1
            aa = g.amino_acids if isinstance(g.amino_acids, str) else "/".join(g.amino_acids)
            a = g.anticodon.upper()
            if g.anticodon_length != 3 or a in ("CCG", "TCG") or "Arg" in aa or "Gly" in aa:
                seen.add((aa, a, g.anticodon_length))
            if g.anticodon_length != 3 or a in ("CCG", "TCG"):
                odd.append((path.split("/")[-1], c, g.begin, aa, a, g.anticodon_length))
    for x in seen:
        ac[x] += 1
print(f"genomes {n_g}; anticodon lengths over all calls {dict(loops)}")
print("genomes with each Arg/Gly/odd tRNA:", *sorted(ac.items()), sep="\n  ")
print("odd calls:", *odd, sep="\n  ")
