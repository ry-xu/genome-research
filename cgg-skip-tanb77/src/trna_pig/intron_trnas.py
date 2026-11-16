"""Identity of the long 'Undet' tRNAscan-SE calls (>= 120 bp, an insertion in the anticodon arm): try every excision
[a, b) with a in 25..45 and a mature length of 73..95, call ARAGORN (codon-hunt's settings) on each mature sequence in
random flanks, keep the best-energy call per tRNA. Also reports group I splice-site signs (5' exon ends in T, intron
ends in G). usage: python -I intron_trnas.py <ss files...>"""
import random
import re
import sys
from collections import Counter

import pyaragorn

random.seed(1)
FL = "".join(random.choices("ACGT", weights=[55, 45, 45, 55], k=800))
f = pyaragorn.RNAFinder(11, tmrna=False, linear=True)
res = Counter()
for path in sys.argv[1:]:
    for b in open(path).read().split("\n\n"):
        m = re.search(r"^(\S+) \((\d+)-(\d+)\)\s+Length: (\d+).*\nType: Undet", b, re.M)
        if not m or int(m[4]) < 120:
            continue
        s = re.search(r"^Seq: (\S+)", b, re.M)[1].upper()
        best = None
        for a in range(25, 46):
            for e in range(a + 1, len(s)):
                L = len(s) - (e - a)
                if not 73 <= L <= 95:
                    continue
                mat = s[:a] + s[e:]
                for g in f.find_rna(FL + mat + FL):
                    if g.strand == 1 and abs(g.begin - len(FL) - 1) <= 3 and (best is None or g.energy > best[0]):
                        best = (g.energy, a, e, g.amino_acids, g.anticodon, mat)
        if best:
            en, a, e, aa, ac, mat = best
            res[(aa if isinstance(aa, str) else "/".join(aa), ac)] += 1
            print(path.split("/")[-1].split(".")[0], m[1], len(s), f"intron {a}-{e} ({e - a} nt)",
                  f"5'exon ends {s[a - 1]}, intron ends {s[e - 1]}", aa, ac, round(en, 1), mat, sep="\t")
        else:
            print(path.split("/")[-1].split(".")[0], m[1], len(s), "no ARAGORN call", sep="\t")
print(res)
