"""T1: RNase III processing stems. RNAduplex(leader 200 nt, trailer 200 nt) per 16S and 23S; null = 200 dinucleotide
shuffles of each. Output out/t1c_helix.tsv."""
import random, sys
import RNA
sys.path.insert(0, 'src')
from common import operons, region, dishuffle
FL, NS = 200, 200
def duplex(a, b):
    d = RNA.duplexfold(a.replace('T', 'U'), b.replace('T', 'U'))
    left, right = d.structure.split('&')
    a0 = d.i - len(left) + 1; b0 = d.j
    ia = [a0 + k for k, ch in enumerate(left) if ch == '(']
    ib = [b0 + k for k, ch in enumerate(right) if ch == ')'][::-1]
    pairs = list(zip(ia, ib))
    best, cur = 0, 1
    for (x1, y1), (x2, y2) in zip(pairs, pairs[1:]):
        cur = cur + 1 if (x2 - x1 <= 2 and y1 - y2 <= 2) else 1
        best = max(best, cur)
    best = max(best, 1 if pairs else 0)
    return d.energy, len(pairs), best, (min(ia), max(ia), min(ib), max(ib)) if pairs else None, d.structure
rng = random.Random(1)
out = open('out/t1c_helix.tsv', 'w')
out.write('genome\tcontig\tcls\trrna\tdG\tbp\thelix\tdist5\tdist3\tp_null\tnull_med\tcall\tstructure\tleader_seg\ttrailer_seg\n')
for r in operons():
    for rr, a, b in (('16S', r['s16a'], r['s16b']), ('23S', r['l23a'], r['l23b'])):
        if a is None or b is None: continue
        L, T = region(r, a - FL, a - 1), region(r, b + 1, b + FL)
        if len(L) < FL or len(T) < FL or 'N' in L + T: 
            out.write(f"{r['genome']}\t{r['contig']}\t{r['cls']}\t{rr}\t\t\t\t\t\t\t\tnot_assessable\t\t\t\n"); continue
        dG, bp, hx, sp, st = duplex(L, T)
        null = []
        for _ in range(NS):
            null.append(duplex(dishuffle(L, rng), dishuffle(T, rng))[2])
        p = (sum(x >= hx for x in null) + 1) / (NS + 1)
        d5, d3 = FL - sp[1], sp[2] - 1
        call = 'stem' if (p <= 0.05 and bp >= 15 and d5 <= 150 and d3 <= 150) else 'no_stem'
        null.sort()
        out.write(f"{r['genome']}\t{r['contig']}\t{r['cls']}\t{rr}\t{dG:.1f}\t{bp}\t{hx}\t{d5}\t{d3}\t{p:.3f}\t{null[NS//2]:.1f}\t{call}\t{st}\t{L[sp[0]-1:sp[1]]}\t{T[sp[2]-1:sp[3]]}\n")
    out.flush()
