"""T1b (RULES A1): RNAfold MFE of leader200 + rRNA + trailer200; count leader-trailer pairs, longest helix, distance of
the pairing to the mature ends. Usage: t1b_fold.py 16S|23S -> out/t1b_<rrna>.tsv"""
import sys
import RNA
sys.path.insert(0, 'src')
from common import operons, region
FL = 200
CACHE = {}
rr = sys.argv[1]
import os
DONE = set()
if os.path.exists(f'out/t1b_{rr}.tsv'):
    for l in open(f'out/t1b_{rr}.tsv'): DONE.add(tuple(l.split('\t')[:2]))
out = open(f'out/t1b_{rr}.tsv', 'a' if DONE else 'w')
if not DONE: out.write('genome\tcontig\tcls\trrna\tlen\tlt_bp\thelix\tdist5\tdist3\tcall\tleader_db\ttrailer_db\n')
for r in operons():
    if (r['genome'], r['contig']) in DONE: continue
    a, b = (r['s16a'], r['s16b']) if rr == '16S' else (r['l23a'], r['l23b'])
    if a is None or b is None: continue
    s = region(r, a - FL, b + FL)
    if len(s) != (b - a + 1) + 2 * FL or 'N' in s:
        out.write(f"{r['genome']}\t{r['contig']}\t{r['cls']}\t{rr}\t\t\t\t\t\tnot_assessable\t\t\n"); out.flush(); continue
    N = len(s)
    if s not in CACHE: CACHE[s] = RNA.fold(s.replace('T', 'U'))
    st, mfe = CACHE[s]
    pt = RNA.ptable(st)
    pairs = sorted((i, pt[i]) for i in range(1, FL + 1) if pt[i] > N - FL)
    hx, cur = (1 if pairs else 0), 1
    for (x1, y1), (x2, y2) in zip(pairs, pairs[1:]):
        cur = cur + 1 if (x2 - x1 <= 2 and y1 - y2 <= 2) else 1
        hx = max(hx, cur)
    if pairs:
        d5 = FL - max(x for x, _ in pairs); d3 = min(y for _, y in pairs) - (N - FL) - 1
    else: d5 = d3 = ''
    call = 'stem' if (len(pairs) >= 15 and hx >= 8 and d5 <= 150 and d3 <= 150) else 'no_stem'
    out.write(f"{r['genome']}\t{r['contig']}\t{r['cls']}\t{rr}\t{N}\t{len(pairs)}\t{hx}\t{d5}\t{d3}\t{call}\t{st[:FL]}\t{st[-FL:]}\n")
    out.flush()
