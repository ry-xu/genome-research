"""T2 conservation: global NW identity of the 300 nt upstream of mature 16S between genera (one locus per genus),
div vs norm; background = same pairs after dinucleotide shuffling (20 per pair)."""
import itertools, random, sys
sys.path.insert(0, 'src')
from common import dishuffle
REP = {'div': {'Aquicultor_A': 'JAACTV010000069.1', 'CADDYS01': 'CADDYS010000032.1', 'DASQQM01': 'DASQQM010000255.1',
               'DASXVE01': 'DASXVE010000023.1', 'JASEIR01': 'JBFMGP010000004.1', 'Kmv33': 'JACUUO010000012.1'},
       'norm': {'Aquicultor_A': 'JAACTV010000080.1', 'CADDYS01': 'DAHWGC010000001.1', 'DASXVE01': 'DASXVE010000010.1',
                'JASEIR01': 'JBFMGP010000055.1', 'Kmv33': 'JACUUO010000074.1'}}
seq, n = {}, None
for l in open('out/t2_upstream.fa'):
    if l[0] == '>': n = l[1:].split('|')[1]
    else: seq[n] = l.strip()
def nw(a, b, M=1, X=-1, G=-2):
    la, lb = len(a), len(b)
    prev = [j * G for j in range(lb + 1)]; tb = [[0] * (lb + 1) for _ in range(la + 1)]
    for j in range(1, lb + 1): tb[0][j] = 2
    for i in range(1, la + 1):
        cur = [i * G] + [0] * lb; tb[i][0] = 1
        for j in range(1, lb + 1):
            d = prev[j - 1] + (M if a[i - 1] == b[j - 1] else X); u = prev[j] + G; l = cur[j - 1] + G
            if d >= u and d >= l: cur[j] = d; tb[i][j] = 0
            elif u >= l: cur[j] = u; tb[i][j] = 1
            else: cur[j] = l; tb[i][j] = 2
        prev = cur
    i, j, m, cols = la, lb, 0, 0
    while i > 0 or j > 0:
        t = tb[i][j]
        if t == 0: m += a[i - 1] == b[j - 1]; i -= 1; j -= 1
        elif t == 1: i -= 1
        else: j -= 1
        cols += 1
    return m / cols
rng = random.Random(3)
for cls, d in REP.items():
    obs, bg = [], []
    for (g1, c1), (g2, c2) in itertools.combinations(d.items(), 2):
        a, b = seq[c1], seq[c2]
        obs.append(nw(a, b)); bg.extend(nw(dishuffle(a, rng), dishuffle(b, rng)) for _ in range(5))
        print(f'{cls}\t{g1}\t{g2}\t{obs[-1]:.3f}')
    print(f'## {cls}: mean pairwise id {sum(obs)/len(obs):.3f} (n={len(obs)}, range {min(obs):.3f}-{max(obs):.3f}); shuffled background mean {sum(bg)/len(bg):.3f} (95th pct {sorted(bg)[int(.95*len(bg))]:.3f})')
