"""T2: promoter boxes, UP element, discriminator, boxA in 600 nt upstream of the mature 16S 5' end; boxA also in the
16S-23S spacer. Output out/t2_upstream.tsv and out/t2_upstream.fa (300 nt upstream, for conservation)."""
import sys
sys.path.insert(0, 'src')
from common import operons, region
UP = 600
BOXA = 'TGCTCTTTAACA'
def hd(a, b): return sum(x != y for x, y in zip(a, b))
def m(a, b): return sum(x == y for x, y in zip(a, b))
def best_boxa(s):
    best = None
    for i in range(len(s) - 12 + 1):
        d = hd(s[i:i + 12], BOXA)
        if best is None or d < best[0]: best = (d, i, s[i:i + 12])
    return best
def best_prom(s):
    best = None
    for i in range(len(s) - 6):
        m35 = m(s[i:i + 6], 'TTGACA')
        for sp in (16, 17, 18):
            j = i + 6 + sp
            if j + 6 + 8 > len(s): continue
            m10 = m(s[j:j + 6], 'TATAAT')
            sc = m35 + m10
            if best is None or sc > best[0]:
                upel = s[max(0, i - 20):i]; at = (upel.count('A') + upel.count('T')) / max(1, len(upel))
                disc = s[j + 6:j + 6 + 7]
                ext = s[j - 3:j - 1] == 'TG'
                best = (sc, i, s[i:i + 6], sp, s[j:j + 6], ext, round(at, 2), disc, disc.count('G') + disc.count('C'), j + 6)
    return best
