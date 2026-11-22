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
o = open('out/t2_upstream.tsv', 'w'); fa = open('out/t2_upstream.fa', 'w')
o.write('genome\tcontig\tcls\tup_len\tboxA_mm\tboxA_pos_from16S\tboxA\tprom_score\tm35\tspacer\tm10\text10\tUP_AT\tdisc7\tdisc_GC\tprom_m10_end_from16S\tspacer_len\tspacer_boxA_mm\tspacer_boxA\n')
seen = set()
for r in operons():
    if r['s16a'] is None: continue
    up = ''
    for L in range(UP, 49, -10):
        up = region(r, r['s16a'] - L, r['s16a'] - 1)
        if up: break
    if not up: continue
    k = (r['genome'], r['contig'], r['s16a'] + r['lo'])
    bx = best_boxa(up); pr = best_prom(up)
    sp = region(r, r['s16b'] + 1, r['l23a'] - 1) if (r['s16b'] and r['l23a']) else ''
    sbx = best_boxa(sp) if len(sp) >= 12 else None
    o.write('\t'.join(map(str, [r['genome'], r['contig'], r['cls'], len(up), bx[0], len(up) - bx[1], bx[2], pr[0], pr[2], pr[3], pr[4], pr[5], pr[6], pr[7], pr[8], len(up) - pr[9],
                                len(sp) if sp else '', sbx[0] if sbx else '', sbx[2] if sbx else ''])) + '\n')
    fa.write(f">{r['genome']}|{r['contig']}|{r['cls']}\n{up[-300:]}\n")
