"""T2 context: nearest predicted protein gene upstream of each 16S (3 kb window, rRNA sense) and downstream of 23S/5S;
its strand relative to the operon and the gap. Also the promoter-score null (200 dinucleotide shuffles of the 600 nt)."""
import random, sys
import pyrodigal
sys.path.insert(0, 'src')
from common import operons, region, dishuffle
import t2_upstream_lib as U
gf = pyrodigal.GeneFinder(meta=True)
rng = random.Random(2)
print('genome\tcontig\tcls\tup_gene_gap\tup_gene_strand\tup_gene_len_aa\tdown_gene_gap\tdown_gene_strand\tprom_score\tprom_p\tboxA_mm\tboxA_p')
for r in operons():
    if r['s16a'] is None: continue
    up = ''
    for L in range(3000, 199, -100):
        up = region(r, r['s16a'] - L, r['s16a'] - 1)
        if up: break
    ug = ('', '', '')
    if len(up) >= 300:
        gs = [g for g in gf.find_genes(up.encode())]
        if gs:
            g = max(gs, key=lambda g: g.end)
            ug = (len(up) - g.end, '+' if g.strand == 1 else '-', len(g.translate()) - 1)
        else: ug = (f'>{len(up)}', '', '')
    endp = r['r5b'] or r['l23b']
    dg = ('', '')
    if endp:
        dn = ''
        for L in range(3000, 199, -100):
            dn = region(r, endp + 1, endp + L)
            if dn: break
        if len(dn) >= 300:
            gs = list(gf.find_genes(dn.encode()))
            if gs:
                g = min(gs, key=lambda g: g.begin); dg = (g.begin - 1, '+' if g.strand == 1 else '-')
            else: dg = (f'>{len(dn)}', '')
    u6 = region(r, r['s16a'] - 600, r['s16a'] - 1)
    if u6:
        ps = U.best_prom(u6)[0]; bx = U.best_boxa(u6)[0]
        nps, nbx = [], []
        for _ in range(200):
            sh = dishuffle(u6, rng); nps.append(U.best_prom(sh)[0]); nbx.append(U.best_boxa(sh)[0])
        pp = (sum(x >= ps for x in nps) + 1) / 201; pb = (sum(x <= bx for x in nbx) + 1) / 201
    else: ps = bx = pp = pb = ''
    print('\t'.join(map(str, [r['genome'], r['contig'], r['cls'], *ug, *dg, ps, f'{pp:.3f}' if pp != '' else '', bx, f'{pb:.3f}' if pb != '' else ''])))
