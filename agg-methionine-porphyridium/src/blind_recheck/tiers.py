"""Apply the frozen calling rule to Codetta + Method 2 (v2 or v3) results."""
import sys, json
T = sys.argv[1]; ver = sys.argv[2]; names = sys.argv[3:]
B='TCAG'; codons=[a+b+c for a in B for b in B for c in B]
aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; STD=dict(zip(codons,aas))
cd_all = json.load(open(T + '/work/codetta_parsed.json'))
out = {}
for g in names:
    m = json.load(open('%s/work/m2/%s.%s.json' % (T, g, 'm2' if ver == 'v2' else 'm2v3')))['pass2']
    cd = cd_all[g]
    rows = []
    for c in codons:
        std = STD[c]; r = m[c]; n = r['n']; ll = r['ll']; X2 = r['best']
        cdl = cd[c]['inf']
        def m2_strong(X):
            if X != r['best'] or X == std: return False
            if std == '*': return n >= 10 and r['llr_2nd'] >= 3
            return n >= 10 and ll[X] - ll[std] >= 10 and r['llr_2nd'] >= 3
        def m2_support(X):
            if X != r['best'] or X == std: return False
            if std == '*': return n >= 5 and r['llr_2nd'] >= 1
            return n >= 5 and ll[X] - ll[std] >= 3
        def contradicted(X):
            if cdl != '?' and cdl != X: return True
            if r['best'] != X and n >= 5 and ll[r['best']] - ll[X] >= 3: return True
            return False
        cands = set()
        if X2 != std and m2_support(X2): cands.add(X2)
        if cdl not in ('?', std): cands.add(cdl)
        tier = None; X = None
        for X in sorted(cands):
            if not contradicted(X) and ((m2_strong(X) and cdl in (X, '?')) or (cdl == X and m2_support(X))):
                tier = 'call'; break
        if tier is None and cands:
            tier = 'weak'; X = '/'.join(sorted(cands))
        if tier is None and r.get('absent_flag'):
            tier = 'weak-absent'; X = 'stop/unassigned?'
        if tier is None and std == '*':
            tier = 'no call' if (cdl == '?' and n < 5) else 'no sense evidence'
            X = '-'
        if tier is None: continue
        cdd = cd[c]
        gs = cdd.get('gap_std'); g2 = cdd.get('gap_2nd')
        rows.append(dict(codon=c, std=std, X=X, tier=tier, m2_n=n, m2_best=r['best'], m2_2nd=r['second'],
                         m2_llr_std=r['llr_std'], m2_llr_2nd=r['llr_2nd'], m2_top=r['top_ref'],
                         cd_inf=cdl, cd_nused=cdd['nused'], cd_best=cdd.get('best'), cd_gap_std=gs, cd_gap_2nd=g2,
                         others_n=r.get('others_n')))
    out[g] = rows
    for x in rows:
        ls = '%.1f' % x['m2_llr_std'] if x['m2_llr_std'] is not None else 'n/a'
        cs = '%.1f' % x['cd_gap_std'] if x['cd_gap_std'] is not None else 'n/a'
        print('%s %s %s->%s [%s] M2 n=%d best=%s LLRstd=%s LLR2nd=%.1f top=(%s) | CD %s Nused=%d best=%s gapstd=%s gap2nd=%.1f' % (
            g, x['codon'], x['std'], x['X'], x['tier'], x['m2_n'], x['m2_best'], ls, x['m2_llr_2nd'], x['m2_top'],
            x['cd_inf'], x['cd_nused'], x['cd_best'], cs, x['cd_gap_2nd'] or 0))
json.dump(out, open('%s/work/tiers_%s.json' % (T, ver), 'w'), indent=1)
