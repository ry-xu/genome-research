import sys, csv
R = 'results/cun_pro/'  # the lab tree's folder; its reads tables are data/reads/ here
S = sys.argv[1]
for a in ['BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1']:
    keys = {tuple(l.split()) for l in open(f'{S}/{a}.census_ctn_p.txt') if l.strip()}
    for d in ('reads', 'reads_mapq20'):
        rows = {}
        for r in csv.DictReader(open(f'{R}{d}/{a}.tsv'), delimiter='\t'):
            k = (r['codon'], r['model_col'])
            if k in keys:
                # several ORF copies can carry the same site; keep the one with most reads
                if k not in rows or int(r['reads']) > int(rows[k]['reads']):
                    rows[k] = r
        cov = [r for r in rows.values() if int(r['reads']) > 0]
        n = sum(int(r['reads']) for r in rows.values()); s = sum(int(r['same']) for r in rows.values())
        sh = sorted((int(r['same']) / int(r['reads']), r['model_col'], r['codon'], r['same'], r['reads']) for r in cov if int(r['reads']) >= 3)
        extra = sorted(keys - set(rows))
        print(a, d, 'census sites', len(keys), 'found', len(rows), 'covered', len(cov), 'reads', n, 'same', s,
              f'{100*s/max(n,1):.2f}%', 'min(>=3)', sh[0] if sh else None, 'missing', extra)
