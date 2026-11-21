import csv, sys
R = 'results/cun_pro/'  # the lab tree's folder; its reads tables are data/reads/ here
def check(tsv, keyf, label):
    keys = {tuple(l.split()) for l in open(keyf) if l.strip()}
    rows = {}
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        k = (r['codon'], r['model_col'])
        if k in keys and (k not in rows or int(r['reads']) > int(rows[k]['reads'])): rows[k] = r
    v = [(int(r['reads']), int(r['same']), k) for k, r in rows.items()]
    low = sorted(v)[:4]
    minority = [x for x in v if x[0] and x[1] / x[0] <= 0.5]
    print(label, 'sites', len(keys), 'found', len(rows), 'reads total', sum(x[0] for x in v), 'same', sum(x[1] for x in v),
          'sites <3 reads', sum(1 for x in v if x[0] < 3), 'zero-read', sum(1 for x in v if x[0] == 0), 'majority-not-matching', minority, 'lowest-depth', low)
for a in ['BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1']:
    check(f'{R}reads/{a}.tsv', f'{R}audit/{a}.census_ctn_p.txt', a)
check(f'{R}reads_isolates/BK074538.1_ERR11532760.tsv', f'{R}audit/BK074538.1.census_ctn_p.txt', 'ERR11532760')
