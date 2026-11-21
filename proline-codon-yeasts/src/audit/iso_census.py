import sys, csv
keys = {tuple(l.split()) for l in open(sys.argv[2]) if l.strip()}
rows = {}
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    k = (r['codon'], r['model_col'])
    if k in keys and (k not in rows or int(r['reads']) > int(rows[k]['reads'])):
        rows[k] = r
print('census sites', len(keys), 'found', len(rows), 'covered', sum(int(r['reads']) > 0 for r in rows.values()))
for c in ('CTA', 'CTT'):
    v = [r for r in rows.values() if r['codon'] == c]
    print(c, len(v), sum(int(r['same']) for r in v), '/', sum(int(r['reads']) for r in v))
n = sum(int(r['reads']) for r in rows.values()); s = sum(int(r['same']) for r in rows.values())
print(f'overall {s}/{n} {100*s/n:.2f}%', 'min', min((int(r['same'])/int(r['reads']), r['model_col'], r['same'], r['reads']) for r in rows.values() if int(r['reads'])))
print('extra key', sorted(keys - {(r['codon'], r['model_col']) for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if r['strict'] == '1' and r['cons'] == 'P'}))
