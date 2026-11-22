import csv, sys, ast
for f in ('BK074411.1', 'BK074391.1', 'BK074449.1', 'BK074121.1'):
    rows = list(csv.DictReader(open(f'lab/work/new-stops-euk/results/aga_ala/reads/{f}.tsv'), delimiter='\t'))
    aga = [r for r in rows if r['codon'] == 'AGA']
    ala_strict = [r for r in aga if r['cons'] == 'A' and r['strict'] == '1']
    ala_any = [r for r in aga if r['cons'] == 'A']
    def tally(rs):
        n = sum(int(r['reads']) for r in rs); s = sum(int(r['same']) for r in rs)
        maj = sum(1 for r in rs if int(r['reads']) > 0 and int(r['same']) > int(r['reads']) / 2)
        cov = sum(1 for r in rs if int(r['reads']) > 0)
        mind = min((int(r['reads']) for r in rs), default=0)
        return f'{len(rs)} sites, covered {cov}, AGA majority {maj}, reads {n}, agree {s} ({s/max(n,1):.4f}), min depth {mind}'
    print(f, 'all AGA:', tally(aga))
    print('   AGA at Ala-consensus (any):', tally(ala_any))
    print('   AGA at strict Ala:', tally(ala_strict))
    others = [(r['codon_index'], r['reads'], r['same'], r['other']) for r in aga if r['other'] not in ('{}', '')]
    print('   sites with other triplets:', len(others), others[:6])
