"""Recount the headline numbers from the tables in data/census/ (no census files needed).
usage: uv run python -I src/check_numbers.py"""
import csv
from collections import Counter, defaultdict
from scipy.stats import fisher_exact

FOUR = ['BK074532.1', 'BK074538.1', 'BK074542.1', 'BK074262.1']


def rows(path, skip=0):
    with open(path) as f:
        for _ in range(skip):
            next(f)
        yield from csv.DictReader(f, delimiter='\t')


# 1. CUU/CUA at strict columns, by the column's consensus residue (src/cun_pro/cun_by_res.py)
by = defaultdict(Counter)
for r in rows('data/census/cun_sites.tsv'):
    if r['strict'] == '1':
        by[r['acc']][r['consensus']] += 1
tot = Counter()
for a in FOUR:
    tot.update(by[a])
    print(a, 'CUU/CUA at strict columns', sum(by[a].values()), dict(by[a].most_common()))
print('all four: at proline', tot['P'], 'of', sum(tot.values()), dict(tot.most_common()))

# 2. Codon counts at strict Pro / Leu columns (src/codes_census/pro_leu.py)
pl = {r['acc']: r for r in rows('data/census/pro_leu_wss.tsv')}
for a in FOUR:
    r = pl[a]
    print(a, 'Pro: CUN', r['P:CTN'], 'CCN', r['P:CCN'], '| Leu: CUN', r['L:CTN'], 'UUA/UUG', r['L:TTR'])
others = [r for a, r in pl.items() if r['genus_species'] == 'Wickerhamiella' and a not in FOUR]
cc = sorted(int(r['P:CCN']) for r in others)
cn = sorted(int(r['P:CTN']) for r in others)
print(f'other Wickerhamiella: {len(others)} genomes, CCN at Pro {cc[0]}-{cc[-1]}, CUN at Pro {cn[0]}-{cn[-1]}')

# 3. One-sided Fisher tests, recomputed from the counts (src/codes_census/codon_contrast.py)
for codon in ('CTA', 'CTT'):
    ps = {}
    for r in rows(f'data/census/cun_{codon}.tsv', skip=1):
        if r['acc'] == 'pooled':
            continue
        t = [[int(r['codon@Y']), int(r['codon@STD'])], [int(r['syn@Y']), int(r['syn@STD'])]]
        ps[r['acc']] = fisher_exact(t, alternative='greater')[1]
    four = [f'{ps[a]:.1e}' for a in FOUR]
    rest = [p for a, p in ps.items() if a not in FOUR]
    print(f'{codon} at Pro vs Leu, one-sided Fisher: the four {four}; {len(rest)} other genomes, lowest p {min(rest):.3f}')

# 4. All codons at strict proline columns, per genome (src/cun_pro/export_sites.py)
for r in rows('data/census/strict_pro_codons.tsv'):
    print(r['acc'], r['species'], 'strict proline positions', r['strict_pro_sites'], ': CUU/CUA', r['CUU_CUA'],
          'CCN', r['CCN'], 'other', r['other'])
