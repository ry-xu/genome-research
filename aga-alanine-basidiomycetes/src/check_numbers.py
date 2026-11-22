"""Recount the README's numbers from the tables in ../data/ (no census needed).
usage (from the folder root): uv run python -I src/check_numbers.py"""
import csv
import os
import re

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def tsv(*p):
    return list(csv.DictReader(open(os.path.join(D, *p)), delimiter='\t'))


# Per-genome AGA placement and tRNA(UCU) acceptor (src/aga_ala/trna_vs_usage.py)
rows = tsv('census', 'trna_vs_usage_fungi.tsv')
ala = [r for r in rows if r['class'] == 'ala']
arg = [r for r in rows if r['class'] == 'arg']
print(f'fungal mitochondrial genomes in the census: {len(rows):,}')
print(f'AGA at >= 1 strict Ala column and 0 strict Arg columns: {len(ala)} genomes, '
      f'{sum(int(r["AGA_at_Ala"]) for r in ala)} AGA-at-Ala sites, '
      f'{sum(int(r["AGA_at_Ala"]) == 1 for r in ala)} on a single site, '
      f'per genome {min(int(r["AGA_at_Ala"]) for r in ala)} to {max(int(r["AGA_at_Ala"]) for r in ala)}')
print(f'  with a G3:U70 tRNA(UCU): {sum(int(r["G3U70"]) > 0 for r in ala)}; '
      f'all of those with A73: {all(int(r["G3U70_A73"]) == int(r["G3U70"]) for r in ala)}')
print(f'AGA at >= 3 strict Arg columns and 0 Ala columns: {len(arg):,} genomes; '
      f'with a tRNA(UCU): {sum(int(r["tRNA_UCU"]) > 0 for r in arg):,}; '
      f'of those with G3:U70: {sum(int(r["G3U70"]) > 0 for r in arg)}')
off = [r for r in rows if r['class'] == 'other' and int(r['G3U70']) > 0 and int(r['AGA_at_Ala']) > 0]
print('G3:U70 carriers off the axis (AGA at Ala and at Arg):',
      ', '.join(f'{r["acc"]} {r["genus"]} {r["AGA_at_Ala"]} Ala / {r["AGA_at_Arg"]} Arg' for r in off))

# Selection-free contrast per order (src/aga_ala/lineage_nosel.py; needs the census, output kept)
print('\nselection-free contrast (distinct strict columns):')
print('\n'.join(l for l in open(os.path.join(D, 'census', 'lineage_nosel.txt')).read().split('\n') if l and not l.startswith('[')))

# All aligned columns (src/aga_ala/all_columns.py; needs the census, output kept)
print('\nall aligned columns:')
for line in open(os.path.join(D, 'census', 'all_columns.txt')):
    if line.startswith(('==', 'AGA', 'CGN', 'GCN')) or 'genomes' in line:
        print(' ', line.rstrip()[:110])

# Raw reads: AGA sites at Ala-consensus columns (src/codes_census/reads_codons.py)
print('\nraw reads, AGA at Ala-consensus columns:')
for acc in ('BK074411.1', 'BK074449.1', 'BK074391.1', 'BK074121.1'):
    s = [r for r in tsv('reads', f'{acc}.tsv') if r['codon'] == 'AGA' and r['cons'] == 'A']
    n, same = sum(int(r['reads']) for r in s), sum(int(r['same']) for r in s)
    maj = sum(int(r['same']) * 2 > int(r['reads']) for r in s)
    print(f'  {acc}: {len(s)} sites, AGA the majority at {maj}, {same:,} of {n:,} reads agree ({same / n:.1%})')

# Proteome tallies (src/prot/tally_aga.py)
print('\nproteome PXD025163, PSMs at spectrum q <= 0.01 (target / decoy):')
tot = {}
for s in ('setA', 'setB', 'setC'):
    t = open(os.path.join(D, 'proteome', s, 'tally.txt')).read()
    c = {m.group(1): (int(m.group(2)), int(m.group(3))) for m in re.finditer(r'^  (\S+)\s+(\d+)\s+(\d+)$', t, re.M)}
    print(f'  {s}: ' + ', '.join(f'{k} {v[0]}/{v[1]}' for k, v in c.items()))
    for k, v in c.items():
        tot[k] = tot.get(k, 0) + v[0]
print(f'  total: CGN = Arg (RA) {tot.get("RA", 0)}, AGA = Ala (A) {tot.get("A", 0)}, '
      f'AGA = Arg (RC) {tot.get("RC", 0)}, CGN = Ala (C) {tot.get("C", 0)}')
