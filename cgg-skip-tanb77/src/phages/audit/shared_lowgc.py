# per shared vFAM, GC<30% phages only: TANB77-host vs R skip/O
import sys, csv, collections
pre, metap = sys.argv[1:3]
meta = {r['id']: r for r in csv.DictReader(open(metap), delimiter='\t')}
grp = {l.split('\t')[0]: l.split('\t')[1] for l in open('phages_work/data/groups/run_ids.tsv') if not l.startswith('id\t')}
O = collections.Counter(); K = collections.Counter()
for r in csv.DictReader(open(pre + '.cgg_sites.tsv'), delimiter='\t'):
    O[r['id']] += 1; K[r['id']] += int(r['skip'])
TAN = {'T1', 'T0', 'Tf', 'Tx'}
c = collections.Counter()
for i, s in grp.items():
    if float(meta[i]['gc']) >= 30 or meta[i]['vfam'].endswith('NA'): continue
    g = 'TAN' if s in TAN else ('R' if s == 'R' else None)
    if g: c[(meta[i]['vfam'], g, 'k')] += K[i]; c[(meta[i]['vfam'], g, 'o')] += O[i]; c[(meta[i]['vfam'], g, 'n')] += 1
fams = sorted({k[0] for k in c})
tot = collections.Counter()
for f in fams:
    if c[(f, 'TAN', 'o')] and c[(f, 'R', 'o')]:
        print(f"{f}: TAN {c[(f,'TAN','k')]}/{c[(f,'TAN','o')]} (n={c[(f,'TAN','n')]}) | R {c[(f,'R','k')]}/{c[(f,'R','o')]} (n={c[(f,'R','n')]})")
        if f != 'vFAM-00025':
            for g in ('TAN', 'R'): tot[(g, 'k')] += c[(f, g, 'k')]; tot[(g, 'o')] += c[(f, g, 'o')]
print('sum over vFAMs with both groups at GC<30, excluding vFAM-00025: TAN %d/%d, R %d/%d' % (tot[('TAN','k')], tot[('TAN','o')], tot[('R','k')], tot[('R','o')]))
