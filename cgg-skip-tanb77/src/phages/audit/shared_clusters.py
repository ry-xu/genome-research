# Within viral clusters (vFAM / vGENUS) that contain both TANB77-host phages and R phages: CGG skip/O by host group
import sys, csv, collections
pre, metap = sys.argv[1:3]
meta = {r['id']: r for r in csv.DictReader(open(metap), delimiter='\t')}
grp = {l.split('\t')[0]: l.split('\t')[1] for l in open('phages_work/data/groups/run_ids.tsv') if not l.startswith('id\t')}
O = collections.Counter(); K = collections.Counter()
for r in csv.DictReader(open(pre + '.cgg_sites.tsv'), delimiter='\t'):
    O[r['id']] += 1; K[r['id']] += int(r['skip'])
TAN = {'T1', 'T0', 'Tf', 'Tx'}
for col in ('vgenus', 'vfam'):
    cl = collections.defaultdict(set)
    for i, s in grp.items():
        key = meta[i][col]
        if key.endswith('NA'): continue
        cl[key].add(s)
    shared = {k for k, v in cl.items() if v & TAN and 'R' in v}
    t = collections.Counter(); totals = collections.Counter()
    for i, s in grp.items():
        if meta[i][col] in shared:
            g = 'TAN' if s in TAN else s
            t[(g, 'k')] += K[i]; t[(g, 'o')] += O[i]; t[(g, 'n')] += 1
            lowgc = float(meta[i]['gc']) < 30
            if lowgc: t[(g, 'k_low')] += K[i]; t[(g, 'o_low')] += O[i]; t[(g, 'n_low')] += 1
    print(f'{col}: shared clusters {len(shared)}')
    for g in ('TAN', 'R'):
        print(f"   {g}: phages {t[(g,'n')]} skip/O {t[(g,'k')]}/{t[(g,'o')]}; GC<30%: phages {t[(g,'n_low')]} skip/O {t[(g,'k_low')]}/{t[(g,'o_low')]}")
    for k in sorted(shared):
        rows = [(grp[i], K[i], O[i], meta[i]['gc']) for i in grp if meta[i][col] == k]
        tk = sum(r[1] for r in rows if r[0] in TAN); to = sum(r[2] for r in rows if r[0] in TAN)
        rk = sum(r[1] for r in rows if r[0] == 'R'); ro = sum(r[2] for r in rows if r[0] == 'R')
        rgc = sorted(float(r[3]) for r in rows if r[0] == 'R')
        if col == 'vfam' or to + ro > 0:
            print(f'     {k}: TAN {tk}/{to} | R {rk}/{ro} (R n={len(rgc)}, R GC {rgc[0]:.1f}-{rgc[-1]:.1f})')
