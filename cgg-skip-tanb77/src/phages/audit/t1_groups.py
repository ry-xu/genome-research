# Task 1: rebuild T1/T0/Tf/R from UHGV metadata + host skip calls; compare with run_ids.tsv
import sys, csv, collections, math
meta_path, calls_path, runids_path, outp = sys.argv[1:5]
calls = collections.defaultdict(list)
for r in csv.DictReader(open(calls_path), delimiter='\t'):
    calls[r['genus']].append(int(r['call']))
skip_genera = {g for g, v in calls.items() if all(x == 1 for x in v)}
mixed = {g for g, v in calls.items() if len(set(v)) > 1}
print('skipper genera', len(skip_genera), 'mixed', mixed)
import subprocess
p = subprocess.Popen(['zstd', '-dcq', meta_path], stdout=subprocess.PIPE, text=True)
rd = csv.DictReader(p.stdout, delimiter='\t')
mine = {}; info = {}
for r in rd:
    q = r['checkv_quality']; gcode = r['genetic_code']
    h = r['crispr_host_lineage_gtdb_r207']
    if q not in ('Complete', 'High-quality', 'Medium-quality', 'High', 'Medium'):
        continue
    if gcode != '11': continue
    if 'c__Clostridia' not in h: continue
    try: L = int(r['genome_length'])
    except: continue
    ranks = dict(x.split('__', 1) for x in h.split(';') if '__' in x)
    gc = float(r['gc_percent']); gbin = int(math.floor(gc / 2) * 2)
    s = None
    if ranks.get('o') == 'TANB77':
        g = ranks.get('g'); f = ranks.get('f')
        if g in skip_genera: s = 'T1'
        elif g in ('CAG-492', 'CAG-465'): s = 'T0'
        elif g is None and f == 'CAG-508': s = 'Tf'
        else: s = 'Tx_other'
    elif ranks.get('o'):
        s = 'R'
    if s:
        mine[r['uhgv_genome']] = s
        info[r['uhgv_genome']] = (gbin, L, q, h)
qs = collections.Counter()
run = {}
for r in csv.DictReader(open(runids_path), delimiter='\t'):
    run[r['id']] = (r['set'], int(r['gc_bin']))
print('my sets (all lengths):', collections.Counter(mine.values()))
print('my sets (>=5kb):', collections.Counter(s for k, s in mine.items() if info[k][1] >= 5000))
# agreement
for S in ['T1', 'T0', 'Tf', 'R']:
    a = {k for k, v in mine.items() if v == S and info[k][1] >= 5000}
    b = {k for k, v in run.items() if v[0] == S}
    print(f'{S}: mine={len(a)} run={len(b)} both={len(a & b)} mine_only={len(a - b)} run_only={len(b - a)}')
    ro = list(b - a)[:5]
    for k in ro: print('   run_only ex', k, mine.get(k), info.get(k))
    mo = list(a - b)[:5]
    for k in mo: print('   mine_only ex', k, run.get(k), info.get(k))
# gc bin agreement
dis = sum(1 for k, v in run.items() if k in info and info[k][0] != v[1])
print('gc_bin disagreements among overlapping ids:', dis, 'of', sum(1 for k in run if k in info))
# R per bin: pool vs run
pool = collections.Counter(info[k][0] for k, v in mine.items() if v == 'R' and info[k][1] >= 5000)
rr = collections.Counter(v[1] for k, v in run.items() if v[0] == 'R')
print('R pool vs run per bin:', sorted((b, pool[b], rr[b]) for b in set(pool) | set(rr)))
with open(outp, 'w') as f:
    f.write('id\tset\tgc_bin\tlength\tquality\thost\n')
    for k, v in mine.items():
        f.write(f'{k}\t{v}\t{info[k][0]}\t{info[k][1]}\t{info[k][2]}\t{info[k][3]}\n')
