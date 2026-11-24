"""RULES.md §5: rank screened runs by spanning templates per GB, then lay out the main list within the 8 GB RNA budget.
usage: plan_main.py  -> meta/screen_rank.tsv, meta/main_list.txt (RUN SPOTS), prints the budget"""
import csv, subprocess, glob, os
os.chdir('/Users/ryan/lab/work/tack-rna')

def info(run):
    z = q = 0
    for l in subprocess.run(['vdb-dump', '--info', run], capture_output=True, text=True).stdout.splitlines():
        k, _, v = l.partition(':')
        if k.strip() == 'size':
            z = int(v.strip().replace(',', ''))
        elif k.strip() == 'SEQ':
            q = int(v.strip().replace(',', ''))
    return z, q

def streamed(d):
    return sum(int(open(f).read().split()[0]) for f in glob.glob(f'{d}/*.bytes') if os.path.getsize(f))

rows = []
for r in csv.DictReader(open('meta/screen_runs.tsv'), delimiter='\t'):
    run = r['run']
    c = f'out/screen/{run}.counts.tsv'
    if not os.path.exists(f'out/screen/{run}.done'):
        rows.append(dict(study=r['study'], run=run, spanning=-1, gb=0, per_gb=-1)); continue
    sp = sum(int(x['joined']) + int(x['through']) + int(x['other']) for x in csv.DictReader(open(c), delimiter='\t'))
    gb = int(open(f'out/screen/{run}.bytes').read()) / 1e9
    rows.append(dict(study=r['study'], run=run, spanning=sp, gb=round(gb, 4), per_gb=round(sp / gb, 1) if gb else 0))
rows.sort(key=lambda x: -x['per_gb'])
with open('meta/screen_rank.tsv', 'w') as f:
    f.write('study\trun\tspanning\tgb\tper_gb\n')
    for x in rows:
        f.write(f"{x['study']}\t{x['run']}\t{x['spanning']}\t{x['gb']}\t{x['per_gb']}\n")
spent = streamed('out/pc') + streamed('out/pc2') + streamed('out/screen')
print(f'streamed so far (pc, pc2, screen): {spent / 1e9:.3f} GB')
main, tot = [], spent
for st, bp in (('S2', 'PRJNA916083'), ('S3', 'PRJNA1163300')):
    for r in csv.DictReader(open(f'meta/runinfo_{bp}.csv')):
        if 'TRANSCRIPT' not in r['LibrarySource'] and r['LibraryStrategy'] != 'RNA-Seq':
            continue
        z, q = info(r['Run'])
        b = z * min(1, 2_000_000 / q)
        main.append((st, r['Run'], 2_000_000, b)); tot += b
for x in rows:
    if x['study'] != 'S1' or x['per_gb'] <= 0:
        continue
    z, q = info(x['run'])
    b = z * min(1, 15_000_000 / q)
    if tot + b > 8e9:
        print(f"stop before {x['run']}: would reach {(tot + b) / 1e9:.2f} GB")
        break
    main.append(('S1', x['run'], 15_000_000, b)); tot += b
with open('meta/main_list.txt', 'w') as f:
    for st, run, sp, b in main:
        f.write(f'{run} {sp}\n')
for st, run, sp, b in main:
    print(st, run, sp, f'{b / 1e9:.3f} GB')
print(f'planned RNA total {tot / 1e9:.3f} GB')
