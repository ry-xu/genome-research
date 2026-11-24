"""DNA check, census part (RULES.md §3): is the joined junction (20 + 20 nt, every slide position) present in any census genome?
usage: dnacheck.py EVENTS_FILE OUT   (EVENTS_FILE: one event id per line)
  -> OUT.queries.tsv (event, 40-mer), OUT.hits.tsv (acc, 40-mer), OUT.tsv (event, n_genomes, genomes)"""
import sys, subprocess, csv
from collections import defaultdict
sys.path.insert(0, '/Users/ryan/lab/work/archaeal-intron-census/src')
from bhbset import slide_set

CEN = '/Users/ryan/lab/work/archaeal-intron-census/data'
ids = {l.strip() for l in open(sys.argv[1]) if l.strip()}
out = sys.argv[2]
seqs, name = {}, None
for l in open('/Users/ryan/lab/work/tack-rna/out/regions.fa'):
    if l[0] == '>':
        name = l[1:].strip()
    else:
        seqs[name] = l.strip().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
Q = defaultdict(set)    # 40-mer -> events
for r in csv.DictReader(open('/Users/ryan/lab/work/tack-rna/out/events.tsv'), delimiter='\t'):
    if r['id'] not in ids:
        continue
    s, x, y = seqs[r['region']], int(r['rx']), int(r['ry'])
    pairs, _ = slide_set(s, x, y)
    for a, b in pairs:
        if a >= 20 and b + 20 <= len(s):
            j = s[a - 20:a] + s[b:b + 20]
            Q[j].add(r['id']); Q[j.translate(comp)[::-1]].add(r['id'])
with open(out + '.queries.tsv', 'w') as f:
    for q, es in Q.items():
        for e in es:
            f.write(f'{e}\t{q}\n')
open(out + '.pat', 'w').write('\n'.join(Q) + '\n')
paths = {}
for fn in ('run1.tsv', 'run2.tsv', 'run3.tsv'):
    for l in open(f'{CEN}/{fn}'):
        a = l.rstrip('\n').split('\t')
        paths.setdefault(a[0], a[1])
hits = defaultdict(set)
with open(out + '.hits.tsv', 'w') as f:
    import os
    sel = sorted(paths.items())[:int(os.environ.get("LIMIT", 10**9))] + [(a, paths[a]) for a in os.environ.get("ALSO", "").split(",") if a]
    for i, (acc, p) in enumerate(sel):
        cmd = f"set -o pipefail; gzcat {p} | sed 's/^>.*/N/' | tr -d '\\n' | rg -o -i -F -f {out}.pat | tr a-z A-Z"  # contigs joined with an N between them
        pr = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True)
        if pr.returncode not in (0, 1) or pr.stderr:   # rg exits 1 on no match
            sys.exit(f'dnacheck failed on {acc}: rc {pr.returncode} {pr.stderr[:200]}')
        res = pr.stdout.split()
        for m in set(res):
            f.write(f'{acc}\t{m}\n')
            for e in Q[m]:
                hits[e].add(acc)
        if i % 500 == 0:
            print(i, len(paths), file=sys.stderr)
with open(out + '.tsv', 'w') as f:
    f.write('event\tn_genomes\tgenomes\n')
    for e in sorted(ids):
        f.write(f"{e}\t{len(hits[e])}\t{','.join(sorted(hits[e]))}\n")
print('events', len(ids), 'with_dna_joined', sum(1 for e in ids if hits[e]), file=sys.stderr)
