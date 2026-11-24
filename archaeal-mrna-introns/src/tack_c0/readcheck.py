"""DNA check, read part (RULES.md amendment A3): can a joined read be explained by a genome that already lacks the insert?
Each joined read (primary SEQ from RUN.sam.gz) is aligned with minimap2 -x sr to every census genome of the same GTDB order as
the event's genome. A read is DNA-explained when one alignment covers the junction +-20 nt of the read, has no gap >= 5 and
identity >= 0.97 (same formula as classify.py). A template is explained when any of its joined reads is.
usage: readcheck.py OUT EVENTS_FILE DIR:RUN,RUN,... [DIR:RUN,...]   (EVENTS_FILE: one event id per line)
  -> OUT.reads.tsv (event, run, read, junction, len, explained, contig), OUT.tsv (event, joined_templates, explained_templates, contigs)
env THREADS (default 1); env IDX: index directory with events.tsv (default out)."""
import sys, os, re, csv, gzip, subprocess
from collections import defaultdict

TK = '/Users/ryan/lab/work/tack-rna'
CEN = '/Users/ryan/lab/work/archaeal-intron-census/data'
out, evfile, *dirs = sys.argv[1:]
ids = {l.strip() for l in open(evfile) if l.strip()}
FL, IDN, BIG = 20, 0.97, 5
cig_re = re.compile(r'(\d+)([MIDNSHP=X])')
EV = {r['id']: r for r in csv.DictReader(open(os.environ.get('IDX', f'{TK}/out') + '/events.tsv'), delimiter='\t') if r['id'] in ids}
order_of = {r['acc']: r['order'] for r in csv.DictReader(open(f'{CEN}/genomes_meta.tsv'), delimiter='\t')}
paths = {}
for fn in ('run1.tsv', 'run2.tsv', 'run3.tsv'):
    for l in open(f'{CEN}/{fn}'):
        a = l.rstrip('\n').split('\t')
        paths.setdefault(a[0], a[1])

def junction(cigar, rpos, rx, L, k):
    """query offset of the joined gap (length L, ref start within rx +- k+1); None if absent"""
    q, p = 0, rpos
    for n, o in ((int(n), o) for n, o in cig_re.findall(cigar)):
        if o in 'SI':
            q += n
        elif o in 'M=X':
            q += n; p += n
        elif o in 'DN':
            if n == L and rx - k - 1 <= p <= rx + k + 1:
                return q
            p += n
    return None

# joined reads per run
want = {}                                   # (run, read) -> [(event, junction)]
for dr in dirs:
    d, runs = dr.split(':', 1)
    for run in runs.split(','):
        if not os.path.exists(f'{d}/{run}.done'):
            sys.exit(f'{d}/{run} not done')
        for x in csv.DictReader(open(f'{d}/{run}.reads.tsv'), delimiter='\t'):
            if x['cls'] != 'joined' or x['event'] not in EV:
                continue
            e = EV[x['event']]
            j = junction(x['cigar'], int(x['rstart']), int(e['rx']), int(e['L']), int(e['k']))
            want.setdefault((run, x['read']), []).append((x['event'], j, e['region']))
seqs = {}                                   # (run, read, region) -> SEQ
for dr in dirs:
    d, runs = dr.split(':', 1)
    for run in runs.split(','):
        with gzip.open(f'{d}/{run}.sam.gz', 'rt') as f:
            for l in f:
                if l[0] == '@':
                    continue
                a = l.split('\t', 11)
                if int(a[1]) & 0x904 or (run, a[0]) not in want:
                    continue
                seqs[(run, a[0], a[2])] = a[9]
by_order = defaultdict(list)                # order -> fasta records
rows = []
for (run, rd), lst in want.items():
    for ev, j, reg in lst:
        s = seqs.get((run, rd, reg))
        if s is None or j is None:
            sys.exit(f'missing SEQ or junction for {run} {rd} {ev}')
        name = f'{ev}|{run}|{rd}|{j}|{len(s)}'
        by_order[order_of[EV[ev]['acc']]].append(f'>{name}\n{s}\n')
        rows.append(name)
expl = {}                                   # name -> contig
os.makedirs(f'{out}.tmp', exist_ok=True)
thr = os.environ.get('THREADS', '1')
for order, recs in by_order.items():
    accs = sorted(a for a, o in order_of.items() if o == order and a in paths)
    ref = f'{out}.tmp/ref.fa'
    with open(ref, 'w') as f:
        for a in accs:
            r = subprocess.run(['bash', '-c', f"set -o pipefail; gzcat {paths[a]} | sed 's/^>\\([^ ]*\\).*/>{a}|\\1/'"],
                               stdout=f, stderr=subprocess.PIPE, text=True)
            if r.returncode:
                sys.exit(f'ref build failed {a}: {r.stderr[:200]}')
    open(f'{out}.tmp/q.fa', 'w').write(''.join(recs))
    mm = subprocess.run(['minimap2', '-ax', 'sr', '--eqx', '-Y', '--secondary=yes', '-N', '20', '-I', '400M', '-t', thr,
                         ref, f'{out}.tmp/q.fa'], capture_output=True, text=True)
    if mm.returncode:
        sys.exit(f'minimap2 failed on {order}: {mm.stderr[-300:]}')
    print(f'{order}: {len(accs)} genomes, {len(recs)} reads', file=sys.stderr)
    for l in mm.stdout.splitlines():
        if l[0] == '@':
            continue
        a = l.split('\t', 11)
        flag = int(a[1])
        if flag & 4:
            continue
        _, _, _, j, ln = a[0].split('|')
        j, ln = int(j), int(ln)
        ops = [(int(n), o) for n, o in cig_re.findall(a[5])]
        mm_ = ins = dsmall = alnd = 0
        big = False
        for n, o in ops:
            if o in '=X':
                alnd += n; mm_ += n if o == 'X' else 0
            elif o == 'M':
                alnd += n
            elif o == 'I':
                ins += n; big |= n >= BIG
            elif o in 'DN':
                dsmall += n; big |= n >= BIG
        if big or not alnd:
            continue
        if 1 - (mm_ + ins + dsmall) / (alnd + dsmall) < IDN:
            continue
        qs = ops[0][0] if ops[0][1] in 'SH' else 0
        qe = qs + sum(n for n, o in ops if o in 'M=XI')
        if flag & 16:
            qs, qe = ln - qe, ln - qs
        if qs <= j - FL and qe >= j + FL:
            expl.setdefault(a[0], a[2])
    os.remove(ref)
with open(out + '.reads.tsv', 'w') as f:
    f.write('event\trun\tread\tjunction\tlen\texplained\tcontig\n')
    for n in rows:
        ev, run, rd, j, ln = n.split('|')
        f.write(f'{ev}\t{run}\t{rd}\t{j}\t{ln}\t{int(n in expl)}\t{expl.get(n, "")}\n')
T = defaultdict(lambda: defaultdict(bool))
C = defaultdict(set)
for n in rows:
    ev, run, rd = n.split('|')[:3]
    t = (run, rd.rsplit('.', 1)[0] if rd[-2:] in ('.1', '.2') else rd)
    T[ev][t] |= n in expl
    if n in expl:
        C[ev].add(expl[n].split('|')[0])
with open(out + '.tsv', 'w') as f:
    f.write('event\tjoined_templates\texplained_templates\texplaining_genomes\n')
    for ev in sorted(ids):
        f.write(f"{ev}\t{len(T[ev])}\t{sum(T[ev].values())}\t{','.join(sorted(C[ev]))}\n")
print(f'events {len(ids)} with joined reads {len(T)} with any explained template {sum(1 for e in T if any(T[e].values()))}', file=sys.stderr)
