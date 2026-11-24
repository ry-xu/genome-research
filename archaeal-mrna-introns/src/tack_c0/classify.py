"""Classify primary alignments against events (RULES.md §2). SAM on stdin.
usage: classify.py PREFIX  -> PREFIX.counts.tsv (event, joined, through, other templates), PREFIX.reads.tsv (every spanning read)
v2 (20 Nov, after the independent review, review/REVIEW.md): identity denominator = aligned + small deletions (§2);
the pre-filter admits a joined read with a 20-nt flank at the far end of the [x-k-1, x+k+1] window; duplicates are
collapsed on the whole template (region, lowest start over both mates, highest end, mate-1 strand), not on one mate.
env IDX: index directory holding events.tsv (default out; out/rep for the A5 replication)."""
import sys, re, csv, os
from collections import defaultdict

pre = sys.argv[1]
EV = defaultdict(list)
for r in csv.DictReader(open(os.environ.get('IDX', '/Users/ryan/lab/work/tack-rna/out') + '/events.tsv'), delimiter='\t'):
    EV[r['region']].append((r['id'], int(r['rx']), int(r['ry']), int(r['L']), int(r['k'])))
FL, IDN, BIG = 20, 0.97, 5
RANK = {'joined': 3, 'through': 2, 'other': 1}
cig_re = re.compile(r'(\d+)([MIDNSHP=X])')

best = {}                                       # (event, template) -> class
ext = {}                                        # (region, template) -> [min start, max end, mate-1 strand]
out_reads = open(pre + '.reads.tsv', 'w')
out_reads.write('event\tcls\tread\trstart\trend\tflag\tident\tcigar\n')
nprim = 0
for line in sys.stdin:
    if line[0] == '@':
        continue
    f = line.split('\t', 11)
    flag = int(f[1])
    if flag & 0x904 or f[2] not in EV:
        continue
    nprim += 1
    pos = int(f[3]) - 1
    ops = [(int(n), o) for n, o in cig_re.findall(f[5])]
    mm = ins = dsmall = alnd = 0
    blocks, gaps = [], []          # aligned ref blocks; big gaps (start, len)
    p = pos; bstart = pos
    for n, o in ops:
        if o in '=X':
            alnd += n
            if o == 'X':
                mm += n
            p += n
        elif o == 'M':
            alnd += n; p += n
        elif o == 'I':
            ins += n
        elif o in 'DN':
            if n >= BIG:
                blocks.append((bstart, p)); gaps.append((p, n)); p += n; bstart = p
            else:
                dsmall += n; p += n
    blocks.append((bstart, p))
    end = p
    den = alnd + dsmall
    ident = 1 - (mm + ins + dsmall) / den if den else 0
    if ident < IDN:
        continue
    tmpl = f[0].rsplit('.', 1)[0] if f[0][-2:] in ('.1', '.2') else f[0]
    m1 = not f[0].endswith('.2')
    e = ext.setdefault((f[2], tmpl), [pos, end, None])
    e[0] = min(e[0], pos); e[1] = max(e[1], end)
    if m1:
        e[2] = (flag & 16) > 0
    elif e[2] is None:
        e[2] = not (flag & 16)     # mate 2 only: mate-1 strand is the opposite

    def aligned_in(a, b):
        return sum(max(0, min(b, y) - max(a, x)) for x, y in blocks)
    for eid, rx, ry, L, k in EV[f[2]]:
        if pos > rx + k + 1 - FL or end < ry - k - 1 + FL:
            continue
        g_in = [(s, n) for s, n in gaps if s < ry + FL and s + n > rx - FL]
        cls = None
        if len(g_in) == 1:
            s, n = g_in[0]
            if n == L and rx - k - 1 <= s <= rx + k + 1 and aligned_in(s - FL, s) >= FL and aligned_in(s + n, s + n + FL) >= FL:
                cls = 'joined'
        if cls is None and not g_in and pos <= rx - FL and end >= ry + FL and aligned_in(rx - FL, ry + FL) == ry - rx + 2 * FL:
            cls = 'through'
        if cls is None:
            if pos <= rx - FL and end >= ry + FL:
                cls = 'other'
            else:
                continue
        key = (eid, f[2], tmpl)
        if key not in best or RANK[cls] > RANK[best[key]]:
            best[key] = cls
        out_reads.write(f'{eid}\t{cls}\t{f[0]}\t{pos}\t{end}\t{flag}\t{ident:.4f}\t{f[5]}\n')
out_reads.close()
cnt = defaultdict(lambda: defaultdict(set))
for (eid, reg, tmpl), cls in best.items():
    s, e, st = ext[(reg, tmpl)]
    cnt[eid][cls].add((s, e, st))
with open(pre + '.counts.tsv', 'w') as o:
    o.write('event\tjoined\tthrough\tother\n')
    for eid in sorted(cnt):
        c = cnt[eid]
        o.write(f"{eid}\t{len(c['joined'])}\t{len(c['through'])}\t{len(c['other'])}\n")
print(f'primary_on_index {nprim} events_spanned {len(cnt)}', file=sys.stderr)
