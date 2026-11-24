"""Expand group ids (as named by words.py: lowest id of events with identical +-100 nt context) to all member ids.
usage: expand_groups.py IDS_FILE > ALL_IDS"""
import sys, csv
from collections import defaultdict
seqs, name = {}, None
for l in open('/Users/ryan/lab/work/tack-rna/out/regions.fa'):
    if l[0] == '>':
        name = l[1:].strip()
    else:
        seqs[name] = l.strip().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
members = defaultdict(list)
for r in csv.DictReader(open('/Users/ryan/lab/work/tack-rna/out/events.tsv'), delimiter='\t'):
    s, x, y = seqs[r['region']], int(r['rx']), int(r['ry'])
    c = s[max(0, x - 100):y + 100]
    members[(r['kind'], min(c, c.translate(comp)[::-1]))].append(r['id'])
grp = {}
for ids in members.values():
    ids.sort()
    for e in ids:
        grp[e] = ids
want = {l.strip() for l in open(sys.argv[1]) if l.strip()}
out = sorted({m for e in want for m in grp[e]})
print('\n'.join(out))
print(f'{len(want)} groups -> {len(out)} events', file=sys.stderr)
