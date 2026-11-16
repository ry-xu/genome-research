# For the best spanning HSP per query: target residues aligned to query positions J-1, J, J+1, J+2 (J+1 = first +1-frame aa, from GGN)
import sys, csv, collections
qtsv, hits_p = sys.argv[1:3]
Q = {r['qid']: r for r in csv.DictReader(open(qtsv), delimiter='\t')}
best = {}
for l in open(hits_p):
    f = l.rstrip('\n').split('\t'); qid = f[0]
    if f[1].split('|')[0] == Q[qid]['id']: continue
    J = int(Q[qid]['junction']); qs, qe = int(f[4]), int(f[5])
    if not (qs <= J - 19 and qe >= J + 20): continue
    if qid not in best or float(f[9]) > float(best[qid][9]): best[qid] = f
cnt = collections.defaultdict(collections.Counter)
for qid, f in sorted(best.items()):
    J = int(Q[qid]['junction']); qp = int(f[4]) - 1; res = {}
    for a, b in zip(f[12], f[13]):
        if a != '-': qp += 1
        if a != '-' and J - 1 <= qp <= J + 2: res[qp - J] = (a, b)
    s = ' '.join(f'{k:+d}:{res[k][0]}/{res[k][1]}' if k in res else f'{k:+d}:?' for k in (-1, 0, 1, 2))
    print(qid, s)
    if 1 in res: cnt[Q[qid]['class']][res[1][1]] += 1
for c in cnt: print(c, 'target residue at J+1 (query Gly from GGN):', dict(cnt[c].most_common()))
