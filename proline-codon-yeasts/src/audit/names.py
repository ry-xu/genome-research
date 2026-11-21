import sys, json
d = json.load(open(sys.argv[1]))['result']
ref = {l.split('\t')[0]: l.split('\t')[1] for l in open(sys.argv[2]) if not l.startswith('acc')}
bad = 0
for u in d['uids']:
    r = d[u]; acc = r['accessionversion']; t = r['title']
    nm = ref.get(acc)
    ok = nm and nm in t
    if not ok: bad += 1; print('DIFF', acc, '| fig:', nm, '| ncbi:', t)
print('checked', len(d['uids']), 'mismatches', bad)
