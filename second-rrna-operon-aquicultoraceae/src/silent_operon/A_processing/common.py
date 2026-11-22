import random
rc = lambda s: s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
_G = None
def genome():
    global _G
    if _G is None:
        _G, n = {}, None
        for l in open('data/all.fna'):
            if l[0] == '>': n = l[1:].split()[0]; _G[n] = []
            else: _G[n].append(l.strip().upper())
        _G = {k: ''.join(v) for k, v in _G.items()}
    return _G
def operons():
    rows = []
    f = open('out/operons.tsv'); hd = f.readline().rstrip('\n').split('\t')
    for l in f:
        r = dict(zip(hd, l.rstrip('\n').split('\t')))
        for k in ('s16a', 's16b', 'l23a', 'l23b', 'r5a', 'r5b'):
            r[k] = None if r[k] == 'None' else int(r[k])
        w = r['win'].split('|'); r['key'] = f"{w[0]}|{w[1]}"; r['st'] = w[3]; r['lo'] = int(w[7]); r['hi'] = int(w[8])
        rows.append(r)
    return rows
def region(r, p, q):
    """window 1-based positions p..q (may extend beyond window) in rRNA sense, from the full contig; '' if off contig."""
    s = genome()[r['key']]
    if r['st'] == '+':
        a, b = r['lo'] + p, r['lo'] + q          # genomic 1-based
        if a < 1 or b > len(s): return ''
        return s[a - 1:b]
    a, b = r['hi'] - q + 1, r['hi'] - p + 1
    if a < 1 or b > len(s): return ''
    return rc(s[a - 1:b])
def dishuffle(s, rng):
    s = s.replace('U', 'T')
    if len(s) < 3: return s
    last = s[-1]; E = {}
    for a, b in zip(s, s[1:]): E.setdefault(a, []).append(b)
    while True:
        lastE = {v: rng.choice(E[v]) for v in E if v != last}
        ok = True
        for v in lastE:
            seen, u = set(), v
            while u != last:
                if u in seen: ok = False; break
                seen.add(u); u = lastE[u]
            if not ok: break
        if ok: break
    L = {}
    for v, es in E.items():
        es = es[:]
        if v in lastE: es.remove(lastE[v]); rng.shuffle(es); es.append(lastE[v])
        else: rng.shuffle(es)
        L[v] = es
    out, u, idx = [s[0]], s[0], {v: 0 for v in L}
    for _ in range(len(s) - 1):
        nx = L[u][idx[u]]; idx[u] += 1; out.append(nx); u = nx
    return ''.join(out)
