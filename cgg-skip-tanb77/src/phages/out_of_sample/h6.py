# H6 (RULES amendment): collapse in-frame CGG sites into conserved 33-nt window clusters (Hamming <= 6, single
# linkage), then compare skip share per (cluster, side). T = T1+Tf, R = R (kill1 site table, bins 24-30); KT/KR from h4.
import sys, collections, math, numpy as np, pyrodigal
CW = 'phages_work'
COMP = str.maketrans('ACGTN', 'TGCAN'); STOP = ('TAA', 'TAG', 'TGA')
vg = {}
for l in open('meta_cols.tsv'):
    f = l.rstrip('\n').split('\t'); vg[f[0]] = next((x for x in f[5].split(';') if x.startswith('vGENUS')), 'NA:' + f[0])
def fasta(p, need):
    n, b, out = None, [], {}
    for line in open(p):
        if line[0] == '>':
            if n in need: out[n] = ''.join(b).upper()
            n = line[1:].split()[0]; b = []
        elif n in need: b.append(line.strip())
    if n in need: out[n] = ''.join(b).upper()
    return out
W = []   # (window, side, id, skip)
rows = []
for l in open(f'{CW}/kill1/out/meta.cgg_sites.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'id' or not (24 <= int(f[2]) <= 30): continue
    side = {'T1': 'T', 'Tf': 'T', 'R': 'R'}.get(f[1])
    if side: rows.append((f, side))
S = fasta(f'{CW}/data/phages.fna', {f[0] for f, _ in rows})
for f, side in rows:
    g = S[f[0]]; L = len(g); b, e, st, k = int(f[3]), int(f[4]), int(f[5]), int(f[6])
    s = g if st == 1 else g.translate(COMP)[::-1]
    p = (b - 1 if st == 1 else L - e) + 3 * k
    assert s[p:p+3] == 'CGG'
    w = s[p-15:p+18]
    if len(w) == 33: W.append((w, side, f[0], int(f[14]), f[1]))
# KT / KR, same rule as src/skip.py
lab = {l.split('\t')[0]: (l.split('\t')[1], int(l.split('\t')[2])) for l in open('h4.ids')}
H = fasta('h4.fna', set(lab)); gf = pyrodigal.GeneFinder(meta=True)
for n, s0 in H.items():
    if len(s0) < 5000 or not (24 <= lab[n][1] <= 30): continue
    r = s0.translate(COMP)[::-1]; L = len(s0)
    for g in gf.find_genes(s0):
        if g.partial_begin or g.partial_end or g.translation_table != 11 or g.end - g.begin + 1 < 300: continue
        t = s0 if g.strand == 1 else r; a = g.begin - 1 if g.strand == 1 else L - g.end
        cds = t[a:a + g.end - g.begin + 1]
        if cds[-3:] not in STOP: continue
        nc = len(cds) // 3 - 1
        for k in range(nc):
            if cds[3*k:3*k+3] != 'CGG': continue
            p = a + 3*k; sk = int(k >= nc - 16 and all(len(t[p+1+3*j:p+4+3*j]) == 3 and t[p+1+3*j:p+4+3*j] not in STOP for j in range(50)))
            w = t[p-15:p+18]
            if len(w) == 33: W.append((w, lab[n][0], n, sk, lab[n][0]))
print('windows', collections.Counter(x[1] for x in W), 'skip', collections.Counter(x[1] for x in W if x[3]))
A = np.frombuffer(''.join(x[0] for x in W).encode(), dtype=np.uint8).reshape(len(W), 33)
par = list(range(len(W)))
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for i in range(0, len(W), 500):
    d = (A[i:i+500, None, :] != A[None, :, :]).sum(2)
    for a, b in zip(*np.nonzero(d <= 6)):
        if a + i < b: par[find(a + i)] = find(b)
cl = collections.defaultdict(list)
for i in range(len(W)): cl[find(i)].append(i)
print('clusters', len(cl))
def fisher(a, b, c, d):   # two-sided, sum of tables with prob <= observed
    n1, n2, m = a + b, c + d, a + c; N = n1 + n2
    pr = lambda x: math.comb(n1, x) * math.comb(n2, m - x) / math.comb(N, m)
    po = pr(a); return sum(pr(x) for x in range(max(0, m - n2), min(n1, m) + 1) if pr(x) <= po * (1 + 1e-9))
def units(single, tsets=('T1', 'Tf')):
    U = collections.Counter()
    for ix in cl.values():
        genera = {vg[W[i][2]] for i in ix}
        if single and len(genera) > 1: continue
        for side in 'TR':
            m = [W[i][3] for i in ix if W[i][1] == side and (side == 'R' or W[i][4] in tsets)]
            if m: U[side, 'skip' if 2 * sum(m) >= len(m) else 'no'] += 1
    return U
for single in (False, True):
    U = units(single); ts, tn, rs, rn = U['T', 'skip'], U['T', 'no'], U['R', 'skip'], U['R', 'no']
    tsh, rsh = ts / (ts + tn), rs / (rs + rn)
    print(f'H6a{" singletons(one vGENUS)" if single else ""}: T skip units {ts}/{ts+tn} = {tsh:.3f}; R {rs}/{rs+rn} = {rsh:.4f}; ratio {tsh/rsh if rsh else float("inf"):.1f}; Fisher p {fisher(ts, tn, rs, rn):.2e}')
U = units(False, ('T1',)); print(f"descr T1-only units: skip {U['T','skip']}/{U['T','skip']+U['T','no']}; R {U['R','skip']}/{U['R','skip']+U['R','no']}")
U = units(True, ('T1',)); print(f"descr T1-only one-vGENUS clusters: skip {U['T','skip']}/{U['T','skip']+U['T','no']}; R {U['R','skip']}/{U['R','skip']+U['R','no']}")
print('T1 CGGs', sum(1 for x in W if x[4]=='T1'), 'skip', sum(x[3] for x in W if x[4]=='T1'), '| R CGGs', sum(1 for x in W if x[1]=='R'), 'skip', sum(x[3] for x in W if x[1]=='R'))
# H6b
mix = [ix for ix in cl.values() if {W[i][1] for i in ix} >= {'T', 'R'}]
ts = sum(W[i][3] for ix in mix for i in ix if W[i][1] == 'T'); to = sum(1 for ix in mix for i in ix if W[i][1] == 'T')
rs = sum(W[i][3] for ix in mix for i in ix if W[i][1] == 'R'); ro = sum(1 for ix in mix for i in ix if W[i][1] == 'R')
tm = rm = 0
for ix in mix:
    t = [W[i][3] for i in ix if W[i][1] == 'T']; r = [W[i][3] for i in ix if W[i][1] == 'R']
    T_, R_ = 2 * sum(t) >= len(t), 2 * sum(r) >= len(r)
    tm += T_ and not R_; rm += R_ and not T_
print(f'H6b mixed clusters {len(mix)}: T skip {ts}/{to}, R skip {rs}/{ro}; T-skip-only {tm}, R-skip-only {rm}')
for ix in sorted(mix, key=len, reverse=True)[:10]:
    t = [W[i][3] for i in ix if W[i][1] == 'T']; r = [W[i][3] for i in ix if W[i][1] == 'R']
    print('  cluster', len(ix), 'T', sum(t), '/', len(t), 'R', sum(r), '/', len(r), W[ix[0]][0][:15].lower() + W[ix[0]][0][15:18] + W[ix[0]][0][18:].lower())
# H6c
Tcl = {find(i) for i in range(len(W)) if W[i][1] == 'T'}
kt = [i for i in range(len(W)) if W[i][1] == 'KT' and W[i][3]]
print(f'H6c KT skip sites {len(kt)}, in a cluster with a T1/Tf CGG {sum(find(i) in Tcl for i in kt)}')
kn = [i for i in range(len(W)) if W[i][1] == 'KT']; kr = [i for i in range(len(W)) if W[i][1] == 'KR']
print(f'    KT CGGs outside T clusters: skip {sum(W[i][3] for i in kn if find(i) not in Tcl)}/{sum(1 for i in kn if find(i) not in Tcl)}; KR outside: {sum(W[i][3] for i in kr if find(i) not in Tcl)}/{sum(1 for i in kr if find(i) not in Tcl)}')
