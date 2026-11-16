# descriptive T1-only variant of h6.py (T1 and R clustered alone, per-side genus) for the phage agent's addendum numbers
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
    side = {'T1': 'T', 'R': 'R'}.get(f[1])
    if side: rows.append((f, side))
S = fasta(f'{CW}/data/phages.fna', {f[0] for f, _ in rows})
for f, side in rows:
    g = S[f[0]]; L = len(g); b, e, st, k = int(f[3]), int(f[4]), int(f[5]), int(f[6])
    s = g if st == 1 else g.translate(COMP)[::-1]
    p = (b - 1 if st == 1 else L - e) + 3 * k
    assert s[p:p+3] == 'CGG'
    w = s[p-15:p+18]
    if len(w) == 33: W.append((w, side, f[0], int(f[14]), f[1]))
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
for side in 'TR':
    cs=[[i for i in ix if W[i][1]==side] for ix in cl.values()]; cs=[c for c in cs if c]
    sk=lambda c: 2*sum(W[i][3] for i in c)>=len(c)
    one=[c for c in cs if len({vg[W[i][2]] for i in c})==1]; many=[c for c in cs if len({vg[W[i][2]] for i in c})>=5]
    print(side,'CGGs',sum(map(len,cs)),'sites',len(cs),'skip sites',sum(map(sk,cs)),'| one-genus sites',len(one),'skip',sum(map(sk,one)),'CGGs skip/O',sum(W[i][3] for c in one for i in c),'/',sum(map(len,one)),'| >=5 genera sites',len(many),'CGGs',sum(map(len,many)))
