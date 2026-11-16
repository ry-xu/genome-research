# H4 readings (RULES.md): per-phage skip table -> pooled skip/O by label, MH over GC bins KT vs KR.
import sys, collections
tab, novel = sys.argv[1], (set(l.split()[0] for l in open(sys.argv[2])) if len(sys.argv) > 2 else None)
D = collections.defaultdict(lambda: [0, 0, 0])   # (label, bin) -> skip, O, phages
for l in open(tab):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'id': continue
    if novel is not None and f[1] == 'KT' and f[0] not in novel: continue
    d = D[f[1], int(f[2])]; d[0] += int(f[4]); d[1] += int(f[3]); d[2] += 1
def pool(lab, lo, hi):
    s = sum(D[k][0] for k in D if k[0] == lab and lo <= k[1] <= hi); o = sum(D[k][1] for k in D if k[0] == lab and lo <= k[1] <= hi)
    n = sum(D[k][2] for k in D if k[0] == lab and lo <= k[1] <= hi); return s, o, n
def mh(lo, hi):
    num = den = 0
    for b in range(lo, hi + 1, 2):
        a, o1, _ = D[('KT', b)]; c, o2, _ = D[('KR', b)]
        if o1 and o2: N = o1 + o2; num += a * o2 / N; den += c * o1 / N
    return num / den if den else float('inf')
for name, lo, hi in (('low 24-30', 24, 30), ('all', 0, 100)):
    kt = pool('KT', lo, hi); kr = pool('KR', lo, hi)
    print(f'{name}: KT skip/O {kt[0]}/{kt[1]} = {kt[0]/max(kt[1],1):.3f} ({kt[2]} phages); KR {kr[0]}/{kr[1]} = {kr[0]/max(kr[1],1):.3f} ({kr[2]} phages); MH {mh(lo, hi):.2f}')
for k in sorted(D): print(*k, *D[k], sep='\t')
