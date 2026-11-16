# H1-H3, H5 of RULES.md: kill1 per-CGG site table joined to UHGV k-mer host and vGENUS.
import sys, collections
sites, meta = sys.argv[1:3]
M = {}
for l in open(meta):
    f = l.rstrip('\n').split('\t')
    vg = next((x for x in f[5].split(';') if x.startswith('vGENUS')), 'NA:' + f[0])
    M[f[0]] = (f[7], vg, f[6])
def kcls(k):
    if 'o__TANB77' in k: return 'kTANB77'
    if k in ('NULL', '', 'NA'): return 'kNone'
    return 'kOther'
S = collections.defaultdict(lambda: [0, 0])   # (set, kcls, lowgc) -> [skip, O]
phg = collections.defaultdict(set)
VG = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))  # vg -> side -> [skip,O]
for l in open(sites):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'id': continue
    pid, st, b, sk = f[0], f[1], int(f[2]), int(f[14])
    k = kcls(M[pid][0]); low = 24 <= b <= 30
    for key in ((st, k, 'low' if low else 'high'), (st, k, 'all')):
        S[key][0] += sk; S[key][1] += 1; phg[key].add(pid)
    side = 'T' if st in ('T1', 'T0', 'Tf') else ('R' if st == 'R' else None)
    if side: VG[M[pid][1]][side][0] += sk; VG[M[pid][1]][side][1] += 1
print('set\tkmer_host\tgc\tphages_with_CGG\tskip\tO\tshare')
for key in sorted(S):
    s, o = S[key]; print(*key, len(phg[key]), s, o, f'{s/o:.3f}', sep='\t')
# H3: k-mer host share among T1 phages by GC (per phage, from the group table, not only phages with CGG)
grp = {}
for l in open('phages_work/data/groups/run_ids.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'id': grp[f[0]] = (f[1], int(f[2]))
c = collections.Counter()
for pid, (st, b) in grp.items():
    if st == 'T1': c[('low' if b <= 30 else ('high32+' if b >= 32 else 'x')), kcls(M[pid][0])] += 1
print('H3 T1 phages by GC and k-mer host:', dict(sorted(c.items())))
# H5: vGENUS with both sides having CGG
num = den = 0; ts = to = rs = ro = 0; n = 0
for vg, d in VG.items():
    if d['T'][1] and d['R'][1]:
        a, o1 = d['T']; cc, o2 = d['R']; N = o1 + o2; n += 1
        ts += a; to += o1; rs += cc; ro += o2
        num += a * (o2) / N; den += cc * (o1) / N
print(f'H5 vGENUS with both: {n}; T skip/O {ts}/{to} = {ts/max(to,1):.3f}; R {rs}/{ro} = {rs/max(ro,1):.3f}; MH ratio {num/den if den else float("inf"):.2f}')
