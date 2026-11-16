# Task 2/4 stats: skip/O per set per GC bin, MH RR vs R with phage-cluster bootstrap, per-codon table, filters.
# usage: stats.py <prefix> <meta_sub.tsv> [filter]   filter in: none, vgenus, vfam, noprophage, far200, hostgenus
import sys, collections, random, csv
pre, metap = sys.argv[1:3]; filt = sys.argv[3] if len(sys.argv) > 3 else 'none'
meta = {r['id']: r for r in csv.DictReader(open(metap), delimiter='\t')}
grp = {}
for l in open('phages_work/data/groups/run_ids.tsv'):
    if l.startswith('id\t'): continue
    f = l.rstrip('\n').split('\t'); grp[f[0]] = (f[1], int(f[2]), f[3])
rng = random.Random(7)
# phage subset
keep = set(grp)
if filt in ('vgenus', 'vfam', 'hostgenus'):
    col = {'vgenus': 'vgenus', 'vfam': 'vfam'}.get(filt)
    clus = collections.defaultdict(list)
    for i in sorted(grp):
        if filt == 'hostgenus':
            h = grp[i][2]; key = h.split(';s__')[0]  # genus-level lineage string (or higher if no genus)
        else:
            key = meta[i][col]
            if key.endswith('NA'): key = 'NA_' + i
        clus[(grp[i][0], key)].append(i)
    keep = set(rng.choice(v) for v in clus.values())
elif filt == 'noprophage':
    keep = {i for i in grp if meta[i]['checkv_trimmed'] != 'Yes' and int(meta[i]['checkv_host_markers'] or 0) == 0}
print(f'filter={filt} phages kept {len(keep)}; per set', dict(collections.Counter(grp[i][0] for i in keep)))
# per phage CGG counts
O = collections.Counter(); K = collections.Counter()
if filt == 'far200':
    for r in csv.DictReader(open(pre + '.cgg_sites.tsv'), delimiter='\t'):
        if int(r['dist_end']) < 200: continue
        O[r['id']] += 1; K[r['id']] += int(r['skip'])
else:
    for r in csv.DictReader(open(pre + '.cgg_sites.tsv'), delimiter='\t'):
        O[r['id']] += 1; K[r['id']] += int(r['skip'])
called = set(l.split('\t')[0] for l in open(pre + '.phage_codon.tsv'))
ph = collections.defaultdict(list)   # (set, bin) -> [(id, skip, O)]
for i in keep:
    if i not in called: continue
    S, B, _ = grp[i]
    ph[(S, B)].append((i, K[i], O[i]))
def tab(S, B): v = ph.get((S, B), []); return sum(x[1] for x in v), sum(x[2] for x in v), len(v)
SETS = ['T1', 'T0', 'Tf', 'Tx', 'R']
print('\nCGG skip/O by set and GC bin (skip/O [n phages])')
for B in range(22, 54, 2):
    cells = []
    for S in SETS:
        k, o, n = tab(S, B)
        cells.append(f'{S} {k}/{o}' + (f'={k/o:.3f}' if o else '') + f' [{n}]')
    if any(tab(S, B)[2] for S in SETS): print(B, ' | '.join(cells))
for bins_name, BINS in [('24-28', [24, 26, 28]), ('24-30', [24, 26, 28, 30])]:
    print(f'\nPooled bins {bins_name}:')
    for S in SETS:
        k = sum(tab(S, B)[0] for B in BINS); o = sum(tab(S, B)[1] for B in BINS)
        print(f'  {S}: {k}/{o} = {k/o if o else float("nan"):.3f}')
    def mh(data, S):
        num = den = 0.0
        for B in BINS:
            a = sum(x[1] for x in data.get((S, B), [])); n1 = sum(x[2] for x in data.get((S, B), []))
            c = sum(x[1] for x in data.get(('R', B), [])); n0 = sum(x[2] for x in data.get(('R', B), []))
            N = n1 + n0
            if n1 == 0 or n0 == 0: continue
            num += a * n0 / N; den += c * n1 / N
        return num / den if den > 0 else float('inf')
    for S in ['T1', 'T0', 'Tf']:
        est = mh(ph, S)
        bs = []
        for rep in range(1000):
            d = {}
            for key in [(S, B) for B in BINS] + [('R', B) for B in BINS]:
                v = ph.get(key, [])
                if v: d[key] = [rng.choice(v) for _ in v]
            bs.append(mh(d, S))
        bs.sort()
        print(f'  MH RR_skip {S} vs R: {est:.2f}  cluster(phage)-bootstrap 95% CI {bs[25]:.2f}-{bs[974]:.2f}')
# per-codon skip/O, bins 24-30, T1 and R
if filt == 'none':
    CO = collections.defaultdict(lambda: [0, 0])
    for l in open(pre + '.phage_codon.tsv'):
        if l.startswith('id\t'): continue
        i, S, B, c, o, k = l.rstrip('\n').split('\t')
        if int(B) not in (24, 26, 28, 30) or S not in ('T1', 'R', 'T0', 'Tf'): continue
        CO[(S, c)][0] += int(k); CO[(S, c)][1] += int(o)
    rows = []
    for (S, c), (k, o) in list(CO.items()):
        if S == 'T1': rows.append((k / o if o else 0, c, k, o, CO[('R', c)][0], CO[('R', c)][1]))
    rows.sort(reverse=True)
    print('\nPer-codon skip/O bins 24-30, T1 vs R (top 12):')
    for fr, c, k, o, kr, orr in rows[:12]:
        print(f'  {c}: T1 {k}/{o}={fr:.3f}  R {kr}/{orr}={kr/orr if orr else 0:.3f}')
    # CGG usage: CGG per 1000 codons and CGG share of Arg codons
    for S in ('T1', 'T0', 'Tf', 'R'):
        tot = sum(v[1] for (s, c), v in list(CO.items()) if s == S)
        arg = sum(CO[(S, c)][1] for c in ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'))
        print(f'  {S}: CGG per 1000 codons {1000*CO[(S,"CGG")][1]/tot:.3f}; CGG/Arg {CO[(S,"CGG")][1]/arg:.4f}; AGG/Arg {CO[(S,"AGG")][1]/arg:.4f}; CGA/Arg {CO[(S,"CGA")][1]/arg:.4f}')
