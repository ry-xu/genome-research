# Task 3 build: (1) target proteins = complete ORFs (pyrodigal meta, not partial, >=60 aa) from phages of chosen sets;
# (2) seeded sample of fused queries: upstream codons 0..k-1 + (+1 frame from the codon's 2nd base to its stop).
import sys, csv, random, collections
fa, genes_p, cgg_p, ctrl_p, outdir = sys.argv[1:6]
grp = {}
for l in open('phages_work/data/groups/run_ids.tsv'):
    if l.startswith('id\t'): continue
    f = l.rstrip('\n').split('\t'); grp[f[0]] = (f[1], int(f[2]))
CODE = {}
bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bases:
    for b in bases:
        for c in bases: CODE[a+b+c] = aas[i]; i += 1
def tr(s): return ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3))
COMP = str.maketrans('ACGTN', 'TGCAN')
seqs = {}; name = None; buf = []
for line in open(fa):
    if line.startswith('>'):
        if name: seqs[name] = ''.join(buf).upper()
        name = line[1:].split()[0]; buf = []
    else: buf.append(line.strip())
seqs[name] = ''.join(buf).upper()
rcs = {}
def strand_seq(i, st):
    if st == 1: return seqs[i]
    if i not in rcs: rcs[i] = seqs[i].translate(COMP)[::-1]
    return rcs[i]
# targets
nt = 0
with open(outdir + '/targets_R.faa', 'w') as fR, open(outdir + '/targets_all.faa', 'w') as fA:
    for l in open(genes_p):
        if l.startswith('id\t'): continue
        i, b, e, st, pl, pr, tt = l.rstrip('\n').split('\t')
        if pl == '1' or pr == '1': continue
        b, e, st = int(b), int(e), int(st)
        s = strand_seq(i, st); L = len(s); s0 = b - 1 if st == 1 else L - e
        p = tr(s[s0:s0 + e - b + 1]).rstrip('*')
        if len(p) < 60 or '*' in p: continue
        rec = f'>{i}|{b}|{e}|{st}\n{p}\n'
        fA.write(rec)
        if grp.get(i, ('',))[0] == 'R': fR.write(rec); nt += 1
print('R target proteins', nt)
rng = random.Random(11)
def load(p): return list(csv.DictReader(open(p), delimiter='\t'))
cgg = load(cgg_p); ctrl = load(ctrl_p)
LOWGC = {24, 26, 28, 30}
t1 = [r for r in cgg if r['set'] == 'T1' and r['skip'] == '1']
rsk = [r for r in cgg if r['set'] == 'R' and r['skip'] == '1' and int(r['gc_bin']) in LOWGC]
genes_with_cggskip = {(r['id'], r['begin']) for r in cgg if r['skip'] == '1'}
byg = collections.defaultdict(list)
for r in ctrl:
    if r['set'] == 'T1' and r['skip'] == '1' and (r['id'], r['begin']) not in genes_with_cggskip:
        byg[(r['id'], r['begin'])].append(r)
t1c = [rng.choice(v) for k, v in sorted(byg.items())]
print('pool sizes: T1 CGG skip', len(t1), '| T1 non-CGG open gene ends (genes)', len(t1c), '| R CGG skip (GC 24-30)', len(rsk))
samp = {'T1skip': rng.sample(t1, 40), 'T1ctrl': rng.sample(t1c, 40), 'Rskip': rng.sample(rsk, min(40, len(rsk)))}
with open(outdir + '/queries.faa', 'w') as fq, open(outdir + '/queries.tsv', 'w') as ft:
    ft.write('qid\tclass\tid\tset\tgc_bin\tcodon\tk\tn\tjunction\tqlen\tnext3\n')
    for cls, rows in samp.items():
        for j, r in enumerate(rows):
            i, b, e, st, k = r['id'], int(r['begin']), int(r['end']), int(r['strand']), int(r['k'])
            s = strand_seq(i, st); L = len(s); s0 = b - 1 if st == 1 else L - e
            up = tr(s[s0:s0 + 3*k])
            p = s0 + 3*k + 1; down = []
            while p + 3 <= L:
                a = CODE.get(s[p:p+3], 'X')
                if a == '*': break
                down.append(a); p += 3
            q = up + ''.join(down)
            qid = f'{cls}_{j:02d}'
            fq.write(f'>{qid}\n{q}\n')
            ft.write(f'{qid}\t{cls}\t{i}\t{r["set"]}\t{r["gc_bin"]}\t{r["codon"]}\t{k}\t{r["n"]}\t{len(up)}\t{len(q)}\t{r["next3"]}\n')
print('queries written', {k: len(v) for k, v in samp.items()})
