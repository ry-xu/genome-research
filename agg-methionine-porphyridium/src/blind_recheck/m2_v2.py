"""Method 2: reference-protein alignment to six-frame (stops->X) translation; per-codon likelihood of meaning."""
import sys, os, subprocess, collections, math, json
genome_fa, refs_fa, outdir, name = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)
B = 'TCAG'
codons = [a+b+c for a in B for b in B for c in B]
aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
STD = dict(zip(codons, aas))
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
comp = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')
seq = ''.join(l.strip() for l in open(genome_fa) if not l.startswith('>')).upper()
rc = seq.translate(comp)[::-1]
frames = {}
for strand, s in (('+', seq), ('-', rc)):
    for f in range(3):
        cl = [s[i:i+3] for i in range(f, len(s) - 2, 3)]
        prot = ''.join(STD[c] if c in STD and STD[c] != '*' else 'X' for c in cl)
        frames['%s%d' % (strand, f)] = (cl, prot)
ffa = os.path.join(outdir, name + '.frames.faa')
with open(ffa, 'w') as fo:
    for k, (cl, p) in frames.items():
        fo.write('>%s\n%s\n' % (k, p))
subprocess.run(['makeblastdb', '-in', ffa, '-dbtype', 'prot', '-out', ffa], check=True, capture_output=True)
bo = os.path.join(outdir, name + '.blastp.tsv')
if not os.path.exists(bo):
    subprocess.run(['blastp', '-query', refs_fa, '-db', ffa, '-evalue', '1e-5', '-num_threads', '1',
                    '-max_target_seqs', '6', '-max_hsps', '20',
                    '-outfmt', '6 qseqid sseqid evalue bitscore qstart qend sstart send qseq sseq', '-out', bo], check=True)
TRIM = 5
best = {}  # (ref, site) -> (bits, aa)
genes_hit = collections.Counter()
for line in open(bo):
    q, s, ev, bits, qs, qe, ss, se, qa, sa = line.rstrip('\n').split('\t')
    bits = float(bits); spos = int(ss) - 1
    pairs = []
    for x, y in zip(qa, sa):
        if y != '-':
            if x != '-':
                pairs.append((spos, x))
            spos += 1
    genes_hit[q.split('|')[1]] += 1
    for sp, x in pairs[TRIM:len(pairs) - TRIM]:
        key = (q, (s, sp))
        if key not in best or best[key][0] < bits:
            best[key] = (bits, x)
site_refs = collections.defaultdict(dict)
for (q, site), (bits, x) in best.items():
    site_refs[site][q] = x
# per site weight vector
codon_sites = collections.defaultdict(list)  # codon -> list of Counter(weights)
for (fr, sp), refs in site_refs.items():
    c = frames[fr][0][sp]
    if c not in STD: continue
    vals = [x for x in refs.values() if x in AA20]
    if not vals: continue
    w = collections.Counter()
    for x in vals: w[x] += 1.0 / len(vals)
    codon_sites[c].append(w)
MINPROF = 30
def profiles(assign, exclude):
    prof = {}
    for a in AA20:
        tot = collections.Counter(); n = 0
        rest = sum(len(codon_sites[c]) for c in codons if c != exclude and assign.get(c) == a)
        for c in codons:
            if assign.get(c) != a: continue
            if c == exclude and rest >= MINPROF: continue  # leave-one-out only if >=30 sites remain (v2)
            for w in codon_sites[c]:
                tot.update(w); n += 1
        prof[a] = {b: (tot[b] + 1.0) / (n + 20.0) for b in AA20}
    return prof
def score(assign):
    res = {}
    for c in codons:
        sites = codon_sites[c]
        n = len(sites)
        prof = profiles(assign, c)
        ll = {a: sum(sum(wb * math.log(prof[a][b]) for b, wb in w.items()) for w in sites) for a in AA20}
        order = sorted(AA20, key=lambda a: -ll[a])
        X, Y = order[0], order[1]
        std = STD[c]
        llr_std = ll[X] - ll[std] if std != '*' else None
        llr_2nd = ll[X] - ll[Y]
        tot = collections.Counter()
        for w in sites: tot.update(w)
        top = ', '.join('%s %.2f' % (b, v / n) for b, v in tot.most_common(3)) if n else ''
        if std != '*':
            strong = X != std and n >= 10 and llr_std >= 10 and llr_2nd >= 3
            support = X != std and n >= 5 and llr_std >= 3
        else:
            strong = n >= 10 and llr_2nd >= 3
            support = n >= 5 and llr_2nd >= 1
        res[c] = dict(codon=c, std=std, n=n, best=X, second=Y, llr_std=llr_std, llr_2nd=llr_2nd,
                      ll_std=(ll[std] if std != '*' else None), top_ref=top, strong=strong, support=support,
                      ll={a: ll[a] for a in AA20})
    return res
assign1 = {c: STD[c] for c in codons}
r1 = score(assign1)
assign2 = dict(assign1)
moved = {c: r1[c]['best'] for c in codons if r1[c]['strong']}
assign2.update(moved)
r2 = score(assign2)
# absent-codon flag
for c in codons:
    a = STD[c]
    if a == '*': r2[c]['absent_flag'] = False; continue
    others = sum(len(codon_sites[d]) for d in codons if d != c and STD[d] == a)
    r2[c]['absent_flag'] = len(codon_sites[c]) <= 1 and others >= 100
    r2[c]['others_n'] = others
out = dict(name=name, genome_len=len(seq), n_sites=sum(len(v) for v in codon_sites.values()),
           genes_hit=dict(genes_hit), moved_pass1=moved, pass1={c: {k: r1[c][k] for k in ('n','best','llr_std','llr_2nd','strong')} for c in codons}, pass2=r2)
json.dump(out, open(os.path.join(outdir, name + '.m2.json'), 'w'), indent=1)
print(name, 'len', len(seq), 'sites', out['n_sites'], 'genes', len(genes_hit))
