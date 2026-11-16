# Recompute skip sites from a gene table. Rule: CDS not partial, >=300 nt (incl stop), ttable 11 (others counted);
# in-frame sense codons; site = codon in the last 16 sense codons whose +1 frame read from the codon's 2nd base
# has no TAA/TAG/TGA in its first 50 codons (all 50 inside the contig).
# Outputs: <pre>.cgg_sites.tsv (every in-frame CGG), <pre>.ctrl_sites.tsv (non-CGG last-16 codons with open +1),
#          <pre>.phage_codon.tsv (per phage per codon: O, skip)
import sys, collections
fa, genes_path, groups_path, pre = sys.argv[1:5]
grp = {}
for l in open(groups_path):
    if l.startswith('id\t'): continue
    f = l.rstrip('\n').split('\t'); grp[f[0]] = (f[1], f[2])
genes = collections.defaultdict(list)
for l in open(genes_path):
    if l.startswith('id\t'): continue
    f = l.rstrip('\n').split('\t')
    genes[f[0]].append((int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6])))
def read_fa(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith('>'):
            if name: yield name, ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: yield name, ''.join(buf)
COMP = str.maketrans('ACGTN', 'TGCAN')
STOPS = {'TAA', 'TAG', 'TGA'}
def plus1(s, p, cap=1000):
    # count codons from p until stop; return (n_open, reached_end)
    n = 0; L = len(s)
    while n < cap:
        if p + 3 > L: return n, True
        if s[p:p+3] in STOPS: return n, False
        n += 1; p += 3
    return n, False
fc = open(pre + '.cgg_sites.tsv', 'w'); fk = open(pre + '.ctrl_sites.tsv', 'w'); fp = open(pre + '.phage_codon.tsv', 'w')
hdr = 'id\tset\tgc_bin\tbegin\tend\tstrand\tk\tn\tcodon\topen\treached_end\tdist_end\tnext3\tlast16\tskip\n'
fc.write(hdr); fk.write(hdr); fp.write('id\tset\tgc_bin\tcodon\tO\tskip\n')
stat = collections.Counter()
for name, seq in read_fa(fa):
    if name not in genes or name not in grp: continue
    S, B = grp[name]
    seq = seq.upper(); L = len(seq); rc = seq.translate(COMP)[::-1]
    O = collections.Counter(); K = collections.Counter()
    for (b, e, st, pl, pr, tt) in genes[name]:
        stat['genes'] += 1
        if pl or pr: stat['partial'] += 1; continue
        if e - b + 1 < 300: stat['short'] += 1; continue
        if tt != 11: stat['ttable_' + str(tt)] += 1; continue
        s = seq if st == 1 else rc
        s0 = b - 1 if st == 1 else L - e
        g = s[s0:s0 + (e - b + 1)]
        if g[-3:] not in STOPS: stat['no_stop_end'] += 1; continue
        n = len(g) // 3 - 1
        stat['used'] += 1
        for k in range(n):
            c = g[3*k:3*k+3]
            O[c] += 1
            last16 = k >= n - 16
            is_cgg = c == 'CGG'
            if not (last16 or is_cgg): continue
            p = s0 + 3*k + 1
            op, re_ = plus1(s, p)
            skip = int(last16 and op >= 50)
            if skip: K[c] += 1
            if is_cgg or skip:
                gpos = (s0 + 3*k) if st == 1 else (L - (s0 + 3*k) - 3)
                dist = min(gpos, L - gpos - 3)
                row = f'{name}\t{S}\t{B}\t{b}\t{e}\t{st}\t{k}\t{n}\t{c}\t{op}\t{int(re_)}\t{dist}\t{s[s0+3*k+3:s0+3*k+6]}\t{int(last16)}\t{skip}\n'
                (fc if is_cgg else fk).write(row)
    for c in O:
        fp.write(f'{name}\t{S}\t{B}\t{c}\t{O[c]}\t{K[c]}\n')
print(dict(stat))
