# The referee's own skip-signature counter (RULES.md H4). usage: skip.py FASTA IDS.tsv(id<TAB>label<TAB>gc_bin) OUT.tsv
# pyrodigal meta; table-11, complete CDS >= 300 nt ending in a stop; per in-frame CGG: is it in the last 16 sense
# codons and is the +1 frame from its 2nd base free of stops for 50 codons (inside the contig)?
import sys, pyrodigal
fa, idp, out = sys.argv[1:4]
lab = {}
for l in open(idp):
    f = l.rstrip('\n').split('\t'); lab[f[0]] = (f[1], f[2])
STOP = ('TAA', 'TAG', 'TGA')
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def seqs(p):
    n, b = None, []
    for line in open(p):
        if line[0] == '>':
            if n in lab: yield n, ''.join(b)
            n = line[1:].split()[0]; b = []
        elif n in lab: b.append(line.strip())
    if n in lab: yield n, ''.join(b)
gf = pyrodigal.GeneFinder(meta=True)
with open(out, 'w') as fo:
    fo.write('id\tlabel\tgc_bin\tcgg\tskip\n')
    for n, s in seqs(fa):
        s = s.upper()
        if len(s) < 5000: continue
        r = rc(s); L = len(s); cg = sk = 0
        for g in gf.find_genes(s):
            if g.partial_begin or g.partial_end or g.translation_table != 11 or g.end - g.begin + 1 < 300: continue
            t = s if g.strand == 1 else r
            a = g.begin - 1 if g.strand == 1 else L - g.end
            cds = t[a:a + g.end - g.begin + 1]
            if cds[-3:] not in STOP: continue
            nc = len(cds) // 3 - 1
            for k in range(nc):
                if cds[3*k:3*k+3] != 'CGG': continue
                cg += 1
                if k < nc - 16: continue
                p = a + 3*k + 1; ok = True
                for j in range(50):
                    c = t[p + 3*j:p + 3*j + 3]
                    if len(c) < 3 or c in STOP: ok = False; break
                sk += ok
        fo.write(f'{n}\t{lab[n][0]}\t{lab[n][1]}\t{cg}\t{sk}\n')
