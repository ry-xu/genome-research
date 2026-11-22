"""Codon at each P. jirovecii TGA site (Arg-consensus columns, from the RNA tally table) in other P. jirovecii
mitogenomes, found by an exact 18-nt flank match (left flank, else right flank; either strand).
usage: python isolates.py NC_020331.gb sites.tsv others.fa out.tsv"""
import sys
from collections import Counter
from Bio import SeqIO

gb, sites, others, out = sys.argv[1:]
rec = SeqIO.read(gb, 'genbank')
G = str(rec.seq).upper()
cds = {}
for f in rec.features:
    if f.type == 'CDS':
        cds.setdefault(f.qualifiers.get('gene', ['?'])[0], f)
rows = [l.rstrip('\n').split('\t') for l in open(sites)][1:]
want = []
for r in rows:
    if r[2] != 'TGA' or r[4] != 'R':
        continue
    f = cds[r[0]]; i = int(r[1])
    s = str(f.extract(rec.seq)).upper()
    assert s[3 * i:3 * i + 3] == 'TGA', (r[0], i)
    ctx = s[3 * i - 18:3 * i + 21]   # 18 + codon + 18 in CDS orientation
    want.append((r[0], i, r[6], ctx))
rc = lambda x: x[::-1].translate(str.maketrans('ACGT', 'TGCA'))
tot, per = Counter(), []
for o in SeqIO.parse(others, 'fasta'):
    q = str(o.seq).upper(); Q = (q, rc(q))
    for g, i, strict, ctx in want:
        L, R = ctx[:18], ctx[21:]
        hit = '-'
        for t in Q:
            if t.count(L) == 1:
                k = t.index(L) + 18; hit = t[k:k + 3]; break
            if t.count(R) == 1:
                k = t.index(R) - 3; hit = t[k:k + 3]; break
        tot[hit] += 1; per.append((o.id, g, i, strict, hit))
with open(out, 'w') as fh:
    fh.write('acc\tgene\tcodon_index\tstrict\tcodon\n')
    for p in per:
        fh.write('\t'.join(map(str, p)) + '\n')
print('sites', len(want), 'strict', sum(w[2] == '1' for w in want), 'genomes', len(per) // len(want))
print('all sites', sorted(tot.items(), key=lambda t: -t[1]))
st = Counter(p[4] for p in per if p[3] == '1')
print('strict sites', sorted(st.items(), key=lambda t: -t[1]))
