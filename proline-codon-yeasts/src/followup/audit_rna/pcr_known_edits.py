"""(b) PC-R hand check: RefSeq NC_037304.1 /translation is the edited protein. Sites where a transcript C->U in the
genomic codon turns the genomic amino acid into the /translation amino acid are 'annotated edit sites'.
Compare with the caller's RDD list, and report raw pileup at a few sites (both gene strands).
usage: pcr_known_edits.py gb fa bam rdd.tsv"""
import sys, csv
sys.path.insert(0, 'src/followup')  # run from the folder root
import rdd
from Bio import SeqIO
from Bio.Seq import Seq
gb, fa, bam, rddtsv = sys.argv[1:5]
rec = SeqIO.read(gb, 'genbank')
G = str(rec.seq).upper()
COMP = str.maketrans('ACGT', 'TGCA')
ann = {}   # pos1 -> (gene, strand, cp)
nmis = 0; ncds = 0
for f in rec.features:
    if f.type != 'CDS' or 'translation' not in f.qualifiers: continue
    ncds += 1
    name = f.qualifiers.get('gene', ['?'])[0]
    tr = f.qualifiers['translation'][0]
    pos = []
    for part in f.location.parts:
        rng = range(int(part.start), int(part.end))
        if part.strand == -1: rng = reversed(rng)
        pos += [(p + 1, part.strand) for p in rng]
    nt = ''.join(G[p - 1] if s == 1 else G[p - 1].translate(COMP) for p, s in pos)
    for i in range(min(len(tr), len(nt) // 3)):
        cod = nt[3 * i:3 * i + 3]
        aa = str(Seq(cod).translate(table=1))
        if aa == tr[i]: continue
        nmis += 1
        # try single and double C->U
        for k in range(3):
            if cod[k] == 'C':
                c2 = cod[:k] + 'T' + cod[k + 1:]
                if str(Seq(c2).translate(table=1)) == tr[i] or (i == 0 and c2 == 'ATG'):
                    p, s = pos[3 * i + k]; ann[p] = (name, s, k + 1)
        if all(str(Seq(cod[:k] + 'T' + cod[k + 1:]).translate(table=1)) != tr[i] for k in range(3) if cod[k] == 'C'):
            cs = [k for k in range(3) if cod[k] == 'C']
            if len(cs) >= 2:
                c2 = ''.join('T' if j in cs else cod[j] for j in range(3))
                if str(Seq(c2).translate(table=1)) == tr[i]:
                    for k in cs:
                        p, s = pos[3 * i + k]; ann[p] = (name, s, k + 1)
print(f'CDS with translation {ncds}; aa differing genomic vs /translation {nmis}; annotated C->U edit nt sites {len(ann)}')
called = {int(r['pos']): r for r in csv.DictReader(open(rddtsv), delimiter='\t')}
pr = rdd.pileup(fa, bam)
dep = {p: sum(c.values()) for p, c in pr.items()}
a20 = [p for p in ann if dep.get(p, 0) >= 20]
hit = [p for p in a20 if p in called]
print(f'annotated sites at depth>=20: {len(a20)}; called as RDD: {len(hit)} ({len(hit)/max(1,len(a20)):.2f})')
print(f'called RDDs that are annotated edit sites: {sum(1 for p in called if p in ann)} of {len(called)}')
# edit fraction distribution at annotated depth>=20 sites
fr = []
for p in a20:
    name, s, cp = ann[p]; c = pr[p]; ref = G[p - 1]
    alt = 'T' if s == 1 else 'A'
    fr.append(c[alt] / sum(c.values()))
fr.sort()
if fr: print('edit fraction at annotated sites: quartiles', [round(fr[int(q * (len(fr) - 1))], 3) for q in (0, .25, .5, .75, 1)])
print('missed annotated sites (depth>=20), first 10: pos gene strand cp ref counts')
for p in [p for p in a20 if p not in called][:10]:
    print(' ', p, *ann[p], G[p - 1], dict(pr[p]))
print('hand check, 3 + strand and 3 - strand called sites: pos gene strand cp type | genome base, pileup counts (+ strand bases)')
for st in ('1', '-1'):
    for p in [p for p, r in called.items() if r['strand'] == st][:3]:
        r = called[p]
        print(' ', p, r['gene'], st, r['codon_pos'], r['type'], '| ref', G[p - 1], dict(pr[p]), 'annotated' if p in ann else 'NOT annotated')
