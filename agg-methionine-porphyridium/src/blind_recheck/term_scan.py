"""For each gene family: best HSP reaching the reference C-terminus (qend >= qlen-20); scan downstream in-frame
up to 60 codons and report the first codon that is a table-1 stop or in the candidate list."""
import sys, collections
genome_fa, blast_tsv, refs_fa = sys.argv[1:4]; cand = set(sys.argv[4].split(',')) if len(sys.argv) > 4 and sys.argv[4] else set()
comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
seq=''.join(l.strip() for l in open(genome_fa) if not l.startswith('>')).upper(); rc=seq.translate(comp)[::-1]
fr={}
for st,s in (('+',seq),('-',rc)):
    for f in range(3): fr['%s%d'%(st,f)]=[s[i:i+3] for i in range(f,len(s)-2,3)]
qlen={}; cur=None
for l in open(refs_fa):
    if l.startswith('>'): cur=l[1:].split()[0]; qlen[cur]=0
    else: qlen[cur]+=len(l.strip())
best={}
for line in open(blast_tsv):
    q,s,ev,bits,qs,qe,ss,se,qa,sa=line.rstrip('\n').split('\t')
    fam=q.split('|')[1]
    if int(qe) < qlen[q]-20: continue
    if fam not in best or float(bits) > best[fam][0]: best[fam]=(float(bits), s, int(se), q)
stops={'TAA','TAG','TGA'}
tally=collections.Counter()
for fam,(bits,s,se,q) in sorted(best.items()):
    cl=fr[s]; hit=None
    for k in range(se, min(len(cl), se+60)):   # se is 1-based last aligned residue -> index se is next codon
        if cl[k] in stops or cl[k] in cand: hit=(k-se+1, cl[k]); break
    tally[hit[1] if hit else 'none']+=1
    print(fam, s, 'alnEnd', se, 'first terminator-like codon after', hit)
print('TALLY', dict(tally))
