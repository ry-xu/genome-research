"""Tally the 4 bases after a 27-nt P. freudenreichii SRL upstream anchor (ending at E. coli C2658) in reads from stdin (fasta)."""
import sys,collections
A='AATCTTGAGAAGGGCTGTCCTTAGTAC'
comp=str.maketrans('ACGTN','TGCAN')
c=collections.Counter(); n=0; gc=0; tot=0
for l in sys.stdin:
    if l[0]=='>': continue
    r=l.strip().upper(); n+=1
    if n<=200000: gc+=r.count('G')+r.count('C'); tot+=len(r)
    for x in (r, r.translate(comp)[::-1]):
        i=x.find(A)
        if i>=0 and i+len(A)+8<=len(x):
            c[x[i+len(A):i+len(A)+4]+'|'+x[i+len(A)+4:i+len(A)+8]]+=1
print(sys.argv[1],'reads',n,'GC%%(first200k) %.1f'%(100*gc/max(tot,1)),'anchor_hits',sum(c.values()),dict(c.most_common(6)))
