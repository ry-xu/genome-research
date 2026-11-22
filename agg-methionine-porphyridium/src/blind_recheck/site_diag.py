"""Diagnose sites of one codon: per site, local identity (+-10 aligned cols, X excluded), refs, genes."""
import sys, collections
genome_fa, blast_tsv, codon = sys.argv[1:4]
B='TCAG'; cods=[a+b+c for a in B for b in B for c in B]
aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; STD=dict(zip(cods,aas))
comp=str.maketrans('ACGTRYKMSWBDHVN','TGCAYRMKSWVHDBN')
seq=''.join(l.strip() for l in open(genome_fa) if not l.startswith('>')).upper(); rc=seq.translate(comp)[::-1]
fr={}
for st,s in (('+',seq),('-',rc)):
    for f in range(3): fr['%s%d'%(st,f)]=[s[i:i+3] for i in range(f,len(s)-2,3)]
best={}
for line in open(blast_tsv):
    q,s,ev,bits,qs,qe,ss,se,qa,sa=line.rstrip('\n').split('\t'); bits=float(bits); sp=int(ss)-1
    cols=[]
    for x,y in zip(qa,sa):
        if y!='-':
            cols.append((sp if x!='-' else None, x, y)); sp+=1
        else: cols.append((None,x,y))
    pairs=[i for i,c in enumerate(cols) if c[0] is not None]
    keep=pairs[5:len(pairs)-5]
    for i in keep:
        spos,x,y=cols[i]
        if fr[s][spos]!=codon: continue
        win=[cols[j] for j in range(max(0,i-10),min(len(cols),i+11)) if j!=i]
        idn=[1 if a==b else 0 for _,a,b in win if a!='-' and b!='-' and b!='X']
        lid=sum(idn)/len(idn) if idn else 0
        key=(q,(s,spos))
        if key not in best or best[key][0]<bits: best[key]=(bits,x,lid,q.split('|')[1])
persite=collections.defaultdict(list)
for (q,site),(bits,x,lid,g) in best.items(): persite[site].append((x,lid,g))
for site,v in sorted(persite.items()):
    lids=[l for _,l,_ in v]; genes=collections.Counter(g for *_,g in v)
    print(site, 'nref',len(v),'medLocalId %.2f'%sorted(lids)[len(lids)//2], 'refs',''.join(sorted(x for x,_,_ in v)), dict(genes))
