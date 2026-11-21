import re,collections
summ={l.split('\t')[0]:l.rstrip('\n').split('\t') for l in list(open('out_pf/summary.tsv'))[1:]}
meta=[l.rstrip('\n').split('\t') for l in open('meta_out.txt') if l.startswith('GCA_')]
chrom=[m for m in meta if int(summ[m[0]][8])>0]
def norm(n): n=n.upper().replace(' ','').replace('-GENOME',''); n=re.sub(r'^PFRJS-?','JS',n); n=re.sub(r'^PFREUDJS1$','JS1',n); return n
names=[norm(m[1]) for m in chrom]
print('chromosome assemblies:',len(chrom),' distinct raw names:',len(set(m[1] for m in chrom)),' distinct normalized names:',len(set(names)))
base=[re.sub(r'-\d$','',n) for n in names]
c=collections.Counter(base); print('names collapsing sub-isolate suffix (-1):',len(set(base)),' repeated:',{k:v for k,v in c.items() if v>1})
# also treat identical-sequence MGYG-HGUT-02401 as P.UF1
base2=['P.UF1' if b=='MGYG-HGUT-02401' else b for b in base]; print('after MGYG=P.UF1:',len(set(base2)))
print('type material flagged:',[(m[0],m[1]) for m in chrom if m[4]])
print('submitters:',dict(collections.Counter(m[8] for m in chrom)))
print('tech:',dict(collections.Counter(m[9] for m in chrom)))
