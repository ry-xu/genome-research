import sys, collections
fa, gp = sys.argv[1:3]
seqs = {}; name=None
for l in open(fa):
    if l[0]=='>': name=l[1:].split()[0]; seqs[name]=[]
    else: seqs[name].append(l.strip())
seqs = {k:''.join(v).upper() for k,v in seqs.items()}
C = str.maketrans('ACGT','TGCA'); c = collections.Counter()
for l in open(gp):
    f=l.split('\t')
    if f[0]=='id' or f[4]=='1' or f[5]=='1': continue
    s=seqs[f[0]]; b,e,st=int(f[1]),int(f[2]),int(f[3])
    g=s[b-1:e] if st==1 else s[b-1:e].translate(C)[::-1]
    c[(g[:3],g[-3:])]+=1
print(c.most_common(8))
