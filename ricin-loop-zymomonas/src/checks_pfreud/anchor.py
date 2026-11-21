import sys,glob
comp=str.maketrans('ACGT','TGCA')
fa=glob.glob('data/dl/GCA_000940845.1/*.fna')[0]
s=''.join(l.strip() for l in open(fa) if l[0]!='>').upper()
rc=s.translate(comp)[::-1]
for strand,x in (('+',s),('-',rc)):
    i=x.find('TTAGTACGCAAGGACC')
    while i>=0:
        print(strand,i,'up30',x[i-22:i+8],'loop',x[i+8:i+12],'down',x[i+12:i+30]); i=x.find('TTAGTACGCAAGGACC',i+1)
