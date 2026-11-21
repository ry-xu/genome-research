import sys
fa=sys.argv[1]
seq=''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
print(len(seq))
# rrnH region per my recollection of NC_000913.3: rrsH ~223771..225312, rrlH ~225759..228662
for name,s,e in [('16S',223771,225312),('23S',225759,228662)]:
    g=seq[s-1:e]
    print(name,len(g),g[:40],'...',g[-30:])
    if name=='23S':
        print('2652-2667:',g[2651:2667])
    open(f'ref_ecoli_{name}.fna','w').write(f'>ecoli_rrnH_{name}_{s}_{e}\n{g}\n')
