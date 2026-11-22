import sys, json
B='TCAG'; codons=[a+b+c for a in B for b in B for c in B]
aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; STD=dict(zip(codons,aas))
def parse(p):
    inf={}; lp={}; hdr=None; sec=None
    for l in open(p):
        if l.startswith('# Codon inferences'): sec='inf'; continue
        if l.startswith('# Log decoding'): sec='lp'; continue
        if l.startswith('# Final'): sec=None; continue
        if l.startswith('# codon') and sec=='lp': hdr=[x[5] for x in l.split()[2:]]; continue
        if l.startswith('#') or not l.strip(): continue
        f=l.split()
        if sec=='inf' and f[0] in STD: inf[f[0]]=dict(inf=f[1], nal=int(f[4]), nused=int(f[5]))
        elif sec=='lp' and f[0] in STD: lp[f[0]]=dict(zip(hdr, map(float, f[1:])))
    out={}
    for c in codons:
        d=dict(inf[c]); v=lp.get(c)
        if v:
            o=sorted(v, key=lambda a:-v[a]); d['best']=o[0]; d['second']=o[1]
            d['gap_2nd']=v[o[0]]-v[o[1]]
            d['gap_std']=(v[o[0]]-v[STD[c]]) if STD[c]!='*' else None
            d['logP']=v
        out[c]=d
    return out
if __name__=='__main__':
    res={}
    for name,p in zip(sys.argv[1::2], sys.argv[2::2]):
        res[name]=parse(p)
        for c in codons:
            d=res[name][c]
            if (d['inf']!='?' and d['inf']!=STD[c]) or STD[c]=='*':
                gs = '%.1f'%d['gap_std'] if d.get('gap_std') is not None else 'n/a'
                print(name, c, 'std', STD[c], 'codetta', d['inf'], 'Nal', d['nal'], 'Nused', d['nused'], 'best', d.get('best'), '2nd', d.get('second'), 'gap2nd %.1f'%d.get('gap_2nd',0), 'gapstd', gs)
    json.dump(res, open(T+'/work/codetta_parsed.json','w') if False else sys.stdout.__class__ and open(sys.argv[0].rsplit('/',2)[0]+'/work/codetta_parsed.json','w'), indent=0)
