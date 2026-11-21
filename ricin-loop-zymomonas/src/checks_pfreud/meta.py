import sys,json,glob,os,hashlib,collections
rows=[]
for d in sorted(glob.glob('data/dl/GCA_*')):
    acc=os.path.basename(d)
    r=json.loads(open(f'{d}/assembly_data_report.jsonl').readline())
    ai=r['assemblyInfo']; bs=ai.get('biosample',{})
    attrs={a['name']:a.get('value','') for a in bs.get('attributes',[])}
    org=r.get('organism',{}); inf=org.get('infraspecificNames') or {}
    name=inf.get('strain') or inf.get('isolate') or attrs.get('strain') or attrs.get('isolate') or ''
    fna=glob.glob(f'{d}/*_genomic.fna')[0]
    seqs={};n=None
    for l in open(fna):
        if l[0]=='>': n=l[1:].split()[0]; seqs[n]=[]
        else: seqs[n].append(l.strip().upper())
    seqs={k:''.join(v) for k,v in seqs.items()}
    big=max(seqs,key=lambda k:len(seqs[k]))
    md5={k:hashlib.md5(v.encode()).hexdigest()[:10] for k,v in seqs.items()}
    rows.append(dict(acc=acc,name=name,strain_attr=attrs.get('strain',''),isolate=attrs.get('isolate',''),
        sub=attrs.get('sub_species',org.get('organismName','')),type=(r.get('typeMaterial') or {}).get('typeLabel',''),
        refcat=ai.get('refseqCategory',''),biosample=bs.get('accession',''),bioproject=ai.get('bioprojectAccession',''),
        submitter=ai.get('submitter',''),tech=ai.get('sequencingTech',''),src=attrs.get('isolation_source',''),
        recs=';'.join(f'{k}:{len(v)}:{md5[k]}' for k,v in seqs.items()),bigmd5=md5[big],biglen=len(seqs[big]),
        seqtitle=bs.get('description',{}).get('title','')[:60],paired=(ai.get('pairedAssembly') or {}).get('accession','')))
for x in rows: print('\t'.join(str(x[k]) for k in ['acc','name','strain_attr','isolate','type','refcat','biosample','bioproject','submitter','tech','src','biglen','bigmd5','recs','paired','seqtitle']))
print()
c=collections.Counter(x['name'].replace(' ','').upper() for x in rows); print('distinct names',len(c)); print({k:v for k,v in c.items() if v>1})
m=collections.defaultdict(list)
for x in rows: m[x['bigmd5']].append(x['acc']+'/'+x['name'])
print('identical largest-record groups:',[v for v in m.values() if len(v)>1])
b=collections.Counter(x['biosample'] for x in rows); print('dup biosamples',{k:v for k,v in b.items() if v>1})
