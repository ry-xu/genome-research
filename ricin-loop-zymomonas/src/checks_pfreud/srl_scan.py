"""Independent SRL check: blastn E. coli rrnH 23S/16S against each genome,
map E. coli 23S positions 2652-2667 through the alignment, read the loop (2659-2662)."""
import sys, os, json, subprocess, glob, collections
REF23='ref_ecoli_23S.fna'; REF16='ref_ecoli_16S.fna'
OUT=sys.argv[1]; dirs=sys.argv[2:]
tmp=os.path.join(OUT,'blastdb'); os.makedirs(tmp,exist_ok=True)

def readfa(p):
    d={};n=None;b=[]
    for l in open(p):
        if l.startswith('>'):
            if n: d[n]=''.join(b).upper()
            n=l[1:].split()[0]; hdr[n]=l[1:].strip(); b=[]
        else: b.append(l.strip())
    if n: d[n]=''.join(b).upper()
    return d
comp=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(comp)[::-1]

def blast(q,db):
    fmt='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq slen'
    r=subprocess.run(['blastn','-task','blastn','-query',q,'-db',db,'-evalue','1e-10','-dust','no',
        '-num_threads','1','-max_target_seqs','500','-outfmt',fmt],capture_output=True,text=True,check=True)
    hits=[]
    for l in r.stdout.splitlines():
        f=l.split('\t')
        hits.append(dict(s=f[1],pid=float(f[2]),qs=int(f[4]),qe=int(f[5]),ss=int(f[6]),se=int(f[7]),
                         bits=float(f[9]),qseq=f[10],sseq=f[11],slen=int(f[12])))
    return hits

def loci(hits,gap=4000):
    # merge HSPs on same contig+strand within gap into loci
    out=[]
    for h in sorted(hits,key=lambda h:(h['s'],h['ss']<=h['se'],min(h['ss'],h['se']))):
        lo,hi=min(h['ss'],h['se']),max(h['ss'],h['se']); st='+' if h['ss']<=h['se'] else '-'
        if out and out[-1]['s']==h['s'] and out[-1]['st']==st and lo<=out[-1]['hi']+gap:
            L=out[-1]; L['hi']=max(L['hi'],hi); L['hsps'].append(h)
        else: out.append(dict(s=h['s'],st=st,lo=lo,hi=hi,hsps=[h],slen=h['slen']))
    for L in out:
        cov=set()
        for h in L['hsps']: cov.update(range(h['qs'],h['qe']+1))
        L['qcov']=len(cov)
    return out

def map_q(h,a,b):
    """return subject letters aligned to query positions a..b (incl. inserts), and subject coords"""
    qp=h['qs']-1; sp=h['ss']; step=1 if h['se']>=h['ss'] else -1; sp-=step
    seg=[];scoords=[]
    for qc,sc in zip(h['qseq'],h['sseq']):
        if qc!='-': qp+=1
        if sc!='-': sp+=step
        if a<=qp<=b and not (qc=='-' and qp==b):
            seg.append(sc if qc!='-' else sc.lower()); 
            if sc!='-': scoords.append(sp)
    return ''.join(seg),scoords

hdr={}
rows=[];summ=[]
for d in dirs:
    acc=os.path.basename(d.rstrip('/'))
    fna=glob.glob(os.path.join(d,'*_genomic.fna'))[0]
    rep=json.loads(open(os.path.join(d,'assembly_data_report.jsonl')).readline())
    seqs=readfa(fna)
    db=os.path.join(tmp,acc)
    subprocess.run(['makeblastdb','-in',fna,'-dbtype','nucl','-out',db],capture_output=True,check=True)
    h23=blast(REF23,db); h16=blast(REF16,db)
    L23=loci(h23); L16=loci(h16)
    ai=rep.get('assemblyInfo',{}); org=rep.get('organism',{})
    strain=(org.get('infraspecificNames') or {}).get('strain') or ai.get('biosample',{}).get('strain') or ''
    tm=(rep.get('typeMaterial') or {}).get('typeLabel','')
    loops=[]
    for i,L in enumerate(L23):
        cand=[h for h in L['hsps'] if h['qs']<=2652 and h['qe']>=2667]
        partial_end = L['lo']<=50 or L['hi']>=L['slen']-50
        if cand:
            h=max(cand,key=lambda h:h['bits'])
            ctx,_=map_q(h,2652,2667); loop,sc=map_q(h,2659,2662)
            # independent genomic re-read: take subject coords of the loop, pull from contig with strand
            s=seqs[h['s']]
            if h['se']>=h['ss']: g=s[min(sc)-1:max(sc)]
            else: g=rc(s[min(sc)-1:max(sc)])
            loop_g=g
        else:
            ctx=loop=loop_g='NA'
        loops.append(loop)
        rows.append([acc,org.get('organismName',''),strain,L['s'],hdr[L['s']][:70],L['st'],L['lo'],L['hi'],L['qcov'],
                     'contig_end' if partial_end else '',len(L['hsps']),loop,loop_g,ctx])
    nc=len(seqs)
    summ.append([acc,org.get('organismName',''),strain,ai.get('assemblyLevel',''),tm,nc,sum(len(v) for v in seqs.values()),
                 len(L16),len(L23),sum(1 for L in L23 if L['qcov']>=2700),dict(collections.Counter(loops)),
                 ai.get('sequencingTech',''),ai.get('biosample',{}).get('accession',''),ai.get('submitter','')])
    print('SUMMARY',*summ[-1],sep='\t',flush=True)
with open(os.path.join(OUT,'copies.tsv'),'w') as f:
    f.write('\t'.join('acc organism strain contig header strand lo hi qcov23 contig_end nhsp loop_aln loop_genomic ctx2652_2667'.split())+'\n')
    for r in rows: f.write('\t'.join(map(str,r))+'\n')
with open(os.path.join(OUT,'summary.tsv'),'w') as f:
    f.write('\t'.join('acc organism strain level type ncontig len n16S_loci n23S_loci n23S_full loops tech biosample submitter'.split())+'\n')
    for r in summ: f.write('\t'.join(map(str,r))+'\n')
