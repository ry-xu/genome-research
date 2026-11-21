"""Missed-copy check: blastn a P. freudenreichii 23S locus (DSM 20271, from out_pf/copies.tsv) at relaxed e-value against every
P. freudenreichii assembly; list every HSP locus, flag any not among the loci already read."""
import sys,glob,os,subprocess,collections
comp=str.maketrans('ACGT','TGCA'); rc=lambda x:x.translate(comp)[::-1]
fa=glob.glob('data/dl/GCA_000940845.1/*.fna')[0]
s=''.join(l.strip() for l in open(fa) if l[0]!='>').upper()
q=rc(s[1503211-1:1506217]); open('ref_pfreud_23S_locus.fna','w').write('>pf_DSM20271_23S_locus\n'+q+'\n')
known=collections.defaultdict(list)
for l in list(open('out_pf/copies.tsv'))[1:]:
    f=l.split('\t'); known[(f[0],f[3])].append((int(f[6]),int(f[7])))
tot=collections.Counter(); extra=[]
for d in sorted(glob.glob('data/dl/GCA_*')):
    acc=os.path.basename(d); db=f'out_pf/blastdb/{acc}'
    r=subprocess.run(['blastn','-task','blastn','-query','ref_pfreud_23S_locus.fna','-db',db,'-evalue','1e-5','-dust','no','-num_threads','1',
        '-outfmt','6 sseqid sstart send length pident qstart qend'],capture_output=True,text=True,check=True).stdout
    n=0
    for l in r.splitlines():
        sid,ss,se,ln,pid,qs,qe=l.split('\t'); lo,hi=sorted((int(ss),int(se))); n+=1
        if not any(lo<=b+50 and hi>=a-50 for a,b in known[(acc,sid)]): extra.append((acc,sid,lo,hi,ln,pid,qs,qe))
    tot[acc]=n
print('HSPs per assembly:',dict(tot))
print('HSPs outside read loci:',len(extra))
for e in extra: print(' ',*e)
