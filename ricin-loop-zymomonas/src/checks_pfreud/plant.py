"""Plant: change the SRL loop of ONE 23S copy in C. acnes and E. coli to GCAA / UUCG (by exact context replacement at first occurrence on either strand)."""
import sys,glob,os,shutil
comp=str.maketrans('ACGT','TGCA'); rc=lambda x:x.translate(comp)[::-1]
for acc,ctx,new in [('GCF_000008345.1','TTAGTACGAGAGGACC','TTAGTACGCAAGGACC'),('GCF_000005845.2','CTAGTACGAGAGGACC','CTAGTACTTCGGGACC')]:
    d=f'data/dl/{acc}'; od=f'data/plant/{acc}'; os.makedirs(od,exist_ok=True)
    shutil.copy(f'{d}/assembly_data_report.jsonl',od)
    fa=glob.glob(f'{d}/*_genomic.fna')[0]
    hdr=open(fa).readline(); s=''.join(l.strip() for l in open(fa) if l[0]!='>').upper()
    n_f=s.count(ctx); n_r=s.count(rc(ctx))
    if n_f: s=s.replace(ctx,new,1); where='+'
    else: s=s.replace(rc(ctx),rc(new),1); where='-'
    open(f'{od}/planted_genomic.fna','w').write(hdr+s+'\n')
    print(acc,'ctx occurrences +',n_f,'-',n_r,'planted one copy on',where,'->',new)
