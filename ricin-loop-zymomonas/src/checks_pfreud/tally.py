import collections
rows=[l.rstrip('\n').split('\t') for l in list(open('out_pf/copies.tsv'))[1:]]
summ=[l.rstrip('\n').split('\t') for l in list(open('out_pf/summary.tsv'))[1:]]
print('P. freudenreichii zips scanned:',len(summ))
with23=[s for s in summ if int(s[8])>0]; no23=[s for s in summ if int(s[8])==0]
print('assemblies with >=1 23S copy:',len(with23),' without:',len(no23),[(s[0],s[6]) for s in no23])
print('23S copies total:',len(rows),' loop counts:',dict(collections.Counter(r[11] for r in rows)),' genomic re-read:',dict(collections.Counter(r[12] for r in rows)))
print('aln==genomic all:',all(r[11]==r[12] for r in rows),' ctx counts:',dict(collections.Counter(r[13] for r in rows)))
print('16S loci == 23S loci in every assembly:',all(s[7]==s[8] for s in summ))
rec=collections.Counter((r[0],r[3]) for r in rows); print('copies per chromosome record:',dict(collections.Counter(rec.values())),' records:',len(rec))
non=[r for r in rows if r[11]!='GCAA']; print('non-GCAA copies:',len(non))
# distinct genomes: drop exact-duplicate sequences (meta_out identical groups) and same-name re-assemblies
dup_drop={'GCA_902386715.1','GCA_059735485.1','GCA_059735495.1','GCA_965286265.1'}
print('chromosome assemblies after removing identical-sequence duplicates (P.UF1=MGYG-HGUT-02401, TL18 x2, PB4 x2) and J117 re-assembly:',len([s for s in with23 if s[0] not in dup_drop]))
copies_dedup=[r for r in rows if r[0] not in dup_drop]; print('copies in deduplicated set:',len(copies_dedup),dict(collections.Counter(r[11] for r in copies_dedup)))
