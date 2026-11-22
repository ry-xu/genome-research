#!/bin/zsh
# Protein-level scan of Tder reads (SRR28034346, first 4M spots) against all Shikimatogenerans proteins (+E. coli RF2).
F=${FACTCHECK_DIR:-factcheck2}
R=$F/dl_reads2
cd $R
diamond makedb --in $F/rfdb.faa -d $R/rfdb > mk.log 2>&1
fastq-dump -X 4000000 --split-spot --skip-technical --stdout SRR28034346 2> fqd.err \
 | awk 'NR%4==1{print ">" substr($1,2) "_" (++n)} NR%4==2{print}' | gzip -1 > reads.fa.gz
diamond blastx -q reads.fa.gz -d rfdb --query-gencode 4 --threads 2 --evalue 1e-5 --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore -o hits.tsv > dmd.log 2>&1
echo "done $(date +%H:%M:%S)" > done.txt
