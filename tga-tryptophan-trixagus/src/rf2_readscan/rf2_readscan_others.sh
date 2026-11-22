#!/bin/zsh
# Same protein-level RF2 read scan for the other two readers (Tser, Ttur), first 4M spots each, run one after the other.
F=${FACTCHECK_DIR:-factcheck2}
for RUN in SRR28034353 SRR28034362; do
  R=$F/dl_reads_$RUN; mkdir -p $R; cd $R
  diamond makedb --in $F/rfdb.faa -d $R/rfdb > mk.log 2>&1
  fastq-dump -X 4000000 --split-spot --skip-technical --stdout $RUN 2> fqd.err \
   | awk 'NR%4==1{print ">" substr($1,2) "_" (++n)} NR%4==2{print}' | gzip -1 > reads.fa.gz
  diamond blastx -q reads.fa.gz -d rfdb --query-gencode 4 --threads 2 --evalue 1e-5 --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore -o hits.tsv > dmd.log 2>&1
  echo "done $(date +%H:%M:%S)" > done.txt
done
