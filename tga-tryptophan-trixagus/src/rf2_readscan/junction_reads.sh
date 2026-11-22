#!/bin/zsh
# Stream Tder reads (SRR28034346), map to Tder genome + Tcar/Tmey RF2 regions; keep mapped records only.
D=${FACTCHECK_DIR:-factcheck2}/dl_reads
cd $D
fastq-dump -X 4000000 --split-spot --skip-technical --stdout SRR28034346 2> fqd.err \
 | minimap2 -t 2 -ax sr -k 15 -w 5 ref.fa - 2> mm2.err \
 | awk '/^@/ || $3!="*"' > mapped.sam
echo "done $(date +%H:%M:%S)" > done.txt
