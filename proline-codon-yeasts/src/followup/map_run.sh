#!/bin/bash
# usage: map_run.sh RUN ref.fa out_prefix [spots=4000000] [threads]
# Streams RUN with fastq-dump (first N spots, split spots), maps with minimap2 -ax sr, keeps MAPQ>=20 primary, sorted BAM.
set -euo pipefail
RUN=$1; REF=$2; OUT=$3; N=${4:-4000000}; T=${5:-1}
date +%H:%M:%S > $OUT.log
fastq-dump -X $N --split-spot --stdout $RUN 2>>$OUT.log | minimap2 -ax sr -t $T $REF - 2>>$OUT.log \
  | samtools view -u -q 20 -F 0x904 - | samtools sort -m 500M -o $OUT.bam - && samtools index $OUT.bam
echo "spots_requested $N" >> $OUT.log; samtools view -c $OUT.bam >> $OUT.log; date +%H:%M:%S >> $OUT.log
