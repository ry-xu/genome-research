#!/bin/bash
# R2: map RNA + same-strain WGS for each relative (sequential, 1 thread, load-gated), species check, call.
cd "$(dirname "$0")/../results/rdd"
R="uv run --no-project python -I ../../src"
while read rna dna ref; do
  for run in $rna $dna; do
    until [ "$(sysctl -n vm.loadavg | awk '{print ($2<12)}')" = 1 ]; do sleep 30; done
    [ -s $run.bam.bai ] || ../../src/map_run.sh $run $ref.mmi $run 4000000 1
    $R/species_check.py $run
  done
  $R/rdd.py call $ref.fa $ref.genes.tsv $rna.bam $dna.bam > r2_$rna.txt
  cat r2_$rna.txt
done <<LIST
SRR5487602 SRR5938400 BK074530.1
SRR5942407 SRR5938403 BK074543.1
SRR5487603 SRR16989166 BK074251.1
LIST
echo done $(date +%H:%M)
