#!/bin/bash
# PC-T5a (RULES_T5.md): 500 simulated pairs from each A. secundus 16S copy, same map/filter/count path as the real run.
set -euo pipefail
D=$1; REF=out/t5/ref.fa
for c in div norm; do
  awk -v c=$c '/^>/{p=($1==">"c)} p' out/t5/sec16.fa > $D/$c.fa
  wgsim -N 500 -1 151 -2 151 -d 300 -s 50 -e 0.01 -r 0 -R 0 -X 0 -S 7 $D/$c.fa $D/$c.1.fq $D/$c.2.fq > /dev/null
  # interleave, both mates named sim_<class>.<n> as fastq-dump --split-spot does
  paste <(paste - - - - < $D/$c.1.fq) <(paste - - - - < $D/$c.2.fq) | awk -F'\t' -v c=$c '{n="@sim_" c "." NR; print n"\n"$2"\n+\n"$4"\n"n"\n"$6"\n+\n"$8}'
done | minimap2 -x sr -t 1 -c --secondary=no $REF - 2> $D/mm2.log | awk -F'\t' '$11>=60 && $10/$11>=0.98' > $D/sim.paf
