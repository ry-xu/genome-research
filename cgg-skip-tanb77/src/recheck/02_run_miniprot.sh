#!/bin/bash
# miniprot, no splicing (-S), standard table (-T 1 default; CGG scored as Arg), 4 threads.
W=recheck_work
for L in G1 G2 G3 G4 G5 G6 G7 G8 G9; do
  [ -s $W/out/$L.paf ] || miniprot -S -t 4 $W/genomes/$L.fna $W/refs/refs_all.faa > $W/out/$L.paf 2> $W/out/$L.log
  # prfB positive control: full-length RF2 (B. subtilis, C. perfringens, E. coli) from Swiss-Prot
  miniprot -S -t 2 --aln $W/genomes/$L.fna $W/refs/rf2_full.faa > $W/out/$L.rf2.paf 2> /dev/null
done
