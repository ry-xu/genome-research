#!/bin/bash
# tRNAscan-SE -B (bacterial models; reports isotype-specific CM, which separates CAU tRNAs into iMet / Met / Ile2) on each genome.
set -u
B=lab/work/organelle-reading/checks/chor_wgs/dbcheck
T=~/dev/tools/bin
export PATH=$T:$PATH
for acc in NC_085326 CAXNAD020000071.1 CAXNAE020000396.1 NC_013473.1 ON552993.1 NC_007684.1; do
  nice -n 10 $T/tRNAscan-SE -B --thread 1 -q -Q --detail -o $B/trna/$acc.tscanB.txt -s $B/trna/$acc.tscanB.iso.txt -f $B/trna/$acc.tscanB.ss $B/data/$acc.fa > $B/trna/$acc.tscanB.log 2>&1
done
echo done
