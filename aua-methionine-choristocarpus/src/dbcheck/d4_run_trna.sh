#!/bin/bash
# Run tRNAscan-SE 2.0.13 (-O organellar, and -G general) and ARAGORN 1.2.41 (-t) on each genome, single thread.
set -u
B=lab/work/organelle-reading/checks/chor_wgs/dbcheck
T=~/dev/tools/bin
export PATH=$T:$PATH
mkdir -p $B/trna
for acc in NC_085326 CAXNAD020000071.1 CAXNAE020000396.1 NC_013473.1 ON552993.1 NC_007684.1; do
  f=$B/data/$acc.fa
  nice -n 10 $T/tRNAscan-SE -O --thread 1 -q -Q --detail -o $B/trna/$acc.tscanO.txt -f $B/trna/$acc.tscanO.ss $f > $B/trna/$acc.tscanO.log 2>&1
  nice -n 10 $T/tRNAscan-SE -G --thread 1 -q -Q --detail -o $B/trna/$acc.tscanG.txt $f > $B/trna/$acc.tscanG.log 2>&1
  nice -n 10 $T/aragorn1.2.41 -t -gc1 -w -o $B/trna/$acc.aragorn.txt $f
  nice -n 10 $T/aragorn1.2.41 -t -gc1 -ps -w -o $B/trna/$acc.aragorn_ps.txt $f
done
echo done
