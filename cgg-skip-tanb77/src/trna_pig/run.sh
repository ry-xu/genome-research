#!/bin/zsh
# tRNAscan-SE 2 (Infernal) on the 24 pig TANB77 representatives: every tRNA with its secondary structure, to look for
# a tRNA-Gly (or any tRNA) with an 8-nt anticodon loop that could read CGGN as Gly (ARAGORN may not call such loops).
cd .
for f in data/proteomics/mgnify_pig/*.fna; do
  g=${f:t:r}
  [[ -s scratch/trnascan_pig/$g.ss ]] && continue
  tRNAscan-SE -B -Q --thread 2 -o scratch/trnascan_pig/$g.out -f scratch/trnascan_pig/$g.ss $f > scratch/trnascan_pig/$g.log 2>&1 || echo "FAILED $g"
  echo "$(date +%H:%M) $g $(grep -c 'Type:' scratch/trnascan_pig/$g.ss)"
done
echo "$(date +%H:%M) all done"
