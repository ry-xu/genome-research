#!/bin/bash
# Fetch the first N spots of one SRA run (mags fetch, shared store), map to the frozen rRNA reference, count (RULES.md).
#   bash src/run_one.sh RUN N
set -euo pipefail
RUN=$1; N=$2; HERE=lab/work/silent-operon/C_rna
REF=lab/work/orthogonal-ribosome/out/t5/ref.fa
export MAGS_AGENT=silent-C_rna
O=$HERE/out/${RUN}_head$N
mkdir $O.lock 2>/dev/null || { echo "skip $RUN (taken)"; exit 0; }
P=$(~/dev/projects/better-mags/infra/mags fetch sra $RUN --spots $N 2> $O.fetch.log | tail -2)
R1=$(echo "$P" | head -1); R2=$(echo "$P" | tail -1)
echo "$R1 $R2" > $O.paths
ls -la $R1 $R2 >> $O.paths
TOT=$(( $(gzip -dc $R1 | wc -l) / 4 ))
seqtk mergepe $R1 $R2 | minimap2 -x sr -t 1 -c --secondary=no $REF - 2> $O.mm2.log \
  | awk -F'\t' '$11>=60 && $10/$11>=0.85' > $O.paf
python3 -I $HERE/src/count.py $O.paf $HERE/ref/classes.tsv $TOT > $O.counts.txt
echo "done $RUN $N pairs $TOT"
