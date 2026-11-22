#!/bin/bash
# Re-map an already fetched head (no download) and count. bash src/map_only.sh RUN N
set -euo pipefail
RUN=$1; N=$2; HERE=lab/work/silent-operon/C_rna
REF=lab/work/orthogonal-ribosome/out/t5/ref.fa
O=$HERE/out/${RUN}_head$N
D=~/dev/projects/better-mags/data/shared/sra/${RUN}_head$N
R1=$D/${RUN}_head${N}_1.fq.gz; R2=$D/${RUN}_head${N}_2.fq.gz
TOT=$(( $(gzip -dc $R1 | wc -l) / 4 ))
seqtk mergepe $R1 $R2 | minimap2 -x sr -t 1 -c --secondary=no $REF - 2> $O.mm2.log \
  | awk -F'\t' '$11>=60 && $10/$11>=0.85' > $O.paf
python3 -I $HERE/src/count.py $O.paf $HERE/ref/classes.tsv $TOT > $O.counts.txt
echo "done $RUN $N pairs $TOT (remap)"
