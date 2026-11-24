#!/bin/bash
# RNA step (RULES.md amendment 3): map HC1 then SC1 (one at a time), then annot.py unchanged per unit -> out/units/<AP acc>.{removals,runs}.tsv
set -e -o pipefail
cd /Users/ryan/lab/work/archaeal-intron-census
mkdir -p out/units
for u in "GCA_055385335.1 DRR620298 AP029002.1" "GCA_055394625.1 DRR620299 AP029247.1"; do
  set -- $u
  echo "$(date '+%T') map $1 $2"
  src/rna_map.sh $1 $2 8000000 2
  echo "$(date '+%T') n=$(cat out/rna/$1__$2.n) mm2: $(tail -1 out/rna/$1__$2.mm2.log)"
  uv run --no-project --with numpy python -I /Users/ryan/lab/work/archaeal-introns/survey/src/annot.py $1 out/rna/$1__$2 \
    --genome data/ref/$1/genomic.fna.gz --gff data/ref/$1/genomic.gff.gz --out out/units/$3.removals.tsv
  echo "$(date '+%T') annot $3 done: $(cat out/units/$3.runs.tsv)"
done
