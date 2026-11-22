#!/bin/bash
# usage: run_codetta.sh R01 R02 ...  (sequential lane)
C=scratchpad/blind_cun/dl.bllj/codetta
T=scratchpad/blindagg.4ceG
export HMMER_NCPU=1
for g in "$@"; do
  mkdir -p $T/work/codetta/$g
  cd $T/work/codetta/$g
  echo "$(date) start $g"
  uv run --no-project --with numpy,scipy python $C/codetta_align.py $T/work/genomes/$g.fa --align_output $T/work/codetta/$g/$g > align.log 2>&1
  uv run --no-project --with numpy,scipy python $C/codetta_summary.py $T/work/codetta/$g/$g > summary.log 2>&1
  echo "$(date) done align $g"
done
