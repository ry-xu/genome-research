#!/bin/bash
T=scratchpad/blindagg.4ceG
g=$1
if [ $g = CTRL ]; then fa=$T/work/control/CTRL.fa; else fa=$T/work/genomes/$g.fa; fi
cd $T/work/m2
for v in v2 v3; do uv run --no-project python -I $T/scripts/m2.py $fa $T/data/uniprot/refs.fa $T/work/m2 $g $v; done
