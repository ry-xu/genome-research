#!/bin/bash
# Pilot (RULES.md): first 1 M spots of each new run, one at a time, 1 thread.
cd lab/work/silent-operon/C_rna
for r in SRR11558439 SRR6261398 SRR6241672 SRR6241673 SRR6261448 SRR33092280 SRR33092281 SRR33092282 SRR7589086; do
  bash src/run_one.sh $r 1000000 || echo "fail $r"
done
