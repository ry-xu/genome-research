#!/bin/bash
cd lab/work/silent-operon/C_rna
for r in SRR33092282 SRR33092281 SRR33092280 SRR7589086 SRR6261448; do bash src/run_one.sh $r 1000000 || echo "fail $r"; done
