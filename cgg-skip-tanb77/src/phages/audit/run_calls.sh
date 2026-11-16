#!/bin/bash
# sequential heavy steps; waits for 1-min load <= 10 before each
cd phages_work/kill1
waitload() { while true; do l=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,); if awk -v l="$l" 'BEGIN{exit !(l<=10)}'; then break; fi; echo "$(date +%H:%M) load $l, waiting"; sleep 60; done; }
for m in meta single orf; do
  waitload; echo "$(date +%H:%M) call $m"
  /usr/bin/time -l uv run --no-project --with pyrodigal python -I call_genes.py ../data/phages.fna ../data/groups/run_ids.tsv $m out/genes_$m.tsv > out/call_$m.log 2>&1 || { echo fail $m; exit 1; }
  waitload; echo "$(date +%H:%M) sites $m"
  uv run --no-project python -I skip_sites.py ../data/phages.fna out/genes_$m.tsv ../data/groups/run_ids.tsv out/$m > out/sites_$m.log 2>&1 || { echo fail sites $m; exit 1; }
done
echo "$(date +%H:%M) all done"
