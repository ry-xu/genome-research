#!/bin/bash
# Cut finished ledger batches whenever the 1-min load is under 8; exit when RUN has N batches and all are cut.
# usage: src/cut_loop.sh out/ledger1 out/cuts1 85
cd /Users/ryan/lab/work/archaeal-intron-census
run=$1; out=$2; n=$3
while true; do
  done_b=$(ls $run/genomes 2>/dev/null | wc -l | tr -d ' ')
  load=$(sysctl -n vm.loadavg | awk '{print $2}')
  if awk -v l=$load 'BEGIN{exit !(l<8)}'; then
    nice -n 10 uv run --project /Users/ryan/lab/work/reading-ledger python src/cuts.py $run $out >> $out.log 2>&1
  fi
  cut_b=$(ls $out 2>/dev/null | grep -c parquet)
  echo "$(date '+%T') load $load batches $done_b cut $cut_b" >> $out.log
  [ "$done_b" -ge "$n" ] && [ "$cut_b" -ge "$n" ] && break
  sleep 120
done
