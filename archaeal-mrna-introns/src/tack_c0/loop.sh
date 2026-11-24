#!/bin/zsh
# Sequential runs from a list (RUN SPOTS per line); 2 threads if load<8, else 1 thread under nice (BRIEF 11:30 priority note:
# biology agents use 1 thread whenever load > 12; better-MAGs has CPU priority, hence nice). usage: src/loop.sh LIST OUTDIR
cd /Users/ryan/lab/work/tack-rna
while read run spots; do
  [[ -e $2/$run.done ]] && continue
  l=$(sysctl -n vm.loadavg | awk '{print $2}')
  if (( l < 8 )); then t=2; else t=1; fi
  df_free=$(df -g /Users/ryan | awk 'NR==2{print $4}')
  (( df_free < 60 )) && { echo "$(date '+%T') disk $df_free GB, stop"; exit 3; }
  echo "$(date '+%T') $run $spots t=$t"
  nice -n 10 src/map_run.sh $run $spots $t $2 || echo "$(date '+%T') FAIL $run"
  echo "$(date '+%T') $run n=$(cat $2/$run.n) $(cat $2/$run.cls.log) bytes=$(cat $2/$run.bytes)"
done < $1
echo "$(date '+%T') loop done"
