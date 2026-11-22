#!/bin/zsh
# Run Sage with a memory guard (kills only its own PID if RSS passes LIMIT_GB, default 7). Copied from work/proteome-code-census/src/run_sage.sh (Nettle), lock moved to this project; after Jade's
# work/tag-proteomics/src/run_sage.sh. usage: src/run_sage.sh <config.json> <fasta> <outdir> <spectra...>
SAGE=~/dev/tools/sage/sage-v0.14.7-aarch64-apple-darwin/sage
# One of this project's Sage jobs at a time (shared, loaded machine): re-run under an exclusive lockf lock.
[[ -z $SAGE_LOCKED ]] && SAGE_LOCKED=1 exec lockf -k lab/work/new-stops-euk/results/aga_ala/prot/sage.lock $0 "$@"
cfg=$1; fasta=$2; out=$3; shift 3
mkdir -p $out
$SAGE --disable-telemetry-i-dont-want-to-improve-sage --batch-size ${BATCH:-4} -f $fasta -o $out $cfg "$@" > $out/sage.log 2>&1 &
pid=$!
echo $pid > $out/sage.pid
peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' ')
  [[ -n $rss ]] && (( rss > peak )) && peak=$rss
  if [[ -n $rss ]] && (( rss > ${LIMIT_GB:-7} * 1048576 )); then
    kill $pid; echo "killed: RSS ${rss} kB over limit" >> $out/sage.log
  fi
  sleep 2
done
wait $pid; st=$?
echo "exit $st peak_rss_GB $(( peak / 1048576.0 ))" >> $out/sage.log
