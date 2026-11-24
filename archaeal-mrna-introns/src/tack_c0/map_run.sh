#!/bin/zsh
# One run through the frozen path (RULES.md §2). usage: [IDX=out/rep] src/map_run.sh RUN SPOTS THREADS OUTDIR
#   IDX: index directory with regions.fa and events.tsv (default out)
#   -> OUTDIR/RUN.{counts.tsv,reads.tsv,n,mm2.log,fqd.log,bytes,done}
set -o pipefail
cd /Users/ryan/lab/work/tack-rna
export IDX=${IDX:-/Users/ryan/lab/work/tack-rna/out}
run=$1; spots=$2; t=${3:-2}; od=${4:-out/runs}; mkdir -p $od; pre=$od/$run
[[ -e $pre.done ]] && { echo "done already $pre"; exit 0; }
export TMPDIR=/Users/ryan/lab/work/tack-rna/out/tmp/$run; mkdir -p $TMPDIR
fastq-dump --stdout --split-spot --readids --skip-technical -X $spots $run 2> $pre.fqd.log \
  | tee >(awk 'END{print NR/4}' > $pre.n) \
  | minimap2 -ax splice:sr -un --eqx --secondary=no -t $t --sam-hit-only $IDX/regions.fa - 2> $pre.mm2.log \
  | tee >(gzip > $pre.sam.gz) \
  | python3 -I src/classify.py $pre 2> $pre.cls.log
rc=$?
vdb-dump --info $run 2>/dev/null | awk -v s=$spots '/^size/{gsub(",","",$3);z=$3} /^SEQ/{gsub(",","",$3);q=$3} END{f=s/q; if(f>1)f=1; printf "%d\n", z*f}' > $pre.bytes
rm -rf /Users/ryan/lab/work/tack-rna/out/tmp/$run
[[ $rc == 0 ]] && date '+%F %T' > $pre.done || rm -f $pre.counts.tsv
exit $rc
