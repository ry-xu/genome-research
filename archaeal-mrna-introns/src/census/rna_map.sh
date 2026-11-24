#!/bin/zsh
# One run against one unit's genome (PREREG_SURVEY.md). Only gapped reads, coverage and stats are kept.
# usage: src/survey_map.sh ASSEMBLY RUN [SPOTS=8000000] [THREADS=2]   -> out/runs/<ASSEMBLY>__<RUN>.{gapreads.tsv,cov.npz,stats.tsv,n,mm2.log}
set -o pipefail
cd /Users/ryan/lab/work/archaeal-intron-census
asm=$1; run=$2; spots=${3:-8000000}; t=${4:-2}
ref=data/ref/$asm/genomic.fna.gz
[[ -e data/ref/$asm/md5ok ]] || { echo "no verified ref for $asm"; exit 2; }
od=${SURVEY_OUT:-out/rna}; mkdir -p $od; pre=$od/${asm}__${run}
[[ -e $pre.done ]] && { echo "done already $pre"; exit 0; }
fastq-dump --stdout --split-spot --readids --skip-technical -X $spots $run 2> $pre.fqd.log \
  | tee >(awk 'END{print NR/4}' > $pre.n) \
  | minimap2 -ax splice:sr -un --eqx -t $t --sam-hit-only $ref - 2> $pre.mm2.log \
  | uv run --no-project --with numpy python -I /Users/ryan/lab/work/archaeal-introns/survey/src/gapcount.py $pre $ref \
  && echo "$(date '+%F %T')" > $pre.done
