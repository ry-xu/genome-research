#!/bin/bash
# Read-level check of an SRL call in a MAG. usage: src/srl_readcheck.sh ACC RUN [RUN...]
# Per run: download the .sra from the SRA AWS bucket into its own directory, stream reads (fastq-dump on the local file)
# into minimap2 -x sr against the MAG's census 23S, keep the PAF, delete the .sra. Then tally the dodecamer
# (src/srl_readcheck.py). Outputs in results/lsu/bac/reads/runs/ACC/.
set -euo pipefail
W=${FROZEN_CORE:-frozen-core}
ACC=$1; shift
D=$W/results/lsu/bac/reads/runs/$ACC; mkdir -p $D
awk -v a=">$ACC" '/^>/{p=($1==a)} p' $W/data/lsu/bac_lsu_all.fna > $D/ref.fa
[ -s $D/ref.fa ] || { echo "no census 23S for $ACC"; exit 1; }
for RUN in "$@"; do
  [ -s $D/$RUN.paf.done ] && continue
  DL=$W/scratch/sra_dl/$RUN; mkdir -p $DL
  URL=https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN
  SIZE=$(curl -sI $URL | awk 'tolower($1)=="content-length:"{print $2+0}')
  echo "$(date +%H:%M) $RUN size $SIZE"
  for i in 1 2 3 4 5 6; do
    curl -s -C - --speed-limit 200000 --speed-time 60 -o $DL/$RUN.sra $URL || true
    [ "$(stat -f %z $DL/$RUN.sra 2>/dev/null || echo 0)" = "$SIZE" ] && break
  done
  [ "$(stat -f %z $DL/$RUN.sra)" = "$SIZE" ] || { echo "download incomplete $RUN"; exit 1; }
  echo "$(date +%H:%M) $RUN downloaded"
  TMPDIR=$DL fastq-dump --stdout --skip-technical --split-spot --readids $DL/$RUN.sra 2> $D/$RUN.fqd.log \
    | minimap2 -x sr -c --cs=long -t 2 --secondary=no $D/ref.fa - 2> $D/$RUN.mm2.log > $D/$RUN.paf
  echo "$(date +%H:%M) $RUN mapped $(wc -l < $D/$RUN.paf) alignments"
  rm -f $DL/$RUN.sra; rmdir $DL
  echo ok > $D/$RUN.paf.done
done
cat $D/[SED]RR*.paf > $D/all.paf.tmp && mv $D/all.paf.tmp $D/all.paf
cd $W && uv run -q python -I src/srl_readcheck.py $ACC $D/all.paf $D/tally.tsv
