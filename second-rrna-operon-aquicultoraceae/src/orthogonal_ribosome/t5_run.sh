#!/bin/bash
# T5 real run (RULES_T5.md): stream SRR11558543 through the frozen map/filter, then count. No FASTQ kept.
set -euo pipefail
SRA=data/sra/SRR11558543/SRR11558543.sra; O=out/t5/SRR11558543
fastq-dump --split-spot --stdout $SRA 2> $O.fqd.log \
  | minimap2 -x sr -t 1 -c --secondary=no out/t5/ref.fa - 2> $O.mm2.log \
  | awk -F'\t' '$11>=60 && $10/$11>=0.98' > $O.paf
uv run python src/t5_count.py $O.paf out/t5/class16.tsv > $O.counts.txt
