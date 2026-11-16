#!/bin/bash
# usage: run_genomes.sh SET.tsv OUT.tsv  -- one row per genome (whole-genome CGG stats, pyrodigal single mode)
set -u
SET=$1; OUT=$2
CACHE_PY='import sys; from codon_atlas.cache import genome_path; print(genome_path(sys.argv[1]))'
echo -e "id\tlen\tgc\tn_cds\tn_codons\tO\tn_fam\tthird_share\tE\tr\tskip_sites\ttrna_ccg\ttrna_tcg\ttrna_all\tO_CGA\tE_CGA\tO_CGT\tE_CGT\tO_CGC\tE_CGC" > $OUT
for gca in $(cut -f2 $SET | tail -n +2); do
  f=$(uv run python -I -c "$CACHE_PY" $gca)
  gzip -dc $f > /tmp/cw_$$.fna 2>/dev/null || { echo "skip $gca"; continue; }
  uv run --no-project --with pyrodigal --with pyaragorn python -I phages_work/src/codon_stats.py - /tmp/cw_$$.fna $OUT --whole $gca --mode single --minlen 1000 --append
done
rm -f /tmp/cw_$$.fna
