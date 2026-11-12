#!/bin/zsh
cd "$(dirname "$0")"
while IFS=$'\t' read lab acc role genus rf2; do
  [ -s genes/$lab.faa ] && continue
  uv run --no-project --with pyrodigal python -I call_genes.py ../../../data/genomes/$acc.fna.gz $lab genes && echo "done $lab $(date +%T)"
done < panel.tsv
