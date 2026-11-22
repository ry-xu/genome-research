#!/bin/bash
# Fetch genome fna for each acc in results/lsu/bac/cca/aq/genomes.tsv (GB_/RS_ prefixes) into aq/fna/. Skips existing.
set -uo pipefail
D=data/external/aq; mkdir -p $D/fna
tail -n +2 $D/genomes.tsv | cut -f1 | while read a; do
  acc=${a#*_}; [ -s $D/fna/$acc.fna ] && continue
  pre=${acc:0:3}; n=${acc:4:9}; d=$pre/${n:0:3}/${n:3:3}/${n:6:3}
  base=$(curl -s --max-time 60 https://ftp.ncbi.nlm.nih.gov/genomes/all/$d/ | grep -o "${acc}_[^/\"]*" | head -1)
  [ -n "$base" ] || { echo "nodir $acc"; continue; }
  curl -s --retry 3 --max-time 300 https://ftp.ncbi.nlm.nih.gov/genomes/all/$d/$base/${base}_genomic.fna.gz | gunzip -c > $D/fna/$acc.fna.tmp && mv $D/fna/$acc.fna.tmp $D/fna/$acc.fna || echo "fail $acc"
done
echo "$(ls $D/fna/*.fna | wc -l) genomes"
