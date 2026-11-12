#!/bin/zsh
# Fetch MGnify genomes (GFF3 with ##FASTA section; the FTP has no separate .fna) one at a time (EBI limit),
# and extract contigs to data/catalogues/mgnify/fna/<id>.fna.gz.
# usage: zsh src/catalogues/mgnify_fetch.sh ID URL [ID URL ...]
set -e
D=data/catalogues/mgnify
mkdir -p $D/gff $D/fna
while [ $# -ge 2 ]; do
  id=$1; url=$2; shift 2
  [ -s $D/fna/$id.fna.gz ] && continue
  url=${url//$'\r'/}
  curl -sS --fail --retry 3 --max-time 900 -o $D/gff/$id.gff.gz.part "$url" && mv $D/gff/$id.gff.gz.part $D/gff/$id.gff.gz
  gzcat $D/gff/$id.gff.gz | awk 'f{print} /^##FASTA/{f=1}' | gzip -c > $D/fna/$id.fna.gz.part && mv $D/fna/$id.fna.gz.part $D/fna/$id.fna.gz
  rm -f $D/gff/$id.gff.gz
  echo "$id $(gzcat $D/fna/$id.fna.gz | grep -c '^>') contigs"
done
