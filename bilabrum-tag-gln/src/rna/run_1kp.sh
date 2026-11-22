#!/bin/bash
# Fetch 1KP runs that carry Bilabrum sp. RNA (Vierock et al. 2022 Sci Adv, Text S1) from the SRA AWS mirror, stream reads
# (fasterq-dump --split-spot to stdout, mates as single reads) into minimap2 (sr preset, PAF with long cs) against the
# Bilabrum genome, read codons at the reviewer's TAG sites. No FASTQ is written; each .sra is deleted after mapping.
# Run from results/bilabrum/rna_map in the lab tree (here: src/rna/ with data/review/bil.sites.tsv as ../review/).
set -u
G=lab/work/stop-sense-euk/data/dl_bil/GCA_965178125.1_puBilCfla1.1_genomic.fna
D=lab/work/stop-sense-euk/data/dl_1kp_sra
T="${TMPDIR:-/tmp}/fqd_tmp"
mkdir -p "$T"
one() {
  R=$1; SIZE=$2
  echo "$R start $(date +%H:%M)"
  mkdir -p "$D/$R"
  if [ "$(stat -f %z "$D/$R/$R.sra" 2>/dev/null)" != "$SIZE" ]; then
    curl -s -S --retry 3 -C - -o "$D/$R/$R.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$R/$R" || { echo "$R curl failed"; return 1; }
  fi
  [ "$(stat -f %z "$D/$R/$R.sra")" = "$SIZE" ] || { echo "$R size mismatch"; return 1; }
  echo "$R fetched $(date +%H:%M)"
  set -o pipefail
  fasterq-dump --split-spot -Z -e 2 -t "$T" "$D/$R/$R.sra" 2> "$R.fqd.log" \
    | minimap2 -x sr -c --cs=long -t 3 "$G" - > "$R.paf" 2> "$R.mm2.log" || { echo "$R map failed"; return 1; }
  set +o pipefail
  echo "$R mapped $(wc -l < "$R.paf") alignments $(date +%H:%M)"
  uv run --no-project python -I site_reads.py "$G" ../review/bil.sites.tsv "$R.paf" "$R"
}
one ERR2041128 1426207087 && rm -f "$D"/ERR2041128/ERR2041128.sra
one ERR2041020 "$(curl -s -I https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR2041020/ERR2041020 | awk 'tolower($1)=="content-length:"{print $2+0}')" \
  && rm -f "$D"/ERR2041020/ERR2041020.sra
one ERR3487419 "$(curl -s -I https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR3487419/ERR3487419 | awk 'tolower($1)=="content-length:"{print $2+0}')" \
  && rm -f "$D"/ERR3487419/ERR3487419.sra
echo "all done $(date +%H:%M)"
