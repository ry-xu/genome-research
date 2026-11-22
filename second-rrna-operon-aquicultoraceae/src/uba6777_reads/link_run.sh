#!/bin/bash
# Physical-link test (ux/NOTES.md "Physical-link test" + amendment): fetch ERR599038 from SRA AWS, pass 1 map vs ref5 for names
# (any CEN/SSU hit or MAG|DKGF01000269.1), pass 2 pick pairs, delete the .sra, metaspades -t 3 -m 8.
# usage: link_run.sh WORKDIR   (WORKDIR holds ref5.fa)
set -u
W=$1; D=$(cd "$(dirname "$0")" && pwd)
cd "${W:?}" || exit 1
export TMPDIR="$W/tmpdir"; mkdir -p "$TMPDIR"
r=ERR599038
url=https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r
want=$(curl -sI --max-time 60 "$url" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{print $2}')
for i in $(seq 1 30); do
  have=$(stat -f %z "$r.sra" 2>/dev/null || echo 0)
  [ "$have" = "$want" ] && break
  curl -s -C - --speed-limit 200000 --speed-time 60 -o "$r.sra" "$url"
done
have=$(stat -f %z "$r.sra" 2>/dev/null || echo 0)
[ "$have" = "$want" ] || { echo "size $have != $want $(date +%H:%M)"; exit 1; }
echo "downloaded $want $(date +%H:%M)"
fastq-dump --split-spot -Z "./$r.sra" 2> fqd1.err | nice -n 10 minimap2 -t 3 -ax sr -k13 -w4 -m20 -s30 --sam-hit-only --secondary=no ref5.fa - 2> mm2.err \
 | awk -F'\t' '!/^@/ && ($3 ~ /^(SSU|CEN)\|/ || $3=="MAG|DKGF01000269.1"){print $1}' | sort -u > names.txt
echo "pass1 ${PIPESTATUS[*]} $(wc -l < names.txt) names $(date +%H:%M)"
tail -1 fqd1.err
fastq-dump --split-spot -Z "./$r.sra" 2> fqd2.err | python3 -I "$D/pick.py" names.txt r12.fq rs.fq
echo "pass2 $(( $(wc -l < r12.fq) / 8 )) pairs $(date +%H:%M)"
tail -1 fqd2.err
rm -f "${W:?}/ERR599038.sra"
nice -n 10 metaspades.py -t 3 -m 8 --12 r12.fq -o spades > spades.log 2>&1
echo "spades $? $(date +%H:%M)"
