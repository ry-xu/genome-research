#!/bin/bash
# Co-abundance mapping (ux/NOTES.md, "Co-abundance test"): for each run, fetch the .sra from SRA AWS (resumable),
# check its size, map all spots to ref4 with the read-test settings, then delete the .sra.
# usage: coab_run.sh WORKDIR RUN...
set -u
W=$1; shift
cd "${W:?}" || exit 1
export TMPDIR="$W/tmpdir"; mkdir -p "$TMPDIR"
for r in "$@"; do
  url=https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r
  want=$(curl -sI --max-time 60 "$url" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{print $2}')
  for i in $(seq 1 30); do
    have=$(stat -f %z "$r.sra" 2>/dev/null || echo 0)
    [ "$have" = "$want" ] && break
    curl -s -C - --speed-limit 200000 --speed-time 60 -o "$r.sra" "$url"
  done
  have=$(stat -f %z "$r.sra" 2>/dev/null || echo 0)
  [ "$have" = "$want" ] || { echo "$r size $have != $want, skipped $(date +%H:%M)"; continue; }
  echo "$r downloaded $(date +%H:%M)"
  fastq-dump --split-spot -Z "./$r.sra" 2> "$r.fqd.err" | minimap2 -t 2 -ax sr -k13 -w4 -m20 -s30 --sam-hit-only --secondary=no ref4.fa - > "$r.ref4.sam" 2> "$r.mm2.err"
  echo "$r mapped ${PIPESTATUS[*]} $(grep -vc '^@' "$r.ref4.sam") lines $(date +%H:%M)"
  tail -1 "$r.fqd.err"
  rm -f "${W:?}/${r:?}.sra"
done
