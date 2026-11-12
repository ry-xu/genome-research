#!/bin/bash
# usage: fetch_list.sh <list.tsv: id<TAB>url> <dest_dir> <suffix e.g. .fa.gz> [parallel=2]
# Downloads each URL to <dest_dir>/<id><suffix>; skips files already present; logs failures.
list="$1"; dest="$2"; suf="$3"; par="${4:-2}"
cd "$dest" || exit 1
fetch_one() {
  id="$1"; url="$2"; suf="$3"
  out="${id}${suf}"
  [ -s "$out" ] && return 0
  if curl -sS -L -f --retry 3 --retry-delay 5 -m 600 -o "${out}.part" "$url"; then
    mv "${out}.part" "$out"
  else
    rm -f "${out}.part"; echo -e "FAIL\t$id\t$url" >> fetch_failures.log
  fi
}
export -f fetch_one
tr '\t' ' ' < "$list" | xargs -P "$par" -L 1 bash -c 'fetch_one "$0" "$1" "'"$suf"'"'
echo "done $(date)"
