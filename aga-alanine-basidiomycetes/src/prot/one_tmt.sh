#!/bin/zsh
# usage: one_tmt.sh URL RAWDIR MZDIR: resumable download of one .raw (curl -C -, no time cap; a stall below 2 KB/s for
# 120 s restarts the attempt from where it stopped), checked against Content-Length; convert to mzML; delete the raw
u=$1; raw=$2; mz=$3; f=${u:t}; n=
[[ -s $mz/${f:r}.mzML ]] && { echo "have ${f:r}"; exit 0; }
for try in {1..40}; do
  [[ -z $n ]] && n=$(curl -sI --max-time 60 $u | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  have=$(stat -f %z $raw/$f 2>/dev/null || echo 0)
  [[ -n $n && $have == $n ]] && break
  curl -s -C - --speed-limit 2000 --speed-time 120 -o $raw/$f $u
  sleep 10
done
have=$(stat -f %z $raw/$f 2>/dev/null || echo 0)
[[ -n $n && $have == $n ]] || { echo "FAIL dl $f $have/$n $(date +%H:%M)"; exit 1; }
export DOTNET_ROOT=/opt/homebrew/opt/dotnet/libexec DOTNET_ROLL_FORWARD=Major
dotnet ~/dev/tools/trfp/ThermoRawFileParser.dll -i $raw/$f -o $mz -f 2 -L 2 > $mz/${f:r}.trfp.log 2>&1
echo "$f $have -> $(stat -f %z $mz/${f:r}.mzML 2>/dev/null) $(date +%H:%M)"
[[ -s $mz/${f:r}.mzML ]] && rm -f $raw/$f
exit 0
