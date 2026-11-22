#!/bin/zsh
# Download one Thermo .raw and convert it to centroided indexed mzML (ThermoRawFileParser via Homebrew dotnet,
# as in work/proteome-code-census/src/trfp.sh). usage: fetch_convert.sh <url> <rawdir> <mzmldir>
url=$1; raw=$2; mz=$3; f=${url:t}
[[ -s $mz/${f:r}.mzML ]] && { echo "have ${f:r}"; exit 0; }
curl -s --retry 3 --max-time 3600 -o $raw/$f "$url" || { echo "FAIL dl $f"; exit 1; }
export DOTNET_ROOT=/opt/homebrew/opt/dotnet/libexec DOTNET_ROLL_FORWARD=Major
dotnet ~/dev/tools/trfp/ThermoRawFileParser.dll -i $raw/$f -o $mz -f 2 -L 2 > $mz/${f:r}.trfp.log 2>&1
echo "$f $(stat -f %z $raw/$f) -> $(stat -f %z $mz/${f:r}.mzML 2>/dev/null) $(date +%H:%M)"
