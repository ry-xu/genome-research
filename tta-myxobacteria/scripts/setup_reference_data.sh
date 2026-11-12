#!/usr/bin/env bash
# Download the public reference data the analysis needs into data/ (all gitignored).
#   - GTDB R232 metadata -> data/gtdb_meta/{bac120,ar53}_metadata.tsv.gz
# The screen reuses code, the Swiss-Prot DIAMOND database (data/sprot.dmnd), the GTDB tree and the marker sets
# (data/reference/markers.faa.gz, data/markers.dmnd) of the TAG project, in the sibling folder ../tag-bacteria of
# this repository. Run its setup first (the original run used the Swiss-Prot release current on 2026-11-04; a newer
# release can shift individual hits slightly):
#   ../tag-bacteria/scripts/setup_reference_data.sh
# Genomes are fetched on demand from NCBI by src/fetch.py (cached in data/genomes/). Needs: curl, diamond.
set -euo pipefail
cd "$(dirname "$0")/.."
D=data
GTDB=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p $D/gtdb_meta
for f in bac120_metadata ar53_metadata; do
  [ -s $D/gtdb_meta/$f.tsv.gz ] || curl -fL -o $D/gtdb_meta/$f.tsv.gz $GTDB/${f}_r232.tsv.gz
done
echo "reference data ready in $D/"
