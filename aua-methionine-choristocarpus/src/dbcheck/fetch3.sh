#!/bin/bash
set -u
D=lab/work/organelle-reading/checks/chor_wgs/dbcheck/data
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
curl -s --max-time 60 "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=116065&lvl=3&lin=f&keep=1&srchmode=1&unlock" -o $D/wwwtax_116065.html; sleep 1
curl -s --max-time 60 "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=2870&lvl=3&lin=f&keep=1&srchmode=1&unlock" -o $D/wwwtax_2870.html; sleep 1
for t in Xanthophyceae Labyrinthulomycetes Sphaerosorus Schizocladiophyceae Bacillariophyta Oomycota; do
  curl -s --max-time 60 "$E/esearch.fcgi?db=taxonomy&term=${t}%5BScientific+Name%5D" -o $D/tax_search_$t.xml; sleep 0.5
done
