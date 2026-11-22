#!/bin/bash
# Fetch NCBI taxonomy + nuccore records. No email/tool parameters.
set -u
D=lab/work/organelle-reading/checks/chor_wgs/dbcheck/data
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
for t in "Choristocarpus+tenellus" "Choristocarpus" "Discosporangium" "Discosporangiales" "Phaeophyceae" "Discosporangium+mesarthrocarpum" "Discosporangiaceae" "Choristocarpaceae"; do
  f=$(echo $t | tr '+' '_')
  curl -s --max-time 60 "$E/esearch.fcgi?db=taxonomy&term=${t}%5BScientific+Name%5D" -o $D/tax_search_$f.xml; sleep 0.5
done
for acc in ON552980 NC_085326; do
  curl -s --max-time 120 "$E/efetch.fcgi?db=nuccore&id=$acc&rettype=gbwithparts&retmode=text" -o $D/$acc.gb; sleep 0.5
  curl -s --max-time 120 "$E/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" -o $D/$acc.fa; sleep 0.5
done
curl -s --max-time 60 "$E/esearch.fcgi?db=nuccore&term=Discosporangium%5BOrganism%5D+AND+mitochondrion%5BFilter%5D&retmax=50" -o $D/disco_mito_search.xml
sleep 0.5
curl -s --max-time 60 "$E/esearch.fcgi?db=nuccore&term=Discosporangium%5BOrganism%5D+AND+mitochondrion+complete+genome&retmax=50" -o $D/disco_mito_search2.xml
