#!/bin/bash
set -u
D=lab/work/organelle-reading/checks/chor_wgs/dbcheck/data
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
curl -s --max-time 60 "$E/efetch.fcgi?db=taxonomy&id=116065,116064,116066,376270,376269,376268,2870&retmode=xml" -o $D/taxonomy.xml; sleep 0.5
curl -s --max-time 60 "$E/esummary.fcgi?db=nuccore&id=886850108,781786540" -o $D/disco_mito_esummary.xml; sleep 0.5
for acc in CAXNAE020000396.1 CAXNAD020000071.1; do
  curl -s --max-time 120 "$E/efetch.fcgi?db=nuccore&id=$acc&rettype=gbwithparts&retmode=text" -o $D/$acc.gb; sleep 0.5
  curl -s --max-time 120 "$E/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" -o $D/$acc.fa; sleep 0.5
done
curl -s --max-time 60 "https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt" -o $D/gc.prt; sleep 0.5
curl -s --max-time 60 "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi" -o $D/wprintgc.html
