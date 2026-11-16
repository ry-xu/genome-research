#!/bin/bash
cd phages_work/kill1/t3
waitload() { while true; do l=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,); if awk -v l="$l" 'BEGIN{exit !(l<=10)}'; then break; fi; echo "$(date +%H:%M) load $l, waiting"; sleep 60; done; }
D=/opt/homebrew/bin/diamond
F="qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq_gapped sseq_gapped"
for db in R all; do
  waitload; echo "$(date +%H:%M) makedb $db"
  /usr/bin/time -l $D makedb -p 1 --in targets_$db.faa -d targets_$db > makedb_$db.log 2>&1 || { echo fail makedb; exit 1; }
  waitload; echo "$(date +%H:%M) blastp $db"
  /usr/bin/time -l $D blastp -p 1 --sensitive -b 0.4 -c 1 --max-hsps 0 -k 200 -e 1e-5 -q queries.faa -d targets_$db -o hits_$db.tsv -f 6 $F > blastp_$db.log 2>&1 || { echo fail blastp; exit 1; }
  grep -E 'maximum resident|real' blastp_$db.log
done
echo "$(date +%H:%M) all done"
