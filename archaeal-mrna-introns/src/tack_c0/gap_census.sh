#!/bin/bash
# Library-level RT-skip control (diagnostic, not in RULES): every deletion/N >= 10 nt in primary reads of the given runs,
# grouped by (window, ref start, length); each gap labelled by the event it matches (exact rx..ry +- k+1 and length L) or 'none'.
# usage: src/gap_census.sh DIR IDXDIR RUN,RUN,...   -> prints gap groups with >= 2 reads, and totals by label
cd /Users/ryan/lab/work/tack-rna
d=$1; idx=$2; IFS=, read -ra runs <<< "$3"
for r in "${runs[@]}"; do gzcat $d/$r.sam.gz; done | awk -F'\t' -v OFS='\t' '
NR==FNR { if (FNR>1) { n[$1]++; ex[$1,n[$1]]=$5; ey[$1,n[$1]]=$6; eL[$1,n[$1]]=$7; ek[$1,n[$1]]=$8; eid[$1,n[$1]]=$9 } next }
/^@/ { next }
int($2/4)%2==0 && int($2/256)%2==0 && int($2/2048)%2==0 {
  p=$4-1; c=$6
  while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
    len=substr(c,1,RLENGTH-1)+0; op=substr(c,RLENGTH,1); c=substr(c,RLENGTH+1)
    if (op ~ /[M=XDN]/) {
      if (op ~ /[DN]/ && len >= 10) {
        lab="none"
        for (i=1;i<=n[$3];i++) if (len==eL[$3,i] && p>=ex[$3,i]-ek[$3,i]-1 && p<=ex[$3,i]+ek[$3,i]+1) { lab=eid[$3,i]; break }
        print $3, p, len, lab
      }
      p+=len
    }
  }
}' $idx/events.tsv - | sort | uniq -c | sort -k1,1nr > /tmp/gc_$$.txt
echo "gap groups with >= 2 reads (reads window start len label):"; awk '$1>=2' /tmp/gc_$$.txt | head -40
echo "totals (reads, distinct gap positions) by label kind:"
awk '{k=($5=="none")?"none":($5 ~ /:t/)?"tRNA intron":($5 ~ /[mp]$/)?"decoy":"cut"; r[k]+=$1; g[k]++} END{for(k in r) print k, r[k], g[k]}' /tmp/gc_$$.txt
rm -f /tmp/gc_$$.txt
