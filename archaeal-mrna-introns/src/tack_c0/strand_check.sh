#!/bin/bash
# Strand check for S2 (diagnostic, not in RULES): orientation of reads relative to the feature strand.
# usage: src/strand_check.sh DIR RUN,RUN,...   -> prints counts; dUTP = mate 1 antisense to the feature.
# Part 1: every event read in reads.tsv vs the event strand (cut: 0 fwd / 1 rev; tRNA: 1 fwd / -1 rev).
# Part 2: primary SAM reads in M. verstraetei windows vs the GenBank CDS strand (meta/CP186953.1.ft.txt).
cd /Users/ryan/lab/work/tack-rna
d=$1; IFS=, read -ra runs <<< "$2"
tmp=$(mktemp -d)
for r in "${runs[@]}"; do tail -n +2 $d/$r.reads.tsv; done > $tmp/reads.tsv
echo "part 1: kind  event  orientation  class  reads"
awk -F'\t' 'NR==FNR{if(FNR>1){st[$9]=$11; kd[$9]=$2}; next}
  {m=substr($3,length($3)); r=int($6/16)%2; r1rev=(m=="1")?r:1-r
   fwd=(kd[$1]=="trna")?(st[$1]==1):(st[$1]==0)
   print kd[$1], ($1=="GCA_049570475.1:c0"?"c0":"other"), (r1rev==fwd)?"sense":"antisense", $2}' out/events.tsv $tmp/reads.tsv | sort | uniq -c
awk -F'\t' '$1 ~ /^[<>]?[0-9]+$/ && $2 ~ /^[<>]?[0-9]+$/ && $3=="CDS"{a=$1;b=$2;gsub(/[<>]/,"",a);gsub(/[<>]/,"",b); print (a<b?a:b), (a<b?b:a), (a<b?"+":"-")}' meta/CP186953.1.ft.txt > $tmp/cds.txt
for r in "${runs[@]}"; do gzcat $d/$r.sam.gz | awk -F'\t' '!/^@/ && $3 ~ /^GCA_049570475.1\|CP186953.1\|/ && int($2/4)%2==0 && int($2/256)%2==0 && int($2/2048)%2==0 {split($3,w,"|"); print w[3]+$4+50, ($1 ~ /\.1$/)?1:2, int($2/16)%2}'; done > $tmp/mv.txt
echo "part 2: M. verstraetei window reads inside one CDS (read start + 50): dUTP-like vs opposite"
awk 'NR==FNR{lo[NR]=$1;hi[NR]=$2;s[NR]=$3;n=NR;next} {g="";for(i=1;i<=n;i++) if($1>=lo[i]&&$1<=hi[i]){g=s[i];break} if(g=="")next; r1rev=($2==1)?$3:1-$3; print (r1rev==(g=="+"))?"dUTP-like":"opposite"}' $tmp/cds.txt $tmp/mv.txt | sort | uniq -c
rm -r "$tmp"
