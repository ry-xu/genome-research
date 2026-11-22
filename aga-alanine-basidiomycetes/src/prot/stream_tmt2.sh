#!/bin/zsh
# usage: stream_tmt2.sh SET [LANES]: one_tmt.sh over the set's 16 .raw files, LANES at a time (default 4, enforced by
# xargs -P; replaces stream_tmt.sh, whose lane throttle didn't hold and whose curl restarted from 0 every hour), then Sage + tally
set -u
P=lab/work/new-stops-euk; S=scratch
X=$1; L=${2:-4}; raw=$S/raw_tasahii_$X; mz=$S/mzml_tasahii_$X; mkdir -p $raw $mz
date
for i in {1..16}; do echo "http://download.iprox.org/IPX0002946000/IPX0002946001/BPI_20056_4plex_TMT_$X$i.raw"; done > $S/raw_tasahii_$X.txt
xargs -P $L -I{} $S/scripts/one_tmt.sh {} $raw $mz < $S/raw_tasahii_$X.txt
date
ls $mz/*.mzML | wc -l
cd $P; mkdir -p $S/sage_set$X
src/prot/run_sage.sh results/aga_ala/prot/sage_tmt.json results/aga_ala/prot/db.fasta $S/sage_set$X $mz/*.mzML
tail -1 $S/sage_set$X/sage.log
$S/upy src/prot/tally_aga.py $S/sage_set$X/results.sage.tsv results/aga_ala/prot/db.fasta > $S/sage_set$X/tally.txt
cat $S/sage_set$X/tally.txt
date
