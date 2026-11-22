#!/bin/sh
# T1 gene trees: cmalign (match columns) -> drop columns with >50% gaps -> IQ-TREE 3 GTR+F+R4, UFBoot 1000, 1 thread
set -e
cd lab/work/silent-operon/B_origin
T=~/dev/tools/bin; R=lab/work/frozen-core/data/rfam
for k in 23S 16S; do
  cm=$R/RF02541.cm; [ $k = 16S ] && cm=$R/RF00177.cm
  $T/cmalign --cpu 1 --matchonly --outformat AFA -o out/t1_$k.afa $cm out/t1_${k}_in.fna > out/t1_${k}_cmalign.log
  uv run -q --with biopython python -I -c "
from Bio import AlignIO
a=AlignIO.read('out/t1_$k.afa','fasta'); n=len(a)
keep=[i for i in range(a.get_alignment_length()) if sum(r.seq[i] in '-.' for r in a)/n<=0.5]
with open('out/t1_$k.trim.fa','w') as fh:
    for r in a: fh.write('>'+r.id+'\n'+''.join(str(r.seq[i]) for i in keep).upper().replace('U','T')+'\n')
print('$k', n, 'seqs', len(keep), 'cols')"
  iqtree3 -s out/t1_$k.trim.fa -m GTR+F+R4 -B 1000 -T 1 --prefix out/t1_$k -redo -quiet
done
echo done > out/t1_tree.done
