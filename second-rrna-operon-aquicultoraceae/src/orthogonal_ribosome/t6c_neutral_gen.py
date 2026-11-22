"""T6-C neutral null sequences (post hoc, reported only): A. secundus div 16S with 15 % random substitutions, 5 seeds.
usage: uv run python src/t6c_neutral_gen.py OUTDIR [TS]   TS = transition share of substitutions (default: uniform, ~1/3)"""
import random, sys
od = sys.argv[1]; ts = float(sys.argv[2]) if len(sys.argv) > 2 else None
TR = {'A': 'G', 'G': 'A', 'C': 'U', 'U': 'C'}
s = open('out/t5/sec16.fa').read().split('>')[1:]
d = {x.split('\n')[0]: ''.join(x.split('\n')[1:]).upper().replace('T', 'U') for x in s}
ec = ''.join(l.split('\t')[1] for l in open('out/t6/ecoli16S_num.tsv'))
with open(f'{od}/in.fa', 'w') as o:
    o.write(f'>ecoli\n{ec}\n>sec_div\n{d["div"]}\n>sec_norm\n{d["norm"]}\n')
    for r in range(5):
        rng = random.Random(100 + r); q = list(d['div'])
        for i in range(len(q)):
            if rng.random() < 0.15:
                if ts is None: q[i] = rng.choice([b for b in 'ACGU' if b != q[i]])
                else: q[i] = TR[q[i]] if rng.random() < ts else rng.choice([b for b in 'ACGU' if b not in (q[i], TR[q[i]])])
        o.write(f'>SIM{r}_div\n{"".join(q)}\n')
