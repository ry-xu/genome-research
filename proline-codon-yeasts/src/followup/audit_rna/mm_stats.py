"""(d) read-level stats from a SAM stream (minimap2 --cs, all reads). Substitution types on the genome + strand
(all BK074543 genes are + strand, so genome T>C = transcript U>C). usage: samtools view x.bam | mm_stats.py"""
import sys, re
from collections import Counter
mq = Counter(); clip = Counter(); per = {}; nread = 0; strand = Counter()
hyper = Counter()   # reads with >=3 mismatches of one type
tot = Counter()
for line in sys.stdin:
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x900: continue
    nread += 1
    if flag & 4: mq['unmapped'] += 1; continue
    q = int(f[4]); b = 'q>=20' if q >= 20 else ('q0' if q == 0 else 'q1-19'); mq[b] += 1
    strand[('rev' if flag & 16 else 'fwd', b)] += 1
    cig = f[5]; s = sum(int(n) for n, op in re.findall(r'(\d+)([SH])', cig)); L = len(f[9])
    clip[(b, 'clip>=10' if s >= 10 else 'clip<10')] += 1
    cs = next((t[5:] for t in f[11:] if t.startswith('cs:Z:')), '')
    ty = Counter(f'{a.upper()}>{c.upper()}' for a, c in re.findall(r'\*([acgtn])([acgtn])', cs))
    for t, n in ty.items():
        tot[(b, t)] += n
        if n >= 3: hyper[(b, t)] += 1
print('reads', nread, dict(mq))
print('strand by MAPQ bin', dict(strand))
print('soft/hard clip >=10 nt by MAPQ bin', dict(clip))
types = sorted({t for _, t in tot})
for b in ('q>=20', 'q1-19', 'q0'):
    print(b, 'total mismatches by genome + strand type:', {t: tot[(b, t)] for t in types})
    print(b, 'reads with >=3 mismatches of one type:', {t: hyper[(b, t)] for t in types if hyper[(b, t)]})
