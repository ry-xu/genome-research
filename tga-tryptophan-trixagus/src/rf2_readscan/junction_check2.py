import sys, re, collections
sam, ref = sys.argv[1], sys.argv[2]
seqs = {}; name = None
for l in open(ref):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l)
seqs = {k: ''.join(v) for k, v in seqs.items()}
g = seqs['CP157895.2']
print('ref 39800-39840:', g[39799:39840])
print('N/ambig in 39333-39846:', sum(1 for c in g[39332:39846] if c.upper() not in 'ACGT'))
L = len(g); depth = [0]*(L+2)
gap_reads = []
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x904 or f[2] != 'CP157895.2': continue
    pos = int(f[3]); cig = f[5]
    p = pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=X':
            for i in range(p, p+n): depth[i] += 1
            p += n
        elif op in 'DN': p += n
    if pos <= 39806 and p - 1 >= 39835:
        sa = [t for t in f[11:] if t.startswith('SA:')]
        gap_reads.append((pos, cig, f[4], sa))
import statistics
d = depth[1:L+1]
print('genome depth median', statistics.median(d), 'mean', round(sum(d)/L, 1))
# 100-bp window medians across the genome, percentile of junction
win = [statistics.median(d[i:i+100]) for i in range(0, L-100, 100)]
jw = [statistics.median(d[i:i+100]) for i in range(39332, 39846-100, 100)]
print('junction 100bp-window medians:', jw)
ws = sorted(win)
for q in (0.01, 0.05, 0.1, 0.25, 0.5):
    print(' genome window q%.2f' % q, ws[int(q*len(ws))])
print('depth at 39700..39900 step 20:', [(i, depth[i]) for i in range(39700, 39901, 20)])
print('flanking: 38800-39300 median', statistics.median(d[38799:39300]), '; 39900-40400 median', statistics.median(d[39899:40400]))
for r in gap_reads: print('gap-spanning', r)
