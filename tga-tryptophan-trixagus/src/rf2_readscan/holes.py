import sys, re, statistics
sam, ref = sys.argv[1], sys.argv[2]
seqs = {}; name = None
for l in open(ref):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l)
g = ''.join(seqs['CP157895.2']); L = len(g); depth = [0]*(L+2)
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x904 or f[2] != 'CP157895.2': continue
    p = int(f[3])
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
        n = int(n)
        if op in 'M=X':
            for i in range(p, p+n): depth[i] += 1
            p += n
        elif op in 'DN': p += n
# runs of depth <= 2 of length >= 150
runs = []; i = 1
while i <= L:
    if depth[i] <= 2:
        j = i
        while j <= L and depth[j] <= 2: j += 1
        if j - i >= 150: runs.append((i, j-1))
        i = j
    else: i += 1
gc = lambda s: sum(c in 'GCgc' for c in s) / max(1, sum(c in 'ACGTacgt' for c in s))
print('genome GC %.3f' % gc(g), 'length', L)
print('low-depth runs (depth<=2, >=150 bp):', len(runs), 'total bp', sum(b-a+1 for a, b in runs))
for a, b in runs:
    s = g[a-1:b]
    print('  %d-%d len %d GC %.3f N %d %s' % (a, b, b-a+1, gc(s), sum(c not in 'ACGTacgt' for c in s), '<-- RF2 junction' if a <= 39846 and b >= 39333 else ''))
