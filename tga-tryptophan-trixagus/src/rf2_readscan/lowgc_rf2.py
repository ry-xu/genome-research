import sys, gzip, statistics
hits, reads = sys.argv[1], sys.argv[2]
best = {}
for l in open(hits):
    f = l.rstrip('\n').split('\t')
    if f[0] not in best or float(f[9]) > best[f[0]][2]: best[f[0]] = (f[1], float(f[2]), float(f[9]))
want = {q for q, v in best.items() if v[0] in ('sp|P07012|RF2_ECOLI', 'S_bostrichidophilus|0_75')}
out = open(sys.argv[3], 'w'); pid = {}
name = None
for l in gzip.open(reads, 'rt'):
    if l[0] == '>': name = l[1:].strip(); continue
    if name in want:
        s = l.strip(); g = (s.count('G') + s.count('C')) / len(s)
        if g < 0.20: out.write('>%s\n%s\n' % (name, s)); pid[name] = best[name]
p = [v[1] for v in pid.values()]
print('low-GC RF2-hit reads:', len(p), 'pident median %.1f range %.1f-%.1f' % (statistics.median(p), min(p), max(p)))
