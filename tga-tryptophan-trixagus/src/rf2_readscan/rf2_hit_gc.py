import sys, gzip, collections, statistics
hits, reads, own = sys.argv[1], sys.argv[2], sys.argv[3]
RF2 = {'Tcar_carinifrons|0_130', 'Tduv_duvalii|0_133', 'Tmey|0_26', 'Aulonothroscus_AUS03|0_70',
       'S_bostrichidophilus|0_75', 'sp|P07012|RF2_ECOLI'}
best = {}
for l in open(hits):
    f = l.rstrip('\n').split('\t')
    if f[0] not in best or float(f[9]) > best[f[0]][2]: best[f[0]] = (f[1], float(f[2]), float(f[9]))
want = {q: v for q, v in best.items() if v[0] in RF2}
ownq = set(q for q, v in best.items() if v[0].startswith(own + '|'))
sbq = set(q for q, v in best.items() if v[0].startswith('S_bostrichidophilus|'))
gc = {}; owngc = []; sbgc = []
name = None
for l in gzip.open(reads, 'rt'):
    if l[0] == '>': name = l[1:].strip(); continue
    s = l.strip(); g = (s.count('G') + s.count('C')) / max(1, len(s))
    if name in want: gc[name] = g
    if name in ownq and len(owngc) < 20000: owngc.append(g)
    if name in sbq and len(sbgc) < 20000: sbgc.append(g)
print('GC of reads best-hitting own genome: median %.3f, 95th pct %.3f' % (statistics.median(owngc), sorted(owngc)[int(.95*len(owngc))]))
print('GC of reads best-hitting S_bost proteins: median %.3f' % statistics.median(sbgc))
by = collections.defaultdict(list)
for q, v in want.items(): by[v[0]].append(gc.get(q, -1))
for s, g in by.items():
    print('RF2 hits to %-28s n=%d GC median %.3f; reads with GC < 0.20: %d' % (s, len(g), statistics.median(g), sum(x < 0.20 for x in g)))
