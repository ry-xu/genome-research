import sys, collections, statistics
hits, own, db = sys.argv[1], sys.argv[2], sys.argv[3]
plen = {}; name = None
for l in open(db):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; plen[name] = 0
    else: plen[name] += len(l)
RF2 = {'Tcar_carinifrons|0_130', 'Tduv_duvalii|0_133', 'Tmey|0_26', 'Aulonothroscus_AUS03|0_70',
       'S_bostrichidophilus|0_75', 'sp|P07012|RF2_ECOLI'}
best = {}
for l in open(hits):
    f = l.rstrip('\n').split('\t')
    q, s, pid, bits = f[0], f[1], float(f[2]), float(f[9])
    if q not in best or bits > best[q][2]: best[q] = (s, pid, bits, int(f[3]))
per_prot = collections.Counter(v[0] for v in best.values())
per_genome = collections.Counter(v[0].split('|')[0] for v in best.values())
print('reads with a hit:', len(best), '; best-hit genome:', dict(per_genome.most_common()))
ownp = [p for p in plen if p.startswith(own + '|')]
dens = sorted(per_prot.get(p, 0) / plen[p] for p in ownp if plen[p] >= 100)
med = statistics.median(dens)
print('own proteins (>=100 aa):', len(dens), 'median reads per aa %.3f -> expected reads for a 345-aa RF2: %.0f; 10th pct %.0f' % (med, med * 345, dens[len(dens)//10] * 345))
r2 = [(q, v) for q, v in best.items() if v[0] in RF2]
print('reads whose best hit is an RF2:', len(r2))
by = collections.defaultdict(list)
for q, v in r2: by[v[0]].append(v[1])
for s, p in by.items():
    print('  %-28s n=%d pident median %.1f range %.1f-%.1f' % (s, len(p), statistics.median(p), min(p), max(p)))
hi = [v for q, v in r2 if v[1] >= 60 and v[0] != 'sp|P07012|RF2_ECOLI']
print('  of these, >=60%% identity to a Shikimatogenerans RF2: %d' % len(hi))
