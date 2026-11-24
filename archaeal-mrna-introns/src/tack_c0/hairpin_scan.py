"""RT-skip-at-hairpin background for c0 (diagnostic, not in RULES; answers Sabik 12:34 (a) in L02's terms).
Altair's scan (agents/altair/kills/bhb_fold/scan.py) adapted to region-relative SAMs: in every index window, each W-nt stretch
with >= MINCOV primary reads (MAPQ >= 20) spanning it +-5 nt unbroken is folded on both strands (ViennaRNA MFE, stretch alone);
removal reads = reads with an N/D gap >= 15 nt lying inside the stretch +-10. Stretches are collapsed to loci (best dG per 2W),
and any locus overlapping a known event (cut or tRNA intron, +-10) is reported apart from the background.
Also folds the stretch +-40 nt (Vega 12:44: c0's region is -48.7) and writes all loci to OUT.tsv if OUT is given.
usage: agents/altair/rnaenv/bin/python -I src/hairpin_scan.py IDXDIR W DIR:RUN,RUN,... [OUT]"""
import sys, re, gzip, csv, numpy as np, RNA
from collections import defaultdict
idx, W, spec = sys.argv[1], int(sys.argv[2]), sys.argv[3]
out = sys.argv[4] if len(sys.argv) > 4 else None
d, runs = spec.split(':')
MINCOV, MAPQ, GAPMIN, PAD, FLANK = 30, 20, 15, 5, 10
CIG = re.compile(r'(\d+)([MIDNSHP=X])'); RC = str.maketrans('ACGT', 'TGCA')
seq, name = {}, None
for l in open(f'{idx}/regions.fa'):
    if l[0] == '>': name = l[1:].strip()
    else: seq[name] = l.strip().upper()
ev = defaultdict(list)
for r in csv.DictReader(open(f'{idx}/events.tsv'), delimiter='\t'):
    if r['kind'] in ('cut', 'trna'):
        ev[r['region']].append((int(r['rx']), int(r['ry']), r['id']))
blocks, gaps = defaultdict(list), defaultdict(list)
for run in runs.split(','):
    for line in gzip.open(f'{d}/{run}.sam.gz', 'rt'):
        if line[0] == '@': continue
        f = line.split('\t', 6)
        if int(f[1]) & 0x904 or int(f[4]) < MAPQ: continue
        rg = f[2]; pos = int(f[3]) - 1; bs = pos
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in 'M=XD' and not (op == 'D' and n >= GAPMIN): pos += n
            elif op in 'ND':
                blocks[rg].append((bs, pos)); gaps[rg].append((pos, n)); pos += n; bs = pos
        blocks[rg].append((bs, pos))
loci = []
S = W + 2 * PAD
for rg, bl in blocks.items():
    s0 = seq[rg]; N = len(s0)
    b = np.array(bl); b = b[(b[:, 1] - b[:, 0]) >= S]
    if not len(b): continue
    dd = np.zeros(N + 2, np.int64); np.add.at(dd, b[:, 0], 1); np.add.at(dd, b[:, 1] - S + 1, -1)
    cov = np.cumsum(dd)[:N]
    g = np.array(gaps[rg]).reshape(-1, 2)
    rows = []
    for x in np.nonzero(cov >= MINCOV)[0]:
        a = x + PAD
        if a + W + FLANK >= N: continue
        st = s0[a:a + W]
        if 'N' in st: continue
        dg = min(RNA.fold(st.replace('T', 'U'))[1], RNA.fold(st.translate(RC)[::-1].replace('T', 'U'))[1])
        j = int(np.sum((g[:, 0] >= a - FLANK) & (g[:, 0] + g[:, 1] <= a + W + FLANK))) if len(g) else 0
        rows.append((a, dg, int(cov[x]), j))
    taken = np.zeros(N, bool)
    for a, dg, c, j in sorted(rows, key=lambda r: r[1]):
        if taken[a]: continue
        taken[max(0, a - W):a + W] = True
        lab = next((i for x0, y0, i in ev[rg] if a < y0 + FLANK and a + W > x0 - FLANK), '')
        r40 = s0[max(0, a - 40):a + W + 40]
        d40 = min(RNA.fold(r40.replace('T', 'U'))[1], RNA.fold(r40.translate(RC)[::-1].replace('T', 'U'))[1]) if len(r40) == W + 80 else float('nan')
        loci.append((rg, a, dg, c, j, lab, d40))
bg = [l for l in loci if not l[5]]
print(f'W{W}: loci {len(loci)}, background (no known event) {len(bg)}')
for cut in (-3, -5, -7, -9):
    m = [l for l in bg if l[2] <= cut]
    print(f'  background dG<={cut}: loci {len(m)}, >=1 removal read {sum(l[4] >= 1 for l in m)}, >=3 {sum(l[4] >= 3 for l in m)}, '
          f'removal reads {sum(l[4] for l in m)} vs spanning {sum(l[3] for l in m)}')
m40 = [l for l in bg if l[6] == l[6]]
print(f'  background +-40 nt fold: loci with full context {len(m40)}; <=-30 {sum(l[6] <= -30 for l in m40)}, <=-40 {sum(l[6] <= -40 for l in m40)}, '
      f'<=-44 {sum(l[6] <= -44 for l in m40)} (removal reads {sum(l[4] for l in m40 if l[6] <= -44)} vs spanning {sum(l[3] for l in m40 if l[6] <= -44)}), '
      f'min {min(l[6] for l in m40):.1f}, median {sorted(l[6] for l in m40)[len(m40) // 2]:.1f}')
if out:
    with open(out, 'w') as f:
        f.write('region\tstart\tdG\tspanning\tremoval\tevent\tdG40\n')
        for l in loci: f.write('\t'.join(map(str, l)) + '\n')
print('  loci at known events (id, dG, spanning, removal reads):')
for l in sorted((l for l in loci if l[5]), key=lambda l: l[5]):
    print(f'    {l[5]}\t{l[2]:.1f}\t{l[3]}\t{l[4]}')
print('  strongest background loci (region, start, dG, spanning, removal):')
for l in sorted(bg, key=lambda l: l[2])[:8]:
    print(f'    {l[0]}\t{l[1]}\t{l[2]:.1f}\t{l[3]}\t{l[4]}')
