"""Build one row per rRNA locus from full-CM hits in windows: mature 16S/23S/5S coordinates in window (1-based),
dedupe by genomic position, classify div/norm by blastn identity to DASXVE01 div vs norm 16S/23S."""
import subprocess, tempfile, os
W = {}
n = None
for l in open('out/windows.fna'):
    if l[0] == '>': n = l[1:].strip(); W[n] = ''
    else: W[n] += l.strip()
def rd(f, emax=1e-10):
    d = {}
    for l in open(f):
        if l[0] == '#': continue
        x = l.split()
        if float(x[15]) > emax or x[9] != '+': continue
        d.setdefault(x[0], []).append(dict(mf=int(x[5]), mt=int(x[6]), a=int(x[7]), b=int(x[8]), sc=float(x[14])))
    return d
ssu, lsu, s5 = rd('out/w_ssu.tbl'), rd('out/w_lsu.tbl'), rd('out/w_5s.tbl', 1e-5)
CL = {'ssu': 1533, 'lsu': 2925, '5s': 119}
rows, seen = [], set()
for w, s in W.items():
    f = w.split('|'); g, c, st, lo, hi = f[0], f[1], f[3], int(f[7]), int(f[8])
    def gpos(p):  # window 1-based -> genomic 1-based
        return lo + p if st == '+' else hi - p + 1
    # 16S whose hit starts at the window's anchor (within 400 nt), else none (LSU windows)
    anchor = int(f[4]) + 1
    s16 = [h for h in ssu.get(w, []) if abs(h['a'] - anchor) < 400] if not f[5].startswith('LSU') else []
    s16 = max(s16, key=lambda h: h['sc']) if s16 else None
    if s16:
        l23 = [h for h in lsu.get(w, []) if 0 < h['a'] - s16['b'] < 2500]
    else:
        l23 = [h for h in lsu.get(w, []) if abs(h['a'] - anchor) < 400]
    l23 = min(l23, key=lambda h: h['a']) if l23 else None
    # 23S may be split by an insert: take the farthest downstream piece within 1500 nt of the first as the 3' end
    l23e = None
    if l23:
        pcs = [h for h in lsu.get(w, []) if h['a'] >= l23['a'] and h['a'] - l23['b'] < 1500]
        l23e = max(pcs, key=lambda h: h['b'])
    r5 = [h for h in s5.get(w, []) if l23e and 0 < h['a'] - l23e['b'] < 600]
    r5 = r5[0] if r5 else None
    key = (g, c, gpos(s16['a']) if s16 else None, gpos(l23['a']) if l23 else None)
    if key in seen: continue
    seen.add(key)
    rows.append(dict(win=w, genome=g, contig=c, wlen=len(s),
                     s16a=s16['a'] if s16 and s16['mf'] == 1 else None, s16b=s16['b'] if s16 and s16['mt'] == CL['ssu'] else None,
                     s16sc=s16['sc'] if s16 else None,
                     l23a=l23['a'] if l23 and l23['mf'] == 1 else None, l23b=l23e['b'] if l23e and l23e['mt'] == CL['lsu'] else None,
                     l23sc=l23['sc'] if l23 else None, l23pieces=len([1]) if l23 else 0,
                     r5a=r5['a'] if r5 else None, r5b=r5['b'] if r5 else None, r5sc=r5['sc'] if r5 else None,
                     r5full=(r5['mf'] == 1 and r5['mt'] >= CL['5s']) if r5 else None,
                     seq16=s[s16['a']-1:s16['b']] if s16 else '', seq23=s[l23['a']-1:l23e['b']] if l23 else ''))
# classify by blastn vs DASXVE01 refs
ref = {r['contig']: r for r in rows if r['genome'] == 'GCA_036504445.1'}
td = tempfile.mkdtemp()
with open(f'{td}/ref.fa', 'w') as o:
    for k, lab in (('DASXVE010000023.1', 'div'), ('DASXVE010000010.1', 'norm')):
        o.write(f'>{lab}_16\n{ref[k]["seq16"]}\n>{lab}_23\n{ref[k]["seq23"]}\n')
with open(f'{td}/q.fa', 'w') as o:
    for i, r in enumerate(rows):
        if r['seq16']: o.write(f'>{i}_16\n{r["seq16"]}\n')
        if r['seq23']: o.write(f'>{i}_23\n{r["seq23"]}\n')
out = subprocess.run(['blastn', '-query', f'{td}/q.fa', '-subject', f'{td}/ref.fa', '-outfmt', '6 qseqid sseqid pident length bitscore',
                      '-task', 'blastn'], capture_output=True, text=True).stdout
best = {}
for l in out.splitlines():
    q, sub, pid, ln, bs = l.split('\t')
    i, gene = q.split('_'); lab, sg = sub.split('_')
    if sg != gene: continue
    k = (int(i), gene, lab); best[k] = max(best.get(k, 0), float(bs))
for i, r in enumerate(rows):
    sc = {lab: sum(best.get((i, gn, lab), 0) for gn in ('16', '23')) for lab in ('div', 'norm')}
    if r['genome'] in ('ECOLI', 'TBIS'): r['cls'] = 'ctrl'
    else: r['cls'] = 'div' if sc['div'] > sc['norm'] else 'norm'
cols = ['genome', 'contig', 'cls', 'wlen', 's16a', 's16b', 's16sc', 'l23a', 'l23b', 'l23sc', 'r5a', 'r5b', 'r5sc', 'r5full', 'win']
with open('out/operons.tsv', 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in rows: o.write('\t'.join(str(r[c]) for c in cols) + '\n')
