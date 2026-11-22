"""Group 16S/23S/5S HMM hits into operons and extract a window per 16S hit (16S start-1500 .. +7500, rRNA sense)."""
import sys
def rd(f, emax):
    out = []
    for l in open(f):
        if l[0] == '#': continue
        x = l.split()
        a, b = int(x[7]), int(x[8])
        if float(x[15]) > emax: continue
        out.append(dict(c=x[0], mf=int(x[5]), mt=int(x[6]), a=a, b=b, st=x[9], sc=float(x[14])))
    return out
seqs, n = {}, None
for l in open('data/all.fna'):
    if l[0] == '>': n = l[1:].split()[0]; seqs[n] = []
    else: seqs[n].append(l.strip().upper())
seqs = {k: ''.join(v) for k, v in seqs.items()}
rc = lambda s: s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
ssu = rd('out/ssu_hmm.tbl', 1e-10)
w = open('out/windows.fna', 'w')
for h in ssu:
    s = seqs[h['c']]; L = len(s)
    if h['st'] == '+':
        lo, hi = max(0, h['a'] - 1 - 1500), min(L, h['a'] - 1 + 7500); sub = s[lo:hi]; off = h['a'] - 1 - lo
    else:
        hi, lo = min(L, h['a'] + 1500), max(0, h['a'] - 7500); sub = rc(s[lo:hi]); off = hi - h['a']
    # header: contig|16Sstart|strand|offset_of_16S_hit_in_window|mf|mt|lo|hi
    w.write(f">{h['c']}|{h['a']}|{h['st']}|{off}|{h['mf']}|{h['mt']}|{lo}|{hi}\n{sub}\n")
# 23S starts (mf<=10) not inside any 16S window: window 23S start-3000 .. +4500
cov = {}
for h in ssu:
    cov.setdefault(h['c'], []).append(h)
for h in rd('out/lsu_hmm.tbl', 1e-10):
    if h['mf'] > 10: continue
    near = any(abs(g['a'] - h['a']) < 6000 and g['st'] == h['st'] for g in cov.get(h['c'], []))
    if near: continue
    s = seqs[h['c']]; L = len(s)
    if h['st'] == '+':
        lo, hi = max(0, h['a'] - 1 - 3000), min(L, h['a'] - 1 + 4500); sub = s[lo:hi]; off = h['a'] - 1 - lo
    else:
        hi, lo = min(L, h['a'] + 3000), max(0, h['a'] - 4500); sub = rc(s[lo:hi]); off = hi - h['a']
    w.write(f">{h['c']}|{h['a']}|{h['st']}|{off}|LSU{h['mf']}|{h['mt']}|{lo}|{hi}\n{sub}\n")
w.close()
