"""(d) count reads rescued by a conversion: set S kept if forward-mapped, set A if reverse-mapped; MAPQ>=20,
aligned >=60 nt, converted-space NM <= 2; then count original X>Y mismatches (+ strand) against the real reference.
usage: rescue_count.py ref.fa X Y S.sam A.sam"""
import sys, re
fa, X, Y, ssam, asam = sys.argv[1:6]
g = ''.join(l.strip() for l in open(fa) if not l.startswith('>')).upper()
C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(s): return ''.join(C.get(c, 'N') for c in reversed(s))
res = []
for path, want_rev in ((ssam, False), (asam, True)):
    for line in open(path):
        if line.startswith('@'): continue
        f = line.split('\t'); flag = int(f[1])
        if flag & 0x904 or bool(flag & 16) != want_rev or int(f[4]) < 20: continue
        nm = next(int(t[5:]) for t in f[11:] if t.startswith('NM:i:'))
        ops = re.findall(r'(\d+)([MIDSH])', f[5])
        alen = sum(int(n) for n, o in ops if o == 'M')
        if alen < 60 or nm > 2: continue
        orig = f[0].split('|')[1]
        q = rc(orig) if flag & 16 else orig   # read as aligned to + strand
        pos = int(f[3]) - 1; qi = 0; xy = other = 0
        for n, o in ops:
            n = int(n)
            if o in 'SI': qi += n
            elif o == 'D': pos += n
            elif o == 'M':
                for k in range(n):
                    r, b = g[pos + k], q[qi + k]
                    if r != b and b != 'N':
                        if r == X and b == Y: xy += 1
                        else: other += 1
                pos += n; qi += n
        res.append((xy, other))
n3 = sum(1 for xy, o in res if xy >= 3)
clean = sum(1 for xy, o in res if xy >= 3 and o <= 1)
print(f'{X}>{Y}: rescued reads {len(res)}; with >=3 {X}>{Y} mismatches {n3}; of those with <=1 other mismatch {clean}')
