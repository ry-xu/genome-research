import sys, re
sam, out = sys.argv[1], sys.argv[2]
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
names = set(); recs = {}
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    recs.setdefault(f[0], []).append(f)
    if flag & 0x904 or f[2] != 'CP157895.2': continue
    pos = int(f[3]); end = pos + sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', f[5])) - 1
    if end >= 39380 and pos <= 39830: names.add(f[0])
w = open(out, 'w'); n = 0
for nm in names:
    for f in recs[nm]:
        flag = int(f[1])
        if flag & 0x900 or f[9] == '*': continue
        seq = f[9]
        if flag & 16: seq = rc(seq)      # back to read orientation
        mate = 'r1' if flag & 64 else 'r2'
        w.write('>%s_%s_%s:%s:%s:%s\n%s\n' % (nm, mate, f[2], f[3], f[5], f[4], seq)); n += 1
print('reads (both mates) of fragments touching 39380-39830:', len(names), 'records', n)
