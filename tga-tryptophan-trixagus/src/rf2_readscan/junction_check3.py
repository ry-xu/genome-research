import sys, re, statistics, collections
sam, outfa = sys.argv[1], sys.argv[2]
def rend(pos, cig):
    return pos + sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', cig)) - 1
fa = open(outfa, 'w'); tl = []; straddle = []; prop = collections.Counter()
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x904: continue
    if f[2].endswith('RF2_region'):
        fa.write('>%s_%s_%s_%s\n%s\n' % (f[0], f[2], f[3], f[5], f[9]))
        continue
    if f[2] != 'CP157895.2' or not (flag & 1) or (flag & 8) or f[6] != '=': continue
    pos, tlen = int(f[3]), int(f[8])
    if tlen <= 0: continue          # count each pair once (leftmost mate)
    if flag & 0x2: tl.append(tlen)
    s, e = pos, pos + tlen - 1      # fragment span on reference
    if s <= 39420 and e >= 39740:
        straddle.append((s, e, tlen, f[4], f[5], f[7], bool(flag & 2)))
tl.sort()
q = lambda p: tl[int(p * (len(tl) - 1))]
print('proper pairs on Tder:', len(tl), 'TLEN q01/q05/q50/q95/q99:', q(.01), q(.05), q(.5), q(.95), q(.99), 'max', tl[-1])
print('fragments spanning 39420..39740 (the depth hole):', len(straddle))
for r in straddle: print('  ', r)
