import sys, re
sam, a, b, step = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
depth = {}; clipR = {}; clipL = {}
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x904 or f[2] != 'CP157895.2': continue
    p = int(f[3]); ops = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
    if p > b + 200 or p < a - 400: continue
    if ops[0][1] == 'S' and int(ops[0][0]) >= 20: clipL[p] = clipL.get(p, 0) + 1
    for n, op in ops:
        n = int(n)
        if op in 'M=X':
            for i in range(p, p+n): depth[i] = depth.get(i, 0) + 1
            p += n
        elif op in 'DN': p += n
    if ops[-1][1] == 'S' and int(ops[-1][0]) >= 20: clipR[p-1] = clipR.get(p-1, 0) + 1
print('depth:', [(i, depth.get(i, 0)) for i in range(a, b+1, step)])
print('reads with >=20bp right clip ending at:', sorted((k, v) for k, v in clipR.items() if a <= k <= b))
print('reads with >=20bp left clip starting at:', sorted((k, v) for k, v in clipL.items() if a <= k <= b))
