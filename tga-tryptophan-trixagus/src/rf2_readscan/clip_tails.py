import sys, re
sam, out = sys.argv[1], sys.argv[2]
w = open(out, 'w')
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t'); flag = int(f[1])
    if flag & 0x904 or f[2] != 'CP157895.2': continue
    p = int(f[3]); ops = re.findall(r'(\d+)([MIDNSHP=X])', f[5])
    end = p + sum(int(n) for n, op in ops if op in 'MDN=X') - 1
    seq = f[9]
    if ops[-1][1] == 'S' and int(ops[-1][0]) >= 20 and 39440 <= end <= 39460:
        k = int(ops[-1][0]); w.write('>R_%s_end%d_mq%s\n%s\n' % (f[0], end, f[4], seq[-k-20:]))   # tail + 20 bp anchor
    if ops[0][1] == 'S' and int(ops[0][0]) >= 20 and 39750 <= p <= 39800:
        k = int(ops[0][0]); w.write('>L_%s_start%d_mq%s\n%s\n' % (f[0], p, f[4], seq[:k+20]))
