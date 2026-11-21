import sys, re
txt = open(sys.argv[1]).read()
blocks = re.findall(r'\{(.*?)\}', txt, re.S)
B = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
b2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
b3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
codes = {}
for bl in blocks:
    m = re.search(r'id (\d+)', bl); n = re.search(r'name "([^"]+)"', bl.replace('\n', ' '))
    a = re.search(r'ncbieaa\s+"([^"]+)"', bl)
    if m and a: codes[int(m.group(1))] = (n.group(1), a.group(1))
print('n tables', len(codes), sorted(codes))
std = codes[1][1]
for i, (n, aa) in sorted(codes.items()):
    pro = [B[k]+b2[k]+b3[k] for k in range(64) if aa[k] == 'P']
    diff = [(B[k]+b2[k]+b3[k], std[k], aa[k]) for k in range(64) if aa[k] != std[k]]
    print(i, n[:40], 'Pro:', pro, 'n_diff_from_std:', len(diff), diff[:8])
