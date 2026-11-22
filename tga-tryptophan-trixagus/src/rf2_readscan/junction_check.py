import sys, re, collections
sam = sys.argv[1]
J0, J1 = 39333, 39846      # RF2 junction (flanking gene ends), 1-based
G0, G1 = 39816, 39825      # assembly gap
def ref_span(pos, cigar):
    L = sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', cigar))
    return pos, pos + L - 1
cnt = collections.Counter()
span_gap = []; overlap = []; rf2hits = collections.defaultdict(list)
per_ref = collections.Counter()
for line in open(sam):
    if line[0] == '@': continue
    f = line.rstrip('\n').split('\t')
    flag = int(f[1]); rname = f[2]; pos = int(f[3]); mapq = int(f[4]); cig = f[5]
    if flag & 4 or cig == '*': continue
    sec = bool(flag & 0x900)
    if not sec: per_ref[rname] += 1
    s, e = ref_span(pos, cig)
    nm = [t for t in f[11:] if t.startswith('NM:i:')]
    nm = int(nm[0][5:]) if nm else -1
    alen = sum(int(n) for n, op in re.findall(r'(\d+)([MI=X])', cig))
    if rname == 'CP157895.2' and not sec and e >= J0 and s <= J1:
        overlap.append((s, e, mapq, cig, nm, f[9]))
        if s <= G0 - 10 and e >= G1 + 10:
            span_gap.append((s, e, mapq, cig, nm, f[9]))
    if rname != 'CP157895.2':
        rf2hits[rname].append((s, e, mapq, cig, nm, sec, f[0]))
print('primary mapped per ref:', dict(per_ref))
print('reads overlapping junction', J0, J1, ':', len(overlap))
import statistics
# depth per 50bp window across junction +/- 500
for w0 in range(J0 - 500, J1 + 500, 50):
    d = sum(1 for s, e, *_ in overlap if s <= w0 + 49 and e >= w0) if J0 <= w0 <= J1 else None
print('reads spanning gap', G0, G1, '(>=10bp each side):', len(span_gap))
for r in span_gap[:30]:
    print('  span', r[0], r[1], 'mapq', r[2], r[3], 'NM', r[4])
print('reads with soft clip ending/starting within 30bp of gap:')
for s, e, mq, cig, nm, seq in overlap:
    if ('S' in cig) and (abs(e - G0) < 30 or abs(s - G1) < 30):
        print('  ', s, e, mq, cig, nm)
for k, v in rf2hits.items():
    print(k, 'alignments:', len(v), 'primary:', sum(1 for x in v if not x[5]))
    for x in v[:40]:
        print('  ', x)
