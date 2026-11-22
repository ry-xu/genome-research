"""C2 (RULES.md): how probe-like is each rRNA copy? For every position of a gene span, the best identity of a 50-nt
query window (centred there) inside any blastn HSP against E. coli rrnB / B. subtilis rrnA. Prints the fraction of
positions with window identity >= 0.90, >= 0.85 and >= 0.80, per copy, gene and probe-target organism.

usage: python3 -I src/probe_windows.py BLAST6_WITH_QSEQ_SSEQ QSTART QEND LABEL
(BLAST outfmt "6 sseqid qstart qend qseq sseq"; QSTART/QEND = gene span in query coordinates)
"""
import sys
W = 50
f, a, b, lab = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
best = {}
for l in open(f):
    s, qs, qe, qa, sa = l.rstrip('\n').split('\t')
    org = 'Ecoli' if 'NC_000913' in s else 'Bsub'
    qs, qe = int(qs), int(qe)
    step = 1 if qe >= qs else -1
    pos, q = [], qs  # query coordinate of each alignment column (None for query gaps)
    for c in qa:
        if c == '-': pos.append(None)
        else: pos.append(q); q += step
    match = [x == y and x != '-' for x, y in zip(qa, sa)]
    # windows of W query bases
    qidx = [i for i, p in enumerate(pos) if p is not None]
    for k in range(0, len(qidx) - W + 1):
        i0, i1 = qidx[k], qidx[k + W - 1]
        cols = i1 - i0 + 1
        ident = sum(match[i0:i1 + 1]) / max(cols, W)
        centre = pos[qidx[k + W // 2]]
        d = best.setdefault(org, {})
        if ident > d.get(centre, 0): d[centre] = ident
n = b - a + 1
for org in ('Ecoli', 'Bsub'):
    d = best.get(org, {})
    fr = [sum(1 for p in range(a, b + 1) if d.get(p, 0) >= t) / n for t in (0.90, 0.85, 0.80)]
    mx = max((d.get(p, 0) for p in range(a, b + 1)), default=0)
    print(f'{lab}\t{org}\t>=0.90 {fr[0]:.3f}\t>=0.85 {fr[1]:.3f}\t>=0.80 {fr[2]:.3f}\tmax50 {mx:.2f}')
