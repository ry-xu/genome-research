"""Read-level check of the DASXVE01 divergent operon from minimap2 PAF (--cs=long) of SRR11794808.
usage: python3 -I pile.py PAF REF.fa
Prints: (1) read depth (id>=0.97, mapq>=20) and mate-pair fragment coverage in 250-bp windows over DASXVE010000023.1:1-9000;
(2) clipped-alignment ends (>=15 nt unaligned read end) per window vs genome-wide rate; (3) base calls at PTC coords from all
alignments (no filters) with their identities."""
import collections, re, sys
paf, ref = sys.argv[1:3]
C = "DASXVE010000023.1"
ptc = {6447: "2062", 6625: "2251", 6626: "2252", 6627: "2253", 6825: "2449", 6827: "2451", 6828: "2452", 6829: "2453",
       6879: "2503", 6880: "2504", 6958: "2583", 6959: "2584"}
depth = collections.Counter(); clip = collections.Counter(); calls = collections.defaultdict(collections.Counter)
mates = collections.defaultdict(list); clip_all = 0; aln_bases = 0; tot_len = 0
for l in open(ref):
    if not l.startswith(">"): tot_len += len(l.strip())
lowid = collections.Counter()
for line in open(paf):
    f = line.rstrip("\n").split("\t")
    qn, ql, qs, qe, st, tn, tl, ts, te, nm, bl, mq = f[:12]
    ql, qs, qe, ts, te, nm, bl, mq = map(int, (ql, qs, qe, ts, te, nm, bl, mq))
    ident = nm / bl
    clipped = qs >= 15 or ql - qe >= 15
    if mq >= 20: clip_all += clipped; aln_bases += 1
    if tn != C: continue
    if ts > 9000: continue
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    cur = ts
    for op, s in re.findall(r"([=*+\-~:])([A-Za-z0-9]+)", cs):
        if op == "=":
            for k, b in enumerate(s.upper()):
                p = cur + k + 1
                if p in ptc: calls[p][(b, round(ident, 3) >= 0.97)] += 1
            cur += len(s)
        elif op == "*":
            if cur + 1 in ptc: calls[cur + 1][(s[1].upper(), round(ident, 3) >= 0.97)] += 1
            cur += 1
        elif op == "-":
            cur += len(s)
    if mq >= 20 and ident >= 0.97:
        for p in range(ts // 250, (te - 1) // 250 + 1): pass
        for p in range(ts + 1, te + 1): depth[p] += 1
        if clipped:
            # position of the clipped end on the reference
            if qs >= 15: clip[(ts if st == "+" else te) // 250] += 1
            if ql - qe >= 15: clip[(te if st == "+" else ts) // 250] += 1
        mates[qn.rsplit(".", 1)[0]].append((ts, te))
    elif mq >= 20:
        lowid[ts // 250] += 1
frag = collections.Counter(); ins = []
for k, v in mates.items():
    if len(v) == 2:
        a = min(x[0] for x in v); b = max(x[1] for x in v)
        if b - a < 1500:
            ins.append(b - a)
            for p in range(a + 1, b + 1): frag[p] += 1
print(f"genome-wide clipped-alignment rate (mapq>=20): {clip_all}/{aln_bases} = {clip_all/aln_bases:.4f}")
ins.sort(); print(f"fragments on contig window: {len(ins)}, median insert {ins[len(ins)//2]}")
print("window\tmean_read_depth\tmin_read_depth\tmin_fragment_cov\tclipped_ends\tlowid_mapq20_alns")
for w in range(0, 36):
    a, b = w * 250 + 1, w * 250 + 250
    d = [depth[p] for p in range(a, b + 1)]; fr = [frag[p] for p in range(a, b + 1)]
    print(f"{a}-{b}\t{sum(d)/len(d):.1f}\t{min(d)}\t{min(fr)}\t{clip[w]}\t{lowid[w]}")
print("PTC base calls (base, id>=0.97): count")
g = "".join(l.strip() for l in open(ref) if True)  # unused
for p, e in ptc.items():
    print(p, e, dict(calls[p]))
