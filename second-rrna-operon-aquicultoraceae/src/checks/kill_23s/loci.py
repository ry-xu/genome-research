"""Merge blastn hits of 4 23S queries (dx_div23, dx_norm23, kmv_div23, cad_div23) on the 65 Aquicultoria genomes
into loci (subject intervals), then label each locus by which query class (div vs norm) gives the higher identity-weighted
bitscore over the locus. Prints per locus: genome, contig, start, end, span, best div pident/aligned, best norm pident/aligned."""
import collections, sys
tax = {}
for l in open("../aq/genomes.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "acc": tax[f[0][3:]] = (f[1], f[2])
hits = collections.defaultdict(list)
for l in open(sys.argv[1]):
    q, s, pid, ln, qs, qe, ss, se, ev, bs = l.split("\t")
    g, c = s.split("|")
    a, b = sorted((int(ss), int(se)))
    hits[(g, c)].append((a, b, q, float(pid), int(ln), float(bs)))
rows = []
for (g, c), hs in hits.items():
    hs.sort()
    groups = []
    for h in hs:
        if groups and h[0] <= groups[-1][1] + 300:
            groups[-1][1] = max(groups[-1][1], h[1]); groups[-1][2].append(h)
        else:
            groups.append([h[0], h[1], [h]])
    for a, b, hh in groups:
        sc = collections.defaultdict(lambda: [0.0, 0, 0.0])
        for _, _, q, pid, ln, bs in hh:
            sc[q][0] += bs; sc[q][1] += ln; sc[q][2] += pid * ln
        def best(cls):
            qs = [q for q in sc if (("div" in q) == (cls == "div"))]
            if not qs: return (0, 0, 0)
            q = max(qs, key=lambda q: sc[q][0])
            return (sc[q][0], sc[q][1], sc[q][2] / sc[q][1])
        d, n = best("div"), best("norm")
        fam, gen = tax.get(g, ("?", "?"))
        lab = "div" if d[0] > n[0] else "norm"
        rows.append((fam, gen, g, c, a, b, b - a + 1, lab, d[1], round(d[2], 1), n[1], round(n[2], 1)))
rows.sort()
print("family\tgenus\tgenome\tcontig\tstart\tend\tspan\tlabel\tdiv_aln\tdiv_pid\tnorm_aln\tnorm_pid")
for r in rows:
    print("\t".join(map(str, r)))
